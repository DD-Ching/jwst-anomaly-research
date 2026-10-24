import struct
import sys
import zlib
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import s1063_wslap as sw  # noqa: E402
import spt0546_inputs as sp  # noqa: E402


def _cd_record(name: bytes, crc, csize, usize, offset, zip64=False):
    extra = b""
    if zip64:
        extra = struct.pack("<HHQQQ", 1, 24, usize, csize, offset)
        csize = usize = offset = 0xFFFFFFFF
    return (
        b"PK\x01\x02"
        + struct.pack("<HHHHHH", 20, 20, 0, 8, 0, 0)
        + struct.pack("<IIIHHH", crc, csize, usize, len(name), len(extra), 0)
        + struct.pack("<HHII", 0, 0, 0, offset)
        + name
        + extra
    )


def test_parse_central_directory_plain_and_zip64():
    blob = _cd_record(b"a.fits", 7, 10, 20, 0) + _cd_record(
        b"dir/b.fits", 9, 5_000_000_000, 6_000_000_000, 4_500_000_000, zip64=True
    )
    d = sw.parse_central_directory(blob)
    assert d["a.fits"] == {"method": 8, "crc": 7, "csize": 10, "usize": 20, "offset": 0}
    assert d["dir/b.fits"]["usize"] == 6_000_000_000
    assert d["dir/b.fits"]["csize"] == 5_000_000_000
    assert d["dir/b.fits"]["offset"] == 4_500_000_000
    with pytest.raises(ValueError):
        sw.parse_central_directory(b"junk")


class _Resp:
    def __init__(self, content, status=206):
        self.content, self.status_code = content, status

    def raise_for_status(self):
        pass

    def close(self):
        pass

    def iter_content(self, n):
        for k in range(0, len(self.content), 7):
            yield self.content[k : k + 7]


class _Session:
    def __init__(self, blob, status=206):
        self.blob, self.status = blob, status

    def get(self, url, headers, **kw):
        a, b = headers["Range"].removeprefix("bytes=").split("-")
        return _Resp(self.blob[int(a) : int(b) + 1], self.status)


def _zip_with(payload: bytes):
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    comp = c.compress(payload) + c.flush()
    local = b"PK\x03\x04" + b"\0" * 22 + struct.pack("<HH", 3, 2) + b"m.f" + b"xx"
    blob = b"PREFIX" + local + comp
    entry = {
        "method": 8,
        "crc": zlib.crc32(payload),
        "csize": len(comp),
        "usize": len(payload),
        "offset": 6,
    }
    return blob, entry


def test_fetch_member_inflates_and_checks_crc(tmp_path):
    payload = bytes(range(256)) * 50
    blob, entry = _zip_with(payload)
    sw.fetch_member(_Session(blob), entry, tmp_path / "m.fits")
    assert (tmp_path / "m.fits").read_bytes() == payload
    bad = dict(entry, crc=entry["crc"] ^ 1)
    with pytest.raises(RuntimeError):
        sw.fetch_member(_Session(blob), bad, tmp_path / "n.fits")
    assert not (tmp_path / "n.fits").exists()


class _HeaderOk(_Session):
    """Answers the 30-byte local-header read with 206 and the member stream with 200."""

    def get(self, url, headers, **kw):
        a, b = headers["Range"].removeprefix("bytes=").split("-")
        return _Resp(self.blob[int(a) : int(b) + 1], 206 if int(b) - int(a) == 29 else 200)


def test_fetch_member_refuses_a_non_range_answer(tmp_path):
    blob, entry = _zip_with(b"abc" * 100)
    with pytest.raises(RuntimeError):
        sw.fetch_member(_HeaderOk(blob), entry, tmp_path / "m.fits")


def test_eocd_location_plain_and_zip64():
    plain = b"xx" + b"PK\x05\x06" + b"\0" * 8 + struct.pack("<IIH", 300, 1000, 0)
    assert sw.eocd_location(plain) == (300, 1000)
    loc = b"PK\x06\x07" + struct.pack("<IQI", 0, 123_456_789_012, 1)
    z64 = b"yy" + loc + b"PK\x05\x06" + b"\0" * 8 + struct.pack("<IIH", 0xFFFFFFFF, 0xFFFFFFFF, 0)
    assert sw.eocd_location(z64) == (None, 123_456_789_012)
    with pytest.raises(ValueError):
        sw.eocd_location(b"no record here")
    with pytest.raises(ValueError):
        sw.eocd_location(z64[22:])  # locator cut off


@pytest.mark.network
def test_zenodo_zip_directory_lists_the_alpha_members():
    import requests

    size, d = sw.zip_directory(requests.Session())
    assert size == 1_016_075_135
    for m, crc in zip(sw.MAPS, (0xCC397319, 0x87568E0D), strict=True):
        e = d[sw.ZIP_DIR + m]
        assert e["method"] == 8 and e["usize"] == 141_122_880 and e["crc"] == crc


def test_unit_check_exactly_one_reading():
    def uc(v):
        return sp.unit_check(v, sw.U_MAX_MEDIAN_ARCSEC)

    ok = uc({"A_arcsec": {"1": 0.4, "2": 0.9}, "P_pixels": {"1": 30.0, "2": 40.0}})
    assert ok["passes"] and ok["reading"] == "A_arcsec"
    assert ok["median_arcsec"]["A_arcsec"] == pytest.approx(0.65)
    both = uc({"A_arcsec": {"1": 0.4}, "P_pixels": {"1": 1.0}})
    assert not both["passes"] and both["reading"] is None
    missing = uc({"A_arcsec": {"1": 0.4, "2": None}, "P_pixels": {"1": 5.0}})
    assert not missing["passes"] and np.isnan(missing["median_arcsec"]["A_arcsec"])


def test_spec_images_are_the_twelve_zspec_systems():
    t = sw.spec_images()
    assert len(t) == 36 and len(set(t["system"])) == 12
    assert bool(np.all(t["z_spec"])) and not bool(np.any(t["dup_position"]))


def test_dja_check_matches_field_name():
    html = '<a href="abells1063-grizli-v7.5-fix">x</a><a href="macs0416-grizli-v7.0">y</a>'
    r = sp.dja_check(html, sw.DJA_NAMES)
    assert r["passes"] and r["matches"] == ["abells1063"] and r["fields_listed"] == 2
    with pytest.raises(ValueError):
        sp.dja_check("<html></html>", sw.DJA_NAMES)

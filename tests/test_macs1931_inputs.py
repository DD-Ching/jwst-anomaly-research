"""Offline checks of the MACS J1931 input-check parsers and the committed image list."""

import struct
import sys
import zlib
from pathlib import Path

import pytest
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import macs1931_inputs as mi  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / mi.OUT

ARCS = """# REFERENCE 0
1.1 292.9556015 -26.5741522 0.5 0.5 0. 1.178 25.
1.2 292.9561018 -26.5774723 0.5 0.5 0. 1.178 25.
11.1 292.9525117 -26.5682155 0.5 0.5 0. 4.9 25.
11.2 292.9542430 -26.5677604 0.5 0.5 0. 17.15 25.
c11.3 292.9585968 -26.5683021 0.5 0.5 0. 0. 25.
9.1 292.9555356 -26.5863210 0.5 0.5 0. 0. 25.
"""
MUSE = """# Id RA Dec a b Theta z mag conf AutoXCorMARZ
7 292.9556015 -26.5741522 0.5 0.5 0. 1.1781 25. 3 9.0
8 292.9600000 -26.5700000 0.5 0.5 0. 0.3500 25. 4 9.0
"""
CHIRES = """      Image     ByEye_RMSi [arcsec] Avg_RMSi [arcsec] Bary_RMSi [arcsec]
       1.1             0.3             0.4            0.0
       1.2             0.4             0.3            0.0
"""


def test_parse_arcs_and_counts():
    t = mi.parse_arcs(ARCS)
    assert list(t["used"]) == [True, True, True, True, False, True]
    assert list(t["z_spec"]) == [True, True, False, False, False, False]
    c = mi.list_counts(t)
    assert (c["used_images"], c["used_systems"], c["spec_systems"]) == (5, 3, 1)
    assert c["photo_systems_mixed_z"] == ["11"] and c["photo_images_z0"] == 1
    assert c["candidates_not_used"] == ["c11.3"]
    with pytest.raises(ValueError):
        mi.parse_arcs("1.1 292.9 -26.5 0.5 0.5 0. 1.178\n")
    with pytest.raises(ValueError):
        mi.parse_arcs("A1.1 292.9 -26.5 0.5 0.5 0. 1.178 25.\n")
    with pytest.raises(ValueError):
        mi.parse_arcs("1.1 292.9 -26.5 0.5 0.5 0. 0. 25.\n")  # spectroscopic family without z
    with pytest.raises(ValueError):
        mi.parse_arcs(ARCS + ARCS.splitlines()[1] + "\n")


def test_muse_check():
    t = mi.parse_arcs(ARCS)
    r = mi.muse_check(t, mi.parse_muse(MUSE))
    assert (r["spec_images"], r["images_matched"], r["systems_matched"]) == (2, 1, ["1"])
    assert r["passes"]
    off = mi.muse_check(t, mi.parse_muse(MUSE.replace("1.1781", "1.2000")))
    assert off["systems_unmatched"] == ["1"] and not off["passes"]
    with pytest.raises(ValueError):
        mi.parse_muse("7 292.9 -26.5 0.5\n")


def test_chires_rms():
    r = mi.chires_rms(CHIRES, ["1.1", "1.2", "2.1"])
    assert r["images"] == 2 and r["missing"] == ["2.1"]
    assert r["ByEye_RMSi"] == pytest.approx(0.3536, abs=1e-4) and r["Bary_RMSi"] == 0.0


def _zip_head(bunit="Pixels", method=8):
    h = fits.ImageHDU(data=None).header
    h["NAXIS"], h["NAXIS1"], h["NAXIS2"] = 2, 3001, 3001
    h["BITPIX"] = -64
    h["CTYPE1"], h["CTYPE2"], h["CDELT1"], h["CDELT2"] = (
        "RA---TAN",
        "DEC--TAN",
        -1.1111111e-05,
        1.1111111e-05,
    )
    h["BUNIT"] = bunit
    raw = fits.PrimaryHDU().header.tostring().encode() + h.tostring().encode()
    data = zlib.compress(raw)[2:-4] if method == 8 else raw
    name = b"Alpha_x.fits"
    head = struct.pack(
        "<IHHHHHIIIHH", 0x04034B50, 20, 0, method, 0, 0, 0, len(data), len(raw), len(name), 0
    )
    return head + name + data


@pytest.mark.parametrize("method", [8, 0])
def test_zip_first_member_header_and_unit(method):
    h = mi.zip_first_member_header(_zip_head(method=method))
    assert (h["member"], h["NAXIS1"], h["BUNIT"], h["pixel_arcsec"]) == (
        "Alpha_x.fits",
        3001,
        "Pixels",
        0.04,
    )
    assert not mi.unit_check([h])["agree"]
    assert mi.unit_check([mi.zip_first_member_header(_zip_head("arcsec", method))])["agree"]
    with pytest.raises(ValueError):
        mi.zip_first_member_header(b"\0" * 64)


def test_dja_check():
    assert mi.dja_check("x abell1689-grizli-v7.4 sunrise-grizli-v7.0")["matches"] == []
    assert mi.dja_check("macs1931-grizli-v7.4")["passes"]
    with pytest.raises(ValueError):
        mi.dja_check("<html></html>")


def test_committed_image_list():
    t = Table.read(OUT / mi.IMAGES)
    c = mi.list_counts(t)
    assert (c["used_images"], c["used_systems"], c["spec_images"], c["spec_systems"]) == (
        51,
        19,
        26,
        9,
    )
    assert t.meta["source"].endswith(mi.ARCS)

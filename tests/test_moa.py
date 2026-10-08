"""Offline tests for the MOA-II 9-year adapter on tiny fixture files (D-062)."""

from __future__ import annotations

import gzip
import hashlib
import io
import tarfile

import numpy as np
import pytest

from jwst_anomaly import moa, schema
from jwst_anomaly.signatures import LightCurveSurvey

LC = """\
|       HJD|         flux|cor_flux|    flux_err|obsID|         JD|   fwhm|included|
|    double|       double|  double|      double| long|     double|  float|    char|
|         d|             |        |            |     |          d|       |        |
|       nan|          nan|     nan|         nan|  nan|        nan|    nan|     nan|
  3824.1196  -2352.056396      nan   179.286194   335 3824.119592 1.73912     True
 3825.12153    164.498596      nan   165.900116   339 3825.121425 1.60333     True
 3825.16063   -398.869659      nan   166.603958   340 3825.160525 1.96697    False
 3826.00000     10.000000      nan   100.000000   341 3826.000000 1.90000     True
"""

LC_COR = (
    LC.replace("      nan   179", "     12.0   179")
    .replace("      nan   165", "     13.0   165")
    .replace("      nan   100", "     14.0   100")
)

META = """\
|field|chip|subframe|     id| tag|      x|      y|          ra_j2000|          dec_j2000|pspl_tE|
| long|long|    long|   long|char| double| double|            double|             double| double|
|     |    |        |       |    |    pix|    pix|               deg|                deg|      d|
| null|null|    null|   null|null|   null|   null|              null|               null|   null|
     1    1        0       1 None   12.57  728.99        266.2195875 -33.704547222222224    null
    22    8        5      98 None   10.00   20.00        270.0000000 -25.000000000000000    null
    22    2        1      7  None   11.00   21.00        270.1000000 -25.100000000000000   42.5
"""


def _tar(tmp_path, members: dict[str, str]):
    path = tmp_path / "gb22.tar"
    with tarfile.open(path, "w") as tar:
        for name, text in members.items():
            data = gzip.compress(text.encode())
            info = tarfile.TarInfo(f"exodata/lcurve/gb22/R/8/{name}.ipac.gz")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


def _meta(tmp_path):
    path = tmp_path / "metadata.ipac.tar.gz"
    with tarfile.open(path, "w:gz") as tar:
        data = META.encode()
        info = tarfile.TarInfo("metadata.ipac")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    return path


@pytest.fixture
def field(tmp_path):
    tar = _tar(tmp_path, {"gb22-R-8-5-98": LC, "gb22-R-2-1-7": LC_COR})
    return moa.MoaField(
        "gb22", derived_dir=tmp_path / "derived", tar_path=tar, metadata_path=_meta(tmp_path)
    )


def test_event_ids_round_trip():
    assert moa.parse_event_id("gb22-R-8-5-98") == (22, 8, 5, 98)
    assert moa.make_event_id(22, 8, 5, 98) == "gb22-R-8-5-98"
    with pytest.raises(ValueError):
        moa.parse_event_id("BLG500.01.1")


def test_zero_point_differs_on_chip_2():
    assert moa.mag_to_counts(20.0, 2) == pytest.approx(691.8)
    assert moa.mag_to_counts(20.0, 5) == pytest.approx(1445.0)
    assert moa.counts_to_mag(moa.mag_to_counts(18.3, 7), 7) == pytest.approx(18.3)


def test_parse_lightcurve_drops_the_nan_row_and_reads_included():
    cols = moa.parse_lightcurve(gzip.compress(LC.encode()))
    assert cols["HJD"].size == 4
    assert cols["included"].tolist() == [True, True, False, True]
    assert np.all(np.isnan(cols["cor_flux"]))


def test_select_flux_prefers_detrended_flux_when_complete():
    t, f, sf, kind = moa.select_flux(moa.parse_lightcurve(LC))
    assert kind == "difference" and f.size == 3  # the excluded epoch is dropped
    assert t[0] == pytest.approx(2453824.1196)  # HJD − 2450000 restored
    t, f, sf, kind = moa.select_flux(moa.parse_lightcurve(LC_COR))
    assert kind == "detrended difference"
    assert f.tolist() == [12.0, 13.0, 14.0]


def test_adapter_satisfies_protocol_and_reads_members_in_place(field):
    assert isinstance(field, LightCurveSurvey)
    assert set(field.index()) == {"gb22-R-8-5-98", "gb22-R-2-1-7"}
    lc = field.light_curve("gb22-R-8-5-98")
    assert lc.colnames == list(schema.LIGHT_CURVE_FLUX_COLUMNS)
    assert lc["flux"][0] == pytest.approx(-2352.056396)  # negative difference flux kept
    assert lc.meta["provenance"] == "observed" and "gb22.tar" in lc.meta["source"]
    assert lc.meta["n_excluded"] == 1 and lc.meta["flux_kind"] == "difference"
    with pytest.raises(KeyError):
        field.light_curve("gb22-R-1-1-1")
    assert dict(field.iter_light_curves()).keys() == field.index().keys()
    assert field.efficiency(30.0) is None


def test_events_keep_only_the_field_and_are_cached(field, tmp_path):
    ev = field.events()
    assert sorted(ev["event_id"]) == ["gb22-R-2-1-7", "gb22-R-8-5-98"]
    assert ev.colnames[:3] == ["event_id", "ra", "dec"]
    assert ev.meta["provenance"] == "observed" and ev.meta["source"]
    row = ev[ev["event_id"] == "gb22-R-2-1-7"][0]
    assert row["pspl_tE"] == pytest.approx(42.5)
    assert np.isnan(ev[ev["event_id"] == "gb22-R-8-5-98"][0]["pspl_tE"])
    assert (tmp_path / "derived" / "metadata_gb22.ecsv").exists()


def test_unpinned_field_is_refused(tmp_path):
    with pytest.raises(KeyError):
        moa.MoaField("gb5", raw_dir=tmp_path).tar_path()


def test_star_count_model_brackets_the_published_ratio():
    est, lo, hi = moa.star_count_estimate(18_599)
    assert lo < est < hi
    assert 100 < est / 18_599 < 300


def test_standard_flux_light_curve_keeps_negative_flux_and_drops_bad_rows():
    from jwst_anomaly.signatures import standard_flux_light_curve

    lc = standard_flux_light_curve(
        [3.0, 1.0, 2.0, 4.0],
        [-5.0, 2.0, np.nan, 1.0],
        [1.0, 1.0, 1.0, 0.0],
        "MOA-Red",
        source="fixture",
        time_system="HJD",
        flux_unit="counts",
    )
    assert lc["time"].tolist() == [1.0, 3.0] and lc["flux"].tolist() == [2.0, -5.0]
    assert lc.meta["n_dropped"] == 2 and lc.meta["flux_kind"] == "difference"


def test_fast_parser_equals_the_line_parser_and_reads_column_subsets():
    for text in (LC, LC_COR):
        fast = moa.parse_lightcurve(text)
        slow = moa._parse_lines(text)
        fixed = moa._parse_fixed(text.encode())
        assert fixed is not None  # the fixed-width path is the one taken
        assert fast.keys() == slow.keys() == fixed.keys()
        for k in slow:
            np.testing.assert_array_equal(fast[k], slow[k])
            assert fixed[k].tobytes() == moa._parse_tokens(text.encode())[k].tobytes()
    shifted = LC.replace(" 3825.12153    164", "3825.12153     164")  # token under a bar
    assert moa._parse_fixed(shifted.encode()) is None
    assert moa.parse_lightcurve(shifted)["HJD"].size == 4
    sub = moa.parse_lightcurve(gzip.compress(LC.encode()), columns=("flux", "included"))
    assert set(sub) == {"HJD", "flux", "included"}
    ragged = LC + "  3827.0  1.0\n"  # a short row: the line parser drops it
    assert moa.parse_lightcurve(ragged)["HJD"].size == 4


def _big_tar(tmp_path, n=40, seed=3):
    """A tar like the field tars: directory entries, then gzipped members of varied size."""
    rng = np.random.default_rng(seed)
    path = tmp_path / "gb21.tar"
    with tarfile.open(path, "w", format=tarfile.USTAR_FORMAT) as tar:
        for d in ("exodata/lcurve/gb21", "exodata/lcurve/gb21/R"):
            info = tarfile.TarInfo(d)
            info.type = tarfile.DIRTYPE
            tar.addfile(info)
        for i in range(n):
            data = gzip.compress(rng.bytes(int(rng.integers(100, 5000))))  # incompressible
            info = tarfile.TarInfo(f"exodata/lcurve/gb21/R/3/gb21-R-3-0-{i}.ipac.gz")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


@pytest.mark.parametrize("seg", [512, 1536, 4096, 10240, 1 << 20])
def test_range_segments_cover_every_member_exactly_once(tmp_path, seg):
    from jwst_anomaly import moa_stream

    path = _big_tar(tmp_path)
    with tarfile.open(path) as tar:
        truth = {
            m.name.rsplit("/", 1)[-1][: -len(".ipac.gz")]: (m.offset_data, m.size)
            for m in tar
            if m.isfile()
        }
    total = path.stat().st_size
    reader = moa_stream.RangeReader(path=path)
    got, shas = [], []
    raw = path.read_bytes()
    for a, b in moa_stream.segments(0, total, seg):
        members, sha = moa_stream.read_segment_hashed(reader, a, b, total)
        got.extend(members)
        shas.append(sha)
        assert sha == hashlib.sha256(raw[a:b]).hexdigest()  # content pin of exactly [a, b)
    assert len(moa_stream.range_digest(shas)) == 64
    ids = [g[0] for g in got]
    assert len(ids) == len(set(ids)) == len(truth)
    raw = path.read_bytes()
    for eid, off, size, data in got:
        assert truth[eid] == (off, size)
        assert data == raw[off : off + size]
    some = [(e, *truth[e]) for e in list(truth)[:5]]
    fetched = moa_stream.fetch_members(reader, some)
    assert all(fetched[e] == raw[o : o + s] for e, o, s in some)


def test_tar_header_rejects_data_blocks():
    assert moa.tar_header(b"\0" * 512) is None
    assert moa.tar_header(bytes(range(256)) * 2) is None


def test_field_urls_and_sizes():
    assert moa.tar_url("gb5").endswith("/bulk/gb5.tar") and moa.tar_url(22) in moa.FILES
    assert moa.object_url("gb22-R-6-1-29").endswith("/MOA/gb22/R/6/gb22-R-6-1-29.ipac")
    assert set(moa.TAR_BYTES) == set(moa.CUT0_PER_FIELD) and moa.TAR_BYTES[22] == 3510138880


def test_one_metadata_pass_caches_several_fields(tmp_path):
    from astropy.table import Table

    out = moa.write_metadata_caches(_meta(tmp_path), tmp_path / "d", fields=[1, 22])
    assert [p.name for p in out] == ["metadata_gb1.ecsv", "metadata_gb22.ecsv"]
    assert len(Table.read(out[1])) == 2 and len(Table.read(out[0])) == 1


def test_fixed_parser_falls_back_on_an_unparseable_token(recwarn):
    bad = LC.replace("   335 ", "  null ")  # same width, not a number
    assert len(bad) == len(LC) and bad != LC
    assert moa._parse_fixed(bad.encode()) is None
    assert not [w for w in recwarn if issubclass(w.category, DeprecationWarning)]


def test_long_name_headers_are_refused(tmp_path):
    path = tmp_path / "gb21.tar"
    with tarfile.open(path, "w", format=tarfile.GNU_FORMAT) as tar:
        data = gzip.compress(b"x")
        info = tarfile.TarInfo("exodata/" + "d" * 120 + "/gb21-R-3-0-1.ipac.gz")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    raw = path.read_bytes()
    with pytest.raises(ValueError, match="name extension"):
        list(moa.walk_members(memoryview(raw), 0, 0, len(raw)))


def test_split_metadata_tables_carry_provenance(tmp_path):
    t = moa.split_metadata(META.splitlines(keepends=True), [22])[22]
    assert t.meta["provenance"] == schema.Provenance.OBSERVED.value
    assert "metadata" in t.meta["source"]

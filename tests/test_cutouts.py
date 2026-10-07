"""Tests for cutouts.make_cutouts on a synthetic i2d-like image (plus live MAST/S3 checks)."""

from __future__ import annotations

from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import MaskedColumn, Table
from astropy.wcs import WCS

from jwst_anomaly import cutouts, schema
from jwst_anomaly.cutouts import make_cutouts

SCALE = 0.05  # arcsec / pixel -> 3" cutouts are 61 x 61 pixels
NY, NX = 200, 300

# 0-based (x, y) pixel positions of synthetic sources and what each one exercises.
SOURCES = {
    "interior": (150.0, 100.0),
    "array_edge": (292.0, 60.0),  # 8 px from the right edge of the array
    "footprint_edge": (40.0, 150.0),  # box reaches the no-coverage strip x < 20
    "low_weight": (220.0, 160.0),  # inside a WHT patch at 0.2 x the typical weight
    "nan_core": (100.0, 40.0),  # SCI NaN / WHT 0 at the core (e.g. saturation)
}


def _wcs() -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.crpix = [NX / 2 + 1, NY / 2 + 1]
    w.wcs.cdelt = [-SCALE / 3600, SCALE / 3600]
    theta = np.deg2rad(30.0)
    w.wcs.pc = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    return w


def _write_i2d(path, *, filt="F200W", pupil="CLEAR", with_wht=True) -> WCS:
    """Synthetic level-3 image: Gaussian sources on noise, a no-coverage strip, a shallow patch."""
    rng = np.random.default_rng(42)
    sci = rng.normal(0.0, 0.01, (NY, NX)).astype(np.float32)
    yy, xx = np.mgrid[0:NY, 0:NX]
    for x, y in SOURCES.values():
        sci += np.float32(5.0) * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 2.0**2))
    wht = np.ones((NY, NX), dtype=np.float32)
    wht[120:, 180:260] = 0.2
    x0, y0 = SOURCES["nan_core"]
    core = (xx - x0) ** 2 + (yy - y0) ** 2 <= 2.0**2
    wht[core] = 0.0
    sci[core] = np.nan
    sci[:, :20] = np.nan
    wht[:, :20] = 0.0
    w = _wcs()
    primary = fits.PrimaryHDU()
    primary.header["TELESCOP"] = "JWST"
    primary.header["INSTRUME"] = "NIRCAM"
    if filt:
        primary.header["FILTER"] = filt
    primary.header["PUPIL"] = pupil
    sci_hdr = w.to_header()
    sci_hdr["BUNIT"] = "MJy/sr"
    hdus = [
        primary,
        fits.ImageHDU(sci, header=sci_hdr, name="SCI"),
        fits.ImageHDU(np.full((NY, NX), 0.01, np.float32), header=w.to_header(), name="ERR"),
    ]
    hdus.append(fits.ImageHDU(np.ones((NY, NX), np.int32), name="CON"))  # int, like i2d CON
    if with_wht:
        hdus.append(fits.ImageHDU(wht, header=w.to_header(), name="WHT"))
    fits.HDUList(hdus).writeto(path)
    return w


def _targets(w: WCS, names=None) -> Table:
    names = list(names or SOURCES)
    sky = w.pixel_to_world([SOURCES[n][0] for n in names], [SOURCES[n][1] for n in names])
    return Table({"source_uid": names, "ra": sky.ra.deg, "dec": sky.dec.deg})


@pytest.fixture
def image(tmp_path):
    path = tmp_path / "jw99999-o001_t001_nircam_clear-f200w_i2d.fits"
    return str(path), _write_i2d(path)


def _row(table: Table, uid: str):
    return table[list(table["source_uid"]).index(uid)]


def test_contract_and_metadata(image, tmp_path):
    path, w = image
    out = make_cutouts(path, _targets(w), out_dir=tmp_path / "out")
    schema.validate(out, schema.CUTOUT_COLUMNS)
    assert out.meta["provenance"] == schema.Provenance.OBSERVED.value
    assert out.meta["source"] == path
    assert list(out["source_uid"]) == list(SOURCES)
    assert set(out["band"]) == {"F200W"}
    assert out.meta["size_pix"] == [61, 61]
    assert out.meta["weight_ref"] == pytest.approx(1.0)
    assert set(out.meta["quality_flags"]) >= {"ok", "outside", "edge", "low_weight"}
    assert out.meta["column_provenance"]["quality_flag"] == schema.Provenance.DERIVED.value
    out.meta["quality_flags"]["ok"] = "edited"  # meta holds a copy, not the module constant
    assert cutouts.QUALITY_FLAGS["ok"] != "edited"


def test_interior_cutout_is_centered_and_written(image, tmp_path):
    path, w = image
    targets = _targets(w, ["interior"])
    row = make_cutouts(path, targets, out_dir=tmp_path / "out")[0]
    assert row["quality_flag"] == "ok"
    assert not row["on_edge"] and row["frac_nan"] == 0.0
    assert row["wht_rel"] == pytest.approx(1.0)
    assert (row["x"], row["y"]) == pytest.approx(SOURCES["interior"], abs=1e-6)
    expected = tmp_path / "out" / "jw99999-o001_t001_nircam_clear-f200w_i2d" / "interior.fits"
    assert row["path"] == expected.as_posix()
    with fits.open(row["path"]) as hdul:
        assert [h.name for h in hdul] == ["PRIMARY", "SCI", "WHT"]
        assert hdul[0].header["SRCUID"] == "interior"
        assert hdul[0].header["FILTER"] == "F200W"
        assert hdul["SCI"].header["BUNIT"] == "MJy/sr"
        sci = hdul["SCI"].data
        assert sci.shape == (61, 61)
        assert np.unravel_index(np.argmax(sci), sci.shape) == (30, 30)
        cx, cy = WCS(hdul["SCI"].header).world_to_pixel_values(targets["ra"][0], targets["dec"][0])
        assert (float(cx), float(cy)) == pytest.approx((30.0, 30.0), abs=1e-6)


@pytest.mark.parametrize(
    "name, flag",
    [
        ("array_edge", "edge"),
        ("footprint_edge", "edge"),
        ("low_weight", "low_weight"),
        ("nan_core", "nan_center"),
    ],
)
def test_quality_flags(image, tmp_path, name, flag):
    path, w = image
    row = make_cutouts(path, _targets(w, [name]), out_dir=tmp_path)[0]
    assert flag in row["quality_flag"].split(",")
    assert row["path"]
    if flag == "edge":
        assert row["on_edge"] and row["frac_nan"] > 0
    else:
        assert not row["on_edge"]
    if flag == "low_weight":
        assert row["wht_rel"] == pytest.approx(0.2)


def test_target_off_image_is_reported_not_written(image, tmp_path):
    path, w = image
    targets = Table({"source_uid": ["far", "nan_pos"], "ra": [151.0, np.nan], "dec": [3.0, 2.0]})
    out = make_cutouts(path, targets, out_dir=tmp_path)
    assert list(out["quality_flag"]) == ["outside", "outside"]
    assert list(out["path"]) == ["", ""]
    assert list(out["frac_nan"]) == [1.0, 1.0]
    assert [p.name for p in tmp_path.rglob("*.fits")] == [Path(path).name]  # nothing written


def test_inside_array_but_no_coverage_is_outside(image, tmp_path):
    path, w = image
    sky = w.pixel_to_world(5.0, 100.0)  # in the NaN strip; the 61 px box still reaches data
    targets = Table({"source_uid": ["strip"], "ra": [sky.ra.deg], "dec": [sky.dec.deg]})
    row = make_cutouts(path, targets, size_arcsec=0.5, out_dir=tmp_path)[0]
    assert row["quality_flag"] == "outside"


def test_extensions_band_override_and_quantity_columns(image, tmp_path):
    path, w = image
    targets = _targets(w, ["interior"])
    targets["ra"].unit = u.deg
    targets["dec"] = (targets["dec"] * u.deg).to(u.arcmin)
    out = make_cutouts(path, targets, out_dir=tmp_path, extensions=("SCI", "ERR"), band="X1")
    assert out["band"][0] == "X1" and out["quality_flag"][0] == "ok"
    with fits.open(out["path"][0]) as hdul:
        assert [h.name for h in hdul] == ["PRIMARY", "SCI", "ERR"]


def test_pupil_filter_names_the_band(tmp_path):
    path = tmp_path / "pupil_i2d.fits"
    w = _write_i2d(path, filt="F444W", pupil="F470N")
    out = make_cutouts(str(path), _targets(w, ["interior"]), out_dir=tmp_path)
    assert out["band"][0] == "F470N"


def test_missing_filter_requires_band(tmp_path):
    path = tmp_path / "nofilter_i2d.fits"
    w = _write_i2d(path, filt="", pupil="")
    with pytest.raises(ValueError, match="FILTER"):
        make_cutouts(str(path), _targets(w, ["interior"]), out_dir=tmp_path)


def test_without_wht_flags_from_sci_only(tmp_path):
    path = tmp_path / "nowht_i2d.fits"
    w = _write_i2d(path, with_wht=False)
    out = make_cutouts(str(path), _targets(w, ["interior", "low_weight"]), out_dir=tmp_path)
    assert list(out["quality_flag"]) == ["ok", "ok"]
    assert np.isnan(out["wht_rel"]).all()


def test_default_out_dir_uses_outputs_env(image, tmp_path, monkeypatch):
    path, w = image
    monkeypatch.setenv("JWST_ANOMALY_OUTPUTS", str(tmp_path / "outputs"))
    row = make_cutouts(path, _targets(w, ["interior"]))[0]
    assert row["path"].startswith((tmp_path / "outputs" / "cutouts").as_posix())


def test_empty_targets(image, tmp_path):
    path, _ = image
    empty = Table({"source_uid": np.array([], str), "ra": [], "dec": []})
    out = make_cutouts(path, empty, out_dir=tmp_path)
    assert len(out) == 0
    schema.validate(out, schema.CUTOUT_COLUMNS)


def test_row_strip_section_matches_data_and_reads_full_rows(image, monkeypatch):
    path, _ = image
    with fits.open(path) as hdul:
        hdu = hdul["SCI"]
        section = cutouts._RowStripSection(hdu)
        reads = []
        original = fits.Section.__getitem__

        def spy(self, key):
            reads.append(key)
            return original(self, key)

        monkeypatch.setattr(fits.Section, "__getitem__", spy)
        got = section[10:20, 30:50]
        np.testing.assert_array_equal(got, hdu.data[10:20, 30:50])
        assert reads == [(slice(10, 20), slice(None))]


def test_remote_code_path_with_fsspec_memory(image, tmp_path):
    fsspec = pytest.importorskip("fsspec")
    path, w = image
    fs = fsspec.filesystem("memory")
    with open(path, "rb") as f:
        fs.pipe("/unit4/jw99999_i2d.fits", f.read())
    out = make_cutouts(
        "memory://unit4/jw99999_i2d.fits", _targets(w, ["interior"]), out_dir=tmp_path
    )
    assert out["quality_flag"][0] == "ok"
    assert out["path"][0].endswith("jw99999_i2d/interior.fits")


@pytest.mark.parametrize(
    "uri, expected",
    [
        ("s3://stpubdata/jwst/public/jw02736/L3/t/o001/jw02736_i2d.fits", "jw02736_i2d"),
        (
            "mast:JWST/product/jw02736-o002_t001_miri_f770w_i2d.fits",
            "jw02736-o002_t001_miri_f770w_i2d",
        ),
        (r"C:\data\x_i2d.fits", "x_i2d"),
        ("/data/run#1/img3_i2d.fits", "img3_i2d"),
    ],
)
def test_image_stem(uri, expected):
    assert cutouts._image_stem(cutouts._resolve_uri(uri)) == expected


def test_requested_extensions_are_validated(image, tmp_path):
    path, w = image
    targets = _targets(w, ["interior"])
    with pytest.raises(ValueError, match="ERRR not found"):
        make_cutouts(path, targets, out_dir=tmp_path, extensions=("SCI", "ERRR"))
    with pytest.raises(ValueError, match="CON must be a 2-D float image"):
        make_cutouts(path, targets, out_dir=tmp_path, extensions=("CON",))


def test_uids_never_share_a_file(image, tmp_path):
    path, w = image
    sky = w.pixel_to_world(*SOURCES["interior"])
    uids = ["a/b", "a_b", "A_B", "a_b", "\u6e90-1"]
    targets = Table({"source_uid": uids, "ra": [sky.ra.deg] * 5, "dec": [sky.dec.deg] * 5})
    out = make_cutouts(path, targets, out_dir=tmp_path)
    assert len({p.lower() for p in out["path"]}) == len(uids)
    for uid, p in zip(uids, out["path"], strict=True):
        expected = uid.encode("ascii", "backslashreplace").decode()
        assert fits.getheader(p, 0)["SRCUID"] == expected


def test_masked_coordinates_are_outside(image, tmp_path):
    path, w = image
    targets = _targets(w, ["interior", "low_weight"])
    targets["ra"] = MaskedColumn(targets["ra"], mask=[True, False], unit=u.deg)
    out = make_cutouts(path, targets, out_dir=tmp_path)
    assert list(out["quality_flag"]) == ["outside", "low_weight"]


@pytest.mark.filterwarnings("error")
def test_non_ascii_long_paths_and_file_uri(tmp_path):
    folder = tmp_path / "\u8cc7\u6599 with a fairly long directory name for header cards"
    folder.mkdir()
    path = folder / "img_i2d.fits"
    w = _write_i2d(path)
    out = make_cutouts(path.as_uri(), _targets(w, ["interior"]), out_dir=tmp_path / "out")
    assert out["quality_flag"][0] == "ok"
    assert fits.getheader(out["path"][0], 0)["ORIGURI"].endswith("img_i2d.fits")
    assert cutouts._resolve_uri(path.as_uri()) == str(path)


def test_weight_reference_is_lazy_and_skips_padding(image, tmp_path, monkeypatch):
    path, _ = image
    far = Table({"source_uid": ["far"], "ra": [151.0], "dec": [3.0]})
    assert np.isnan(make_cutouts(path, far, out_dir=tmp_path).meta["weight_ref"])
    wht = np.ones((10, 4))
    wht[0] = wht[-1] = 0.0  # zero-weight padding rows at the array edges
    wht[2] = 5.0
    assert cutouts._reference_weight(wht, 2) == 1.0  # small array: read whole
    monkeypatch.setattr(cutouts, "_WHOLE_WHT_BYTES", 0)
    assert cutouts._reference_weight(wht, 2) == 3.0  # rows 2 and 7, never the padding


# Live checks. Position: label 2596 of jw02736-o001_t001_nircam_clear-f200w_cat.ecsv
# (sky_centroid, MAST, retrieved 2026-10-07); it lies in both the F200W and F770W images.
LIVE_TARGET = Table(
    {"source_uid": ["2736-f200w-2596"], "ra": [110.85801683111302], "dec": [-73.46425234249979]}
)


@pytest.mark.network
def test_live_s3_nircam_byte_range_cutout(tmp_path):
    pytest.importorskip("s3fs")
    uri = (
        "s3://stpubdata/jwst/public/jw02736/L3/t/o001/jw02736-o001_t001_nircam_clear-f200w_i2d.fits"
    )
    out = make_cutouts(uri, LIVE_TARGET, out_dir=tmp_path, n_weight_rows=2)
    row = out[0]
    assert row["band"] == "F200W"
    assert row["quality_flag"] != "outside" and row["frac_nan"] < 0.5
    assert out.meta["fetch"]["bytes"] < 20e6  # the file is 1.76 GB
    sci = fits.getdata(row["path"], extname="SCI")
    peak = np.unravel_index(np.nanargmax(sci), sci.shape)
    assert abs(peak[0] - sci.shape[0] // 2) <= 2 and abs(peak[1] - sci.shape[1] // 2) <= 2


@pytest.mark.network
def test_live_mast_https_miri_cutout(tmp_path):
    pytest.importorskip("aiohttp")
    uri = "mast:JWST/product/jw02736-o002_t001_miri_f770w_i2d.fits"
    out = make_cutouts(uri, LIVE_TARGET, out_dir=tmp_path)
    assert out["band"][0] == "F770W"
    assert out["quality_flag"][0] != "outside"
    assert out.meta["fetch"]["bytes"] < 36_083_520  # less than the whole file


def _spiky_star(n=97, rot_deg=0.0, nan_core=False):
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    c = (n - 1) / 2
    img = 50.0 * np.exp(-((xx - c) ** 2 + (yy - c) ** 2) / (2 * 1.5**2))
    r = np.hypot(xx - c, yy - c)
    ang = np.arctan2(yy - c, xx - c)
    for k in range(6):  # six rays 60 deg apart, fading with radius
        a = np.deg2rad(rot_deg + 60.0 * k)
        d_perp = np.abs(-(xx - c) * np.sin(a) + (yy - c) * np.cos(a))
        along = np.cos(ang - a) > 0
        img += np.where(along, 2.0 * np.exp(-(d_perp**2) / 2.0) / (1 + r / 10), 0.0)
    img += np.random.default_rng(0).normal(0, 0.01, img.shape)
    if nan_core:
        img[r <= 2.5] = np.nan
    return img, c


def _elongated_galaxy(n=97, q=0.3, pa_deg=25.0):
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    c = (n - 1) / 2
    a = np.deg2rad(pa_deg)
    u = (xx - c) * np.cos(a) + (yy - c) * np.sin(a)
    v = -(xx - c) * np.sin(a) + (yy - c) * np.cos(a)
    img = 20.0 * np.exp(-(u**2 / (2 * 8.0**2) + v**2 / (2 * (8.0 * q) ** 2)))
    return img + np.random.default_rng(1).normal(0, 0.01, img.shape), c


def test_spike_statistic_separates_spiky_stars_from_elongated_galaxies():
    for rot in (0.0, 17.0, 41.0):  # orientation-free
        star, c = _spiky_star(rot_deg=rot)
        assert cutouts.spike_statistic(star, c, c, 6, 25) > 5
    gal, c = _elongated_galaxy()
    assert cutouts.spike_statistic(gal, c, c, 6, 25) < 2
    assert np.isnan(cutouts.spike_statistic(gal, c, c, 6, 7))  # fewer than three rings


def test_peak_near_finds_a_saturated_core_and_an_offset_star():
    star, c = _spiky_star(nan_core=True)
    assert cutouts._peak_near(star, c + 4, c - 3, 10) == pytest.approx((c, c), abs=1.0)
    s6 = cutouts.spike_statistic(star, *cutouts._peak_near(star, c + 4, c - 3, 10), 6, 25)
    assert s6 > 5


def test_spike_flag_is_nircam_only_and_recorded(tmp_path, monkeypatch):
    monkeypatch.setattr(cutouts, "spike_statistic", lambda *a, **k: 5.0)
    path = tmp_path / "jw99999-o001_t001_nircam_clear-f200w_i2d.fits"
    w = _write_i2d(path)
    t = cutouts.make_cutouts(
        str(path), _targets(w, ["interior"]), out_dir=tmp_path / "c", spike_radii_arcsec=(0.2, 0.8)
    )
    assert t["spike_s6"][0] == 5.0 and t["quality_flag"][0] == "spikes"
    assert t.meta["spike"]["applied"] and t.meta["spike"]["provenance"] == "assumption"
    off = cutouts.make_cutouts(str(path), _targets(w, ["interior"]), out_dir=tmp_path / "d")
    assert np.isnan(off["spike_s6"][0]) and off.meta["spike"] is None
    with fits.open(path, mode="update") as h:
        h[0].header["INSTRUME"] = "MIRI"
    miri = cutouts.make_cutouts(
        str(path), _targets(w, ["interior"]), out_dir=tmp_path / "e", spike_radii_arcsec=(0.2, 0.8)
    )
    assert np.isnan(miri["spike_s6"][0]) and "spikes" not in miri["quality_flag"][0]
    assert miri.meta["spike"]["applied"] is False


def test_host_ratio_separates_a_bare_star_from_a_nucleus_in_a_galaxy():
    star, c = _spiky_star()
    galaxy, _ = _elongated_galaxy(q=0.8)
    nucleus = galaxy + star  # a point-like nucleus inside a galaxy
    bare = cutouts.host_ratio(star, c, c, 10, 19)
    hosted = cutouts.host_ratio(nucleus, c, c, 10, 19)
    assert bare < 0.01 < hosted
    assert np.isnan(cutouts.host_ratio(np.full((40, 40), np.nan), 20, 20, 5, 10))

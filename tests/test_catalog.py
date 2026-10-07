"""Tests for catalog ingestion and cross-band merging (bootstrap unit 2).

Fixtures ``tests/data/catalog_jw02736-o001_t001_nircam_clear-{f200w,f444w}_cat.ecsv`` are
truncated real JWST pipeline catalogs: the full ECSV header of the MAST products
``jw02736-o001_t001_nircam_clear-f200w_cat.ecsv`` / ``...-f444w_cat.ecsv`` (program 2736,
jwst 2.0.1, retrieved 2026-10-07) plus only the rows within 5 arcsec of
RA 110.685719, Dec -73.470661 (10 and 8 rows). Other tests use synthetic tables.
"""

from __future__ import annotations

import hashlib
import urllib.request
import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Table
from astropy.wcs import WCS

from jwst_anomaly import catalog, schema

DATA = Path(__file__).parent / "data"
F200W_FIXTURE = DATA / "catalog_jw02736-o001_t001_nircam_clear-f200w_cat.ecsv"
F444W_FIXTURE = DATA / "catalog_jw02736-o001_t001_nircam_clear-f444w_cat.ecsv"
BASE = SkyCoord(150.0 * u.deg, 2.0 * u.deg)
FIELD = "jw09999-o001_t001_nircam"


def _pos(dx: float, dy: float) -> SkyCoord:
    """Position offset by (dx, dy) arcsec (east, north) from BASE."""
    return BASE.spherical_offsets_by(dx * u.arcsec, dy * u.arcsec)


def _band(band: str, offsets, labels=None, *, jwst: str = "2.0.1") -> Table:
    """Synthetic single-band catalog shaped like ``load_pipeline_catalog`` output."""
    coords = [_pos(dx, dy) for dx, dy in offsets]
    n = len(coords)
    labels = list(labels) if labels is not None else list(range(1, n + 1))
    t = Table(
        {
            "label": np.array(labels, dtype=np.int32),
            "ra": [c.ra.deg for c in coords],
            "dec": [c.dec.deg for c in coords],
            "aper50_abmag": np.linspace(25.0, 26.0, n),
            "is_extended": np.ones(n, dtype=bool),
            "roundness": np.zeros(n, dtype=np.float32),
        }
    )
    t["ra"].unit = t["dec"].unit = u.deg
    t.meta.update(
        provenance="observed",
        source=f"jw09999-o001_t001_nircam_clear-{band.lower()}_cat.ecsv",
        band=band,
        obs_id=f"{FIELD}_clear-{band.lower()}",
        jwst_version=jwst,
        photutils_version="2.3.0",
    )
    return t


def _write_pipeline_like_ecsv(path: Path, scale_arcsec: float = 0.03) -> Table:
    """Write a small ECSV with the pipeline's structure (SkyCoord mixin, units, version meta)."""
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [BASE.ra.deg, BASE.dec.deg]
    wcs.wcs.crpix = [500.0, 500.0]
    wcs.wcs.cdelt = [-scale_arcsec / 3600, scale_arcsec / 3600]
    x = np.array([10.0, 900.0, 450.5, 120.0, 800.0])
    y = np.array([20.0, 50.0, 700.0, 950.0, 610.0])
    t = Table()
    t["label"] = np.arange(1, 6, dtype=np.int32)
    t["xcentroid"], t["ycentroid"] = x, y
    t["sky_centroid"] = wcs.pixel_to_world(x, y)
    t["aper50_flux"] = np.array([1.0, 2.0, -0.5, 3.0, 4.0]) * 1e-8 * u.Jy
    t["aper50_abmag"] = [28.9, 28.1, np.nan, 27.7, 27.4]
    t["is_extended"] = [False, True, False, True, True]
    t["nn_dist"] = np.ones(5) * u.pix
    t.meta["version"] = {"jwst": "9.9.9", "photutils": "8.8.8", "astropy": "7.2.0"}
    t.meta["abvega_offset"] = 1.6880634
    t.write(path, format="ascii.ecsv")
    return t


# --- load_pipeline_catalog -------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, obs_id, band",
    [
        (
            "jw02736-o001_t001_nircam_clear-f200w_cat.ecsv",
            "jw02736-o001_t001_nircam_clear-f200w",
            "F200W",
        ),
        ("jw02736-o002_t001_miri_f770w_cat.ecsv", "jw02736-o002_t001_miri_f770w", "F770W"),
        (
            "jw01345-c1001_t021_nircam_f150w2-clear_cat.ecsv",
            "jw01345-c1001_t021_nircam_f150w2-clear",
            "F150W2",
        ),
        (
            "catalog_jw02736-o001_t001_nircam_clear-f444w_cat.ecsv",
            "jw02736-o001_t001_nircam_clear-f444w",
            "F444W",
        ),
        (
            "JW02736-O001_T001_NIRCAM_CLEAR-F200W_CAT.ECSV",
            "jw02736-o001_t001_nircam_clear-f200w",
            "F200W",
        ),
        (
            "jw01234-o001_t001_nircam_f444w-f470n_cat.ecsv",
            "jw01234-o001_t001_nircam_f444w-f470n",
            None,
        ),
        ("my_catalog.ecsv", "", None),
    ],
)
def test_parse_catalog_name(name, obs_id, band):
    assert catalog._parse_catalog_name(name) == (obs_id, band)


def test_load_pipeline_like_ecsv(tmp_path):
    path = tmp_path / "jw09999-o001_t001_nircam_clear-f150w_cat.ecsv"
    written = _write_pipeline_like_ecsv(path, scale_arcsec=0.03)
    t = catalog.load_pipeline_catalog(path)

    assert schema.validate(t, schema.BAND_CATALOG_COLUMNS) is t
    assert t.colnames[:3] == ["label", "ra", "dec"]
    assert t["ra"].unit == u.deg
    np.testing.assert_allclose(t["ra"], written["sky_centroid"].ra.deg)
    np.testing.assert_allclose(t["dec"], written["sky_centroid"].dec.deg)
    assert isinstance(t["sky_centroid"], SkyCoord)  # pipeline columns are kept
    assert np.isnan(t["aper50_abmag"][2])
    m = t.meta
    assert m["provenance"] == "observed"
    assert m["source"] == m["input_file"] == path.name
    assert (m["band"], m["obs_id"]) == ("F150W", "jw09999-o001_t001_nircam_clear-f150w")
    assert (m["jwst_version"], m["photutils_version"]) == ("9.9.9", "8.8.8")
    assert m["input_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert m["pixel_scale_arcsec"] == pytest.approx(0.03, rel=1e-4)
    assert m["abvega_offset"] == pytest.approx(1.6880634)  # pipeline meta kept


def test_load_band_override_and_errors(tmp_path):
    path = tmp_path / "renamed.ecsv"
    _write_pipeline_like_ecsv(path)
    with pytest.raises(ValueError, match="band="):
        catalog.load_pipeline_catalog(path)
    t = catalog.load_pipeline_catalog(path, band=" f200w ")
    assert (t.meta["band"], t.meta["obs_id"]) == ("F200W", "")
    with pytest.raises(ValueError, match="letters and digits"):
        catalog.load_pipeline_catalog(path, band="F444W/F470N")

    named = tmp_path / "jw09999-o001_t001_nircam_clear-f200w_cat.ecsv"
    _write_pipeline_like_ecsv(named)
    assert catalog.load_pipeline_catalog(named, band="F200W").meta["band"] == "F200W"
    with pytest.raises(ValueError, match="file name says F200W"):
        catalog.load_pipeline_catalog(named, band="F444W")

    Table({"label": [1], "x": [0.0]}).write(named, format="ascii.ecsv", overwrite=True)
    with pytest.raises(ValueError, match="sky_centroid"):
        catalog.load_pipeline_catalog(named)


def test_load_real_fixture():
    t = catalog.load_pipeline_catalog(F200W_FIXTURE)
    assert len(t) == 10
    assert len(t.colnames) == 56  # 54 pipeline columns (5 SkyCoord mixins) + ra, dec
    assert t.meta["band"] == "F200W"
    assert t.meta["obs_id"] == "jw02736-o001_t001_nircam_clear-f200w"
    assert (t.meta["jwst_version"], t.meta["photutils_version"]) == ("2.0.1", "2.3.0")
    assert t.meta["pixel_scale_arcsec"] == pytest.approx(0.03123, abs=1e-4)  # NIRCam SW mosaic
    assert list(t.meta["aperture_params"]["aperture_ee"]) == [30, 50, 70]
    np.testing.assert_allclose(t["ra"], t["sky_centroid"].ra.deg)
    assert t["aper50_flux"].unit == u.Jy
    lw = catalog.load_pipeline_catalog(F444W_FIXTURE)
    assert lw.meta["pixel_scale_arcsec"] == pytest.approx(0.06291, abs=1e-4)  # NIRCam LW mosaic


# --- merge_bands ------------------------------------------------------------------------------


def _three_band_catalogs() -> dict[str, Table]:
    return {
        "F444W": _band("F444W", [(0.02, 0.0), (4.03, 0.0), (10.0, 0.0)], labels=[11, 12, 13]),
        "F200W": _band("F200W", [(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)]),
        "F090W": _band("F090W", [(0.01, 0.01), (2.01, 0.0), (10.02, 0.0)], labels=[6, 7, 8]),
    }


def test_merge_union():
    m = catalog.merge_bands(_three_band_catalogs(), ref_band="F200W", radius_arcsec=0.1)

    assert schema.validate(m, schema.SOURCE_COLUMNS) is m
    assert m.meta["provenance"] == "derived"
    assert m.colnames[:7] == [
        "source_uid",
        "ra",
        "dec",
        "ref_band",
        "n_bands",
        "anchor_band",
        "ref_nn_sep_arcsec",
    ]
    assert list(m["source_uid"]) == [
        f"{FIELD}_f200w_1",
        f"{FIELD}_f200w_2",
        f"{FIELD}_f200w_3",
        f"{FIELD}_f090w_8",  # no F200W counterpart: anchored on the first band that has it
    ]
    assert list(m["ref_band"]) == ["F200W"] * 4
    assert list(m["anchor_band"]) == ["F200W", "F200W", "F200W", "F090W"]
    assert list(m["n_bands"]) == [3, 2, 2, 2]
    assert list(m["f444w_label"]) == [11, -1, 12, 13]  # F444W #13 matched the F090W-anchored row
    assert list(m["f090w_label"]) == [6, 7, -1, 8]
    assert list(m["f444w_detected"]) == [True, False, True, True]
    assert list(m["f200w_detected"]) == [True, True, True, False]

    # Row position = anchoring detection; separations are measured from it.
    f090w_8 = _pos(10.02, 0.0)
    assert (m["ra"][3], m["dec"][3]) == pytest.approx((f090w_8.ra.deg, f090w_8.dec.deg))
    np.testing.assert_allclose(m["f200w_sep_arcsec"][:3], 0.0)
    np.testing.assert_allclose(m["f444w_sep_arcsec"][[0, 2, 3]], [0.02, 0.03, 0.02], rtol=1e-3)
    assert np.isnan(m["f444w_sep_arcsec"][1]) and np.isnan(m["f200w_sep_arcsec"][3])

    # Non-detections: NaN floats, -1 ints, False bools; units and dtypes preserved.
    assert np.isnan(m["f444w_aper50_abmag"][1]) and np.isnan(m["f444w_ra"][1])
    assert not m["f444w_is_extended"][1]
    assert m["f444w_roundness"].dtype == np.float32
    assert m["f444w_ra"].unit == u.deg and m["f444w_sep_arcsec"].unit == u.arcsec

    # Nearest *other* F200W detection: 2" for the first ref row, ~6" for the F090W-only row.
    assert m["ref_nn_sep_arcsec"][0] == pytest.approx(2.0, rel=1e-3)
    assert m["ref_nn_sep_arcsec"][3] == pytest.approx(6.02, rel=1e-3)

    stats = m.meta["match_stats"]
    assert m.meta["match"]["band_order"] == ["F200W", "F090W", "F444W"]
    assert m.meta["bands"] == ["F090W", "F200W", "F444W"]
    assert (stats["F090W"]["n_matched"], stats["F090W"]["n_unmatched"]) == (2, 1)
    assert (stats["F444W"]["n_matched"], stats["F444W"]["n_matched_to_ref"]) == (3, 2)
    # Medians use only matches to F200W rows (0.02" and 0.03"), not the F090W-anchored one.
    assert stats["F444W"]["median_sep_arcsec"] == pytest.approx(0.025, rel=1e-3)
    assert set(stats["F200W"]) == set(stats["F444W"])
    assert set(m.meta["inputs"]) == {"F090W", "F200W", "F444W"}
    assert m.meta["inputs"]["F444W"]["provenance"] == "observed"
    for name in ("f200w", "f090w", "f444w"):
        assert f"clear-{name}_cat.ecsv" in m.meta["source"]


def test_merge_ref_only():
    m = catalog.merge_bands(_three_band_catalogs(), "F200W", include_unmatched=False)
    assert len(m) == 3 and set(m["anchor_band"]) == {"F200W"}
    assert m.meta["match_stats"]["F090W"]["n_unmatched"] == 1
    assert m.meta["match_stats"]["F444W"]["n_unmatched"] == 1
    assert list(m["n_bands"]) == [3, 2, 2]


def test_merge_is_one_to_one_closest_first():
    cats = {
        "F200W": _band("F200W", [(0.0, 0.0), (0.12, 0.0)]),
        # #1 lies between both ref sources; #2 is closest to ref #1 and wins it.
        "F444W": _band("F444W", [(0.05, 0.0), (0.02, 0.0)]),
    }
    m = catalog.merge_bands(cats, "F200W", radius_arcsec=0.1)
    assert list(m["f444w_label"]) == [2, 1]
    np.testing.assert_allclose(m["f444w_sep_arcsec"], [0.02, 0.07], rtol=1e-3)
    assert m.meta["match_stats"]["F444W"]["n_contested"] == 2

    cats["F200W"] = _band("F200W", [(0.0, 0.0)])
    m = catalog.merge_bands(cats, "F200W", radius_arcsec=0.1)
    assert list(m["source_uid"]) == [f"{FIELD}_f200w_1", f"{FIELD}_f444w_1"]
    assert list(m["f444w_label"]) == [2, 1]
    assert np.isnan(m["ref_nn_sep_arcsec"][0])  # only one F200W detection, and it is its own
    assert m["ref_nn_sep_arcsec"][1] == pytest.approx(0.05, rel=1e-3)  # split-detection flag


@pytest.mark.parametrize("radius, matched", [(0.1, False), (0.2, True)])
def test_merge_radius(radius, matched):
    cats = {"F200W": _band("F200W", [(0.0, 0.0)]), "F150W": _band("F150W", [(0.0, 0.15)])}
    m = catalog.merge_bands(cats, "F200W", radius_arcsec=radius)
    assert len(m) == (1 if matched else 2)
    assert bool(m["f150w_detected"][0]) is matched
    assert m.meta["match"]["radius_arcsec"] == radius


def test_merge_is_deterministic_and_case_insensitive():
    cats = _three_band_catalogs()
    a = catalog.merge_bands(cats, "F200W")
    b = catalog.merge_bands({k.lower(): cats[k] for k in reversed(list(cats))}, "f200w")
    assert a.colnames == b.colnames
    for name in a.colnames:
        assert np.array_equal(
            np.asarray(a[name]), np.asarray(b[name]), equal_nan=a[name].dtype.kind == "f"
        )


def test_merge_reports_and_warns_offsets():
    offsets = [(3.0 * i, 2.0 * (i % 3)) for i in range(12)]
    shifted = [(dx + 0.01, dy + 0.03) for dx, dy in offsets]
    cats = {"F200W": _band("F200W", offsets), "F356W": _band("F356W", shifted)}
    with pytest.warns(UserWarning, match="median offset"):
        m = catalog.merge_bands(cats, "F200W", radius_arcsec=0.1)
    s = m.meta["match_stats"]["F356W"]
    assert s["n_matched_to_ref"] == 12
    assert s["median_dra_arcsec"] == pytest.approx(0.01, abs=1e-4)
    assert s["median_ddec_arcsec"] == pytest.approx(0.03, abs=1e-4)
    assert s["median_sep_arcsec"] == pytest.approx(np.hypot(0.01, 0.03), rel=1e-3)
    # Positions are not corrected.
    np.testing.assert_allclose(m["ra"], cats["F200W"]["ra"])

    # Too few pairs for a meaningful median: reported, but no warning.
    few = {band: cat[:5] for band, cat in cats.items()}
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        m = catalog.merge_bands(few, "F200W", radius_arcsec=0.1)
    assert m.meta["match_stats"]["F356W"]["median_ddec_arcsec"] == pytest.approx(0.03, abs=1e-4)


def test_merge_warns_on_mixed_pipeline_versions():
    cats = {
        "F200W": _band("F200W", [(0.0, 0.0)]),
        "F444W": _band("F444W", [(0.0, 0.0)], jwst="3.0.0"),
    }
    with pytest.warns(UserWarning, match="different jwst/photutils versions"):
        catalog.merge_bands(cats, "F200W")


def test_merge_warns_when_matching_collapses():
    # A misregistered band (offset > radius) matches almost nothing; the offset median cannot
    # be measured, so the low match fraction is what must be reported.
    offsets = [(3.0 * i, 2.0 * (i % 3)) for i in range(20)]
    cats = {
        "F200W": _band("F200W", offsets),
        "F444W": _band("F444W", [(dx + 0.3, dy) for dx, dy in offsets]),
    }
    with pytest.warns(UserWarning, match="only 0 of 20 possible matches"):
        m = catalog.merge_bands(cats, "F200W", radius_arcsec=0.1)
    assert len(m) == 40


def test_merge_provenance_follows_inputs():
    cats = _three_band_catalogs()
    cats["F444W"].meta["provenance"] = "simulated"  # e.g. injection-recovery
    m = catalog.merge_bands(cats, "F200W")
    assert m.meta["provenance"] == "simulated"
    assert m.meta["inputs"]["F444W"]["provenance"] == "simulated"
    cats["F444W"].meta["provenance"] = "model_prediction"
    with pytest.raises(ValueError, match="observed, derived or simulated"):
        catalog.merge_bands(cats, "F200W")


def test_merge_normalizes_units_and_drops_unpositioned_rows():
    cats = _three_band_catalogs()
    q = QTable(cats["F444W"])
    q["ra"], q["dec"] = (q["ra"]).to(u.rad), (q["dec"]).to(u.rad)
    q["dec"][1] = np.nan * u.rad  # F444W #12 has no position
    cats["F444W"] = q
    with pytest.warns(UserWarning, match="F444W: dropping 1 rows with non-finite ra/dec"):
        m = catalog.merge_bands(cats, "F200W")
    assert list(m["f444w_label"]) == [11, -1, -1, 13]
    f444w_13 = _pos(10.0, 0.0)
    assert m["f444w_ra"][3] == pytest.approx(f444w_13.ra.deg)
    assert m["f444w_ra"].unit == u.deg
    assert m.meta["match_stats"]["F444W"]["n_no_position"] == 1
    assert cats["F444W"]["ra"].unit == u.rad  # input untouched


def test_merge_keeps_masked_input_values_masked():
    cats = _three_band_catalogs()
    t = Table(cats["F444W"], masked=True)
    t["is_extended"].mask[0] = True  # F444W #11, matched to F200W #1
    cats["F444W"] = t
    m = catalog.merge_bands(cats, "F200W")
    col = m["f444w_is_extended"]
    assert list(col.mask) == [True, False, False, False]
    assert not col[1]  # not detected: fill value, unmasked; f444w_detected tells them apart


def test_merge_columns_and_field_id():
    m = catalog.merge_bands(
        _three_band_catalogs(), "F200W", field_id="smacs0723_nircam", columns=["aper50_abmag"]
    )
    assert m["source_uid"][0] == "smacs0723_nircam_f200w_1"
    expected = ["detected", "sep_arcsec", "label", "ra", "dec", "aper50_abmag"]
    assert [c for c in m.colnames if c.startswith("f444w_")] == [f"f444w_{q}" for q in expected]
    with pytest.raises(TypeError, match="not a string"):
        catalog.merge_bands(_three_band_catalogs(), "F200W", columns="aper50_abmag")
    with pytest.raises(ValueError, match="aper50_abmagg"):
        catalog.merge_bands(_three_band_catalogs(), "F200W", columns=["aper50_abmagg"])


def test_merge_handles_empty_catalogs():
    empty = _band("F200W", [(0.0, 0.0)])[:0]
    other = _band("F444W", [(0.0, 0.0), (5.0, 0.0)])
    m = catalog.merge_bands({"F200W": empty, "F444W": other}, "F200W")
    assert list(m["source_uid"]) == [f"{FIELD}_f444w_1", f"{FIELD}_f444w_2"]
    assert np.all(np.isnan(m["ref_nn_sep_arcsec"]))
    m = catalog.merge_bands({"F200W": empty, "F444W": other[:0]}, "F200W")
    assert len(m) == 0 and m["source_uid"].dtype.kind == "U"


def _with(t: Table, **changes) -> Table:
    t = t.copy()
    for key, value in changes.items():
        if key == "meta":
            t.meta.clear()
            t.meta.update(value)
        else:
            t[key] = value
    return t


@pytest.mark.parametrize(
    "cats, kwargs, match",
    [
        ({"F200W": _band("F200W", [(0, 0)])}, {"ref_band": "F444W"}, "not in catalogs"),
        (
            {"F200W": _band("F200W", [(0, 0)]), "f200w": _band("F200W", [(0, 0)])},
            {"ref_band": "F200W"},
            "given twice",
        ),
        (
            {"F200W": _band("F200W", [(0, 0), (1, 1)], labels=[3, 3])},
            {"ref_band": "F200W"},
            "unique",
        ),
        (
            {
                "F200W": _with(
                    _band("F200W", [(0, 0)]), meta={"provenance": "observed", "source": "x"}
                )
            },
            {"ref_band": "F200W"},
            "field_id",
        ),
        ({"F200W": _band("F200W", [(0, 0)])}, {"ref_band": "F200W", "field_id": "a:b"}, "field_id"),
        (
            {"F200W": _band("F200W", [(0, 0)])},
            {"ref_band": "F200W", "field_id": "ab\n"},
            "field_id",
        ),
        ({"F200W": _band("F200W", [(0, 0)])}, {"ref_band": "F200W", "radius_arcsec": 0}, "radius"),
        ({"F200W": _with(_band("F200W", [(0, 0)]), meta={})}, {"ref_band": "F200W"}, "provenance"),
        ({}, {"ref_band": "F200W"}, "empty"),
        (
            {"F200W": _band("F200W", [(0, 0)]), "F444W/F470N": _band("F444W", [(0, 0)])},
            {"ref_band": "F200W"},
            "letters and digits",
        ),
        (
            {"F200W": _band("F200W", [(0, 0)]), "F444W": _band("F090W", [(0, 0)])},
            {"ref_band": "F200W"},
            "meta\\['band'\\]",
        ),
    ],
)
def test_merge_rejects_bad_input(cats, kwargs, match):
    with pytest.raises(ValueError, match=match):
        catalog.merge_bands(cats, **kwargs)


def test_merge_real_fixture_and_ecsv_round_trip(tmp_path):
    cats = {
        "F200W": catalog.load_pipeline_catalog(F200W_FIXTURE),
        "F444W": catalog.load_pipeline_catalog(F444W_FIXTURE),
    }
    m = catalog.merge_bands(cats, "F200W")
    assert len(m) == 13  # 10 F200W rows + 3 F444W detections without an F200W counterpart
    assert list(m["n_bands"]).count(2) == 5
    s = m.meta["match_stats"]["F444W"]
    assert s["n_matched"] == 5 and s["median_sep_arcsec"] < 0.02
    lw_only = m[m["source_uid"] == "jw02736-o001_t001_nircam_f444w_544"]
    assert len(lw_only) == 1 and not lw_only["f200w_detected"][0]
    assert lw_only["ref_nn_sep_arcsec"][0] > 1.0
    assert (
        m.meta["inputs"]["F444W"]["input_sha256"]
        == hashlib.sha256(F444W_FIXTURE.read_bytes()).hexdigest()
    )

    path = tmp_path / "merged.ecsv"
    m.write(path, format="ascii.ecsv")
    back = Table.read(path, format="ascii.ecsv")
    assert back.colnames == m.colnames
    assert back.meta["provenance"] == "derived"
    assert back.meta["match_stats"] == m.meta["match_stats"]
    assert list(back["source_uid"]) == list(m["source_uid"])


@pytest.mark.network
def test_merge_live_mast_catalogs(tmp_path):
    # Row counts were 3145 (F200W) and 1877 (F444W) with jwst 2.0.1 on 2026-10-07; MAST may
    # reprocess, so only loose checks here.
    url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/{}"
    cats = {}
    for band in ("F200W", "F444W"):
        name = f"jw02736-o001_t001_nircam_clear-{band.lower()}_cat.ecsv"
        with urllib.request.urlopen(url.format(name), timeout=120) as response:
            (tmp_path / name).write_bytes(response.read())
        cats[band] = catalog.load_pipeline_catalog(tmp_path / name)
        assert len(cats[band]) > 1000
        assert cats[band].meta["jwst_version"]
    m = catalog.merge_bands(cats, "F200W")
    s = m.meta["match_stats"]["F444W"]
    assert s["n_matched_to_ref"] > 1400
    assert s["median_sep_arcsec"] < 0.03
    assert abs(s["median_dra_arcsec"]) < 0.01 and abs(s["median_ddec_arcsec"]) < 0.01
    assert len(set(m["source_uid"])) == len(m)

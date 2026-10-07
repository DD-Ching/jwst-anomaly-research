"""Tests for cutouts.sample_weight_map and quality.assess_sources (D-011)."""

from __future__ import annotations

import numpy as np
import pytest
from astropy.table import Table
from test_cutouts import NX, NY, SCALE, SOURCES, _write_i2d

from jwst_anomaly import cutouts, quality, schema


@pytest.fixture
def image(tmp_path):
    path = tmp_path / "jw99999-o001_t001_nircam_clear-f200w_i2d.fits"
    return str(path), _write_i2d(path)


def _sources(w, positions, ci, snr_err=None) -> Table:
    snr_err = snr_err or {}
    names = list(positions)
    sky = w.pixel_to_world([positions[n][0] for n in names], [positions[n][1] for n in names])
    t = Table(
        {
            "source_uid": names,
            "ra": sky.ra.deg,
            "dec": sky.dec.deg,
            "ref_band": ["F200W"] * len(names),
            "n_bands": [1] * len(names),
            "f200w_CI_50_30": [ci.get(n, 1.9) for n in names],
            "f200w_aper50_abmag_err": [snr_err.get(n, 0.01) for n in names],
            "f200w_detected": [True] * len(names),
        }
    )
    t.meta.update(provenance=schema.Provenance.DERIVED.value, source="unit test")
    return t


def test_sample_weight_map_whole_read(image):
    path, _ = image
    wm = cutouts.sample_weight_map(path, step=10)
    assert wm.band == "F200W"
    assert wm.shape == (NY, NX)
    assert wm.values.shape == (NY // 10, NX // 10)
    assert wm.reference_weight == pytest.approx(1.0)
    assert wm.pixel_scale_arcsec == pytest.approx(SCALE, rel=1e-6)


def test_sample_weight_map_row_reads_match_whole_read(image, monkeypatch):
    path, _ = image
    whole = cutouts.sample_weight_map(path, step=10)
    monkeypatch.setattr(cutouts, "_WHOLE_WHT_BYTES", 0)  # force one read per sampled row
    rows = cutouts.sample_weight_map(path, step=10)
    np.testing.assert_array_equal(rows.values, whole.values)


def test_sample_weight_map_remote_path(image):
    fsspec = pytest.importorskip("fsspec")
    path, _ = image
    fs = fsspec.filesystem("memory")
    with open(path, "rb") as f:
        fs.pipe("/quality/jw99999_i2d.fits", f.read())
    wm = cutouts.sample_weight_map("memory://quality/jw99999_i2d.fits", step=10)
    assert wm.values.shape == (NY // 10, NX // 10)


def test_sample_weight_map_needs_wht(tmp_path):
    path = tmp_path / "nowht_i2d.fits"
    _write_i2d(path, with_wht=False)
    with pytest.raises(ValueError, match="WHT"):
        cutouts.sample_weight_map(str(path), step=10)


def test_weight_map_at(image):
    path, w = image
    wm = cutouts.sample_weight_map(path, step=10)
    names = ["interior", "low_weight", "footprint_edge"]
    sky = w.pixel_to_world([SOURCES[n][0] for n in names], [SOURCES[n][1] for n in names])
    rel, edge, covered = wm.at(sky.ra.deg, sky.dec.deg)
    assert covered.all()
    assert rel[0] == pytest.approx(1.0) and edge[0] > 2.0
    assert rel[1] == pytest.approx(0.2)
    assert edge[2] == pytest.approx(20 * SCALE, abs=0.5 * 10 * SCALE)  # 20 px from the strip
    off = w.pixel_to_world(NX + 50, NY + 50)
    rel_off, edge_off, cov_off = wm.at(off.ra.deg, off.dec.deg)
    assert rel_off[0] == 0 and edge_off[0] == 0 and not cov_off[0]


def test_edge_distance_follows_position_not_cells(image):
    path, w = image
    wm = cutouts.sample_weight_map(path, step=10)
    xs = [24.0, 30.0, 36.0, 42.0, 48.0]  # strip of zero weight at x < 20
    sky = w.pixel_to_world(xs, [100.0] * len(xs))
    _, edge, covered = wm.at(sky.ra.deg, sky.dec.deg)
    assert covered.all()
    assert np.all(np.diff(edge) > 0)  # strictly increasing within and across cells
    np.testing.assert_allclose(edge, (np.array(xs) - 20) * SCALE, atol=0.5 * 10 * SCALE)


def test_grid_arcsec_sets_step_from_pixel_scale(image):
    path, _ = image
    assert cutouts.sample_weight_map(path, grid_arcsec=1.0).step == round(1.0 / SCALE)
    assert cutouts.sample_weight_map(path, grid_arcsec=0.5).step == round(0.5 / SCALE)


def test_sample_weight_map_rejects_zero_weight(tmp_path):
    from astropy.io import fits

    path = tmp_path / "zero_i2d.fits"
    _write_i2d(path)
    with fits.open(path, mode="update") as hdul:
        hdul["WHT"].data[:] = 0
    with pytest.raises(ValueError, match="no positive weight"):
        cutouts.sample_weight_map(str(path), step=10)


def test_assess_sources_flags(image):
    path, w = image
    wm = cutouts.sample_weight_map(path, step=10)
    positions = {
        "interior": SOURCES["interior"],
        "low_weight": SOURCES["low_weight"],
        "footprint_edge": SOURCES["footprint_edge"],
        "off_image": (NX + 50.0, NY + 50.0),
        "hot_pixel": (160.0, 110.0),
    }
    q = quality.assess_sources(
        _sources(w, positions, ci={"hot_pixel": 1.2}), wm, min_edge_arcsec=1.5
    )
    reason = dict(zip(q["source_uid"], q["quality_reason"], strict=True))
    assert reason == {
        "interior": "",
        "low_weight": "low_weight",
        "footprint_edge": "edge",
        "off_image": "no_coverage",
        "hot_pixel": "sharper_than_psf",
    }
    assert list(q["quality_ok"]) == [True, False, False, False, False]
    assert q.meta["provenance"] == "derived"
    assert q.meta["thresholds"]["provenance"] == "assumption"
    counts = quality.summarize(q)
    assert counts == {
        "n_sources": 5,
        "n_ok": 1,
        "low_weight": 1,
        "edge": 1,
        "no_coverage": 1,
        "sharper_than_psf": 1,
    }


def test_assess_sources_without_weight_map_only_tests_ci(image):
    _, w = image
    q = quality.assess_sources(
        _sources(w, {"a": SOURCES["interior"], "b": (10.0, 10.0)}, ci={"b": 1.0}), None
    )
    assert list(q["quality_reason"]) == ["", "sharper_than_psf"]
    assert np.all(np.isnan(q["rel_weight"]))


def test_assess_sources_undetected_ref_band_is_not_sharp(image):
    _, w = image
    src = _sources(w, {"a": SOURCES["interior"]}, ci={"a": float("nan")})
    src["f200w_detected"] = [False]
    q = quality.assess_sources(src, None)
    assert bool(q["quality_ok"][0])


def test_sharpness_ignores_noise_and_nonpositive_ci(image):
    _, w = image
    positions = {n: SOURCES["interior"] for n in ("faint", "negative", "zero", "bright_hot")}
    src = _sources(
        w,
        positions,
        ci={"faint": 1.1, "negative": -3.0, "zero": 0.0, "bright_hot": 1.1},
        snr_err={"faint": 0.5},  # S/N ~2 < 3
    )
    q = quality.assess_sources(src, None)
    assert dict(zip(q["source_uid"], q["sharper_than_psf"], strict=True)) == {
        "faint": False,
        "negative": False,
        "zero": False,
        "bright_hot": True,
    }


def test_band_mismatch_and_wrong_image_raise(image):
    path, w = image
    wm = cutouts.sample_weight_map(path, step=10)
    src = _sources(w, {"a": SOURCES["interior"]}, ci={})
    with pytest.raises(ValueError, match="ref band"):
        quality.assess_sources(src, wm, ref_band="F444W")
    far = _sources(w, {f"s{i}": (NX + 100.0 + i, NY + 100.0) for i in range(10)}, ci={})
    with pytest.raises(ValueError, match="covers only"):
        quality.assess_sources(far, wm)


def test_detection_confirmation_flags(image):
    """D-014: low_snr and single_band (unless confirmed by an independent detection)."""
    _, w = image
    positions = {n: SOURCES["interior"] for n in ("good", "faint", "single", "single_conf")}
    src = _sources(w, positions, ci={}, snr_err={"faint": 0.5})  # S/N ~2.2
    src["n_bands"] = [3, 3, 1, 1]
    src["dja05_match_sep_arcsec"] = [0.01, 0.02, np.nan, 0.03]
    q = quality.assess_sources(
        src,
        None,
        min_detection_snr=5,
        require_multiband=True,
        confirm_column="dja05_match_sep_arcsec",
    )
    assert dict(zip(q["source_uid"], q["quality_reason"], strict=True)) == {
        "good": "",
        "faint": "low_snr",
        "single": "single_band",
        "single_conf": "",
    }
    assert q.meta["thresholds"]["min_detection_snr"] == 5
    unconfirmed = quality.assess_sources(src, None, require_multiband=True)
    assert list(unconfirmed["quality_reason"])[2:] == ["single_band", "single_band"]
    off = quality.assess_sources(src, None)
    assert all(r == "" for r in off["quality_reason"])


def test_detection_confirmation_input_errors(image):
    _, w = image
    src = _sources(w, {"a": SOURCES["interior"]}, ci={})
    with pytest.raises(ValueError, match="confirm_column"):
        quality.assess_sources(src, None, require_multiband=True, confirm_column="nope")
    src.remove_column("f200w_aper50_abmag_err")
    with pytest.raises(ValueError, match="min_detection_snr"):
        quality.assess_sources(src, None, min_detection_snr=5)


def test_low_snr_uses_the_best_detected_band_so_dropouts_survive(image):
    _, w = image
    src = _sources(w, {"dropout": SOURCES["interior"], "faint": SOURCES["interior"]}, ci={})
    src["f200w_detected"] = [False, True]  # the dropout is not detected in the reference band
    src["f200w_aper50_abmag_err"] = [np.nan, 0.5]
    src["f444w_detected"] = [True, True]
    src["f444w_aper50_abmag_err"] = [0.003, 0.6]  # S/N ~360 vs ~1.8
    q = quality.assess_sources(src, None, min_detection_snr=5)
    assert list(q["quality_reason"]) == ["", "low_snr"]

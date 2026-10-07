"""Offline tests for features.build_features on synthetic merged source tables."""

import numpy as np
import pytest
from astropy.table import MaskedColumn, Table

from jwst_anomaly import schema
from jwst_anomaly.features import NAN_POLICY, build_features, discover_bands
from jwst_anomaly.rank import score_anomalies

BANDS = ("F090W", "F150W", "F200W", "F277W")


def make_sources(n=120, bands=BANDS, ref_band="F200W", seed=0, p_missing=0.1):
    """Merged table following the unit-2 contract; ``p_missing`` non-detections per band."""
    rng = np.random.default_rng(seed)
    t = Table()
    t["source_uid"] = [f"s{i:04d}" for i in range(n)]
    t["ra"] = 110.8 + rng.normal(0, 0.01, n)
    t["dec"] = -73.45 + rng.normal(0, 0.01, n)
    t["ref_band"] = [ref_band] * n
    base = rng.normal(27.0, 1.0, n)
    detected_all = []
    for i, band in enumerate(bands):
        det = rng.random(n) >= p_missing
        if band == ref_band:
            det[:] = True
        detected_all.append(det)

        def col(values, det=det):
            return np.where(det, values, np.nan)

        b = band.lower()
        t[f"{b}_aper50_abmag"] = col(base - 0.1 * i + rng.normal(0, 0.05, n))
        t[f"{b}_aper50_abmag_err"] = col(np.full(n, 0.05))
        t[f"{b}_isophotal_area"] = col(rng.lognormal(4, 0.5, n))
        t[f"{b}_ellipticity"] = col(rng.uniform(0, 0.6, n))
        t[f"{b}_CI_50_30"] = col(rng.normal(1.6, 0.1, n))
        t[f"{b}_CI_70_50"] = col(rng.normal(1.3, 0.05, n))
        t[f"{b}_sharpness"] = col(rng.normal(0.5, 0.1, n))
        t[f"{b}_roundness"] = col(rng.normal(0.0, 0.2, n))
        t[f"{b}_nn_dist"] = col(rng.uniform(10, 100, n))
        t[f"{b}_is_extended"] = det & (rng.random(n) > 0.5)
        t[f"{b}_detected"] = det
        t[f"{b}_sep_arcsec"] = col(rng.uniform(0, 0.05, n))
    t["n_bands"] = np.sum(detected_all, axis=0)
    t.meta.update(provenance="derived", source="synthetic merged table")
    return t


def test_contract_and_feature_spec():
    src = make_sources()
    f = build_features(src)
    schema.validate(f, schema.FEATURE_ID_COLUMNS)
    assert f.meta["provenance"] == schema.Provenance.DERIVED.value
    assert "synthetic merged table" in f.meta["source"]
    assert f.meta["nan_policy"] == NAN_POLICY
    spec = f.meta["feature_spec"]
    assert list(spec) == [c for c in f.colnames if c != "source_uid"]
    assert all(isinstance(v, str) and v for v in spec.values())
    assert list(f["source_uid"]) == list(src["source_uid"])
    assert f.meta["bands"] == ["f090w", "f150w", "f200w", "f277w"]
    assert f.meta["ref_band"] == "F200W"
    for name in (
        "ref_mag",
        "ref_log_isophotal_area",
        "ref_ellipticity",
        "ref_log_ci_50_30",
        "ref_log_ci_70_50",
        "ref_sharpness",
        "ref_roundness",
        "ref_log_nn_dist",
        "blue_dropout",
        "red_dropout",
        "n_gaps",
    ):
        assert name in spec
    assert f.meta["min_snr_skipped_bands"] == []
    assert [c for c in spec if c.startswith("color_")] == [
        "color_f090w_f150w",
        "color_f150w_f200w",
        "color_f200w_f277w",
    ]


def test_discover_bands_in_wavelength_order():
    t = Table(
        {
            "f444w_aper50_abmag": [1.0],
            "f1000w_detected": [True],
            "f090w_detected": [True],
            "f150w2_aper50_abmag": [1.0],
            "foo_bar": [1.0],
            "f200w_ellipticity": [0.1],
        }
    )
    # f200w has no detected/abmag column, so it is not a band of this table.
    assert discover_bands(t) == ["f090w", "f150w2", "f444w", "f1000w"]


def test_colors_use_adjacent_aper50_magnitudes():
    src = make_sources()
    f = build_features(src)
    expected = np.asarray(src["f090w_aper50_abmag"]) - np.asarray(src["f150w_aper50_abmag"])
    np.testing.assert_allclose(f["color_f090w_f150w"], expected, equal_nan=True)
    undetected = ~np.asarray(src["f090w_detected"])
    assert undetected.any() and np.all(np.isnan(f["color_f090w_f150w"][undetected]))
    np.testing.assert_allclose(f["ref_mag"], src["f200w_aper50_abmag"])


def test_undefined_values_stay_nan():
    src = make_sources(n=10)
    src["f200w_CI_50_30"][0] = -2.0  # negative aperture flux
    src["f200w_isophotal_area"][1] = 0.0
    src["f200w_aper50_abmag"][2] = np.nan  # pipeline NaN magnitude
    f = build_features(src)
    assert np.isnan(f["ref_log_ci_50_30"][0])
    assert np.isnan(f["ref_log_isophotal_area"][1])
    assert np.isnan(f["ref_mag"][2]) and np.isnan(f["color_f150w_f200w"][2])
    np.testing.assert_allclose(f["ref_log_ci_50_30"][3], np.log10(src["f200w_CI_50_30"][3]))


def test_min_snr_masks_noise_dominated_measurements():
    src = make_sources(n=10)
    src["f200w_aper50_abmag_err"][0] = 1.0  # S/N ~ 1.1
    f = build_features(src)
    for name in ("ref_mag", "ref_sharpness", "ref_log_ci_50_30", "color_f150w_f200w"):
        assert np.isnan(f[name][0]), name
    g = build_features(src, min_snr=None)
    # the detection pattern ignores the S/N floor
    assert f["red_dropout"][0] == g["red_dropout"][0] and f["n_gaps"][0] == g["n_gaps"][0]
    assert np.isfinite(f["ref_mag"][1])
    assert np.isfinite(g["ref_mag"][0]) and np.isfinite(g["ref_sharpness"][0])
    assert f.meta["min_snr"] == 3.0 and g.meta["min_snr"] is None


def test_min_snr_without_error_column_warns_and_is_recorded():
    src = make_sources(n=10)
    del src["f090w_aper50_abmag_err"]
    with pytest.warns(UserWarning, match="f090w"):
        f = build_features(src)
    assert f.meta["min_snr_skipped_bands"] == ["f090w"]


def test_detection_pattern_counts_each_missing_band_once():
    src = make_sources(n=5)
    patterns = [(0, 1, 1, 1), (1, 0, 0, 0), (0, 1, 0, 1), (1, 1, 1, 1), (1, 0, 0, 1)]
    for b, band in enumerate(BANDS):
        src[f"{band.lower()}_detected"] = [bool(p[b]) for p in patterns]
    f = build_features(src)
    assert list(f["blue_dropout"]) == [1, 0, 1, 0, 0]
    assert list(f["red_dropout"]) == [0, 3, 0, 0, 0]
    assert list(f["n_gaps"]) == [0, 0, 1, 0, 2]
    n_missing = np.asarray(f["blue_dropout"]) + f["red_dropout"] + f["n_gaps"]
    assert list(n_missing) == [4 - sum(p) for p in patterns]


def test_float_detected_flag_treats_nan_as_not_detected():
    src = make_sources(n=3, p_missing=0.0)
    src["f090w_detected"] = [1.0, np.nan, 0.0]
    f = build_features(src)
    assert list(f["blue_dropout"]) == [0, 1, 1]


def test_fallback_detection_ignores_integer_sentinels():
    src = make_sources(n=3, p_missing=0.0)
    del src["f090w_detected"]
    for c in [c for c in src.colnames if c.startswith("f090w_") and src[c].dtype.kind == "f"]:
        src[c][1] = np.nan
    src["f090w_label"] = [5, -1, 7]  # merge sentinel for "no match"
    f = build_features(src)
    assert list(f["blue_dropout"]) == [0, 1, 0]


def test_inputs_are_not_modified():
    src = make_sources(n=5)
    src["f200w_sharpness"][0] = np.inf
    before = {c: np.array(src[c]) for c in src.colnames}
    build_features(src)
    for c, values in before.items():
        assert np.array_equal(np.asarray(src[c]), values, equal_nan=values.dtype.kind == "f"), c


def test_unknown_bands_are_rejected():
    src = make_sources(n=5)
    with pytest.raises(ValueError, match="f20w"):
        build_features(src, ref_band="F20W")
    with pytest.raises(ValueError, match="f2000w"):
        build_features(src, bands=["F200W", "F2000W"])
    f = build_features(src, bands="F200W")  # a single band given as a string
    assert f.meta["bands"] == ["f200w"]


def test_masked_ref_band_values_are_ignored():
    src = make_sources(n=5)
    src["ref_band"] = MaskedColumn(["F150W", "", "", "", "F150W"], mask=[0, 1, 1, 1, 0])
    assert build_features(src).meta["ref_band"] == "F150W"


def test_subset_of_bands_and_quantities():
    t = Table()
    t["source_uid"] = ["a", "b", "c"]
    t["ra"] = [1.0, 2.0, 3.0]
    t["dec"] = [0.0, 0.0, 0.0]
    t["ref_band"] = ["F200W"] * 3
    t["n_bands"] = [1, 1, 1]
    t["f200w_aper50_abmag"] = [25.0, 26.0, 27.0]
    t["f200w_detected"] = [True, True, True]
    t.meta.update(provenance="derived", source="minimal")
    f = build_features(t, min_snr=None)  # no error column, so no S/N floor
    assert f.meta["bands"] == ["f200w"]
    assert set(f.meta["feature_spec"]) == {"ref_mag", "blue_dropout", "red_dropout", "n_gaps"}
    assert "ref_sharpness" in f.meta["dropped_features"]
    assert "ref_mag" in f.colnames and "ref_sharpness" not in f.colnames


def test_ref_band_default_and_override():
    src = make_sources(n=6)
    src["ref_band"] = ["F150W", "F200W", "F200W", "F200W", "F150W", "F200W"]
    assert build_features(src).meta["ref_band"] == "F200W"
    f = build_features(src, ref_band="f150w")
    assert f.meta["ref_band"] == "F150W"
    np.testing.assert_allclose(f["ref_mag"], src["f150w_aper50_abmag"], equal_nan=True)


def test_masked_columns_and_missing_detected_flag():
    src = make_sources(n=5, p_missing=0.0)
    del src["f090w_detected"]  # fall back to "any finite value in that band"
    for c in [c for c in src.colnames if c.startswith("f090w_") and src[c].dtype.kind == "f"]:
        src[c] = MaskedColumn(src[c], mask=[False, True, False, False, False])
    f = build_features(src)
    assert np.isnan(f["color_f090w_f150w"][1]) and np.isfinite(f["color_f090w_f150w"][0])
    assert list(f["blue_dropout"]) == [0, 1, 0, 0, 0]


def test_rejects_input_without_provenance():
    src = make_sources(n=3)
    src.meta.clear()
    with pytest.raises(ValueError, match="provenance"):
        build_features(src)


def test_features_feed_score_anomalies():
    f = build_features(make_sources(n=150))
    scores = score_anomalies(f, random_state=0)
    schema.validate(scores, schema.SCORE_COLUMNS)
    assert len(scores) == 150
    assert scores.meta["feature_names"] == list(f.meta["feature_spec"])


def test_daofind_stats_only_for_point_like_sources():
    from jwst_anomaly.features import build_features

    t = Table(
        {
            "source_uid": ["point", "extended", "unknown_ci"],
            "ra": [1.0, 1.0, 1.0],
            "dec": [2.0, 2.0, 2.0],
            "ref_band": ["F200W"] * 3,
            "n_bands": [1, 1, 1],
            "f200w_detected": [True, True, True],
            "f200w_aper50_abmag": [24.0, 24.0, 24.0],
            "f200w_aper50_abmag_err": [0.01, 0.01, 0.01],
            "f200w_CI_50_30": [1.65, 2.3, np.nan],
            "f200w_sharpness": [0.6, 1.8, 0.7],
            "f200w_roundness": [0.1, 0.5, 0.2],
            "f200w_isophotal_area": [20.0, 900.0, 30.0],
        }
    )
    t.meta.update(provenance="derived", source="test")
    plain = build_features(t)
    gated = build_features(t, daofind_max_ci=1.8)
    assert np.isfinite(plain["ref_sharpness"]).all()
    assert gated["ref_sharpness"][0] == pytest.approx(0.6)
    assert np.isnan(gated["ref_sharpness"][1]) and np.isnan(gated["ref_roundness"][1])
    assert np.isnan(gated["ref_sharpness"][2])  # CI unknown -> not point-like
    assert "point-like" in gated.meta["feature_spec"]["ref_sharpness"]
    assert gated.meta["daofind_max_ci"] == 1.8


def test_matched_aperture_colours_keep_pipeline_morphology_for_unmatched_sources():
    """Regression (PR #14 review): morphology is gated on the pipeline aperture, not on DJA S/N."""
    from jwst_anomaly.features import build_features

    nan = np.nan
    t = Table(
        {
            "source_uid": ["matched", "unmatched"],
            "ra": [1.0, 1.0],
            "dec": [2.0, 2.0],
            "ref_band": ["F200W", "F200W"],
            "n_bands": [2, 2],
            "f200w_detected": [True, True],
            "f277w_detected": [True, True],
            "f200w_aper50_abmag": [24.0, 24.0],
            "f200w_aper50_abmag_err": [0.01, 0.01],
            "f277w_aper50_abmag": [23.0, 23.0],
            "f277w_aper50_abmag_err": [0.01, 0.01],
            "f200w_dja05_abmag": [24.2, nan],
            "f200w_dja05_abmag_err": [0.02, nan],
            "f277w_dja05_abmag": [24.3, nan],
            "f277w_dja05_abmag_err": [0.02, nan],
            "f200w_isophotal_area": [500.0, 800.0],
            "f200w_ellipticity": [0.3, 0.6],
            "f200w_CI_50_30": [2.0, 2.2],
        }
    )
    t.meta.update(provenance="derived", source="test")
    out = build_features(t, aperture="dja05")
    assert out["color_f200w_f277w"][0] == pytest.approx(-0.1)
    assert np.isnan(out["color_f200w_f277w"][1])  # no DJA match -> no colour
    assert np.isfinite(out["ref_log_isophotal_area"][1])  # morphology kept
    assert out["ref_ellipticity"][1] == pytest.approx(0.6)
    assert out["ref_mag"][1] == pytest.approx(24.0)  # pipeline aperture
    assert out.meta["morph_aperture"] == "aper50"

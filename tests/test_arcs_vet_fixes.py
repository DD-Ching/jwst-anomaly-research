"""Tests for the PLAN_fixes.md fixes (synthetic data; no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


avs = _load("arcs_vet_spectra")
vd = _load("arcs_vet_deblend")
vce = _load("arcs_vet_census")


def _spectrum(n=1200, seed=1):
    rng = np.random.default_rng(seed)
    return rng.normal(0, 1, n), np.ones(n), np.ones(n, bool)


def test_invalid_channels_excluded_and_cannot_detect():
    f, var, valid = _spectrum()
    f[600:606] += 50.0  # a strong "line" ...
    valid[600:606] = False  # ... on channels with missing spaxels
    sp = avs.spectrum_snr(f, var, valid)
    assert not sp["usable"][600:606].any()
    assert np.all(np.isnan(sp["snr"][600:606]))
    assert np.nanmax(sp["snr"]) < 5  # the line on invalid channels cannot be detected


def test_aperture_spectrum_marks_partial_channels_invalid():
    data = np.ones((5, 3, 3))
    stat = np.ones((5, 3, 3))
    data[2, 1, 1] = np.nan  # one spaxel missing in channel 2
    stat[3, 0, 0] = 0.0  # zero variance in channel 3
    ap = np.ones((3, 3), bool)
    f, var, valid = avs.aperture_spectrum(data, stat, ap)
    assert list(valid) == [True, True, False, False, True]
    assert np.isnan(f[2]) and np.isnan(var[3]) and f[0] == 9.0


def test_sky_line_channels_do_not_inflate_the_rescale():
    f, var, valid = _spectrum(seed=2)
    sky = np.arange(50, 1150, 37)
    for c in sky:  # sky residuals: large flux excursions with raised variance
        f[c] += 40.0 * (-1) ** c
        var[c] = 4.0
    sp = avs.spectrum_snr(f, var, valid)
    assert sp["resc"] < 1.3
    for c in sky:
        assert not sp["usable"][c]


def test_exceedance_runs_counts_separate_peaks():
    snr = np.array([0, 6, 6, 0, 0, 7, 0, 6], float)
    usable = np.ones(8, bool)
    usable[7] = False
    assert avs.exceedance_runs(snr, usable, 5.0) == 2


def test_flux_scale_for_negative_filtered_sum():
    D = np.full(50, -0.03)
    D[:5] = [0.1, 0.2, 0.3, 0.2, 0.1]  # positive bump on an over-subtracted background
    start, up = vd.flux_scale(D)
    assert D.sum() < 0
    assert start > 0.5 and up >= 10 * start
    s, s_pos = vd.flux_scale(np.array([1.0, 2.0, -0.5]))
    assert s == 2.5 and s_pos == 30.0
    tiny = np.r_[np.full(45, -0.0222222), [0.1, 0.2, 0.3, 0.2, 0.1]]  # sum ≈ +1e-6, bump 0.9
    start, _ = vd.flux_scale(tiny)
    assert start > 0.5  # a near-zero positive sum does not set a pinned start


def test_fit_not_pinned_for_faint_source_with_negative_sum():
    yy, xx = np.indices((31, 31), dtype=float)
    img = vd._gauss(xx, yy, 15.0, 15.0, 3.0, 1.2) - 0.03  # sum within r=8 px is negative
    err = np.full_like(img, 0.01)
    m = (xx - 15) ** 2 + (yy - 15) ** 2 <= 64
    assert img[m].sum() < 0
    f = vd.fit_models(img, err, 15.0, 15.0, 8.0)
    assert f["m1"][2] > 1.0  # recovered flux is not pinned near 1e-5
    assert f["dbic"] != 0.0


def test_p_tail_drops_row_from_null():
    probs = {"a": 0.5, "b": 0.5, "c": 0.9}
    labels = ["a", "b", "c"]
    full = vce.p_tail(probs, labels, 2)
    dropped = vce.p_tail(probs, labels, 2, drop=("c",))
    assert np.isclose(dropped, 0.25)  # both remaining rows anti
    assert full > dropped

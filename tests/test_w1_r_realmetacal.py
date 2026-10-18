"""Tests for scripts/w1_r_realmetacal.py (frozen rule, matching, summary, metacal sign)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_r_realmetacal as m  # noqa: E402


def test_rule_outcome_threshold():
    assert m.rule_outcome(0.414, 0.0) == "conservative"
    assert m.rule_outcome(0.40, 0.01) == "conservative"  # within 2σ below
    assert m.rule_outcome(0.39, 0.01) == "optimistic: re-run lower"
    assert m.rule_outcome(None, None) == "undetermined"


def test_match_follows_sheared_positions():
    ref = {"x": np.array([100.0, 900.0]), "y": np.array([500.0, 800.0])}
    c = (511.5, 511.5)
    g1, g2 = 0.01, 0.0
    moved = {
        "x": c[0] + (1 + g1) * (ref["x"] - c[0]),
        "y": c[1] + (1 - g1) * (ref["y"] - c[1]),
    }
    assert list(m._match(ref, moved, g1, g2, c)) == [0, 1]
    # without the shear mapping the far source (≈ 4 px shift) is lost
    assert m._match(ref, moved, 0.0, 0.0, c)[0] == -1


def test_summary_mean_and_components():
    s = m.summary(np.array([1.0, 2.0]), np.array([3.0, 4.0]), n_boot=50)
    assert s["n"] == 2 and s["R"] == pytest.approx(2.5)
    assert s["R11"] == pytest.approx(1.5) and s["R22"] == pytest.approx(3.5)
    assert m.summary(np.zeros(0), np.zeros(0))["R"] is None


def test_good_tiles_skips_gaps():
    shape = (2048, 2048)

    def rows(y0, x0):
        sci = np.ones((m.TILE, m.TILE))
        if (y0, x0) == (0, 0):
            sci[:, :100] = np.nan  # 9.8 % missing
        return sci, np.ones_like(sci)

    assert m.good_tiles(rows, shape, n=2) == [(0, 1024), (1024, 0)]


def test_metacal_shear_sign():
    galsim = pytest.importorskip("galsim")
    pytest.importorskip("photutils")
    gal = galsim.Convolve(galsim.Gaussian(sigma=4 * m.PIX_ARCSEC), m.psf_d())
    arr = gal.drawImage(nx=64, ny=64, scale=m.PIX_ARCSEC).array.astype(float)
    p, q = m.metacal(arr, 0.05, 0.0), m.metacal(arr, -0.05, 0.0)

    def e1(im):
        y, x = np.mgrid[: im.shape[0], : im.shape[1]] - 31.5
        return ((x * x - y * y) * im).sum() / ((x * x + y * y) * im).sum()

    assert e1(p) > 0.03 and e1(q) < -0.03


def test_posthoc_excludes_tile_and_bootstraps(tmp_path):
    tiles = [
        {"y0": 0, "x0": 1024, "n": 10, "R": 0.5, "err": 0.1},
        {"y0": 1024, "x0": 0, "n": 30, "R": 1.0, "err": 0.2},
    ]
    p = tmp_path / "r.json"
    p.write_text('{"variants": {"plan": {"per_tile": ' + str(tiles).replace("'", '"') + "}}}")
    out = m.posthoc(p, tmp_path / "o.json", n_boot=100)["variants"]["plan"]
    assert out["all_tiles"]["R"] == pytest.approx((10 * 0.5 + 30 * 1.0) / 40)
    assert out["excluded"]["R"] == pytest.approx(1.0) and out["excluded"]["n_tiles"] == 1
    assert out["excluded"]["err_tile_bootstrap"] == pytest.approx(0.0)


# ---- v2 (PLAN_v2.md) ----


def test_match_unique_drops_shared_targets():
    ref = {"x": np.array([100.0, 100.8, 400.0]), "y": np.array([100.0, 100.0, 400.0])}
    cat = {"x": np.array([100.4, 400.0]), "y": np.array([100.0, 400.0])}
    c = (511.5, 511.5)
    assert list(m._match(ref, cat, 0.0, 0.0, c)) == [0, 0, 1]
    assert list(m._match(ref, cat, 0.0, 0.0, c, unique=True)) == [-1, -1, 1]


def test_prep_tile_v2_fills_holes_with_noise_and_pads():
    arr = np.full((128, 128), 5.0) + np.random.default_rng(0).normal(0, 1.0, (128, 128))
    arr[60:64, 60:64] = np.nan
    out, near, sky, sig, bad_frac = m.prep_tile_v2(arr, seed=3)
    assert bad_frac == pytest.approx(16 / 128**2)
    assert np.isfinite(out).all()
    assert sky == pytest.approx(5.0, abs=0.1) and sig == pytest.approx(1.0, rel=0.1)
    hole = out[58:66, 58:66]  # hole + 2-px dilation: noise around 0, not zeros
    assert abs(hole.mean()) < 0.6 and hole.std() > 0.5
    assert near[60 - m.MASK_PAD_PX, 61] and not near[60 - m.MASK_PAD_PX - 2, 61]
    assert near[60 - 11, 60 - 11]  # diagonal, Euclidean 15.6 px: inside the pad


def test_tile_gate_and_rule_v2():
    f = m.PSF_SIGMA_FLOOR_PX
    tiles = [
        {"n_detected": 2000, "psf_sigma_raw_px": 1.5},
        {"n_detected": 1800, "psf_sigma_raw_px": 1.3},
        {"n_detected": 94, "psf_sigma_raw_px": 1.4},  # too few detections
        {"n_detected": 1900, "psf_sigma_raw_px": f - 0.01},  # below the floor
        {"n_detected": None, "psf_sigma_raw_px": None},
    ]
    assert m.tile_gate(tiles) == [True, True, False, False, False]
    assert m.rule_outcome_v2(1.0, 0.09) == "adoptable"
    assert m.rule_outcome_v2(1.0, 0.11) == "conservative"
    assert m.rule_outcome_v2(0.3, 0.02) == "optimistic: re-run lower"
    assert m.rule_outcome_v2(None, None) == "undetermined"


def test_tile_bootstrap_zero_for_identical_tiles():
    r = [np.array([1.0, 2.0]), np.array([1.0, 2.0])]
    assert m.tile_bootstrap(r, n_boot=50) == pytest.approx(0.0)
    assert m.tile_bootstrap([np.array([1.0])]) is None

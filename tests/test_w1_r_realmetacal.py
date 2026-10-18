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

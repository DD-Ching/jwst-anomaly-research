"""Offline tests of the 700.1 vetting helpers (scripts/a2744_s700_vet.py)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
s7v = pytest.importorskip("a2744_s700_vet")


def test_pair_chi2_equal_colours_any_ratio():
    fb = np.array([1.0, 2.0, 3.0, 4.0])
    r = s7v.pair_chi2(0.5 * fb, np.full(4, 0.01), fb, np.full(4, 0.01))
    assert r["ratio"] == pytest.approx(0.5, rel=1e-4)
    assert r["chi2"] < 1e-6 and r["p"] > 0.99 and r["n_bands"] == 4


def test_pair_chi2_different_colours_and_nan_bands():
    fa = np.array([1.0, 1.0, 1.0, np.nan])
    fb = np.array([1.0, 2.0, 4.0, 1.0])
    r = s7v.pair_chi2(fa, np.full(4, 0.05), fb, np.full(4, 0.05))
    assert r["n_bands"] == 3
    assert s7v.differs(r["p"]) is True
    assert s7v.pair_chi2([1.0], [0.1], [1.0], [0.1])["p"] is None
    assert s7v.differs(None) is None


def test_decide_frozen_rule_and_post_hoc_branch():
    assert s7v.decide(False, True, True).startswith("undetermined (colour test not validated;")
    assert "no member perturber" in s7v.decide(False, False, False)
    assert s7v.decide(True, True, True) == "one image (H1)"
    assert s7v.decide(True, False, True) == "ordinary perturber possible (H2)"
    assert s7v.decide(True, False, False) == "undetermined (H3)"
    assert s7v.decide(True, None, True) == "undetermined (target colours not measurable)"


def test_pair_chi2_ratio_outside_the_sum_ratio():
    # sum(fb) nearly cancels, so r0 = sum(fa) / sum(fb) ≈ 302 and the old [r0/10, 10 r0] bracket
    # excluded the best ratio near 2
    fb = np.array([1.0, -0.999, 0.5, -0.5])
    fa = 2 * fb + np.array([0.3, 0.0, 0.0, 0.0])
    r = s7v.pair_chi2(fa, np.full(4, 0.1), fb, np.full(4, 0.1))
    brute = np.linspace(0.5, 5.0, 45001)
    c = [np.sum((fa - x * fb) ** 2 / (0.01 + x**2 * 0.01)) for x in brute]
    assert r["ratio"] == pytest.approx(brute[int(np.argmin(c))], abs=2e-4)
    assert r["chi2"] == pytest.approx(min(c), rel=1e-6)
    assert not r["ratio_at_grid_edge"]
    # an image a without flux is not measurable, as an image b without flux
    assert s7v.pair_chi2(-np.abs(fb), np.full(4, 0.1), np.abs(fb), np.full(4, 0.1))["p"] is None


def test_recentre_and_aperture_noise():
    yy, xx = np.mgrid[-50:51, -50:51] * 0.04  # ±2", as the real stamps
    hp = np.exp(-((xx - 0.04) ** 2 + yy**2) / 0.002)
    assert s7v.recentre(hp, xx, yy) == pytest.approx((0.04, 0.0))
    rng = np.random.default_rng(1)
    noise = rng.normal(0, 1.0, xx.shape)
    s = s7v.aperture_noise(noise, xx, yy, np.random.default_rng(2))
    n_pix = np.sum(np.hypot(xx, yy) <= s7v.R_AP)
    assert 0.6 * np.sqrt(n_pix) < s < 1.4 * np.sqrt(n_pix)


def test_perturber_flags_tri_state():
    assert s7v.classify_perturber(False, 18.0) is False
    assert s7v.classify_perturber(True, float("nan")) is None
    assert s7v.classify_perturber(True, 19.7) is True
    assert s7v.classify_perturber(True, 23.0) is False
    assert s7v.tri_any([False, None, True]) is True
    assert s7v.tri_any([False, None]) is None
    assert s7v.tri_any([False]) is False
    assert s7v.decide(False, False, None).endswith("member perturber not assessed)")
    assert s7v.decide(True, False, None) == "undetermined (member perturber not assessed)"
    assert s7v.classify_perturber(None, 18.0) is None
    assert s7v.classify_perturber(None, 27.0) is False  # too faint whatever its z


def test_grid_edge_is_not_measurable():
    fb = np.array([1.0, 2.0, 3.0])
    r = s7v.pair_chi2(1e4 * fb, np.full(3, 0.01), fb, np.full(3, 0.01))
    assert r["ratio_at_grid_edge"] is True
    assert s7v.differs(r["p"], r["ratio_at_grid_edge"]) is None
    assert s7v.pair_chi2([1.0], [0.1], [1.0], [0.1])["ratio_at_grid_edge"] is None

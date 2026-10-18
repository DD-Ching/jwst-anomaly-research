"""Tests for scripts/w1_r_sim.py (plan checks C1, C3 and the frozen rule)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_r_sim as w  # noqa: E402


def test_rule_outcome_threshold():
    assert w.rule_outcome(0.414, 0.0) == "conservative"
    assert w.rule_outcome(0.40, 0.01) == "conservative"  # within 2σ below
    assert w.rule_outcome(0.39, 0.01) == "optimistic"
    assert w.rule_outcome(1.0, 0.01, r_borrowed=1.05) == "optimistic"


def test_weighted_r_nearest_cell():
    cells = [
        {"semimajor_px_median": 2.0, "snr_median": 20.0, "R": {"k": {"R": 0.5, "err": 0.1}}},
        {"semimajor_px_median": 8.0, "snr_median": 200.0, "R": {"k": {"R": 1.0, "err": 0.1}}},
        {"semimajor_px_median": 4.0, "snr_median": 50.0, "R": {"k": {"R": None, "err": None}}},
    ]
    size = np.array([2.1, 1.9, 7.0, 9.0])
    snr = np.array([18.0, 25.0, 150.0, 300.0])
    r, err, counts = w.weighted_r(cells, size, snr, "k")
    assert counts == [2, 2]  # the cell without R is skipped
    assert r == pytest.approx(0.75)
    assert err == pytest.approx(np.sqrt(2 * 0.5**2 * 0.1**2))


def test_c1_noiseless_response_is_one():
    pytest.importorskip("galsim")
    pytest.importorskip("photutils")
    c = w.cell((0.5, 8.0, np.inf, 10, 1, 0.0, 0.0, {}, False))
    assert c["psf_sigma_px"]["stamp"] == 0.0
    assert c["R"]["stamp"]["R"] == pytest.approx(1.0, abs=0.03)


def test_c3_sign_and_psf_subtraction():
    galsim = pytest.importorskip("galsim")
    pytest.importorskip("photutils")
    the_psf = w.psf()
    ps = w.psf_sigma_px(the_psf)
    gal = galsim.Sersic(n=1, half_light_radius=6 * w.PIX_ARCSEC, flux=1.0)
    im, conv = w.draw(gal, the_psf, 0.0, None)
    m = w.measure(im, conv, 0.0)
    assert abs(w.to_eps(m[0], m[1], m[2], ps)) < 0.02  # C2: a round source stays round
    for g1, g2, part in ((0.05, 0.0, "real"), (0.0, 0.05, "imag")):
        im, conv = w.draw(gal.shear(g1=g1, g2=g2), the_psf, 0.0, None)
        m = w.measure(im, conv, 0.0)
        assert getattr(w.to_eps(m[0], m[1], m[2], ps), part) > 0.02

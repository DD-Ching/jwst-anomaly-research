"""Tests for scripts/w1_r_sim.py (plan checks C1, C3 and the frozen rule)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_r_sim as w  # noqa: E402


def test_rule_outcome_threshold():
    assert w.rule_outcome(0.414, 0.0, 0.414) == "conservative"
    assert w.rule_outcome(0.40, 0.01, 0.414) == "conservative"  # within 2σ below
    assert w.rule_outcome(0.39, 0.01, 0.414) == "optimistic"
    assert w.rule_outcome(None, None, 0.414) == "undetermined"
    assert w.rule_outcome(1.0, 0.01, r_borrowed=1.05) == "optimistic"


def _cell(size, snr, r, n_used, n_pairs=100):
    return {
        "semimajor_px_median": size,
        "snr_median": snr,
        "n_pairs": n_pairs,
        "R": {"k": {"R": r, "err": 0.1 if n_used > 1 else None, "n_pairs_used": n_used}},
    }


def test_weighted_r_nearest_cell_and_unrepresented_rows():
    cells = [_cell(2.0, 20.0, 0.5, 50), _cell(8.0, 200.0, 1.0, 100), _cell(4.0, 50.0, None, 0)]
    size = np.array([2.1, 1.9, 7.0, 9.0, 4.0])
    snr = np.array([18.0, 25.0, 150.0, 300.0, 50.0])
    res = w.weighted_r(cells, size, snr, "k")
    assert res["rows_per_cell"] == [2, 2, 1]  # the row nearest the R-less cell stays there
    assert res["n_rows_used"] == 4 and res["n_rows_unrepresented"] == 1
    assert res["R"] == pytest.approx(0.75)
    assert res["err"] == pytest.approx(np.sqrt(2 * 0.5**2 * 0.1**2))
    assert res["n_rows_extrapolated"] == 3  # size 1.9, size 9 / S/N 300 and S/N 18 lie outside
    res = w.weighted_r(cells, size, snr, "k", min_used_frac=0.9)  # cell 0 used 50 of 100
    assert res["n_rows_used"] == 2 and res["R"] == pytest.approx(1.0)
    assert res["n_rows_extrapolated"] == 1  # only the 9 px row counts; cell 0 rows are unused
    empty = w.weighted_r([_cell(None, None, None, 0)], size, snr, "k")
    assert empty["R"] is None and empty["n_rows_unrepresented"] == len(size)


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

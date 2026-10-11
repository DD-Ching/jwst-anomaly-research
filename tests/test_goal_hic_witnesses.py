"""Offline tests for scripts/goal_hic_witnesses.py (Bancal 2012 vs Li et al. 2026 LC4/LC5)."""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("goal_hic_witnesses", _DIR / "goal_hic_witnesses.py")
hw = importlib.util.module_from_spec(_spec)
sys.modules["goal_hic_witnesses"] = hw
_spec.loader.exec_module(hw)


@pytest.mark.parametrize("make", hw.WITNESSES)
def test_positive_controls_reproduce_the_papers(make):
    w = make()
    rho = hw.density(w)
    assert hw.witness_value(rho, w) == pytest.approx(w.quoted["quantum"], abs=1e-4)
    vis = hw.bisect(lambda v: hw.witness_value(hw.white(rho, v), w) - w.bound, 0.0, 1.0)
    assert vis == pytest.approx(w.quoted["white_threshold"], abs=1e-4)
    # the Experiment-v1 estimator is unbiased
    mom = hw.v1_estimator_moments(hw.distribution(rho, w), w)
    assert mom["mean"] == pytest.approx(w.quoted["quantum"], abs=1e-4)


def test_symmetric_readout_scales_k_body_terms():
    # Li et al. Eq. S100-S101: a k-body correlator is multiplied by (1 - 2 p)^k
    w = hw.lc4_witness()
    rho = hw.density(w)
    p = 0.03
    got = hw.witness_value(rho, w, hw.readout_affine(w, flip=p))
    two, three = (1 - 2 * p) ** 2, (1 - 2 * p) ** 3
    # ideal contributions: A0B0, A0B1 = 1/sqrt2 and 2 C0D0 = 2 (two-body); A1B0D0, -A1B1D0 = 1/sqrt2
    # and 2 A0C1D1 = 2 (three-body)
    s = 1 / math.sqrt(2)
    want = two * (s + s + 2) + three * (s + s + 2)
    assert got == pytest.approx(want, abs=1e-12)


def test_generic_lp_reproduces_round_2_signal_bound():
    w = hw.bancal_witness()
    res = hw.hic_lp(hw.distribution(hw.density(w), w), w)
    assert res["min_max_signal"] == pytest.approx(0.0062934355, abs=1e-8)


def test_lc4_witness_is_lp_optimal_and_signal_larger_than_bancal():
    w = hw.lc4_witness()
    assert hw.lp_white_threshold(w) == pytest.approx(w.quoted["white_threshold"], abs=1e-3)
    sig = hw.hic_lp(hw.distribution(hw.density(w), w), w)["min_max_signal"]
    assert sig == pytest.approx((4 + 2 * math.sqrt(2) - 6) / 32, abs=1e-8)


def test_noise_ranking_lc4_beats_lc5_except_white_noise():
    t4 = hw.noise_thresholds(hw.lc4_witness())
    t5 = hw.noise_thresholds(hw.lc5_witness())
    assert t5["white_visibility_min"] < t4["white_visibility_min"]
    assert t4["dephasing_p_max"] > t5["dephasing_p_max"]
    assert t4["readout_flip_max"] > t5["readout_flip_max"]
    assert t4["detection_efficiency_min"]["eta_min"] < t5["detection_efficiency_min"]["eta_min"]


def test_detection_model_is_ideal_at_unit_efficiency():
    w = hw.lc4_witness()
    rho = hw.density(w)
    assign = {(p, s): -1 for p in w.parties for s in (0, 1)}
    ro = hw.readout_affine(w, eta=1.0, assign=assign)
    assert hw.witness_value(rho, w, ro) == pytest.approx(hw.witness_value(rho, w))

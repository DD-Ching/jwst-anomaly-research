"""Offline tests for scripts/goal_bancal_channel.py (Bancal et al. 2012 timing and signalling)."""

from __future__ import annotations

import importlib.util
import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location(
    "goal_bancal_channel", _DIR / "goal_bancal_channel.py"
)
bc = importlib.util.module_from_spec(_spec)
sys.modules["goal_bancal_channel"] = bc
_spec.loader.exec_module(bc)


@pytest.fixture(scope="module")
def pq():
    return bc.quantum_distribution()


def test_positive_control_reproduces_the_paper(pq):
    assert np.linalg.norm(bc.quantum_state()) == pytest.approx(1.0)
    assert bc.bell_value(pq) == pytest.approx(7.2014, abs=1e-4)  # App. B
    # QM is no-signalling: the BCD marginal does not depend on x, the ABC marginal not on w
    for y, z, w in itertools.product((0, 1), repeat=3):
        assert np.allclose(pq[:, :, :, :, 0, y, z, w].sum(0), pq[:, :, :, :, 1, y, z, w].sum(0))
    for x, y, z in itertools.product((0, 1), repeat=3):
        assert np.allclose(pq[:, :, :, :, x, y, z, 0].sum(3), pq[:, :, :, :, x, y, z, 1].sum(3))
    assert bc.quantum_is_no_signalling(pq)
    signalling = pq.copy()
    signalling[0, 0, 0, 0, 0] += 0.01  # B, C, D outcome now depends on x
    signalling[1, 0, 0, 0, 0] -= 0.02
    assert not bc.quantum_is_no_signalling(signalling)


def test_lp_reproduces_lemma_1_and_bounds_the_signal(pq):
    res = bc.signalling_lp(pq)
    assert not res.ns_feasible  # Lemma 1: S > 7 with local BC|AD forces signalling
    assert res.min_bcd_on_x == pytest.approx(0.0, abs=1e-9)  # either direction alone can be silent
    assert res.min_abc_on_w == pytest.approx(0.0, abs=1e-9)
    assert res.min_of_max > 1e-3
    assert res.min_of_sum >= res.min_of_max - 1e-12
    assert res.min_bcd_on_x_if_abc_silent >= res.min_of_sum - 1e-9
    # fixing the best setting pair can only help the receiver
    assert res.min_max_per_setting >= res.min_of_max - 1e-12


@pytest.mark.parametrize("r", [2.0, 10.0, 1e4, 1e9])
def test_fig3_layout_meets_the_conditions(r):
    lay = bc.fig3_layout(r)
    assert all(bc.conditions(lay, r).values())
    ch = bc.channel_timing(lay)
    assert ch["both_superluminal_iii"]


def test_fig3_at_20_ly_needs_15_years_for_large_v():
    ch = bc.channel_timing(bc.fig3_layout(1e9))
    assert ch["A_to_Dprime"]["latency_from_sender"] * 20 == pytest.approx(15.0, abs=1e-6)
    assert ch["D_to_Aprime"]["latency_from_sender"] * 20 == pytest.approx(15.0, abs=1e-6)


def test_no_guaranteed_layout_beats_half_the_light_time():
    g = bc.best_guaranteed_latency(1e9, n=101)
    assert g["analytic_lower_bound_L_over_c"] <= g["worst_latency_L_over_c"] < 0.52


def test_clustered_layout_is_fast_one_way_only():
    r = 1e9
    lay = bc.clustered_layout(r)
    assert all(bc.conditions(lay, r).values())
    ch = bc.channel_timing(lay)
    assert ch["A_to_Dprime"]["latency_from_sender"] < 0.01
    assert ch["D_to_Aprime"]["latency_from_sender"] > 0.99


def test_statistics_helpers():
    assert bc.trials_for_error(0.0, 1e-3) == float("inf")
    assert bc.trials_for_error(0.1, np.exp(-1)) == pytest.approx(200.0)
    assert bc.bc_sync_window_s(1.0, 1.0) == pytest.approx(365.25 * 86400)


def test_r_fraction_matches_closed_forms():
    # rest frame: P(|t_B - t_C| < d / r) = 1 - (1 - w / T)^2 for uniform times on [0, T]
    r, T = 10.0, 0.5
    w = 1.0 / r
    assert bc.r_fraction(r, 0.0, 1.0, T) == pytest.approx(1 - (1 - w / T) ** 2, rel=1e-9)
    # worst case over frames with speed <= beta_max: 1 / (2 r beta_max) at leading order
    row = bc.r_fraction_table((5e6,), (bc.BETA_CMB,))[0]
    assert row["fraction"] == pytest.approx(1 / (2 * 5e6 * bc.BETA_CMB), rel=1e-2)
    # no positive guarantee without a bound on the frame speed: it falls as 1 - beta
    near_c = bc.worst_r_fraction(5e6, 0.999, 2 * 0.999, n_angle=21)["fraction"]
    assert near_c < 3e-10


def test_worst_r_fraction_grid_is_the_minimum():
    rng = np.random.default_rng(1)
    for bm in (bc.BETA_CMB, 0.5):
        win = 2 * bm
        worst = bc.worst_r_fraction(1e3, bm, win, n_angle=41)["fraction"]
        for _ in range(40):
            f = bc.r_fraction(1e3, rng.uniform(0, bm), rng.uniform(-1, 1), win)
            assert f >= worst * (1 - 1e-9)

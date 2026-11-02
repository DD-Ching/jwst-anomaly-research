"""Offline tests for scripts/fgm_oneway_bit.py (one-way bit comparison, a few seconds)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import fgm_oneway_bit as ob  # noqa: E402


def test_readout_thresholds():
    n = ob.helstrom_threshold(0.1)
    assert math.isclose(n, 0.2554, abs_tol=5e-5)
    assert math.isclose(ob.helstrom_error(n), 0.1, rel_tol=1e-12)
    assert math.isclose(ob.pilot_threshold(0.1), math.log(5) / 2, rel_tol=1e-12)
    assert math.isclose(ob.thermal_requirement_min(0.0), n, rel_tol=1e-12)


def test_hankel_tail_envelope_and_halo_formula():
    for q0 in (150.0, 310.0, 600.0):
        q = np.linspace(q0 - 6, q0 + 6, 1201)
        ratio = np.max(np.abs(ob.hankel_bump(q))) / ob.hankel_envelope(q0)
        assert 0.9 < ratio < 1.3
    q = np.arange(0.0, 400.0, 0.01)
    parseval = np.trapezoid(ob.hankel_bump(q) ** 2 * q, q) / ob.aperture_norm()
    assert math.isclose(parseval, 1.0, rel_tol=1e-4)
    qc = 200.0
    q = np.arange(qc, qc + 700.0, 0.02)
    direct = np.trapezoid(ob.hankel_bump(q) ** 2 * q, q) / ob.aperture_norm()
    assert math.isclose(direct, math.exp(ob.halo_fraction_log(qc)), rel_tol=0.15)


def test_emergence_map_matches_fgm_window_form():
    k = 0.25
    # v_out = -kappa^-1 ln(e^{kappa v_max} - e^{kappa v}) with v_max = 0 here
    for dv in (-100.0, -10.0, -2.0):
        assert math.isclose(ob.emergence_offset(dv, k), -math.log(1 - math.exp(k * dv)) / k)
    # equal energy at entry and exit at V = |DeltaV| / 2
    assert math.isclose(math.expm1(k * math.log(2) / k), 1.0)


def test_time_weight_and_no_overlap():
    w = ob.time_weight(-795.0, 0.25)
    assert w > math.exp(0.25 * 795.0)  # the early half of the profile dominates
    res = ob.run()
    iv = res["intervals"]
    assert iv["log10_N_allowed_max_tau_iso"] < math.log10(iv["N_required_min"])
    assert iv["log10_gap_tau_iso_over_critical"] > 50
    assert iv["log10_spectral_term_over_halo_term"] < -50
    assert res["crossover"]["ratio_to_fixed_carrier"] > 5
    assert res["paths"]["wormhole_minus_external_at_B"] > 1000

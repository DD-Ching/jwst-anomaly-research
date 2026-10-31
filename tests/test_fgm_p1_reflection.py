"""Offline tests for scripts/fgm_p1_reflection.py (P1 reflection bound, a few seconds)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import fgm_p1_reflection as pr  # noqa: E402
import fgm_window as fw  # noqa: E402

STORED = json.loads(pr.WINDOW_JSON.read_text())


def test_reproduces_critical_fraction_from_stored_json():
    crit = STORED["closure_into_B"]["dE_crit_centre"]
    e_low = STORED["packets"]["P1_gaussian"]["E_below_1"]
    assert math.isclose(crit / e_low, 1.903e-79, rel_tol=1e-3)
    # the stored threshold follows from FGM (2.20) and (3.3) at the stored arrival label
    assert math.isclose(
        fw.de_crit(STORED["closure_into_B"]["u_B_minus_x0"], 20.0), crit, rel_tol=1e-9
    )
    assert math.isclose(pr.band_energy(0.0, 1.0), e_low, rel_tol=1e-6)


def test_tail_flux_matches_endpoint_asymptote():
    # |G(dv)|^2 with G = (2 pi)^-1/2 int sqrt(w) f(w) e^{-i w dv} dw, f = sqrt(density)
    w = np.linspace(0.0, 1500.0, 1_500_001)
    dens = np.exp(-((w - pr.W0) ** 2) / (2 * pr.SW**2)) / (
        math.sqrt(2 * math.pi) * pr.SW * pr.gauss_norm()
    )
    amp = np.sqrt(w * dens)
    for dv in (60.0, 150.0):
        g = np.trapezoid(amp * np.exp(-1j * w * dv), w) / math.sqrt(2 * math.pi)
        assert math.isclose(abs(g) ** 2, pr.tail_flux_coeff() / dv**3, rel_tol=0.1)


def test_angular_weights():
    a = pr.angular_factors(math.pi)
    assert a["H_decreasing_in_separation"]
    assert math.isclose(a["any_deposit_min"], 1.0)  # the string sits at the kernel minimum
    assert math.isclose(a["uniform_deposit"], 1 / (4 * math.pi * fw.kappa() * a["H_over_r2"]))
    assert a["uniform_deposit"] > a["equatorial_deposit"] > 1.0


def test_setting_s_bound_and_vortex_estimate():
    res = pr.run()
    s = res["setting_S"]
    assert s["ratio_with_mean_delay"] > 1e40 and s["break_even_mean_delay"] > 600
    assert s["future_tail_fraction_after_emergence"] < 1e-4
    rows = res["vortex_completion"]["rows"]
    logs = [r["log10_transfer_max"] for r in rows]
    assert all(x > y for x, y in zip(logs, logs[1:], strict=False))
    assert res["vortex_completion"]["smallest_listed_m_below_threshold"] == 13
    assert res["single_run_S"]["relative_fluctuation_of_E_B"] > 1e7


def test_break_even_lead_is_a_root():
    crit = STORED["closure_into_B"]["dE_crit_centre"]
    t = 10 ** pr.vortex_transfer_log10(13)
    x = pr.break_even_lead(t, crit, 0.25, pr.tail_flux_coeff())
    assert abs(pr.log10_tail_ratio(x, t, crit, 0.25, pr.tail_flux_coeff())) < 1e-6

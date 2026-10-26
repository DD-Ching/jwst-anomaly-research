"""Offline tests of scripts/w3_kepler_physical.py (results/w3_kepler/physical/PLAN.md)."""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import w3_kepler_physical as wp  # noqa: E402

from jwst_anomaly import exotic_sim as es  # noqa: E402


def test_theta_star():
    assert abs(wp.theta_star_uas(1.0, 1000.0) - 4.6505) < 1e-3  # 1 R_sun at 1 kpc


def test_cell_map_is_consistent_with_the_model_einstein_radius():
    r, d, rho, te, v = 1.1, 600.0, 0.1, 1.0, 30.0
    m = wp.cell_map(r, d, rho, te, v)
    assert abs(m["t_star_d"] - rho * te) < 1e-12
    assert 0 < m["x"] < 1
    pc = 3.0857e16
    d_l, d_s = m["D_L_pc"] * pc, d * pc
    th_model = es.einstein_radius(es.eps_bar_point_mass(m["M_msun"]), 1, d_l, d_s, d_s - d_l)
    th_model_uas = th_model * 180 / math.pi * 3600e6
    assert abs(th_model_uas / m["theta_E_uas"] - 1) < 2e-3  # |M| and D_L reproduce theta_E
    v_check = m["mu_rel_mas_yr"] * m["D_L_pc"] / 1000 * 4.74047
    assert abs(v_check / v - 1) < 2e-3  # mu_rel * D_L = v


def test_unreachable_cell_has_no_mass():
    m = wp.cell_map(1.0, 600.0, 0.1, 10.0, 220.0)  # t_* = 1 d needs x = v t_* / R_* >> 1
    assert m["x"] > 1 and "M_msun" not in m


def test_detections_per_density_analytic():
    n = wp.detections_per_density([1000.0], [1.0], 1.0, 365.25)
    tau = wp.FOUR_PI_G_C2 * 1000.0**2 / 6
    assert abs(n / (4 / math.pi * tau) - 1) < 1e-9
    assert abs(wp.detections_per_density([1000.0], [1.0], 0.5, 36.525) / n - 5) < 1e-9

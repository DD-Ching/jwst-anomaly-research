"""Offline tests for scripts/goal_krasnikov_g2.py (gap G2: ANEC on the Krasnikov tube's rays)."""

import importlib.util
import sys
from dataclasses import replace
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("goal_krasnikov_g2", _DIR / "goal_krasnikov_g2.py")
g2 = importlib.util.module_from_spec(_spec)
sys.modules["goal_krasnikov_g2"] = g2
_spec.loader.exec_module(g2)


def test_outgoing_lines_are_null_pregeodesics_with_zero_null_ricci():
    out = g2.outgoing_lines()
    assert out["g_xx"] == "0"
    assert out["pregeodesic"]
    assert out["R_xx"] == "0"  # for an arbitrary k(u, x, rho): front and end caps included


def test_region_with_k_of_x_only_is_flat():
    assert g2.flat_where_k_depends_on_x_only()["nonzero_riemann_components"] == []


def test_outgoing_speed_never_exceeds_light():
    for k in (1.0, 0.3, 1e-3):
        assert g2.max_outgoing_speed(k, n=4000) <= 1.0


def test_incoming_core_rays_stay_in_the_flat_region():
    p = replace(g2.P, d_len=60.0, rho_max=10.0, n_rays=8)
    res = g2.incoming_rays(p)
    assert res["max_dk_drho_or_du_along_rays"] == 0.0
    # they are the shortcut: transit far below the exterior light time
    assert res["transit_time_max_eps"] < 0.5 * res["exterior_light_time_eps"]
    # positive control: a ray sent before completion crosses the front (k varies with u there)
    early = g2.er_k(30.0 + 0.5, 30.0, 0.0, p) - g2.er_k(30.0 + 0.6, 30.0, 0.0, p)
    assert early != 0.0

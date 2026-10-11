"""Offline tests for scripts/goal_bh_wormhole.py (beyond-Horndeski wormhole stability audit)."""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("goal_bh_wormhole", _DIR / "goal_bh_wormhole.py")
bh = importlib.util.module_from_spec(_spec)
sys.modules["goal_bh_wormhole"] = bh
_spec.loader.exec_module(bh)

XS = bh.grid(401)


def test_positive_controls_reproduce_the_paper():
    pc = bh.positive_controls(XS)
    assert pc["F_min"] > 1.0  # 2212.05969 below Eq. (72): F > H = 1
    # Fig. 2 of 2212.05969: K11, det K, G11, det G, M11^(l2), det M^(l2) > 0 everywhere
    for k in ("K11_min_l10", "detK_min_l10", "G11_min_l10", "detG_min_l10", "m11_min_l10"):
        assert pc[k] > 0
    assert pc["detM_l2_min"] > 0
    assert pc["no_ghost_2P1_minus_F_min"] > 0 and pc["N_min"] > 0
    # the leading order in l of Eq. (33) equals Appendix B (M12, M22 tightly; M11 up to points
    # where it vanishes)
    assert pc["eq33_vs_appendixB_max_rel_diff_l1e4"][1] < 1e-5
    assert pc["eq33_vs_appendixB_max_rel_diff_l1e4"][2] < 1e-3


def test_angular_mode_is_superluminal_and_tensor_mode_luminal_far_away():
    sp = bh.angular_speeds(XS)
    assert sp["max_c2"] > 1.5  # 2212.05969 Sec. 4.2 reports superluminal angular propagation
    assert sp["edge_values_c2"][1] == pytest.approx(1.0, abs=1e-3)  # flat far zone: luminal


def test_results_do_not_depend_on_the_throat_scale():
    a = bh.angular_speeds(XS)
    b = bh.angular_speeds(XS, tau_value=10**4)
    assert a["max_c2"] == pytest.approx(b["max_c2"], rel=1e-9)


def test_low_l_sufficient_conditions_fail_and_odd_sector_passes():
    rows = {r["l"]: r for r in bh.low_l_conditions(XS, ells=(2, 50))}
    lo, hi, _ = rows[2]["M22_le_0"]
    assert lo < 0 < hi < 0.6  # (39b) fails around the throat for l = 2
    assert rows[2]["cond39c_fail_fraction_inner_|r|<3tau"] > 0.1
    assert rows[2]["odd_tachyon_le_0"] is None and rows[50]["odd_tachyon_le_0"] is None


def test_crossing_time_integral():
    t = bh.crossing_time()
    assert t["throat_radius_planck"] == pytest.approx(100 * math.log(3))
    assert t["I"] == pytest.approx(1.10661, abs=1e-4)

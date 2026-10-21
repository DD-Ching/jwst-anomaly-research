"""Offline tests for scripts/x4_s1g.py (X4 S1g): pass rule, outcome rule, validation generator and
budget change. No lalsuite needed (lalsuite parts of the script are not exercised here)."""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pytest

_S = Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _S / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


g = _load("x4_s1g")


# ------------------------------------------------------------------ pass rule
@pytest.mark.parametrize(
    ("eps", "ff_case", "want"),
    [
        (4.9e-4, None, "pass"),
        (5e-4, None, "pass"),  # <= is inclusive
        (5.01e-4, None, "fail"),
        (1.5e-3, None, "fail"),  # S1e's V9 seed 11
        (2e-4, 0.999, "fail"),  # 0.1 * (1 - 0.999) = 1e-4 is the stricter limit
        (9e-5, 0.999, "pass"),
        (4e-4, 0.945, "pass"),  # limit 5.5e-3; 5e-4 binds
        (6e-4, 0.945, "fail"),
        (None, None, "unfinished"),
        (None, 0.95, "unfinished"),
        (float("nan"), None, "fail"),
    ],
)
def test_fit_status(eps, ff_case, want):
    assert g.fit_status(eps, ff_case) == want


def test_fit_status_matches_s1e_record():
    """The S1e rows, scored by the function: V9 seed 11 fails, the other three P fits pass."""
    import json

    rows = json.loads((g.OUT / "s1e_validate.json").read_text())["rows"]
    st = {(r["name"], r["seed"]): g.fit_status(r["eps_final"]) for r in rows if r["model"] == "P"}
    assert st == {("V5", 11): "pass", ("V5", 12): "pass", ("V9", 11): "fail", ("V9", 12): "pass"}
    assert g.tool_verdict(st.values()) == "fail"


def test_tool_verdict():
    assert g.tool_verdict(["pass"] * 6) == "pass"
    assert g.tool_verdict(["pass", "fail", "unfinished"]) == "fail"
    assert g.tool_verdict(["pass", "unfinished"]) == "unfinished"
    assert g.tool_verdict([]) == "unfinished"


def test_seeds_agree():
    assert g.seeds_agree([0.95, 0.9499, 0.90])
    assert g.seeds_agree([0.95, 0.9498])  # exactly 2e-4 (inclusive, up to rounding)
    assert not g.seeds_agree([0.95, 0.9497])
    assert not g.seeds_agree([0.95])
    assert not g.seeds_agree([0.95, None])
    # S1f's unlensed seeds: 0.91351 vs 0.89229 disagree
    assert not g.seeds_agree([0.8922899445342309, 0.8922911658790608, 0.9135105524586078])


# ------------------------------------------------------------------ R1 outcome rule
U_S1F = [0.8922899445342309, 0.8922911658790608, 0.9135105524586078]


def test_r1_outcome_counterexample_needs_no_tool():
    o = g.r1_outcome("fail", [0.996], U_S1F, 0.996, True)
    assert o["outcome"] == 1 and o["rho_req"] > 30


def test_r1_outcome_u_reliability_blocks_outcome_2():
    o = g.r1_outcome("pass", [0.95, 0.9499, 0.94], U_S1F, 0.95, False)
    assert o["outcome"] == 3
    assert o["causes"] == [
        "U search reliability (S1f, not re-run): two best seeds differ by > 2e-4"
    ]


def test_r1_outcome_2_and_causes():
    u_ok = [0.9135, 0.9134]
    assert g.r1_outcome("pass", [0.95, 0.9499], u_ok, 0.95, False)["outcome"] == 2
    o = g.r1_outcome("fail", [0.95, 0.94], u_ok, 0.95, False, p_unfinished=1)
    assert o["outcome"] == 3 and len(o["causes"]) == 3
    # unverified FF above the screen is not outcome 1
    o = g.r1_outcome("pass", [0.996, 0.996], u_ok, 0.996, False)
    assert o["outcome"] == 3
    # rho_req <= 30 but intrinsic SNR > 100: reported as the feasibility condition
    o = g.r1_outcome("pass", [0.99, 0.99], u_ok, 0.99, False)
    assert o["outcome"] == "feasibility" and o["intrinsic_snr_needed"] > 100


def test_screen_threshold():
    assert g.FF_SCREEN == pytest.approx(0.994987, abs=1e-6)
    assert g.rho_req(g.FF_SCREEN) == pytest.approx(30.0, rel=1e-9)


# ------------------------------------------------------------------ validation generator
def test_validation_draws_deterministic_and_in_region():
    a, b2 = g.draw_validation(6), g.draw_validation(6)
    assert a == b2
    assert g.draw_validation(2) == a[:2]  # the first draws do not depend on how many are drawn
    r = g.VAL_REGION
    for z in a:
        x = np.asarray(z["x_box"])
        assert g.in_box(x, z["lens"])
        assert r["dln_mc"][0] <= x[0] - np.log(g.m.MC0) <= r["dln_mc"][1]
        assert r["q"][0] <= x[1] <= r["q"][1]
        assert 0.0 <= x[2] <= 0.99 and 0.0 <= x[5] <= 0.99
        assert 0.0 <= x[8] <= np.pi
        assert 300.0 <= z["lens"][0] <= 3000.0 and 0.1 <= z["lens"][1] <= 1.0
        assert 0.0 <= z["psi"] <= np.pi
        # the box mapping of the truth round-trips (the formal path, without lalsuite)
        (mc, q, s1, s2, inc, phi, _psi), _lens = g.spec_from_draw(z)
        assert np.allclose(g.e.from_cart(mc, q, s1, s2, inc, phi), x, atol=1e-12)


def test_validation_draws_differ_from_development_truths():
    v9 = [1373.090166312427, 0.22464965696121997]
    for z in g.draw_validation(2):
        assert abs(np.log(z["lens"][0] / 1000.0)) > 0.05  # V5 lens 1,000 / 0.2
        assert abs(np.log(z["lens"][0] / v9[0])) > 0.05


def test_in_box_rejects_outside():
    x = g.draw_validation(1)[0]["x_box"]
    assert not g.in_box(x, [2e4, 0.2])
    assert not g.in_box(x, [500.0, 6.0])
    bad = list(x)
    bad[1] = 16.0
    assert not g.in_box(bad, None)


# ------------------------------------------------------------------ budget change only
def test_apply_budget_changes_only_the_two_budgets():
    j = g.j
    src = (_S / "x4_s1e_joint.py").read_text()
    assert re.search(r'^DE_P = \{"popsize": 12, "maxiter": 400\}$', src, re.M)
    assert re.search(r"^MAX_ROUNDS = 10$", src, re.M)
    bounds_before = list(j.BOUNDS_P)
    g.apply_budget()
    assert j.DE_P == {"popsize": 12, "maxiter": 1000}
    assert j.MAX_ROUNDS == 20
    assert j.DE_U == {"popsize": 12, "maxiter": 600}
    assert j.TOL == 1e-5 and j.BOUNDS_P == bounds_before
    g.apply_budget()  # idempotent
    assert j.DE_P["maxiter"] == 1000


def test_frozen_schedule():
    names = [n for rnd in g.ROUNDS for n, _ in rnd]
    assert [len(rnd) for rnd in g.ROUNDS] == [3, 3]
    assert ("V9", 11) in g.ROUNDS[0]  # S1e's failing fit runs in round 1
    assert sorted(set(names)) == ["V5", "V9", "VA1", "VA2"]
    assert {g.LABEL[n] for n in ("V5", "V9")} == {"development"}
    assert {g.LABEL[n] for n in ("VA1", "VA2")} == {"validation"}
    assert g.R1_SEEDS == (21, 22, 23)
    assert g.N_PROC <= 3

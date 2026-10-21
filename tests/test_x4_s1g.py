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
    assert g.tool_verdict(["pass", "error", "unfinished"]) == "error"
    assert g.tool_verdict(["fail", "error"]) == "fail"
    assert g.tool_verdict([]) == "unfinished"


def test_row_status():
    assert g.row_status({"status": "finished", "eps_final": 1e-5}) == "pass"
    assert g.row_status({"status": "finished", "eps_final": 1e-5}, 0.9999) == "fail"
    assert g.row_status({"status": "error", "error": "x"}) == "error"
    assert g.row_status({"status": "unfinished", "FF_pad64": 0.99999}) == "unfinished"


class _Raises:
    def get(self):
        raise RuntimeError("boom")


class _Returns:
    def get(self):
        return {"status": "finished", "name": "V5"}


def test_collect_turns_an_exception_into_an_error_row():
    job = ("valid", "VA1", 11, None, "run-x")
    row = g.collect(job, _Raises())
    assert row["status"] == "error" and "boom" in row["error"]
    assert (row["name"], row["seed"], row["set"], row["run_id"]) == (
        "VA1",
        11,
        "validation",
        "run-x",
    )
    assert g.collect(job, _Returns())["status"] == "finished"


def test_stale_state_is_ignored(tmp_path, monkeypatch):
    monkeypatch.setenv("X4_S1G_STATE", str(tmp_path))
    monkeypatch.setattr(g, "_RUN", {"key": "valid_V9_s11", "run_id": "run-old"})
    g._write_state({"stage": "de", "nit": 11, "f": 0.1, "x": [0.0] * 12})
    assert g.read_state("valid_V9_s11", "run-old")["nit"] == 11
    assert g.read_state("valid_V9_s11", "run-new") is None  # an earlier run's file is not used
    assert g.read_state("valid_V5_s11", "run-old") is None


def test_refine_state_keeps_the_best_exact_lens_point():
    x, lens = [1.0, 2.0], [1400.0, 0.2]
    better = g.refine_state({"f": 1e-4, "x": [1.1, 2.1]}, x, 2e-4, lens, 3)
    assert better["improves_exact"] and better["f"] == 1e-4 and better["x"] == [1.1, 2.1]
    same = g.refine_state({"f": 2e-4, "x": [1.1, 2.1]}, x, 2e-4, lens, 3)
    assert not same["improves_exact"] and same["f"] == 2e-4 and same["x"] == x
    assert same["lens"] == lens and same["round"] == 3


def test_r1_fresh_run_ignores_and_resume_checks_the_stored_record(tmp_path, monkeypatch):
    # a fresh run must not reuse committed R1 rows; resume refuses a record made with another budget
    import json

    rec = {"budget": {"DE_P_maxiter": 1, "MAX_ROUNDS": 1}, "rows": [], "saved_rescored": []}
    (tmp_path / "s1g_r1.json").write_text(json.dumps(rec))
    monkeypatch.setattr(g, "OUT", tmp_path)
    with pytest.raises(SystemExit, match="refused to resume"):
        g.main_r1("run-x", resume=True)

    class Fresh(Exception):
        pass

    def saved_best():
        raise Fresh

    monkeypatch.setattr(g.j, "saved_best", saved_best)
    with pytest.raises(Fresh):  # the fresh path builds a new record instead of reading the file
        g.main_r1("run-x")


def test_r1_resume_keeps_finished_rows():
    rows = [{"seed": 21, "status": "finished"}, {"seed": 22, "status": "unfinished"}]
    assert g.r1_seeds_to_run(rows) == [22, 23]
    assert g.r1_seeds_to_run([{"seed": s, "status": "finished"} for s in (21, 22, 23)]) == []


def test_formal_path_accepts_the_committed_validation_truths():
    import json

    rec = json.loads((g.OUT / "s1g_injections.json").read_text())
    for z in rec["accepted"]:
        assert z["round_trip_max_diff"] <= g.ROUND_TRIP_TOL
        assert g.formal_path_ok(z["in_box"], z["round_trip_max_diff"], z["objective_at_truth"])
    assert [z["name"] for z in rec["accepted"]] == ["VA1", "VA2"] and rec["rejected"] == []
    assert not g.formal_path_ok(True, 1e-8, 0.0)
    assert not g.formal_path_ok(True, 0.0, 1e-5)
    assert not g.formal_path_ok(False, 0.0, 0.0)


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
    # rho_req <= 30 but intrinsic SNR > 100, not verified: undetermined, cause named
    o = g.r1_outcome("pass", [0.99, 0.99], u_ok, 0.99, False)
    assert o["outcome"] == 3 and o["causes"] == ["intrinsic SNR needed > 100 but not verified"]


def test_r1_outcome_feasibility_is_constructive_like_outcome_1():
    """Frozen order: outcome 1, then the verified feasibility case, whatever causes exist."""
    o = g.r1_outcome("fail", [0.99, 0.95], U_S1F, 0.99, False, 1, feasibility_verified=True)
    assert o["outcome"] == "feasibility" and o["intrinsic_snr_needed"] > 100
    o = g.r1_outcome("fail", [0.996], U_S1F, 0.996, True, feasibility_verified=False)
    assert o["outcome"] == 1


def test_screen_flags():
    assert g.FF_FEASIBLE == pytest.approx(0.98790, abs=1e-5)
    assert g.rho_req(g.FF_FEASIBLE) / g.SNR_RATIO_R1 == pytest.approx(100.0, rel=1e-9)
    f = g.screen_flags
    assert f(0.996, 1e-7, True) == {
        "counterexample_verified": True,
        "feasibility_fail_verified": False,
    }
    assert f(0.990, 1e-7, True) == {
        "counterexample_verified": False,
        "feasibility_fail_verified": True,
    }
    assert f(0.9476, 0.0, True) == {
        "counterexample_verified": False,
        "feasibility_fail_verified": False,
    }
    assert not any(f(0.999, 0.0, False).values())  # illegal parameters never count
    assert f(0.99499, 1e-5, True)["counterexample_verified"] is False  # spread is subtracted
    assert g.is_legal(True, [0.99, 0.5], 2.8, None)
    assert not g.is_legal(True, [0.995, 0.5], 2.8, None)
    assert not g.is_legal(True, [0.5, 0.5], 15.5, None)
    assert not g.is_legal(True, [0.5, 0.5], 2.8, "RuntimeError()")


def test_stored_r1_summary_inputs_give_outcome_3():
    """R1's committed best (rho_req 9.4, intrinsic SNR 48.6) enters neither constructive branch."""
    import json

    r1 = json.loads((g.OUT / "s1g_r1.json").read_text())
    ver = r1["verify_best_P"]
    legal = g.is_legal(ver["in_box"], ver["chi_mags"], ver["q"], ver["waveform_error"])
    flags = g.screen_flags(ver["FF_conservative"], ver["score_spread"], legal)
    assert flags == {"counterexample_verified": False, "feasibility_fail_verified": False}
    p_ffs = [r["FF_pad64"] for r in r1["rows"] if r["status"] == "finished"]
    o = g.r1_outcome("pass", p_ffs, U_S1F, max(p_ffs), False, 0, feasibility_verified=False)
    assert o["outcome"] == 3 and len(o["causes"]) == 2
    assert o["rho_req"] == pytest.approx(9.394, abs=1e-3)
    assert o["intrinsic_snr_needed"] == pytest.approx(48.58, abs=1e-2)


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
def test_apply_budget_changes_only_the_two_budgets(monkeypatch):
    j = g.j
    # restore the shared x4_s1e_joint module after this test
    monkeypatch.setattr(j, "DE_P", dict(j.DE_P))
    monkeypatch.setattr(j, "MAX_ROUNDS", j.MAX_ROUNDS)
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

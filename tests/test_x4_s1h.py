"""Offline tests for scripts/x4_s1h.py (X4 S1h): frozen schedule, initialisation, pooling and the
decision rule. No lalsuite needed."""

from __future__ import annotations

import importlib.util
import json
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


h = _load("x4_s1h")
S1G_P = [0.9458756195092743, 0.9458833800179149, 0.9476377694069567]
S1F_U = [0.8922899445342309, 0.8922911658790608, 0.9135105524586078]


def test_frozen_schedule():
    assert h.SEEDS == (24, 25, 26)
    assert h.jobs_order() == [("P", 24), ("P", 25), ("P", 26), ("U", 24), ("U", 25), ("U", 26)]
    assert h.N_PROC == 4
    assert h.DEADLINE_UTC == "2026-10-21T08:15:00+00:00"


def test_initialisation_is_the_original_saved_best_not_the_s1g_s1f_bests():
    saved = json.loads((h.OUT / "s1g_r1.json").read_text())["saved_rescored"]
    bp, x0p = h.x0_for("P", saved)
    bu, x0u = h.x0_for("U", saved)
    assert (bp["start"], bu["start"]) == ("generic4", "generic5")
    assert bp["FF_pad64"] == pytest.approx(0.945393, abs=1e-6)
    assert bu["FF_pad64"] == pytest.approx(0.874989, abs=1e-6)
    assert x0p.size == 12 and x0u.size == 10
    assert np.allclose(np.exp(x0p[10:]), bp["lens"])
    # not the S1g P best (seed 21) or the S1f U best (seed 23)
    p21 = next(
        r for r in json.loads((h.OUT / "s1g_r1.json").read_text())["rows"] if r["seed"] == 21
    )
    u23 = next(
        r for r in json.loads((h.OUT / "s1f_case_u.json").read_text())["rows"] if r["seed"] == 23
    )
    assert np.max(np.abs(x0p[:10] - np.asarray(p21["x"]))) > 0.1
    assert np.max(np.abs(x0u - np.asarray(u23["x"]))) > 0.1


def test_seed_ffs_counts_finished_seeds_only():
    rows = [
        {"seed": 24, "status": "finished", "FF_pad64": 0.95},
        {"seed": 25, "status": "unfinished", "FF_pad64": 0.96},
        {"seed": 26, "status": "error"},
    ]
    assert h.seed_ffs(rows) == {24: 0.95}


def test_best_of_never_drops_a_better_old_result():
    c = [{"src": "old", "FF_pad64": 0.9476}, {"src": "new", "FF_pad64": 0.9470}, {"FF_pad64": None}]
    assert h.best_of(c)["src"] == "old"
    assert h.best_of([]) is None


def test_validation_recomputed_with_ff_case():
    eps = {"V9 s11": 8.06e-6, "VA2 s11": 2.79e-5}
    assert h.validation_statuses(eps, 0.9476) == {"V9 s11": "pass", "VA2 s11": "pass"}
    assert h.validation_statuses(eps, 0.9999) == {"V9 s11": "pass", "VA2 s11": "fail"}


def test_outcome_order_and_wording():
    o = h.s1h_outcome("fail", "pass", S1G_P, S1F_U, 0.996, True, False, 2)
    assert o["outcome"] == 1 and o["wording"] == "R1 does not pass the screen"
    o = h.s1h_outcome("fail", "pass", S1G_P, S1F_U, 0.99, False, True, 2)
    assert o["outcome"] == "feasibility"
    ok_p, ok_u = [0.95, 0.94995, 0.94], [0.92, 0.91995]
    o = h.s1h_outcome("pass", "pass", ok_p, ok_u, 0.95, False, False, 0)
    assert o["outcome"] == 2 and o["wording"] == h.ALLOWED_PASS_WORDING
    assert "wormhole" not in o["wording"]
    o = h.s1h_outcome("pass", "pass", ok_p, ok_u, 0.95, False, False, 1)
    assert o["outcome"] == 3 and o["wording"] == "tool/optimization undetermined"
    o = h.s1h_outcome("pass", "pass", S1G_P, S1F_U, 0.9476, False, False, 0)
    assert o["outcome"] == 3 and len(o["causes"]) == 2
    o = h.s1h_outcome("pass", "unfinished", ok_p, [0.92], 0.95, False, False, 0)
    assert o["outcome"] == 3
    assert any("fewer than 2 finished seeds" in c for c in o["causes"])
    assert any(c.startswith("U tool validation") for c in o["causes"])


def test_top_two_gap():
    assert h.top_two_gap(S1G_P)[2] == pytest.approx(1.7544e-3, abs=1e-7)
    assert h.top_two_gap([0.9]) == (0.9, None, None)


# ------------------------------------------------------------------ review fixes (2026-10-21)
def test_eps_limits_report_the_binding_term():
    lim = h.eps_limits(0.9476377694069567)
    assert lim["abs"] == 5e-4 and lim["case"] == pytest.approx(5.236e-3, abs=1e-6)
    assert lim["effective"] == 5e-4
    assert h.eps_limits(0.9999)["effective"] == pytest.approx(1e-5)


def test_row_statuses_keep_the_row_status_guard():
    rows = [
        {"name": "V9", "seed": 11, "status": "finished", "eps_final": 8e-6},
        {"name": "VA1", "seed": 11, "status": "error", "error": "x"},
        {"name": "VA2", "seed": 11, "status": "unfinished", "FF_pad64": 0.9999},
        {"name": "V6", "seed": 11, "eps_final": 2e-10},  # S1e rows carry no status
        {"name": "V7", "seed": 12},  # no final eps: never graded as finished
    ]
    assert h.row_statuses(rows, 0.95) == [
        ("V9 s11", "pass"),
        ("VA1 s11", "error"),
        ("VA2 s11", "unfinished"),
        ("V6 s11", "pass"),
        ("V7 s12", "unfinished"),
    ]


def test_tool_status_duplicates_counts_and_fail_precedence():
    exp = ["A s1", "B s1"]
    assert h.tool_status([("A s1", "pass"), ("B s1", "pass")], exp)["verdict"] == "pass"
    t = h.tool_status([("A s1", "pass"), ("A s1", "pass"), ("B s1", "pass")], exp)
    assert t["verdict"] == "inconsistent" and t["duplicates"] == ["A s1"]
    t = h.tool_status([("A s1", "pass")], exp)
    assert t["verdict"] == "inconsistent" and t["missing"] == ["B s1"]
    t = h.tool_status([("A s1", "fail"), ("A s1", "pass"), ("C s1", "pass")], exp)
    assert t["verdict"] == "fail" and t["extra"] == ["C s1"]  # a real fail is never hidden
    assert h.tool_status([("A s1", "error"), ("B s1", "pass")], exp)["verdict"] == "error"
    assert len(h.EXPECTED_P) == 6 and len(h.EXPECTED_U) == 4


def test_outcome_names_errors_and_unverified_best():
    ok_p, ok_u = [0.95, 0.94995], [0.92, 0.91995]
    o = h.s1h_outcome("pass", "pass", ok_p, ok_u, 0.95, False, False, 0, n_error=1)
    assert o["outcome"] == 3 and o["causes"] == ["fits that raised an error: 1"]
    o = h.s1h_outcome("pass", "pass", ok_p, ok_u, 0.95, False, False, 0, best_verified=False)
    assert o["causes"] == ["best ordinary solution not verified"]


def test_resume_guard_and_per_run_state_keys():
    rec = json.loads((h.OUT / "s1h_r1.json").read_text())
    x0 = rec["x0"]
    assert h.resume_compatible(rec, h.BUDGET, x0) == (True, "")
    bad = {**h.BUDGET, "P": {"DE_maxiter": 2000, "refine_round_cap": 20}}
    assert h.resume_compatible(rec, bad, x0)[0] is False
    assert (
        h.resume_compatible(rec, h.BUDGET, {**x0, "U": {"start": "x", "FF_pad64": 0.9}})[0] is False
    )
    k1 = h._key("P", 24, "2026-10-21T06:11:00Z_971")
    assert k1 != h._key("P", 24, "2026-10-21T07:00:00Z_5") and ":" not in k1


class _Ready:
    def __init__(self, ok):
        self.ok = ok

    def ready(self):
        return self.ok

    def get(self):
        return {"status": "finished", "seed": 24}


def test_finish_pending_collects_results_that_arrived():
    rows = []
    pending = {
        ("P", 24): (("P", 24, None, "r"), _Ready(True)),
        ("P", 25): (("P", 25, None, "r"), _Ready(False)),
    }
    h.finish_pending(
        pending, rows.append, partial=lambda jb: {"status": "unfinished", "seed": jb[1]}
    )
    assert rows == [{"status": "finished", "seed": 24}, {"status": "unfinished", "seed": 25}]


def test_u_verification_scores_equal_s1f_record():
    """The U branch of verify() reproduces S1f's own verify_u scores for the same solution."""
    d = json.loads((h.OUT / "s1h_decision.json").read_text())["verify"]["U"]
    s1f = json.loads((h.OUT / "s1f_verify.json").read_text())
    assert d["label"] == "S1f seed 23" and s1f["best_label"] == "seed23"
    assert d["scores"] == s1f["scores"]
    assert d["FF_conservative"] == s1f["FF_conservative"]

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

"""Offline tests for scripts/x4_s1k.py (X4 S1k): the graded polish runs all three stages, keeps the
best point, wraps periodic angles, builds legal simplices and stops on its budgets; the decision. No
lalsuite or dynesty needed."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

_S = Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _S / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load("x4_s1k")
LO, HI = m.s.LO, m.s.HI
FAR = "2099-01-01T00:00:00+00:00"


def bowl(c):
    def f(z):
        z = np.asarray(z)
        if np.any(z < LO) or np.any(z > HI):
            return 1.0
        return float(np.sum((z - c) ** 2))

    return f


def test_all_three_stages_run_and_never_worse():
    c = (LO + HI) / 2
    r = m.graded(bowl(c), c, deadline=FAR)  # already optimal: no stage can gain, all still run
    assert [st["scale"] for st in r["stages"]] == [1.0, 0.1, 0.01] and r["f"] == 0.0
    r2 = m.graded(bowl(c), c + 0.03, deadline=FAR)
    assert r2["f"] <= r2["f0"] and r2["stop"] == "completed" and r2["nfev"] <= 10_001


def test_simplex_vertices_stay_legal_at_the_bounds():
    x = HI.copy()
    x[m.PERIODIC] = 2 * np.pi - 1e-3
    for sc in m.SCALES:
        v = np.array([m.wrap(z) for z in m.simplex(x, sc)])
        assert np.all(v >= LO) and np.all(v <= HI)
        assert len({tuple(z) for z in v}) == m.s.NDIM + 1  # no repeated vertex


def test_periodic_wrap_lets_the_search_cross_two_pi():
    c = (LO + HI) / 2
    c[7] = 0.01  # optimum just past the 2 pi seam
    x0 = c.copy()
    x0[7] = 2 * np.pi - 0.05
    r = m.graded(bowl(c), x0, deadline=FAR)
    assert r["f"] < 1e-4


def test_budget_and_time_stops(monkeypatch):
    rng = np.random.default_rng(0)
    r = m.graded(lambda z: float(rng.random()), (LO + HI) / 2, stage_fev=(50, 50, 50), deadline=FAR)
    assert r["nfev"] <= 151
    r = m.graded(lambda z: 1.0, (LO + HI) / 2, case_s=-1, deadline=FAR)
    assert r["stop"] == "time_cap"


def test_decide(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "OUT", tmp_path)
    monkeypatch.setattr(m.s, "verify", lambda x: {"FF_conservative": x[0], "spread": 0.0})

    def write(name, rows):
        (tmp_path / name).write_text(json.dumps({"rows": rows}))

    st = {"stop": "maxcall"}
    write(
        "s1k_validate.json",
        [
            {"name": n, "eps_final": 3e-4, "s1i_stage": st, "stop": "completed"}
            for n in ("C1", "C2")
        ],
    )

    def ref(seed, ff, conv=True):
        return {
            "seed": seed,
            "start_FF_pad64_s1i": 0.91,
            "FF_pad64": ff,
            "stop": "completed",
            "last_stage_converged": conv,
            "x": [ff],
        }

    write("s1k_refine.json", [ref(31, 0.9150), ref(32, 0.9151)])
    assert m.decide()["outcome"].startswith("U cross-start consistency passed")
    write("s1k_refine.json", [ref(31, 0.9146), ref(32, 0.9151)])
    assert "gap" in m.decide()["outcome"]
    write("s1k_refine.json", [ref(31, 0.9150, False), ref(32, 0.9151)])
    assert "budget-limited" in m.decide()["outcome"]
    write(
        "s1k_validate.json",
        [{"name": "C1", "eps_final": 6e-4, "s1i_stage": st, "stop": "completed"}],
    )
    d = m.decide()
    assert (
        d["validation"]["C1"]["status"] == "fail"
        and d["validation"]["C2"]["status"] == "unfinished"
    )

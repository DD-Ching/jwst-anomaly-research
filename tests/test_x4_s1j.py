"""Offline tests for scripts/x4_s1j.py (X4 S1j): the extended polish never returns a worse point and
obeys its budget, and the frozen decision. No lalsuite or dynesty needed."""

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


m = _load("x4_s1j")


@pytest.fixture(autouse=True)
def _future_deadline(monkeypatch):
    # DEADLINE_UTC is the frozen run's wall-clock stop; keep it out of the unit tests.
    monkeypatch.setattr(m, "DEADLINE_UTC", "9999-01-01T00:00:00+00:00")


def test_extend_never_worse_and_converges_on_a_bowl(monkeypatch):
    c = (m.s.LO + m.s.HI) / 2
    monkeypatch.setattr(m.e, "objective", lambda z, t: float(np.sum((np.asarray(z) - c) ** 2)))
    x0 = c + 0.05
    f0 = float(np.sum((x0 - c) ** 2))
    r = m.extend(None, x0, f0)
    assert r["f"] <= f0 and r["stop"] == "converged" and r["extra_fev"] <= m.EXTRA_FEV
    r2 = m.extend(None, c, 0.0)  # already optimal: first restart gains nothing
    assert r2["f"] == 0.0 and r2["stop"] == "converged" and len(r2["trajectory"]) == 1


def test_extend_stops_at_the_evaluation_budget(monkeypatch):
    rng = np.random.default_rng(0)
    monkeypatch.setattr(
        m.e, "objective", lambda z, t: float(rng.random())
    )  # noise: always "improves"
    monkeypatch.setattr(m, "EXTRA_FEV", 400)
    r = m.extend(None, (m.s.LO + m.s.HI) / 2, 1.0)
    assert r["extra_fev"] <= 400 and r["stop"] in ("budget_fev", "converged")


def test_decide(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "OUT", tmp_path)
    monkeypatch.setattr(m.s, "verify", lambda x: {"FF_conservative": x[0], "spread": 0.0})

    def write(name, rows):
        (tmp_path / name).write_text(json.dumps({"rows": rows}))

    st = {"stop": "maxcall"}
    v = [
        {"name": n, "eps_final": 3e-4, "s1i_stage": st, "stop": "converged", "extra_fev": 10}
        for n in ("B1", "B2")
    ]
    write("s1j_validate.json", v)

    def ref(seed, ff, stop="converged"):
        return {
            "seed": seed,
            "start_FF_pad64_s1i": 0.9,
            "FF_pad64": ff,
            "stop": stop,
            "extra_fev": 1,
            "bounds": [],
            "x": [ff],
        }

    write("s1j_refine.json", [ref(31, 0.9150), ref(32, 0.9151)])
    assert m.decide()["outcome"].startswith("U cross-start consistency passed")
    write("s1j_refine.json", [ref(31, 0.9146), ref(32, 0.9151)])
    assert "gap" in m.decide()["outcome"]
    write("s1j_refine.json", [ref(31, 0.9150, "budget_fev"), ref(32, 0.9151)])
    assert "did not converge" in m.decide()["outcome"]
    v[1]["eps_final"] = 6e-4
    write("s1j_validate.json", v)
    d = m.decide()
    assert d["validation"]["B2"]["status"] == "fail" and not d["validation_pass"]

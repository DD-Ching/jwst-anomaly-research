"""Offline tests for scripts/w2_j0740_gaia_delay.py (frozen rules and stored result)."""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "w2j0740", ROOT / "scripts" / "w2_j0740_gaia_delay.py"
)
w = importlib.util.module_from_spec(spec)
sys.modules["w2j0740"] = w
spec.loader.exec_module(w)
STORED = ROOT / "results" / "w2_j0740_gaia_delay" / "result.json"


def test_groups_and_error_scale():
    t = np.array([0.0, 0.07, 0.3, 5.0, 5.2, 30.0])
    assert [list(g) for g in w.groups(t)] == [[0, 1, 2], [3, 4], [5]]
    c = {"t": t, "m": np.array([1.0, 1.2, 0.8, 2.0, 2.2, 3.0]), "e": np.full(6, 0.1)}
    q = w.error_scale([c])
    # group 1: mean 1.0, chi2 = 4 + 4 = 8 (2 dof); group 2: mean 2.1, chi2 = 1 + 1 = 2 (1 dof)
    assert q["dof"] == 3 and q["chi2_nu"] == pytest.approx(10 / 3)
    assert q["s"] == pytest.approx(np.sqrt(10 / 3))
    c["m"] = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 3.0])
    assert w.error_scale([c])["s"] == 1.0  # never shrinks errors


def test_intervals():
    g = np.arange(-10.0, 11.0, 5.0)
    assert w.intervals(np.array([1, 1, 0, 1, 0], bool), g) == [[-10.0, -5.0], [5.0, 5.0]]
    assert w.intervals(np.array([0, 0, 0, 1, 1], bool), g) == [[5.0, 10.0]]
    assert w.intervals(np.zeros(5, bool), g) == []


def test_sign_convention_recovers_injected_lag():
    """B repeating A 60 d later must come out as dt = +60 (PyCS3 spline, dense noiseless data)."""
    w._quiet()
    t = np.arange(0.0, 1000.0, 7.0)

    def src(x):
        return 0.15 * np.sin(2 * np.pi * x / 420.0) + 0.05 * np.sin(2 * np.pi * x / 170.0)

    e = np.full(len(t), 0.003)
    a = {"t": t, "m": 18.0 + src(t), "e": e}
    b = {"t": t, "m": 19.3 + src(t - 60.0), "e": e}
    assert w._dt(w._lcs(a, b, 75.0)) == pytest.approx(75.0)
    assert w.estimate(a, b)["dt_hat"] == pytest.approx(60.0, abs=3.0)


@pytest.mark.skipif(not STORED.exists(), reason="no stored result")
def test_stored_verdict_follows_frozen_rules():
    r = json.loads(STORED.read_text())
    t95, t99 = r["thresholds"]["T95"], r["thresholds"]["T99"]
    assert t95 <= t99
    assert r["C2"]["pass"] == (r["C2"]["fpr"] <= 0.05)
    assert r["controls_pass"]["C1_any"] == any(v["fraction"] >= 0.68 for v in r["C1"].values())
    tg = r["target"]
    excl = np.array(tg["scan_D"]) > t99
    assert w.intervals(excl, np.array(tg["scan_grid"])) == tg["excluded_intervals"]
    supports = (
        all(r["controls_pass"].values())
        and tg["D_hat"] is not None
        and tg["D_hat"] <= t95
        and w.ZERO_LAG <= abs(tg["dt_hat"]) <= w.DT_MAX
        and tg["regdiff"]["dt_hat"] is not None
        and abs(tg["regdiff"]["dt_hat"] - tg["dt_hat"]) <= w.TOL
    )
    if supports:
        assert r["verdict"].startswith("supports")
    elif excl.any():
        assert r["verdict"].startswith("excludes")
    else:
        assert r["verdict"] == "data insufficient"

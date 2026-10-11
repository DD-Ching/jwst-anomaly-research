"""Offline tests for scripts/goal_nist_nosignal.py (goal round 10); synthetic counts only."""

import importlib.util
import math
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("goal_nist_nosignal", _DIR / "goal_nist_nosignal.py")
ns = importlib.util.module_from_spec(_spec)
sys.modules["goal_nist_nosignal"] = ns
_spec.loader.exec_module(ns)


def _table(cells):
    return {
        p: dict(zip(("++", "+0", "0+", "00"), c, strict=True)) | {"N": sum(c)}
        for p, c in zip(ns.PAIRS, cells, strict=True)
    }


def test_comparison_mapping():
    # Alice's '+' counts depend on Bob's setting only in the (a, b) vs (a, b') comparison
    n = _table([(10, 90, 0, 900), (20, 80, 0, 900), (5, 5, 0, 990), (5, 5, 0, 990)])
    c = ns.comparisons(n)
    assert c["A|a"] == (100, 1000, 100, 1000)  # ++ and +0 both count for Alice
    assert c["B|b"] == (10, 1000, 5, 1000)  # ++ and 0+ count for Bob; Alice's setting varies


def test_g_test_matches_pooled_z_for_large_counts():
    k1, n1, k2, n2 = 14311, 44449651, 13787, 44454006
    g, p = ns.g_test(k1, n1, k2, n2)
    z = ns.pooled_z(k1, n1, k2, n2)
    assert abs(g - z**2) / z**2 < 1e-3  # likelihood ratio and score tests agree asymptotically
    assert 0 < p < 0.01


def test_mutual_information_small_shift():
    # I ~= delta^2 / (8 p (1 - p) ln 2) for a uniform binary input and a small shift
    p, d = 3e-4, 6e-6
    approx = d**2 / (8 * p * (1 - p) * math.log(2))
    assert abs(ns.bern_mi(p, p + d) / approx - 1) < 0.05
    assert ns.bern_mi(p, p) == 0.0

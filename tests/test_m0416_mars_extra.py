import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m0416_mars_extra as me  # noqa: E402


def test_classify():
    assert me.classify(0.4, True) == "both"
    assert me.classify(3.0, True) == "both"
    assert me.classify(3.1, True) == "CANUCS only"
    assert me.classify(np.inf, True) == "CANUCS only"
    assert me.classify(np.nan, True) == "CANUCS only"
    assert me.classify(0.1, False) == "not testable"


def test_match_observed_one_to_one_greedy():
    # predictions at x = 0, 1, 10; observed at x = 0.2 and 0.6: the closest pair (0.2 -> 0) goes
    # first, then 0.6 takes 1 (0.4 away); 10 is beyond 3" and stays free
    used = me.match_observed([0.0, 1.0, 10.0], [0, 0, 0], [0.2, 0.6], [0, 0])
    assert used.tolist() == [True, True, False]


def test_match_observed_radius_and_empty():
    assert me.match_observed([5.0], [0.0], [0.0], [0.0]).tolist() == [False]
    assert me.match_observed([], [], [0.0], [0.0]).tolist() == []
    assert me.match_observed([1.0], [0.0], [], []).tolist() == [False]

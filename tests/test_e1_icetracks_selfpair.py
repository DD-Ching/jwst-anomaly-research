import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import e1_icetracks_selfpair as sp  # noqa: E402


def _z(a, b, c):
    return {"100s-1h": a, "1h-1d": b, "1d-7d": c}


def test_outcome_explained():
    assert sp.outcome(_z(0.5, -1.9, 1.0), -6.0) == "explained"


def test_outcome_not_explained_either_long_cell():
    assert sp.outcome(_z(0.0, -3.0, 0.0), -6.0) == "not explained"
    assert sp.outcome(_z(0.0, 0.0, -3.5), -6.0) == "not explained"


def test_outcome_undetermined():
    # control does not keep the deficit
    assert sp.outcome(_z(0.0, 0.0, 0.0), -2.9) == "undetermined"
    # a cell between the thresholds
    assert sp.outcome(_z(-2.5, 0.0, 0.0), -6.0) == "undetermined"
    assert sp.outcome(_z(0.0, 2.0, 0.0), -6.0) == "undetermined"


def test_shared_mask():
    day = 86400.0
    ic = np.array([100.0, 200.0, 300.0])
    d = np.array([100.0, 100.0 + 9 / day, 200.0 + 11 / day, 250.0, 300.0 - 5 / day, 50.0, 400.0])
    assert sp.shared_mask(d, ic).tolist() == [True, True, False, False, True, False, False]

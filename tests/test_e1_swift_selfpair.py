import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import e1_swift_selfpair as sp  # noqa: E402


def _z(a, b, c):
    return {"100s-1h": a, "1h-1d": b, "1d-7d": c}


def test_outcome_explained():
    assert sp.outcome(_z(0.5, -1.9, 1.0), -4.0) == "explained"


def test_outcome_not_explained_either_long_cell():
    assert sp.outcome(_z(0.0, -3.0, 0.0), -6.0) == "not explained"
    assert sp.outcome(_z(0.0, 0.0, -3.5), -6.0) == "not explained"


def test_outcome_undetermined():
    # control (1 h-1 d) does not keep the deficit
    assert sp.outcome(_z(0.0, 0.0, 0.0), -2.9) == "undetermined"
    # a cell between the thresholds
    assert sp.outcome(_z(-2.5, 0.0, 0.0), -6.0) == "undetermined"
    assert sp.outcome(_z(0.0, 2.0, 0.0), -6.0) == "undetermined"


def test_shared_window_is_control_window():
    day = 86400.0
    g = np.array([100.0])
    b = np.array([100.0 + 99 / day, 100.0 - 101 / day])
    assert sp.sp.shared_mask(b, g, sp.SHARED_S).tolist() == [True, False]

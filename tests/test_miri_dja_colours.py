import importlib.util
from pathlib import Path

import numpy as np

_spec = importlib.util.spec_from_file_location(
    "miri_dja_colours", Path(__file__).resolve().parents[1] / "scripts" / "miri_dja_colours.py"
)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)


def test_rule_v():
    assert m.rule_v(np.array([0.0, 0.05, -0.05, 0.1, 0.02]))
    assert not m.rule_v(np.array([0.0, 0.5, -0.5, 1.0, -1.0]))
    assert not m.rule_v(np.array([0.0, 0.01]))  # too few rows
    assert not m.rule_v(np.array([np.nan] * 6))


def test_rule_c():
    assert m.rule_c(80, 100)
    assert not m.rule_c(79, 100)
    assert not m.rule_c(0, 0)


def test_rule_b():
    assert m.rule_b(0.05, 0.3) == "no_bias"
    assert m.rule_b(-0.30, 0.19) == "pass"
    assert m.rule_b(0.30, -0.25) == "fail"
    assert m.rule_b(np.nan, 0.0) == "fail"
    assert m.rule_b(0.02, np.nan) == "no_bias"
    assert m.rule_b(0.30, np.nan) == "fail"


def test_decide_and_jaccard():
    assert m.decide(False, True, "pass") == "stopped_v_failed"
    assert m.decide(True, True, "pass") == "adopt"
    assert m.decide(True, False, "pass") == "not_adopted"
    assert m.decide(True, True, "no_bias") == "not_adopted"
    assert m.jaccard(["a", "b"], ["b", "c"]) == 1 / 3

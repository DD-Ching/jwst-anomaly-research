"""Frozen rules of scripts/m0417_extra.py (flag and gate C)."""

import importlib.util
from pathlib import Path

_P = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("m0417_extra", _P / "m0417_extra.py")


def _mod():
    import sys

    sys.path.insert(0, str(_P))
    m = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(m)
    return m


def test_flag_only_absent_unmatched_magnified():
    m = _mod()
    assert m.flag("missing", "absent")
    assert m.flag("candidate", "absent")
    assert not m.flag("observed", "absent")
    assert not m.flag("demagnified", "absent")
    assert not m.flag("outside", "absent")
    assert not m.flag("missing", "ambiguous")
    assert not m.flag("missing", "undetectable")


def test_gate_c_counts_testable_only():
    m = _mod()
    g = m.gate_c(["recovered"] * 4 + ["absent", "undetectable", "no_reference"])
    assert g["n_testable"] == 5 and g["n_recovered"] == 4
    assert g["passes"]
    g = m.gate_c(["recovered"] * 3 + ["absent", "ambiguous"])
    assert not g["passes"]
    assert not m.gate_c(["undetectable"])["passes"]

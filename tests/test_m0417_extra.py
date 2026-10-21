"""Frozen rules of scripts/m0417_extra.py (flag and gate C)."""

import importlib.util
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_SCRIPTS))
_spec = importlib.util.spec_from_file_location("m0417_extra", _SCRIPTS / "m0417_extra.py")
mx = importlib.util.module_from_spec(_spec)
sys.modules["m0417_extra"] = mx
_spec.loader.exec_module(mx)


def test_flag_only_absent_unmatched_magnified():
    assert mx.flag("missing", "absent")
    assert mx.flag("candidate", "absent")
    assert not mx.flag("observed", "absent")
    assert not mx.flag("demagnified", "absent")
    assert not mx.flag("outside", "absent")
    assert not mx.flag("missing", "ambiguous")
    assert not mx.flag("missing", "undetectable")


def test_gate_c_counts_testable_only():
    g = mx.gate_c(["recovered"] * 4 + ["absent", "undetectable", "no_reference"])
    assert g["n_testable"] == 5 and g["n_recovered"] == 4
    assert g["passes"]
    assert not mx.gate_c(["recovered"] * 3 + ["absent", "ambiguous"])["passes"]
    g = mx.gate_c(["undetectable"])
    assert not g["passes"] and g["fraction"] is None

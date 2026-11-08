import importlib.util
from pathlib import Path

import pytest

_P = Path(__file__).resolve().parents[1] / "scripts" / "arcs_shear_error.py"
_spec = importlib.util.spec_from_file_location("arcs_shear_error", _P)
se = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(se)


def test_map_gate():
    assert se.map_gate(0.05) and se.map_gate(1e-5)
    assert not se.map_gate(0.051) and not se.map_gate(float("nan"))


def test_control_gate_and_retention():
    assert se.retention(["aligned", "aligned", "mixed"]) == pytest.approx(2 / 3)
    assert se.control_gate(0.70) and not se.control_gate(0.69)
    assert not se.control_gate(se.retention([]))


def test_models_label():
    assert se.models_label(["anti", "anti"]) == "robust"
    assert se.models_label(["anti", "mixed"]) == "model_sensitive"
    assert se.models_label(["anti", "anti", "mixed"]) == "partly"
    assert se.models_label([]) == "undetermined"


def test_range_label():
    assert se.range_label(0.68, True) == "robust"
    assert se.range_label(0.32, True) == "model_sensitive"
    assert se.range_label(0.5, True) == "partly"
    assert se.range_label(0.9, False) == "undetermined"


def test_outcome():
    assert se.outcome({"a": "robust", "b": "robust"}).startswith("model error does not")
    assert se.outcome({"a": "robust", "b": "model_sensitive"}) == "model error explains b"
    assert se.outcome({"a": "undetermined", "b": "robust"}) == "undetermined"
    assert se.outcome({"a": "partly", "b": "robust"}) == "mixed"

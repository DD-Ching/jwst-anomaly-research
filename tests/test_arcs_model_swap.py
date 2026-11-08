import importlib.util
from pathlib import Path

from astropy.table import Table

_p = Path(__file__).resolve().parents[1] / "scripts" / "arcs_model_swap.py"
_s = importlib.util.spec_from_file_location("arcs_model_swap", _p)
ms = importlib.util.module_from_spec(_s)
_s.loader.exec_module(ms)


def _second():
    return Table(
        {
            "label": [1, 2, 3],
            "ra": [10.0, 10.0, 10.0],
            "dec": [0.0, 0.0, 0.0],
            "strong_shear": [True, False, True],
            "population": ["background", "background", "member_or_foreground"],
            "orientation_class": ["anti", "aligned", "aligned"],
        }
    )


def test_second_class():
    s = _second()
    row = {"label": 1, "ra": 10.0, "dec": 0.0}
    assert ms.second_class(row, s) == "anti"
    assert ms.second_class({"label": 2, "ra": 10.0, "dec": 0.0}, s) == "not_strong_background"
    assert ms.second_class({"label": 3, "ra": 10.0, "dec": 0.0}, s) == "not_strong_background"
    assert ms.second_class({"label": 9, "ra": 10.0, "dec": 0.0}, s) == "unmatched"
    assert ms.second_class({"label": 1, "ra": 10.0, "dec": 0.2 / 3600}, s) == "unmatched"


def test_control_gate():
    g = ms.control_gate(
        ["aligned"] * 7 + ["anti", "intermediate", "unmatched", "not_strong_background"]
    )
    assert g["n_matched"] == 10 and g["aligned"] == 7 and g["passes"]
    assert not ms.control_gate(["aligned"] * 6 + ["anti"] * 4)["passes"]
    assert not ms.control_gate(["unmatched"])["passes"]


def test_outcome():
    assert ms.outcome(True, 6, 6) == "persists"
    assert ms.outcome(True, 6, 5) == "persists"
    assert ms.outcome(True, 6, 4) == "mixed"
    assert ms.outcome(True, 6, 3) == "model-dependent"
    assert ms.outcome(False, 6, 6) == "undetermined"
    assert ms.outcome(True, 3, 3) == "undetermined"

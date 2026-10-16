import importlib.util
from pathlib import Path

from astropy.table import Table

_P = Path(__file__).resolve().parents[1] / "scripts" / "arcs_macs1423_census.py"
_spec = importlib.util.spec_from_file_location("arcs_macs1423_census", _P)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)


def test_anti_sets_use_screen_photoz_class():
    rows = Table(
        {
            "label": [1758, 1818, 2197, 10],
            "orientation_class": ["aligned", "anti", "anti", "aligned"],
            "class_bgsub": ["anti", "mixed", "anti", "aligned"],
        }
    )
    s = m.anti_sets(rows)
    assert s["screen_anti"] == ["1818", "2197"]
    assert s["screen_anti_staying_bgsub"] == ["2197"]  # 1758 (anti at z = 2 only) not counted


def test_final_outcome_gate_fraction():
    o = "excess remains under the clump-aware null"
    assert m.final_outcome(o, 28, 55) == o
    assert m.final_outcome(o, 27, 55).startswith("undetermined")
    assert m.final_outcome(o, 0, 0).startswith("undetermined")

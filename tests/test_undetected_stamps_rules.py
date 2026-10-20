"""Frozen rules of results/dimming/undetected_stamps/PLAN.md (stamps.reason, score.decide)."""

import importlib.util
from pathlib import Path

import numpy as np

_R = Path(__file__).resolve().parents[1] / "results/dimming/undetected_stamps"


def _load(name):
    spec = importlib.util.spec_from_file_location(f"undetected_{name}", _R / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


st = _load("stamps")
sc = _load("score")


def _sci(n=31):
    return np.ones((n, n))


def test_reason_edge_uses_measure_pixel_rule():
    assert st.reason(_sci(), 2.5, 500.0, 2048, 2048, 3.0, 0, 485) == "edge"
    assert st.reason(_sci(), 1000.0, 2045.4, 2048, 2048, 3.0, 985, 2030) == "edge"


def test_reason_aperture_no_data_and_annulus():
    s = _sci()
    s[15, 16] = np.nan  # one pixel inside r = 3 of the centre (15, 15) of a stamp at (985, 985)
    assert st.reason(s, 1000.0, 1000.0, 2048, 2048, 3.0, 985, 985) == "aperture_no_data"
    s = np.full((31, 31), np.nan)
    s[12:19, 12:19] = 1.0  # aperture finite, annulus empty
    assert st.reason(s, 1000.0, 1000.0, 2048, 2048, 3.0, 985, 985) == "annulus"
    assert st.reason(_sci(), 1000.0, 1000.0, 2048, 2048, 3.0, 985, 985) == "not_reproduced"


def test_decide():
    assert sc.decide(0, 0) == "not tested"
    assert sc.decide(40, 41) == "R1"  # Wilson lower 0.87
    assert sc.decide(33, 60) == "R2"

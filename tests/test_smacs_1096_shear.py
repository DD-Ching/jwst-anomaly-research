import importlib.util
from pathlib import Path

import numpy as np

_spec = importlib.util.spec_from_file_location(
    "smacs_1096_shear", Path(__file__).resolve().parents[1] / "scripts" / "smacs_1096_shear.py"
)
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)


def test_decide_frozen_rule():
    assert mod.decide(0.04, 0.001) == "no_power"
    assert mod.decide(np.nan, 0.5) == "no_power"
    assert mod.decide(0.13, 0.009) == "inconsistent"
    assert mod.decide(0.13, 0.01) == "consistent"
    assert mod.decide(0.05, 0.47) == "consistent"


def test_offset_is_axial():
    assert mod.offset_deg(10.0, 170.0) == 20.0
    assert abs(mod.offset_deg(97.38, 122.77) - 25.39) < 1e-9

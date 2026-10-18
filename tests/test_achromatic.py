from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import achromatic_test as at  # noqa: E402


def test_achromatic_chi2_constant_and_outlier():
    c, dof, p, mean = at.achromatic_chi2([0.9, 0.9, 0.9], [0.05, 0.05, 0.05])
    assert c == 0 and dof == 2 and p == 1.0 and abs(mean - 0.9) < 1e-12
    c, dof, p, _ = at.achromatic_chi2([0.5, 1.0, 1.0], [0.05, 0.05, 0.05])
    assert p < 0.01


def test_rule_needs_controls_and_rejects_chromatic():
    base = {"n_controls": 25, "control_median": 1.0, "control_scatter": 0.02}
    cal = {
        "F090W": {**base, "target": [100.0, 1.0, 50.0, 1.0]},
        "F150W": {**base, "target": [100.0, 1.0, 100.0, 1.0]},
        "F200W": {**base, "target": [100.0, 1.0, 100.0, 1.0]},
    }
    assert at.rule(cal)["p"] < at.P_REJECT
    cal["F090W"]["target"] = [100.0, 1.0, 100.0, 1.0]
    assert at.rule(cal)["p"] > at.P_REJECT
    cal["F090W"]["n_controls"] = 19
    assert at.rule(cal)["decision"].startswith("not decided")
    assert np.isclose(at.ratio_err(2.0, 0.0, 1.0, 0.1)[1], 0.05)

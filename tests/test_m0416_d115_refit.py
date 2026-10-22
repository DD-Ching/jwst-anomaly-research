import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import m0416_d115_refit as m  # noqa: E402


def test_family_chi2_one_to_one_and_shared():
    ox, oy, sig = np.array([0.0, 1.0]), np.array([0.0, 0.0]), np.array([0.5, 0.5])
    c, dt, shared = m.family_chi2(ox, oy, sig, Table({"x": [1.1, 0.0, 5.0], "y": [0.0] * 3}))
    assert not shared and np.allclose(dt, [0.0, 0.1]) and np.isclose(c, 0.04)
    c, dt, shared = m.family_chi2(ox, oy, sig, Table({"x": [0.5], "y": [0.0]}))
    assert shared and np.allclose(dt, [0.5, 0.5]) and np.isclose(c, 2.0)
    c, _, shared = m.family_chi2(ox, oy, sig, Table({"x": [], "y": []}))
    assert shared and c == np.inf


def test_verdict():
    assert m.verdict(5.9, False, True) == "ordinary"
    assert m.verdict(6.1, False, True) == "inconclusive"
    assert m.verdict(1.0, True, True) == "inconclusive"
    assert m.verdict(1.0, False, False) == "inconclusive"

"""Frozen rules of ``scripts/m0417_f10_colour.py`` (MACS0417 family-10 colour test)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m0417_f10_colour as fc  # noqa: E402


def test_pair_chi2_scaled_copy_is_zero():
    fb = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    eb = np.full(6, 0.01)
    chi2, dof, r = fc.pair_chi2(2.5 * fb, eb, fb, eb)
    assert chi2 < 1e-8 and dof == 5 and abs(r - 2.5) < 1e-4


def test_pair_chi2_drops_bad_bands():
    f = np.array([1.0, 2.0, np.nan, 4.0, 5.0, 6.0])
    e = np.array([0.1, 0.1, 0.1, 0.0, 0.1, np.nan])
    _, dof, _ = fc.pair_chi2(f, e, f, np.full(6, 0.1))
    assert dof == 2  # bands 0, 1, 4 remain


def test_classify_and_gate_and_decide():
    assert fc.classify(1.0, 3) == "untestable"  # 4 bands < 5
    assert fc.classify(3.0, 7) == "consistent"
    assert fc.classify(40.0, 7) == "differ"
    g = fc.gate(["consistent"] * 8 + ["differ"] * 2 + ["untestable"])
    assert g["n_testable"] == 10 and g["passes"]
    assert not fc.gate(["consistent"] * 7)["passes"]  # too few pairs
    assert not fc.gate(["consistent"] * 7 + ["differ"] * 3)["passes"]  # 70 %
    assert fc.decide(False, ["differ"]) == "not validated"
    assert fc.decide(True, ["consistent", "differ"]).startswith("H2")
    assert fc.decide(True, ["consistent", "consistent"]) == "H1 not excluded"
    assert fc.decide(True, ["consistent", "untestable"]) == "untestable"
    assert fc.decide(True, []) == "untestable"

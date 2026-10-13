"""Fast offline tests for scripts/x2_s5_radial.py (X2 S5, radial perturbations)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import x2_s5_radial as s5  # noqa: E402


def test_ellis_bronnikov_growth_matches_cps_master_equation():
    CF = s5.load_coefficients()
    S = s5.setup("EB", {"b": 1.0}, 15.0, 0.1, CF)
    top = s5.dense_top(S, k=1)[0]
    # -d^2/dx^2 - 3 / (x^2 + 1)^2 on the same grid
    from scipy.linalg import eigh_tridiagonal

    x = S["y"][1:-1]
    e0 = eigh_tridiagonal(
        2 / 0.1**2 - 3 / (x**2 + 1) ** 2,
        -np.ones(x.size - 1) / 0.1**2,
        select="i",
        select_range=(0, 0),
    )[0][0]
    assert abs(top["lam2_re"] + e0) < 1e-3 * abs(e0)
    assert abs(math.sqrt(top["lam2_re"]) - 1.182) < 0.01


def test_c_from_a_satisfies_the_constraint():
    CF = s5.load_coefficients()
    S = s5.setup("SV", {"M": 1.0, "a": 2.01}, 60.0, 0.05, CF)
    A = np.exp(-(S["y"] ** 2) / 50.0)
    C = s5.c_from_a(A, S)
    sig = s5.dx1(C, S) - S["k1"] * C - S["k2"] * A
    m = slice(5, -5)
    assert np.max(np.abs(sig[m])) < 1e-3 * np.max(np.abs(S["k2"][m] * A[m]))
    assert abs(C[-1]) < 1e-6 * np.max(np.abs(C))  # even A: C -> 0 at both ends


def test_sv_background_matter_functions():
    vals, f = s5.background("SV", {"M": 1.0, "a": 2.3}, np.array([0.0, 1.0, 4.0]))
    x = np.array([0.0, 1.0, 4.0])
    R = np.sqrt(x**2 + 2.3**2)
    assert np.allclose(vals[13] ** 2, 2 * 2.3**2 / R**4)  # Phi'^2 (phantom, kappa = 1)
    assert np.allclose(vals[16], 4 * 2.3**2 / (5 * R**5))  # V
    assert np.allclose(vals[20], 6 * 2.3**2 / (5 * R**5))  # U (magnetic NED energy density)
    assert np.allclose(f, 1 - 2 / R)

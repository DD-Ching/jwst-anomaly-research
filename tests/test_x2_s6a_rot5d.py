"""Fast offline tests for scripts/x2_s6a_rot5d.py (X2 S6a, 5D rotating wormhole modes)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import x2_s6a_rot5d as s6a  # noqa: E402


def test_complex_eigenpair_residual_uses_the_full_complex_pair():
    # known complex eigenpairs: a 2D rotation by 90 degrees has lambda = +-i, v = (1, -+i)/sqrt 2
    K = np.array([[0.0, -1.0], [1.0, 0.0]])
    ev, vr = eig(K)
    i = int(np.argmax(ev.imag))
    assert abs(ev[i] - 1j) < 1e-12
    assert s6a.eigpair_residual(K, ev[i], vr[:, i]) < 1e-14
    # the S5 shortcut (real parts only) must fail on the same pair
    assert s6a.eigpair_residual(K, ev[i].real, vr[:, i].real) > 0.1


def test_complex_eigenvalues_of_a_coupled_ode_are_recovered():
    # -u'' + C u = lam u on [0, pi], Dirichlet; C non-symmetric with eigenvalues 2 +- 3i
    # => lam = k^2 + 2 +- 3i for k = 1, 2, ...
    n = 150
    x = np.linspace(0, np.pi, n + 2)[1:-1]
    dx = x[1] - x[0]
    off = np.ones(n - 1)
    D = (np.diag(-2 * np.ones(n)) + np.diag(off, 1) + np.diag(off, -1)) / dx**2
    C = np.array([[2.0, -3.0], [3.0, 2.0]])
    K = np.kron(np.eye(2), -D) + np.kron(C, np.eye(n))
    ev, vr = eig(K)
    i = int(np.argmin(np.abs(ev - (3 + 3j))))
    assert abs(ev[i] - (3 + 3j)) < 1e-3  # O(dx^2) in the real part
    assert s6a.eigpair_residual(K, ev[i], vr[:, i]) < 1e-10


def test_background_gates():
    sol = s6a.solve_background(0.0)
    assert np.max(np.abs(sol.y[[0, 2]])) < 1e-12  # static symmetric: a = q = 0
    sol = s6a.solve_background(0.5, guess=sol)
    P = s6a.background_props(sol, 0.5)
    assert abs(P["smarr_rel_residual"]) < 1e-5
    assert P["Q2_rel_variation"] < 1e-6
    # paper Fig. 1, digitised: M/A^(2/3) = 0.586 at J/A = 0.142
    assert abs(P["J_over_A"] - 0.1420) < 1e-3
    assert abs(P["M_over_A23"] - 0.586) < 0.005


def test_static_mode_and_rotating_pair():
    V, _ = s6a.load_equations()
    sol = s6a.solve_background(0.0)
    area = s6a.background_props(sol, 0.0)["A"]
    ev, _, _, _ = s6a.spectrum(V, sol, 0.0, 5.0, 200)
    assert abs(np.min(ev.real) * area ** (2 / 3) + 21.80) < 0.05  # paper Fig. 5: -21.80
    sol = s6a.background_at(0.5)
    area = s6a.background_props(sol, 0.5)["A"]
    ev, _, _, _ = s6a.spectrum(V, sol, 0.5, 5.0, 300)
    z = ev * area ** (2 / 3)
    i = int(np.argmin(np.abs(z - (-2.319 + 7.462j))))
    assert abs(z[i] - (-2.319 + 7.462j)) < 0.02
    # independent Evans shooting finds the same complex root
    r, info = s6a.Shooter(V, sol, 0.5, S=3.5).root(z[i] / area ** (2 / 3))
    assert info["min_step"] < 1e-6 and info["last_step"] < 1e-6
    assert abs(r * area ** (2 / 3) - (-2.3187 + 7.4620j)) < 1e-3

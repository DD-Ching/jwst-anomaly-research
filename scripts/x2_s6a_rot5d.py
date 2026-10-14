"""X2 S6a: follow the two unstable radial modes of 5D equal-J rotating Ellis wormholes.

Model: Dzhunushaliev et al., arXiv:1309.2448 (PRD 88, 124028), symmetric wormholes (a = q = 0
at both ends), throat parameter r0 = 1. Reported numbers are scale-free: Omega^2 A^(2/3), J/A and
M/A^(2/3). Units for J and M: 4 pi G = 1, fixed by the extremal Myers-Perry point of the paper's
Fig. 1 (J/A = 0.5, M/A^(2/3) = 3.218; see docs/state/units/X2_s6a_5d.md).

Equations: the G0, G1, G2 system *derived* by scripts/x2_s6a_derive5d.py
(results/x2/s6a_equations.json), not the printed V_ij (the printed V10 is a typo). Time dependence
exp(i Omega t): a mode grows iff Im Omega < 0, i.e. iff Omega^2 is not on [0, inf).

Methods:
  * dense: finite differences in s (l = sinh s), even or odd parity at the throat, Dirichlet at
    s = S; all complex eigenpairs kept. A finite box turns the continuum into many spurious complex
    pairs (the mass matrix B is a Jordan block); a mode is accepted only through the checks in
    CRITERIA (box, grid, constraints, shooting). Modes are followed in J by eigenvector overlap,
    not by eigenvalue order.
  * independent residual: the corotating-frame R_tl and the R_ll equations (not used to build the
    G system) evaluated on the complex eigenvector.
  * independent cross-check: Evans-function shooting (decaying data at l = sinh S, integrated in to
    the throat; zero of the parity determinant found by secant iteration in complex Omega^2).

  python scripts/x2_s6a_rot5d.py   (about 30 minutes, 4 processes) -> results/x2/s6a_rot5d.json
"""

from __future__ import annotations

import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_bvp, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import eig

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x2"
EQ_FILE = OUT / "s6a_equations.json"
G_NEWTON = 1 / (4 * np.pi)
BMAT = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
BINV = np.array([[1, -1, 0], [0, 1, 0], [0, 0, 1]])

# Digitised from the gnuplot EPS of 1309.2448 (vector paths, axis ticks as calibration). Data, not
# measurements: about 0.2 % of the axis range. Fig. 5, a_-inf = 0: (|J|/A, Omega^2 A^(2/3)).
PAPER_FIG5 = {
    "static": -21.80,
    "merger_JA": 0.10063,  # quoted in the text
    "branch_static": [
        (0.04957, -18.563),
        (0.07454, -14.487),
        (0.09054, -10.727),
        (0.09704, -8.539),
    ],
    "branch_zero": [(0.04957, -0.360), (0.07454, -1.212), (0.09054, -2.521), (0.09704, -3.663)],
}
# Fig. 1, a_-inf = 0: digitised curve read at our c_omega = 0.2, 0.5, 1, 2, 4; (J/A, M/A^(2/3)).
PAPER_FIG1 = [
    (0.0673, 0.1435),
    (0.1420, 0.5862),
    (0.2144, 1.1776),
    (0.2853, 1.7935),
    (0.3462, 2.2888),
]


def p_of(ell):
    return (ell**2 + 0.5) / (ell**2 + 1)


# ---------------------------------------------------------------- background (Sec. II-IV)
def bg_rhs(ell, y, cw):
    a, ap, q, qp = y
    h = ell**2 + 1
    p = p_of(ell)
    k = (ell / h**2) / p + 3 * ell / h  # xi'/xi, xi = p h^(3/2)
    app = -k * ap + cw**2 * np.exp(4 * a - q) / (2 * p**3 * h**4)
    qpp = -k * qp + 4 * (np.exp(2 * (q - a)) - 1) / h
    return np.vstack([ap, app, qp, qpp])


def solve_background(cw, L=400.0, guess=None, tol=1e-10):
    """Symmetric solution on l in [0, L]: a'(0) = q'(0) = 0, a ~ l^-2, q - a ~ l^-4 at L."""
    grid = np.concatenate([np.linspace(0, 5, 200), np.geomspace(5.05, L, 300)])
    y0 = np.zeros((4, grid.size)) if guess is None else guess.sol(grid)

    def bc(ya, yb):
        return np.array(
            [ya[1], ya[3], L * yb[1] + 2 * yb[0], L * (yb[3] - yb[1]) + 4 * (yb[2] - yb[0])]
        )

    sol = solve_bvp(lambda x, y: bg_rhs(x, y, cw), bc, grid, y0, tol=tol, max_nodes=300000)
    if sol.status != 0:
        raise RuntimeError(f"background did not converge at c_omega = {cw}: {sol.message}")
    return sol


def background_at(cw, prev=None):
    """Continuation in c_omega from a nearby solution (or from the static one)."""
    if prev is not None:
        return solve_background(cw, guess=prev)
    sol = solve_background(0.0)
    for c in np.geomspace(0.05, cw, 12) if cw > 0.05 else [cw]:
        sol = solve_background(c, guess=sol)
    return sol


def background_props(sol, cw):
    L = sol.x[-1]
    a0, q0 = sol.y[0, 0], sol.y[2, 0]
    area = 2 * np.pi**2 * np.sqrt(1 / 8) * np.exp(-(2 * a0 + q0) / 2)  # Eq. (surface) at l = 0
    J = -np.pi * cw / (16 * G_NEWTON)

    def integrand(ell):
        a, _, q, _ = sol.sol(ell)
        return np.exp(4 * a - q) / ((ell**2 + 1) ** 2.5 * p_of(ell) ** 2)

    om0 = -cw * (quad(integrand, 0, L, limit=400, epsabs=1e-14, epsrel=1e-12)[0] + 1 / (4 * L**4))
    M = 3 * np.pi * p_of(L) * (L**2 + 1) ** 1.5 * sol.sol(L)[1] / (8 * G_NEWTON)
    smarr = (2 / 3 * M - 2 * om0 * J) / max(abs(2 / 3 * M), 1e-300)
    ell = np.linspace(0, 30, 3001)
    q2 = q_squared(sol, cw, ell)
    return {
        "c_omega": cw,
        "A": area,
        "J_over_A": abs(J) / area,
        "M_over_A23": M / area ** (2 / 3),
        "omega0": om0,
        "smarr_rel_residual": smarr,
        "Q2": float(np.mean(q2)),
        "Q2_rel_variation": float((q2.max() - q2.min()) / abs(np.mean(q2))),
    }


def q_squared(sol, cw, ell):
    """8 pi G Q^2 = -h^3 p^2 R_ll (Eqs. ric22, Rll-expr); must be constant in l."""
    y = sol.sol(ell)
    a, ap, q, qp = y
    d = bg_rhs(ell, y, cw)
    qpp = d[3]  # R_ll contains q'' but not a''
    h = ell**2 + 1
    p = p_of(ell)
    hp, hpp, pp, ppp = 2 * ell, 2.0, ell / h**2, (1 - 3 * ell**2) / h**3
    wp = cw / (h**2.5 * p**2 * np.exp(q - 4 * a))
    rll = (
        -3 * ppp / (2 * p)
        + qpp / 2
        - 3 * hpp / (2 * h)
        + 3 * hp**2 / (4 * h**2)
        + 3 * pp**2 / (2 * p**2)
        - 2 * ap**2
        - qp**2 / 2
        + ap * qp
        - pp / p * (qp / 2 - ap)
        - hp / (4 * h) * (qp + 3 * pp / p - 4 * ap)
        + np.exp(q - 4 * a) * h * p * wp**2 / 2
    )
    return -(h**3) * p**2 * rll


# ---------------------------------------------------------------- derived perturbation equations
_SYM = {n: sp.Symbol(n) for n in ("l", "r0", "c_omega", "A", "Ap", "Q", "Qp", "Omega2")}
_SYM.update({f"{k}{i}": sp.Symbol(f"{k}{i}") for k in ("g", "gp", "gpp") for i in range(3)})


def load_equations(path=EQ_FILE):
    rec = json.loads(Path(path).read_text())
    args = [_SYM[n] for n in ("l", "A", "Q", "c_omega")]
    V = {
        k: sp.lambdify(args, sp.sympify(v, locals=_SYM).subs(_SYM["r0"], 1), "numpy")
        for k, v in rec["G_system"]["V"].items()
    }
    cargs = [_SYM[n] for n in ("l", "A", "Ap", "Q", "Qp", "c_omega", "Omega2")]
    cargs += [_SYM[f"{k}{i}"] for k in ("g", "gp", "gpp") for i in range(3)]
    cons = {}
    for name, expr in rec["constraints"].items():
        e = sp.expand(sp.sympify(expr, locals=_SYM).subs(_SYM["r0"], 1))
        # keep the individual terms so the residual can be normalised by their size
        cons[name] = [sp.lambdify(cargs, term, "numpy") for term in sp.Add.make_args(e)]
    return V, cons


def vmatrix(V, ell, a, q, cw):
    out = np.zeros((3, 3, np.size(ell)))
    for k, f in V.items():
        out[int(k[0]), int(k[1])] = f(ell, a, q, cw) * np.ones(np.size(ell))
    return out


def operator(V, sol, cw, S, n, parity=+1):
    """K G = Omega^2 G on the s-grid (l = sinh s); returns K, the l-grid and w."""
    s = np.linspace(0, S, n + 1)
    ds = s[1]
    off = 1 if parity < 0 else 0
    si = s[off:-1]
    m = si.size
    ell, ls = np.sinh(si), np.cosh(si)
    sh = s[:-1] + ds / 2
    xih = p_of(np.sinh(sh)) * (np.sinh(sh) ** 2 + 1) ** 1.5 / np.cosh(sh)
    xi = p_of(ell) * (ell**2 + 1) ** 1.5
    D = np.zeros((m, m))
    for i in range(m):
        k = i + off
        cR = xih[k]
        cL = xih[k - 1] if k > 0 else xih[0]
        D[i, i] = -(cR + cL)
        if i + 1 < m:
            D[i, i + 1] = cR
        if k > 0 and i > 0:
            D[i, i - 1] = cL
        if k == 0:
            D[i, i + 1] += cL  # even parity: ghost G(-ds) = G(ds)
        D[i] /= xi[i] * ls[i] * ds**2
    a, _, q, _ = sol.sol(ell)
    Vm = vmatrix(V, ell, a, q, cw)
    w = p_of(ell) * np.exp(-2 * a - q)
    # (1/xi)(xi G')' - V G = -Omega^2 w B G  ->  K = -(1/w) B^-1 (D - V)
    K = np.zeros((3 * m, 3 * m))
    for i in range(3):
        for j in range(3):
            blk = (D if i == j else 0) - np.diag(Vm[i, j])
            for k2 in range(3):
                if BINV[k2, i]:
                    K[k2 * m : (k2 + 1) * m, j * m : (j + 1) * m] -= BINV[k2, i] * blk / w[:, None]
    return K, ell, w


def spectrum(V, sol, cw, S, n, parity=+1):
    K, ell, w = operator(V, sol, cw, S, n, parity)
    ev, vr = eig(K)
    ok = np.isfinite(ev)
    return ev[ok], vr[:, ok], ell, K


def eigpair_residual(K, lam, v):
    """Relative residual of a COMPLEX eigenpair (the full complex lam and v, never .real)."""
    return float(
        np.linalg.norm(K @ v - lam * v) / (np.linalg.norm(K @ v) + abs(lam) * np.linalg.norm(v))
    )


def constraint_residual(cons, sol, cw, lam, v, ell, S, parity=+1):
    """Relative residual of R_tl (corotating frame) and R_ll on the mode, 0 <= l <= sinh(S)/3.

    Normalised as max_l |sum of terms| / max_l sum |terms| (each equation separately)."""
    m = ell.size
    G = v.reshape(3, m)
    s = np.arcsinh(ell)
    ds = s[1] - s[0]
    # derivatives in s with the parity ghost at the throat, then d/dl = (1/cosh s) d/ds
    sign = 1 if parity > 0 else -1
    if parity > 0:
        ext = np.concatenate([sign * G[:, 1:2], G, np.zeros((3, 1))], axis=1)
    else:
        ext = np.concatenate([np.zeros((3, 1)), G, np.zeros((3, 1))], axis=1)
    Gs = (ext[:, 2:] - ext[:, :-2]) / (2 * ds)
    Gss = (ext[:, 2:] - 2 * ext[:, 1:-1] + ext[:, :-2]) / ds**2
    ls, lss = np.cosh(s), np.sinh(s)
    Gp = Gs / ls
    Gpp = (Gss - lss / ls * Gs) / ls**2
    a, ap, q, qp = sol.sol(ell)
    inner = ell <= np.sinh(S) / 3
    out = {}
    for name, terms in cons.items():
        vals = [t(ell, a, ap, q, qp, cw, lam, *G, *Gp, *Gpp) * np.ones(m) for t in terms]
        tot = np.abs(np.sum(vals, axis=0))[inner]
        scale = np.sum(np.abs(vals), axis=0)[inner]
        out[name] = float(np.max(tot) / np.max(scale))
    return out


# ---------------------------------------------------------------- independent: Evans shooting
class Shooter:
    """Decaying solutions integrated from l = sinh(S) to the throat; parity determinant."""

    def __init__(self, V, sol, cw, S=4.5):
        self.S = S
        s = np.linspace(0, S, 6001)
        ell = np.sinh(s)
        a, _, q, _ = sol.sol(ell)
        Vm = vmatrix(V, ell, a, q, cw)
        self.V = CubicSpline(s, np.moveaxis(Vm, 2, 0))
        self.w = CubicSpline(s, p_of(ell) * np.exp(-2 * a - q))

    def det(self, z2, parity=+1):
        S = self.S
        kap = np.sqrt(-z2 + 0j)
        kap = kap if kap.real > 0 else -kap
        L = np.sinh(S)
        xi = p_of(L) * (L * L + 1) ** 1.5
        Bc = BMAT.astype(complex)

        def f(s, Y):
            ell = np.sinh(s)
            xis = p_of(ell) * (ell * ell + 1) ** 1.5
            Y = Y.reshape(6, 3)
            G, P = Y[:3], Y[3:]
            dP = xis * ((self.V(s) - z2 * self.w(s) * Bc) @ G)
            return (np.vstack([P / xis, dP]) * np.cosh(s)).ravel()

        G0 = np.eye(3, dtype=complex)
        Y0 = np.vstack([G0, xi * (-kap - 1.5 / L) * G0]).ravel()
        y = solve_ivp(f, (S, 0.0), Y0, method="DOP853", rtol=1e-10, atol=1e-30).y[:, -1]
        Y = y.reshape(6, 3) * np.exp(-kap * L)  # analytic rescaling of the inward growth
        rows = Y[3:] if parity > 0 else Y[:3]
        return np.linalg.det(rows)

    def root(self, z0, parity=+1, tol=1e-9, floor=1e-6, maxit=30):
        """Secant iteration in complex Omega^2.

        Returns (root, info). The determinant loses digits to cancellation for large S, so the
        iteration can stagnate instead of reaching `tol`. The returned root is the iterate that
        followed the smallest relative step; acceptance (the rule applied since the verification
        run) is min_step <= floor. info also records the last step and every step, so the
        acceptance can be audited against a last-step rule."""
        x0, x1 = z0, z0 * (1 + 1e-4) + 1e-6
        f0, f1 = self.det(x0, parity), self.det(x1, parity)
        best, steps = x1, []
        for _ in range(maxit):
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
            st = abs(x2 - x1) / max(1, abs(x2))
            if not steps or st < min(steps):
                best = x2
            steps.append(st)
            if st < tol:
                break
            x0, f0, x1, f1 = x1, f1, x2, self.det(x2, parity)
        info = {
            "min_step": min(steps),
            "last_step": steps[-1],
            "iterations": len(steps),
            "reached_tol": bool(steps[-1] < tol),
            "steps": [float(x) for x in steps],
        }
        if info["min_step"] > floor:
            raise RuntimeError(f"secant iteration did not converge (min step {min(steps):.1e})")
        return best, info


# ---------------------------------------------------------------- tasks
# Pass criteria for a verified mode, fixed before the verification run (2026-10-14):
CRITERIA = {
    "grid": "dense Richardson (n, 2n) and shooting agree to <= 1e-5 relative",
    "box": "S = 6 and S = 7.5 at the same grid step agree to <= 1e-6 relative",
    "shooting": "outer radius sinh(3.5) vs sinh(4.0) agree to <= 1e-6 relative",
    "eigenpair": "complex eigenpair residual <= 1e-8",
    "constraints": "R_tl and R_ll residuals <= 1e-4 and falling >= 3x when the grid step halves",
    "localised": "eigenvector norm fraction beyond sinh(S)/3 <= 1e-6",
}


def pick(z, target):
    return int(np.argmin(np.abs(z - target)))


def track(cws, start_z, S=6.0, n=480, V=None, cons=None):
    """Follow one mode in c_omega by eigenvector overlap; returns one row per c_omega."""
    if V is None:
        V, cons = load_equations()
    rows, prev, vprev = [], None, None
    for cw in cws:
        sol = background_at(cw, prev)
        prev = sol
        P = background_props(sol, cw)
        ev, vr, ell, K = spectrum(V, sol, cw, S, n)
        z = ev * P["A"] ** (2 / 3)
        if vprev is None:
            i = pick(z, start_z)
            ov, ov2 = 1.0, None  # first point: no previous vector to compare with
        else:
            o = np.abs(vr.conj().T @ vprev) / np.linalg.norm(vr, axis=0) / np.linalg.norm(vprev)
            i = int(np.argmax(o))
            ov, ov2 = float(o[i]), float(np.sort(o)[-2])
        v = vr[:, i]
        root = np.sqrt(z[i] + 0j)
        rows.append(
            {
                "c_omega": cw,
                "J_over_A": P["J_over_A"],
                "z": [z[i].real, z[i].imag],
                "growth_A13": abs(root.imag),  # |Im Omega| A^(1/3)
                "freq_A13": abs(root.real),  # |Re Omega| A^(1/3)
                "overlap": ov,
                "overlap_next_best": ov2,
                "constraint_residual": constraint_residual(cons, sol, cw, ev[i], v, ell, S),
            }
        )
        vprev = v
    return rows


def _track_job(args):
    return track(*args)


def merger_point(V, n, lo=0.29, hi=0.33, iters=10, S=6.0):
    """Bisection in c_omega for the onset of a non-real pair among modes with Re z < -1."""

    def complex_pair(cw):
        sol = background_at(cw)
        P = background_props(sol, cw)
        ev, _, _, _ = spectrum(V, sol, cw, S, n)
        z = ev * P["A"] ** (2 / 3)
        z = np.sort_complex(z[z.real < -1])
        return bool(np.any(np.abs(z.imag) > 1e-8)), P["J_over_A"], z

    for _ in range(iters):
        mid = (lo + hi) / 2
        if complex_pair(mid)[0]:
            hi = mid
        else:
            lo = mid
    _, ja_lo, z_lo = complex_pair(lo)
    _, ja_hi, z_hi = complex_pair(hi)
    return {
        "n": n,
        "c_omega_bracket": [lo, hi],
        "J_over_A_bracket": [ja_lo, ja_hi],
        "z_re_below_minus1_at_lo": [[z.real, z.imag] for z in z_lo],
        "z_re_below_minus1_at_hi": [[z.real, z.imag] for z in z_hi],
    }


def _merger_job(n):
    V, _ = load_equations()
    return merger_point(V, n)


def shooting_and_checks(V, sol, cw, a23, dense):
    """Shooting at two radii and the CRITERIA checks, from the stored dense results."""

    def zc(key):
        return complex(*dense[key]["z"])

    rich = zc("S6.0_n960") + (zc("S6.0_n960") - zc("S6.0_n480")) / 3
    shoot = {}
    for S in (3.5, 4.0):
        r, info = Shooter(V, sol, cw, S=S).root(rich / a23)
        shoot[f"S{S}"] = {"z": [(r * a23).real, (r * a23).imag], **info}
    zs = complex(*shoot["S4.0"]["z"])
    rel = lambda x, y: abs(x - y) / abs(y)  # noqa: E731
    c_lo, c_hi = (
        dense["S6.0_n480"]["constraint_residual"],
        dense["S6.0_n960"]["constraint_residual"],
    )
    checks = {
        "grid": rel(rich, zs),
        "box": rel(zc("S7.5_n600"), zc("S6.0_n480")),
        "shooting": rel(complex(*shoot["S3.5"]["z"]), zs),
        "eigenpair": max(d["eigpair_residual"] for d in dense.values()),
        "constraints": {
            "max_at_n960": max(c_hi.values()),
            "min_ratio_n480_over_n960": min(c_lo[k] / c_hi[k] for k in c_lo),
        },
        "localised": max(d["outer_norm_fraction"] for k, d in dense.items() if k != "S5.0_n400"),
    }
    passed = {
        "grid": checks["grid"] <= 1e-5,
        "box": checks["box"] <= 1e-6,
        "shooting": checks["shooting"] <= 1e-6,
        "eigenpair": checks["eigenpair"] <= 1e-8,
        "constraints": checks["constraints"]["max_at_n960"] <= 1e-4
        and checks["constraints"]["min_ratio_n480_over_n960"] >= 3,
        "localised": checks["localised"] <= 1e-6,
    }
    return {
        "shooting": shoot,
        "checks": checks,
        "passed": passed,
        "status": "verified mode" if all(passed.values()) else "candidate",
    }


def verify_point(cw, start_z):
    """All checks of CRITERIA at one c_omega for the mode nearest start_z (scaled)."""
    V, cons = load_equations()
    sol = background_at(cw)
    P = background_props(sol, cw)
    a23 = P["A"] ** (2 / 3)
    out = {"c_omega": cw, "J_over_A": P["J_over_A"], "background": P}
    dense = {}
    for S, n in ((5.0, 400), (6.0, 480), (7.5, 600), (6.0, 960)):
        ev, vr, ell, K = spectrum(V, sol, cw, S, n)
        z = ev * a23
        i = pick(z, start_z)
        v = vr[:, i]
        G = v.reshape(3, ell.size)
        dense[f"S{S}_n{n}"] = {
            "z": [z[i].real, z[i].imag],
            "eigpair_residual": eigpair_residual(K, ev[i], v),
            "constraint_residual": constraint_residual(cons, sol, cw, ev[i], v, ell, S),
            "outer_norm_fraction": float(
                np.linalg.norm(G[:, ell > np.sinh(S) / 3]) / np.linalg.norm(G)
            ),
        }
    out["dense"] = dense

    def zc(key):
        return complex(*dense[key]["z"])

    out.update(shooting_and_checks(V, sol, cw, a23, dense))
    zs = complex(*out["shooting"]["S4.0"]["z"])
    root = np.sqrt(zs)
    out["result"] = {
        "z": [zs.real, zs.imag],
        "growth_A13": abs(root.imag),
        "freq_A13": abs(root.real),
        "growth_over_freq": abs(root.imag) / abs(root.real) if root.real else None,
    }
    # odd parity screen: unstable eigenvalues present in two boxes with different grid steps
    zz = []
    for S, n in ((6.0, 480), (7.5, 900)):
        ev, _, _, _ = spectrum(V, sol, cw, S, n, parity=-1)
        zz.append(ev * a23)
    out["odd_parity_unstable"] = [
        [z.real, z.imag]
        for z in zz[0]
        if (z.real < 0 or abs(z.imag) > 1e-6)
        and 0.3 < abs(z) < 60
        and np.min(np.abs(zz[1] - z)) < 2e-3 * max(1, abs(z))
    ]
    return out


def _verify_job(args):
    return verify_point(*args)


def static_gate():
    V, _ = load_equations()
    sol0 = solve_background(0.0)
    P = background_props(sol0, 0.0)
    a23 = P["A"] ** (2 / 3)
    st = []
    for n in (480, 960):
        ev, _, _, _ = spectrum(V, sol0, 0.0, 6.0, n)
        st.append(float(np.min(ev.real) * a23))
    r, info = Shooter(V, sol0, 0.0, S=4.0).root(st[-1] / a23)
    return {
        "max_abs_a_q": float(np.max(np.abs(sol0.y[[0, 2]]))),
        "dense_n480_n960": st,
        "richardson": st[1] + (st[1] - st[0]) / 3,
        "shooting": (r * a23).real,
        "shooting_imag": (r * a23).imag,
        "shooting_info": info,
        "paper_fig5_digitised": PAPER_FIG5["static"],
    }


def printed_v10_static(n=480, S=6.0):
    """Static spectrum with the printed (typo) V10 = derived + 16 e^{2q-2a}/h, for the record."""
    V, _ = load_equations()
    derived = V["10"]
    V["10"] = lambda ell, a, q, cw: (
        derived(ell, a, q, cw) + 16 * np.exp(2 * q - 2 * a) / (ell**2 + 1)
    )
    sol = solve_background(0.0)
    a23 = background_props(sol, 0.0)["A"] ** (2 / 3)
    ev, _, _, _ = spectrum(V, sol, 0.0, S, n)
    z = np.sort(ev.real[(ev.real < -1) & (np.abs(ev.imag) < 1e-8)] * a23)
    return {"n": n, "S": S, "z_real_below_minus1": z.tolist()}


def branch_vs_fig5(rows, paper):
    """Tracked branch, linearly interpolated at the digitised J/A (inside the tracked range)."""
    x = np.array([r["J_over_A"] for r in rows])
    y = np.array([r["z"][0] for r in rows])
    return [
        {"J_over_A": ja, "paper_digitised": zp, "ours_linear_interp": float(np.interp(ja, x, y))}
        for ja, zp in paper
        if x[0] <= ja <= x[-1]
    ]


def main() -> int:
    t0 = time.time()
    rec = {
        "label": "derived; linear radial-sector modes of 1309.2448 symmetric 5D equal-J wormholes",
        "conventions": "exp(i Omega t); z = Omega^2 A^(2/3); growth rate |Im Omega|; 4 pi G = 1",
        "criteria": CRITERIA,
    }
    gates, prev = [], None
    for cw in (0.0, 0.2, 0.5, 1.0, 2.0, 4.0):
        prev = solve_background(cw, guess=prev)
        gates.append(background_props(prev, cw))
    rec["background_gates"] = gates
    rec["fig1_comparison"] = [
        {"J_over_A": g["J_over_A"], "M_over_A23": g["M_over_A23"], "paper_digitised": m}
        for g, (_, m) in zip(gates[1:], PAPER_FIG1, strict=True)
    ]
    below = np.linspace(0.14, 0.30, 9)
    reps = [(0.4, -4.1046 + 5.5275j), (0.5, -2.3186 + 7.4618j), (1.0, 2.4 + 8.5j)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        f_static = ex.submit(static_gate)
        f_v10 = ex.submit(printed_v10_static)
        f_track = [ex.submit(_track_job, (below, z)) for z in (-18.7, -0.33)]
        f_merge = [ex.submit(_merger_job, n) for n in (480, 960)]
        f_rep = [ex.submit(_verify_job, r) for r in reps]
        rec["merger"] = [f.result() for f in f_merge]
        # the pair is tracked from the upper end of the merger bracket (no hand-chosen start)
        m0 = rec["merger"][0]
        start = max((complex(*z) for z in m0["z_re_below_minus1_at_hi"]), key=lambda z: z.imag)
        above = np.concatenate([[m0["c_omega_bracket"][1]], np.linspace(0.33, 1.0, 12)])
        f_pair = ex.submit(_track_job, (above, start))
        rec["static_mode"] = f_static.result()
        rec["printed_V10_static"] = f_v10.result()
        rec["branch_static"], rec["branch_zero"] = [f.result() for f in f_track]
        rec["pair"] = f_pair.result()
        rec["representative_points"] = [f.result() for f in f_rep]
    rec["branch_vs_fig5"] = {
        k: branch_vs_fig5(rec[k], PAPER_FIG5[k]) for k in ("branch_static", "branch_zero")
    }
    rec["paper_fig5_digitised"] = PAPER_FIG5
    rec["runtime_s"] = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    text = json.dumps(rec, indent=1, default=float, allow_nan=False)
    (OUT / "s6a_rot5d.json").write_text(text + "\n")
    for r in rec["representative_points"]:
        print(r["c_omega"], round(r["J_over_A"], 5), r["status"], r["result"], r["passed"])
    print("merger", rec["merger"], "static", rec["static_mode"])
    return 0


def redo_shooting(path=OUT / "s6a_rot5d.json"):
    """Re-run only the shooting (and the checks that use it) on the stored dense results."""
    rec = json.loads(Path(path).read_text())
    V, _ = load_equations()
    for pt in rec["representative_points"]:
        cw = pt["c_omega"]
        sol = background_at(cw)
        a23 = background_props(sol, cw)["A"] ** (2 / 3)
        pt["shooting_before_rerun"] = pt["shooting"]
        pt.update(shooting_and_checks(V, sol, cw, a23, pt["dense"]))
        root = np.sqrt(complex(*pt["shooting"]["S4.0"]["z"]))
        pt["result"]["z"] = pt["shooting"]["S4.0"]["z"]
        pt["result"]["growth_A13"] = abs(root.imag)
        pt["result"]["freq_A13"] = abs(root.real)
    rec["static_mode"] = static_gate()
    for key in ("branch_static", "branch_zero", "pair"):
        for row in rec[key]:
            if row["overlap_next_best"] != row["overlap_next_best"]:  # NaN from the first point
                row["overlap_next_best"] = None
    rec["shooting_rerun"] = "2026-10-14: shooting re-run with full step records (redo_shooting)"
    Path(path).write_text(json.dumps(rec, indent=1, default=float, allow_nan=False) + "\n")
    for pt in rec["representative_points"]:
        sh = pt["shooting"]
        print(
            pt["c_omega"], pt["status"], {k: (v["min_step"], v["last_step"]) for k, v in sh.items()}
        )
    return 0


if __name__ == "__main__":
    import sys

    raise SystemExit(redo_shooting() if "--shooting-only" in sys.argv else main())

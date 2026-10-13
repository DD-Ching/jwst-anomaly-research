"""X2 S5: radial (l = 0) growing modes of the S4 background (SV, phantom + magnetic NED), one cell.

Equations: results/x2/s5_radial_equations.json from scripts/x2_s5_radial_derive.py (CPS framework,
delta Phi = 0 gauge, gauge-invariant A = delta gamma / gamma, C = delta r / r). A evolves with its
wave equation; C is rebuilt at every stage from the constraint sigma = 0,
  C' = (alpha'/alpha - r'/r) C + (r'/r) A   ->   C = (alpha / r) int (r' / alpha) A dx,
so the evolution never leaves the constraint surface (evolving C on its own admits a spurious,
constraint-violating pointwise mode). The C (xx) equation is not used; it is an independent check.
Grid: uniform in the tortoise coordinate y (dy = dx / f), Dirichlet ends far from the mode.

Methods:
  * time evolution (RK4) from localised and noisy initial data; growth rate from log-amplitude fits;
  * dense spectrum of the discretised operator (largest real eigenvalues lambda^2);
  * reference: reflection-symmetric Ellis-Bronnikov wormhole, against CPS's master equation.

  python scripts/x2_s5_radial.py   -> results/x2/s5_radial.json   (about 15 minutes on 4 processes)
"""

from __future__ import annotations

import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
OUT = ROOT / "results" / "x2"
NAMES = ("al0", "ga0", "r0", "P0", "Vx")
ARGS = [sp.Symbol(f"{n}_{k}") for n in NAMES for k in range(4)] + [
    sp.Symbol(f"U_{k}") for k in range(3)
]
T_SUN = 4.925490947e-6
M_DET = 68.8


def load_coefficients(path=OUT / "s5_radial_equations.json"):
    rec = json.loads(Path(path).read_text())
    loc = {str(s): s for s in ARGS}
    return {
        eq: {k: sp.lambdify(ARGS, sp.sympify(v, locals=loc), "numpy") for k, v in d.items()}
        for eq, d in rec["coefficients"].items()
        if eq in ("Att", "Ctt")
    }


def background(kind, par, x):
    X = sp.Symbol("X", real=True)
    if kind == "EB":
        b = par["b"]
        R = sp.sqrt(X**2 + b**2)
        fns = {
            "al0": sp.Integer(1),
            "ga0": sp.Integer(1),
            "r0": R,
            "P0": sp.sqrt(2) * sp.atan(X / b),
            "Vx": sp.Integer(0),
        }
        Ur = [lambda r: 0 * r] * 3
        f = sp.Integer(1)
        rr = np.sqrt(x**2 + b**2)
    else:
        M, a = par["M"], par["a"]
        R = sp.sqrt(X**2 + a**2)
        f = 1 - 2 * M / R
        fns = {
            "al0": sp.sqrt(f),
            "ga0": 1 / sp.sqrt(f),
            "r0": R,
            "P0": sp.sqrt(2) * sp.atan(X / a),
            "Vx": 4 * M * a**2 / (5 * R**5),
        }
        Ur = [
            lambda r: 6 * M * a**2 / (5 * r**5),
            lambda r: -6 * M * a**2 / r**6,
            lambda r: 36 * M * a**2 / r**7,
        ]
        rr = np.sqrt(x**2 + a**2)
    vals = []
    for n in NAMES:
        for k in range(4):
            fn = sp.lambdify(X, sp.diff(fns[n], X, k) if k else fns[n], "numpy")
            vals.append(np.broadcast_to(np.asarray(fn(x), float), x.shape).copy())
    vals += [np.asarray(u(rr), float) * np.ones_like(x) for u in Ur]
    speed = np.asarray(sp.lambdify(X, f, "numpy")(x), float) * np.ones_like(x)
    return vals, speed


def setup(kind, par, Y, dy, CF):
    y = np.arange(-Y, Y + dy / 2, dy)
    if kind == "EB":
        x = y.copy()
    else:
        import x2_scatter as xs

        g = xs.Grid(par["a"], 2, 1, xmax=Y + 1, h=dy / 4)
        x = np.sign(y) * np.interp(np.abs(y), g.x, g.r)
    vals, f = background(kind, par, x)
    c = {
        eq: {k: np.asarray(fn(*vals), float) * np.ones_like(x) for k, fn in d.items()}
        for eq, d in CF.items()
    }
    al, al1, al2 = vals[0], vals[1], vals[2]
    r, r1, r2 = vals[8], vals[9], vals[10]
    return dict(
        y=y,
        x=x,
        f=f,
        c=c,
        dy=dy,
        r=r,
        al=al,
        k1=al1 / al - r1 / r,
        k2=r1 / r,
        k1p=al2 / al - (al1 / al) ** 2 - r2 / r + (r1 / r) ** 2,
        k2p=r2 / r - (r1 / r) ** 2,
    )


def dx1(F, S):
    d = np.zeros_like(F)
    d[1:-1] = (F[2:] - F[:-2]) / (2 * S["dy"])
    return d / S["f"]


def dx2(F, S):
    dy, f = S["dy"], S["f"]
    flux = (F[1:] - F[:-1]) / dy / (0.5 * (f[1:] + f[:-1]))
    d = np.zeros_like(F)
    d[1:-1] = (flux[1:] - flux[:-1]) / dy / f[1:-1]
    return d


def c_from_a(A, S):
    integrand = S["k2"] * S["r"] / S["al"] * A  # (r' / alpha) A
    cum = np.concatenate(
        [[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(S["x"]))]
    )
    return S["al"] / S["r"] * cum


def _terms(A, S):
    C = c_from_a(A, S)
    Ap = dx1(A, S)
    Cp = S["k1"] * C + S["k2"] * A
    Cpp = S["k1p"] * C + S["k1"] * Cp + S["k2p"] * A + S["k2"] * Ap
    return C, Ap, Cp, Cpp


def rhs(A, S):
    C, Ap, Cp, Cpp = _terms(A, S)
    g = S["c"]["Att"].get
    out = (
        g("Derivative(A(t, x), (x, 2))", 0.0) * dx2(A, S)
        + g("Derivative(C(t, x), (x, 2))", 0.0) * Cpp
        + g("Derivative(A(t, x), x)", 0.0) * Ap
        + g("Derivative(C(t, x), x)", 0.0) * Cp
        + g("A(t, x)", 0.0) * A
        + g("C(t, x)", 0.0) * C
    )
    out[0] = out[-1] = 0.0
    return out


def c_equation_residual(A, lam2, S):
    """Independent check: lambda^2 C = (C equation RHS) for a mode A ~ e^{lambda t}."""
    C, Ap, Cp, _ = _terms(A, S)
    g = S["c"]["Ctt"].get
    rc = (
        g("C(t, x)", 0.0) * C
        + g("Derivative(C(t, x), x)", 0.0) * Cp
        + g("A(t, x)", 0.0) * A
        + g("Derivative(A(t, x), x)", 0.0) * Ap
    )
    m = slice(20, -20)
    return float(np.max(np.abs(lam2 * C[m] - rc[m])) / np.max(np.abs(lam2 * C[m]) + np.abs(rc[m])))


def evolve(S, T, A0):
    dt = 0.5 * S["dy"]
    A, V = A0.copy(), np.zeros_like(A0)
    every = max(1, int(round(1.0 / dt)))
    ts, amp = [], []
    for i in range(int(T / dt)):
        k1a, k1v = V, rhs(A, S)
        k2a, k2v = V + 0.5 * dt * k1v, rhs(A + 0.5 * dt * k1a, S)
        k3a, k3v = V + 0.5 * dt * k2v, rhs(A + 0.5 * dt * k2a, S)
        k4a, k4v = V + dt * k3v, rhs(A + dt * k3a, S)
        A = A + dt / 6 * (k1a + 2 * k2a + 2 * k3a + k4a)
        V = V + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        if i % every == 0:
            ts.append((i + 1) * dt)
            amp.append(float(np.max(np.abs(A[10:-10]))))
            nrm = amp[-1]
            if nrm > 1e100:
                A, V = A / nrm, V / nrm
                amp = list(np.array(amp) / nrm)  # keep log-slope; rescale history
    ts, amp = np.array(ts), np.array(amp)
    return ts, amp, A


def fit(ts, amp, frac):
    m = ts > ts[-1] * frac
    return float(np.polyfit(ts[m], np.log(amp[m]), 1)[0])


def dense_top(S, k=3):
    n = S["y"].size
    L = np.zeros((n - 2, n - 2))
    e = np.zeros(n)
    for j in range(1, n - 1):
        e[:] = 0.0
        e[j] = 1.0
        L[:, j - 1] = rhs(e, S)[1:-1]
    w, v = np.linalg.eig(L)
    order = np.argsort(-w.real)[:k]
    y = S["y"][1:-1]
    return [
        {
            "lam2_re": float(w[i].real),
            "lam2_im": float(w[i].imag),
            "y_peak": float(abs(y[np.argmax(np.abs(v[:, i]))])),
            # independent check (C / xx equation): physical modes pass it, spurious ones fail
            "c_equation_residual": c_equation_residual(
                np.concatenate([[0.0], v[:, i].real, [0.0]]), float(w[i].real), S
            ),
        }
        for i in order
    ]


def sv_a():
    import x2_axial_ned as an
    import x2_echo_delay as xe

    return 2 + xe.a_for_delay(0.45 / (an.M_DET * an.T_SUN))


def job(spec):
    CF = load_coefficients()
    name, kind, Y, dy, mode, arg = spec
    par = {"b": 1.0} if kind == "EB" else {"M": 1.0, "a": sv_a()}
    S = setup(kind, par, Y, dy, CF)
    if mode == "dense":
        return name, {"Y": Y, "dy": dy, "top": dense_top(S)}
    y = S["y"]
    if mode == "noise":
        rng = np.random.default_rng(1)
        A0 = np.convolve(rng.normal(size=y.size), np.ones(5) / 5, "same") * np.exp(
            -((y / (0.8 * Y)) ** 8)
        )
    else:
        A0 = np.exp(-((y - arg) ** 2) / 18.0) + np.exp(-((y + arg) ** 2) / 18.0)
    T = 40.0 if kind == "EB" else 300.0
    ts, amp, A = evolve(S, T, A0)
    lam = fit(ts, amp, 0.75)
    m = np.abs(A[10:-10])
    return name, {
        "Y": Y,
        "dy": dy,
        "T": T,
        "init": mode if mode == "noise" else f"gaussians at y = +-{arg}",
        "growth_last_half": fit(ts, amp, 0.5),
        "growth_last_quarter": lam,
        "y_peak_end": float(abs(y[10:-10][np.argmax(m)])),
        "x_peak_end": float(abs(S["x"][10:-10][np.argmax(m)])),
        "c_equation_residual": c_equation_residual(A, lam * lam, S),
    }


def main() -> int:
    jobs = [
        ("EB_time", "EB", 40.0, 0.05, "gauss", 0.0),
        ("EB_dense", "EB", 30.0, 0.05, "dense", None),
        ("SV_time_dy0.1", "SV", 600.0, 0.1, "gauss", 315.0),
        ("SV_time_dy0.05", "SV", 600.0, 0.05, "gauss", 315.0),
        ("SV_time_from200", "SV", 600.0, 0.1, "gauss", 200.0),
        ("SV_time_noise", "SV", 600.0, 0.1, "noise", None),
        ("SV_dense_Y1000_dy0.5", "SV", 1000.0, 0.5, "dense", None),
        ("SV_dense_Y2000_dy1", "SV", 2000.0, 1.0, "dense", None),
    ]
    with ProcessPoolExecutor(max_workers=min(len(jobs), os.cpu_count() or 1)) as ex:
        res = dict(ex.map(job, jobs))
    # Ellis-Bronnikov reference: lowest eigenvalue of -d^2/dx^2 - 3 b^2 / (x^2 + b^2)^2
    # (CPS Eq. MasterEllis)
    from scipy.linalg import eigh_tridiagonal

    xg = np.arange(-30.0, 30.0 + 0.025, 0.05)
    e0 = eigh_tridiagonal(
        2 / 0.05**2 - 3 / (xg**2 + 1) ** 2,
        -np.ones(xg.size - 1) / 0.05**2,
        select="i",
        select_range=(0, 0),
    )[0][0]
    a = sv_a()
    lam = res["SV_time_dy0.05"]["growth_last_quarter"]
    s4 = json.loads((OUT / "s4_axial_ned.json").read_text())["comparison_gw150914_cell"]["ned"][
        "dx0.1"
    ]
    tm = M_DET * T_SUN
    t1 = s4["pulses"][0]["arrival_after_free_flight_M"]
    xp = res["SV_time_dy0.05"]["x_peak_end"]
    Rp = math.hypot(xp, a)
    out = {
        "label": "derived; radial l = 0; SV + phantom + magnetic NED (S4 source); one cell",
        "a_minus_2M": a - 2,
        "runs": res,
        "ellis_bronnikov_reference_growth": math.sqrt(-e0),
        "growth_rate_per_M": lam,
        "e_folding_time_M": 1 / lam,
        "e_folding_time_s_GW150914": tm / lam,
        "mode_location": {
            "r_star_M": res["SV_time_dy0.05"]["y_peak_end"],
            "areal_radius_M": Rp,
            "f_there": 1 - 2 / Rp,
            "local_proper_rate_per_M": lam / math.sqrt(1 - 2 / Rp),
        },
        "clock": "coordinate time t of the static background = asymptotic observer time "
        "(alpha -> 1), "
        "the clock of S4's echo times",
        "e_folds": {
            f"per_echo_spacing (S4: {s4['spacing_M']:.1f} M)": lam * s4["spacing_M"],
            f"to_first_echo (S4: {t1:.1f} M after free flight)": lam * t1,
            "to_PMW_start_0.2s": lam * 0.2 / tm,
            "to_PMW_end_1.2s": lam * 1.2 / tm,
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s5_radial.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

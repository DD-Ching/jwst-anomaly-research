"""X2 S6a: re-derive the linear perturbation equations of 5D equal-J rotating Ellis wormholes.

Source: Dzhunushaliev et al., arXiv:1309.2448 (PRD 88, 124028), Sec. V, Eq. (pertlineel).
Field equations R_mu_nu = -8 pi G d_mu Phi d_nu Phi (massless phantom); only the l-derivative
of the scalar is non-zero in the background, so with H_phi = 0 every linearised component is
delta R_mu_nu = 0.

Angular part rewritten exactly in Hopf form (theta~ = 2 theta, chi = psi - phi, eta = psi + phi):
  A h [dtheta^2 + s^2 c^2 (dpsi - dphi)^2] + B h [s^2 (dphi - W dt) + c^2 (dpsi - W dt)]^2
  = A h/4 (dtheta~^2 + sin^2 theta~ dchi^2) + B h/4 (deta + cos theta~ dchi - 2 W dt)^2,
which uses s^2 dphi^2 + c^2 dpsi^2 = (s^2 dphi + c^2 dpsi)^2 + s^2 c^2 (dpsi - dphi)^2.

Checks done here (each prints its residual; all must be 0):
  1. background: every Ricci component except R_ll vanishes on Eqs. (eqa), (eqq), (omprime), (solp);
     R_ll equals the paper's Eq. (Rll-expr);
  2. gauge: a pure-gauge vector field shows that H1 - H3 + H2 can be set to zero, and the scalar
     equation then forces it to a constant when H_phi = 0 (holds for complex Omega);
  3. the G0, G1, G2 system (time dependence exp(i Omega t)) against the printed V_ij.

Writes results/x2/s6a_equations.json: the derived G-system matrix, H_omega', and the two constraint
equations (R_tl in the corotating frame and R_ll) used later as independent residual checks.

  python scripts/x2_s6a_derive5d.py   (about 5 minutes; symbolic)
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x2"

t, ell, th, ch, et = sp.symbols("t l theta chi eta", real=True)
X = [t, ell, th, ch, et]
r0, cw = sp.symbols("r0 c_omega", positive=True)
Om, Om2, eps = sp.symbols("Omega Omega2 epsilon")
hfun, pfun = sp.Function("h")(ell), sp.Function("p")(ell)
a, q, w = [sp.Function(n)(ell) for n in ("a", "q", "w")]
H_EXPL = ell**2 + r0**2
P_EXPL = (ell**2 + r0**2 / 2) / (ell**2 + r0**2)
XI = P_EXPL * H_EXPL ** sp.Rational(3, 2)


def metric(F0, F1, F2, F3, Fw, p=pfun, h=hfun):
    """Eq. (pertlineel) in Hopf coordinates (t, ell, theta~, chi, eta)."""
    c, s = sp.cos(th), sp.sin(th)
    pre = p * sp.exp(-F0)
    N = sp.exp(2 * F0)
    L = pre * sp.exp(F0 - F2 + F1)
    S = pre * sp.exp(F0 - F2) * h / 4
    K = pre * sp.exp(F3 - F0) * h / 4
    one = [-2 * Fw, 0, 0, c, 1]
    g = sp.Matrix(5, 5, lambda i, j: K * one[i] * one[j])
    g[0, 0] += -N
    g[1, 1] += L
    g[2, 2] += S
    g[3, 3] += S * s**2
    gi = sp.zeros(5, 5)
    gi[0, 0] = -1 / N
    gi[0, 4] = gi[4, 0] = -2 * Fw / N
    gi[1, 1] = 1 / L
    gi[2, 2] = 1 / S
    gi[3, 3] = 1 / (S * s**2)
    gi[3, 4] = gi[4, 3] = -c / (S * s**2)
    gi[4, 4] = 1 / K + c**2 / (S * s**2) - 4 * Fw**2 / N
    return g, gi


def ricci(g, gi, comps):
    dg = [[[sp.diff(g[i, j], X[k]) for k in range(5)] for j in range(5)] for i in range(5)]
    Gam = [
        [
            [
                sp.expand(
                    sum(gi[i, m] * (dg[m][j][k] + dg[m][k][j] - dg[j][k][m]) for m in range(5)) / 2
                )
                for k in range(5)
            ]
            for j in range(5)
        ]
        for i in range(5)
    ]

    def ric(j, k):
        return sum(
            sp.diff(Gam[i][j][k], X[i])
            - sp.diff(Gam[i][j][i], X[k])
            + sum(Gam[i][i][m] * Gam[m][j][k] - Gam[i][k][m] * Gam[m][j][i] for m in range(5))
            for i in range(5)
        )

    return {n: ric(j, k) for n, (j, k) in comps.items()}


COMPS = {
    "tt": (0, 0),
    "tl": (0, 1),
    "ll": (1, 1),
    "thth": (2, 2),
    "etet": (4, 4),
    "tet": (0, 4),
    "let": (1, 4),
}


def background_reducers():
    """a'', q'' from (eqa), (eqq); omega' from (omprime); explicit p, h from (solp)."""
    wp = cw / (H_EXPL ** sp.Rational(5, 2) * P_EXPL**2 * sp.exp(q - 4 * a))
    app = sp.solve(
        sp.diff(XI * sp.diff(a, ell), ell) / XI
        - cw**2 * sp.exp(4 * a - q) / (2 * P_EXPL**3 * H_EXPL**4),
        sp.diff(a, ell, 2),
    )[0]
    qpp = sp.solve(
        sp.diff(XI * sp.diff(q, ell), ell) / XI - 4 * (sp.exp(2 * (q - a)) - 1) / H_EXPL,
        sp.diff(q, ell, 2),
    )[0]

    def red(e):
        e = e.subs({pfun: P_EXPL, hfun: H_EXPL}).doit()
        e = e.subs(sp.Derivative(w, (ell, 2)), sp.diff(wp, ell)).subs(sp.Derivative(w, ell), wp)
        e = e.subs(sp.Derivative(a, (ell, 3)), sp.diff(app, ell))
        e = e.subs(sp.Derivative(q, (ell, 3)), sp.diff(qpp, ell))
        return e.subs(sp.Derivative(a, (ell, 2)), app).subs(sp.Derivative(q, (ell, 2)), qpp)

    return red


def check_background(R, F, red):
    F0, F1, F2, F3, Fw = F
    bg = {F0: a, F1: 0, F2: q, F3: q, Fw: w}
    res = {n: str(sp.simplify(red(R[n].subs(bg).doit()))) for n in R if n != "ll"}
    ap, qp, pp, hp = [sp.diff(f, ell) for f in (a, q, pfun, hfun)]
    p, h = pfun, hfun
    rll_paper = (
        -3 * sp.diff(p, ell, 2) / (2 * p)
        + sp.diff(q, ell, 2) / 2
        - 3 * sp.diff(h, ell, 2) / (2 * h)
        + 3 * hp**2 / (4 * h**2)
        + 3 * pp**2 / (2 * p**2)
        - 2 * ap**2
        - qp**2 / 2
        + ap * qp
        - pp / p * (qp / 2 - ap)
        - hp / (4 * h) * (qp + 3 * pp / p - 4 * ap)
        + sp.exp(q - 4 * a) * h * p * sp.diff(w, ell) ** 2 / 2
    )
    res["ll_minus_paper_Rll"] = str(sp.simplify(R["ll"].subs(bg).doit() - rll_paper))
    return res


def check_gauge():
    """Lie derivative of the background along xi = exp(i Omega t) (T, Z, 0, 0, Y)."""
    T, Z, Y = [sp.Function(n)(ell) for n in ("T", "Z", "Y")]
    g0, _ = metric(a, 0, q, q, w)
    E = sp.exp(sp.I * Om * t)
    xi = [T * E, Z * E, 0, 0, Y * E]
    Lg = sp.Matrix(
        5,
        5,
        lambda i, j: sum(
            xi[k] * sp.diff(g0[i, j], X[k])
            + g0[k, j] * sp.diff(xi[k], X[i])
            + g0[i, k] * sp.diff(xi[k], X[j])
            for k in range(5)
        ),
    )
    H = sp.symbols("H0:4")
    Hw = sp.Symbol("Hw")
    e = sp.Symbol("e")
    gp, _ = metric(
        a + e * H[0] * E, e * H[1] * E, q + e * H[2] * E, q + e * H[3] * E, w + e * Hw * E
    )
    d = gp.diff(e).subs(e, 0)
    pairs = [(0, 0), (1, 1), (2, 2), (4, 4), (0, 4), (3, 4)]
    eqs = [sp.simplify((Lg[i, j] - d[i, j]) / E) for i, j in pairs]
    sol = sp.solve(eqs, list(H) + [Hw], dict=True)[0]
    # ansatz has g_tl = g_l,eta = 0: fixes Y' = 2 w T' and T' = i Omega p e^{-2a-q} Z
    Tp = sp.solve(
        sp.simplify(Lg[0, 1] / E).subs(sp.Derivative(Y, ell), 2 * w * sp.Derivative(T, ell)),
        sp.Derivative(T, ell),
    )[0]
    shift = sp.simplify(sol[H[1]] - sol[H[3]] + sol[H[2]])
    return {
        "g_l_eta": str(sp.factor(sp.simplify(Lg[1, 4] / E))),
        "T_prime": str(sp.simplify(Tp)),
        "shift_H1_minus_H3_plus_H2": str(shift),
        "delta_H_phi": "Z * phi'",
        "check_33_component": str(sp.simplify((Lg[3, 3] - d[3, 3]).subs(sol))),
    }


def linearise(R, F, red):
    F0, F1, F2, F3, Fw = F
    G = [sp.Function(n)(ell) for n in ("G0", "G1", "G2")]
    Hw = sp.Function("Hw")(ell)
    E = sp.exp(sp.I * Om * t)
    H0, H2, H3 = G[2] + G[0], -G[2], 2 * G[0] + G[1]
    sub = {
        F0: a + eps * H0 * E,
        F1: eps * (H3 - H2) * E,  # gauge (gauge): H1 = H3 - H2
        F2: q + eps * H2 * E,
        F3: q + eps * H3 * E,
        Fw: w + eps * Hw * E,
    }
    Lin = {}
    for n, e in R.items():
        e = sp.diff(e.subs(sub).doit(), eps).subs(eps, 0)
        Lin[n] = sp.expand(red(sp.expand(e / E)).subs(Om**2, Om2))
    # corotating frame e_t = d_t + 2 w d_eta (background R_t_eta = 0, so these are exact)
    C = {
        "Ett": Lin["tt"] + 4 * w * Lin["tet"] + 4 * w**2 * Lin["etet"],
        "Etl": Lin["tl"] + 2 * w * Lin["let"],
        "ll": Lin["ll"],
        "thth": Lin["thth"],
        "etet": Lin["etet"],
        "leta": Lin["let"],
    }
    A, Ap, Q, Qp = sp.symbols("A Ap Q Qp")
    g, gp, gpp = sp.symbols("g0:3"), sp.symbols("gp0:3"), sp.symbols("gpp0:3")
    hwp = sp.Symbol("hwp")

    def jet(e):
        for i in range(3):
            e = e.subs(sp.Derivative(G[i], (ell, 2)), gpp[i]).subs(sp.Derivative(G[i], ell), gp[i])
            e = e.subs(G[i], g[i])
        e = e.subs(sp.Derivative(Hw, ell), hwp)
        return e.subs(sp.Derivative(a, ell), Ap).subs(sp.Derivative(q, ell), Qp).subs({a: A, q: Q})

    J = {n: sp.simplify(jet(e)) for n, e in C.items()}
    hw = sp.simplify(sp.solve(J["leta"], hwp)[0])
    eqs = [sp.simplify(J[n].subs(hwp, hw)) for n in ("Ett", "etet", "thth")]
    sol = sp.solve(eqs, gpp, dict=True)[0]
    M = {}
    firstder = {}
    for i in range(3):
        e = sp.expand(sp.simplify(sol[gpp[i]]))
        for j in range(3):
            firstder[f"{i}{j}"] = str(
                sp.simplify(e.coeff(gp[j]) + (sp.diff(XI, ell) / XI if i == j else 0))
            )
        e0 = e.subs({gp[0]: 0, gp[1]: 0, gp[2]: 0})
        for j in range(3):
            M[(i, j)] = sp.simplify(e0.coeff(g[j]))
    cons = {
        "Etl_over_Omega": sp.simplify(sp.expand(J["Etl"].subs(hwp, hw)).coeff(Om, 1)),
        "ll": sp.simplify(J["ll"].subs(hwp, hw)),
    }
    return M, firstder, hw, cons, (A, Ap, Q, Qp, g, gp, gpp)


def printed_V(A, Q):
    """1309.2448 Eqs. (eqv00)-(eqv22) as printed."""
    h, p = H_EXPL, P_EXPL
    E2a, E2q, Eq, c2, E6a = sp.exp(2 * A), sp.exp(2 * Q), sp.exp(Q), cw**2, sp.exp(6 * A)
    k = sp.exp(-2 * A - Q) / (h**4 * p**3)
    hp3 = h**3 * p**3
    return {
        (0, 0): -(8 * (2 * E2a - E2q) * Eq * hp3 - c2 * E6a) * k,
        (0, 1): -(16 * (E2a - E2q) * Eq * hp3 + c2 * E6a) * k / 2,
        (0, 2): -2 * (4 * (E2a + E2q) * Eq * hp3 - c2 * E6a) * k,
        (1, 0): 2 * (4 * (4 * E2a + E2q) * Eq * hp3 - c2 * E6a) * k,
        (1, 1): (8 * (2 * E2a - E2q) * Eq * hp3 + c2 * E6a) * k,
        (1, 2): 4 * (2 * (2 * E2a + E2q) * Eq * hp3 - c2 * E6a) * k,
        (2, 0): 8 * (2 * E2a - E2q) / h * sp.exp(-2 * A),
        (2, 1): 8 * (E2a - E2q) / h * sp.exp(-2 * A),
        (2, 2): 8 * (E2a + E2q) / h * sp.exp(-2 * A),
    }


def main() -> int:
    F = [sp.Function(n)(t, ell) for n in ("F0", "F1", "F2", "F3", "Fw")]
    g, gi = metric(*F)
    assert sp.simplify(gi * g - sp.eye(5)) == sp.zeros(5, 5)
    R = ricci(g, gi, COMPS)
    red = background_reducers()
    bgres = check_background(R, F, red)
    print("background residuals:", bgres)
    gauge = check_gauge()
    print("gauge:", gauge)
    M, firstder, hw, cons, (A, Ap, Q, Qp, gs, gps, gpps) = linearise(R, F, red)
    wt = sp.exp(-2 * A - Q) * P_EXPL
    B = [[1, 1, 0], [0, 1, 0], [0, 0, 1]]
    V = {k: sp.simplify(v + Om2 * wt * B[k[0]][k[1]]) for k, v in M.items()}
    pv = printed_V(A, Q)
    diff = {f"V{i}{j}": str(sp.simplify(V[(i, j)] - pv[(i, j)])) for (i, j) in V}
    print("first-derivative structure (0 = paper form):", firstder)
    print("derived V - printed V:", diff)
    # printed Eq. (hop): 8 H0 + H1 + H2 - 3 H3 with H1 = H3 - H2
    h0, h2, h3 = gs[2] + gs[0], -gs[2], 2 * gs[0] + gs[1]
    hw_paper_bracket = 8 * h0 + (h3 - h2) + h2 - 3 * h3
    hw_paper = cw * sp.exp(4 * A - Q) / (P_EXPL**2 * H_EXPL ** sp.Rational(5, 2)) * hw_paper_bracket
    rec = {
        "label": "derived; linear perturbations of 1309.2448 5D equal-J rotating Ellis wormholes",
        "time_convention": (
            "exp(i Omega t); growth iff Im Omega < 0; unstable iff Omega^2 not in [0, inf)"
        ),
        "units": "J/A and M/A^(2/3) use 4 pi G = 1 (1309.2448 figures; see unit notes)",
        "symbols": "A=a(ell), Ap=a'(ell), Q=q(ell), Qp=q'(ell), g_i=G_i, gp_i=G_i', gpp_i=G_i''",
        "background_residuals": bgres,
        "gauge": gauge,
        "G_system": {
            "form": "G_i'' = -(xi'/xi) G_i' + sum_j (V_ij - Omega2 w B_ij) G_j, w = p e^{-2a-q}",
            "first_derivative_residual": firstder,
            "V": {f"{i}{j}": str(v) for (i, j), v in V.items()},
            "V_derived_minus_printed": diff,
        },
        "H_omega_prime": str(hw),
        "H_omega_prime_over_printed_eq_hop": str(sp.simplify(hw / hw_paper)),
        "constraints": {k: str(v) for k, v in cons.items()},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s6a_equations.json").write_text(json.dumps(rec, indent=1) + "\n")
    print("H_omega' derived / printed:", rec["H_omega_prime_over_printed_eq_hop"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

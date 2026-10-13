"""X2 S5: derive the radial (l = 0) perturbation equations of SV with phantom scalar + magnetic NED.

Framework of Cremona, Pizzocchero & Sarbach (arXiv:1911.13103, CPS): metric
  g = -alpha^2 dt^2 + gamma^2 dx^2 + r^2 dOmega^2   (beta = 0),
field equations G^mu_nu = kappa T^mu_nu (kappa = 1 here; Phi is rescaled accordingly) with
  T[Phi]_{mu nu} = -d_mu Phi d_nu Phi + g_{mu nu} (1/2 (dPhi)^2 - V(Phi))   (phantom; CPS action)
  T[NED]^mu_nu   = diag(-U, -U, -U - r U_r / 2, same)                             (magnetic NED)
The magnetic field F = q sin(theta) dtheta ^ dphi is fixed by the Bianchi identity in the radial
sector, so the NED stress is a FIXED function U of the (perturbed) areal radius r: its (t,x)
block is isotropic for a purely magnetic field; the angular part follows from conservation.

Theory constants (fixed under perturbation): the functions V(Phi) and U(r) (i.e. L(F) and q).
Perturbed fields: alpha, gamma, r, Phi. Gauge: delta Phi = 0 (possible because Phi' != 0 everywhere,
also at the throat), beta = 0. delta r is NOT fixed, so the throat radius may move.

Writes results/x2/s5_radial_equations.json (coefficients as sympy strings) and checks:
  * the SV background with this matter: Phi'^2 > 0 (phantom), V(x), U(r) (Bronnikov & Walia forms);
  * U = V = 0, Ellis-Bronnikov: the equations reproduce CPS Eq. (MasterEllis),
    chi_tt - chi_xx - 3 b^2 chi / r^4 = 0.

  python scripts/x2_s5_radial_derive.py   (about 15 minutes; symbolic)
"""

from __future__ import annotations

from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x2"

t, x, th, ph, eps = sp.symbols("t x theta phi epsilon", real=True)
M, a = sp.symbols("M a", positive=True)


def einstein_mixed(g, coords):
    gi = g.inv()
    n = 4
    Gam = [
        [
            [
                sum(
                    gi[i, m]
                    * (
                        sp.diff(g[m, j], coords[k])
                        + sp.diff(g[m, k], coords[j])
                        - sp.diff(g[j, k], coords[m])
                    )
                    for m in range(n)
                )
                / 2
                for k in range(n)
            ]
            for j in range(n)
        ]
        for i in range(n)
    ]

    def ric(j, k):
        return sum(
            sp.diff(Gam[i][j][k], coords[i])
            - sp.diff(Gam[i][j][i], coords[k])
            + sum(Gam[i][i][p] * Gam[p][j][k] - Gam[i][k][p] * Gam[p][j][i] for p in range(n))
            for i in range(n)
        )

    R = sp.Matrix(n, n, lambda j, k: ric(j, k))
    Rs = sum(gi[i, j] * R[i, j] for i in range(n) for j in range(n))
    return gi * R - sp.eye(n) * Rs / 2


def field_equations(al, ga, rr, Phi, Vfun, Ufun):
    """E^mu_nu = G^mu_nu - T^mu_nu for the diagonal metric (functions of t, x)."""
    coords = [t, x, th, ph]
    g = sp.diag(-(al**2), ga**2, rr**2, rr**2 * sp.sin(th) ** 2)
    G = einstein_mixed(g, coords)
    gi = g.inv()
    dP = [sp.diff(Phi, c) for c in coords]
    kin = sum(gi[i, i] * dP[i] ** 2 for i in range(4))
    TPhi = sp.Matrix(4, 4, lambda m, n: -sum(gi[m, k] * dP[k] for k in range(4)) * dP[n])
    TPhi += sp.eye(4) * (kin / 2 - Vfun)
    Ur, rs = Ufun
    U = Ur.subs(rs, rr)
    dU = sp.diff(Ur, rs).subs(rs, rr)
    TN = sp.diag(-U, -U, -U - rr * dU / 2, -U - rr * dU / 2)
    E = G - TPhi - TN
    # Klein-Gordon for the CPS action: box Phi + V'(Phi) = 0 (V' supplied separately by the caller)
    sqrtg = al * ga * rr**2 * sp.sin(th)
    box = sum(sp.diff(sqrtg * gi[i, i] * dP[i], coords[i]) for i in range(4)) / sqrtg
    return E, box


def linearise():
    """Generic static background (alpha0, gamma0, r0, Phi0, V(x), U(r)); delta Phi = 0 gauge."""
    al0, ga0, r0, P0, Vx = [sp.Function(n)(x) for n in ("al0", "ga0", "r0", "P0", "Vx")]
    Uf = sp.Function("Uf")
    rs = sp.symbols("rs", positive=True)
    D, A, C = [sp.Function(n)(t, x) for n in ("D", "A", "C")]
    E, box = field_equations(
        al0 * (1 + eps * D), ga0 * (1 + eps * A), r0 * (1 + eps * C), P0, Vx, (Uf(rs), rs)
    )
    KG = box + sp.diff(Vx, x) / sp.diff(P0, x)  # V'(Phi) at the background (delta Phi = 0)
    lin = {}
    for name, expr in (
        ("tt", E[0, 0]),
        ("xx", E[1, 1]),
        ("thth", E[2, 2]),
        ("tx", E[0, 1]),
        ("KG", KG),
    ):
        lin[name] = sp.expand(sp.simplify(sp.diff(expr, eps).subs(eps, 0)))
    Dp = sp.solve(lin["KG"], sp.Derivative(D, x))[0]
    Ctt = sp.solve(lin["xx"].subs(sp.Derivative(D, x), Dp), sp.Derivative(C, (t, 2)))[0]
    th = lin["thth"].subs(sp.Derivative(D, (x, 2)), sp.diff(Dp, x)).subs(sp.Derivative(D, x), Dp)
    Att = sp.solve(th.subs(sp.Derivative(C, (t, 2)), Ctt), sp.Derivative(A, (t, 2)))[0]
    return lin, Ctt, Att, (al0, ga0, r0, P0, Vx, Uf, D, A, C)


def ellis_bronnikov_check(Ctt, Att, fns):
    """U = V = 0, alpha = gamma = 1, r = sqrt(x^2 + b^2): must give CPS Eq. (MasterEllis)."""
    al0, ga0, r0, P0, Vx, Uf, D, A, C = fns
    b = sp.symbols("b", positive=True)
    R = sp.sqrt(x**2 + b**2)

    def eb(e):
        e = e.replace(Uf, lambda arg: sp.Integer(0))
        return sp.simplify(
            e.subs({al0: 1, ga0: 1, r0: R, P0: sp.sqrt(2) * sp.atan(x / b), Vx: 0}).doit()
        )

    Cp = (-sp.diff(R, x) / R) * C + (sp.diff(R, x) / R) * A

    def red(e):
        for _ in range(3):
            e = e.subs(sp.Derivative(C, (x, 2)), sp.diff(Cp, x)).subs(sp.Derivative(C, x), Cp)
        return sp.simplify(e)

    chi = (A - C) / R
    return sp.simplify(red((eb(Att) - eb(Ctt)) / R - sp.diff(chi, x, 2) - 3 * b**2 / R**4 * chi))


def sv_background():
    """SV metric with this matter: Phi'^2, U(r) and V(x) from the field equations."""
    R = sp.sqrt(x**2 + a**2)
    f = 1 - 2 * M / R
    G = einstein_mixed(sp.diag(-f, 1 / f, R**2, R**2 * sp.sin(th) ** 2), [t, x, th, ph])
    G00, G11, G22 = [sp.simplify(G[i, i]) for i in range(3)]
    P2 = sp.simplify((G00 - G11) / f)
    S = sp.simplify(-(G00 - f * P2 / 2))  # V + U
    r = sp.symbols("r", positive=True)
    Ur = sp.simplify((-2 * (G22 - f * P2 / 2 + S) / R).subs(x, sp.sqrt(r**2 - a**2)))
    U = sp.simplify(-sp.integrate(Ur, (r, r, sp.oo)))
    V = sp.simplify(S - U.subs(r, R))
    return sp.factor(P2), sp.factor(U), sp.factor(V)


def coefficients(lin, Ctt, Att, fns):
    al0, ga0, r0, P0, Vx, Uf, D, A, C = fns
    F = {"al0": al0, "ga0": ga0, "r0": r0, "P0": P0, "Vx": Vx}
    U = [sp.Symbol(f"U_{k}") for k in range(3)]

    def to_sym(e):
        e = e.subs({sp.Derivative(Uf(r0), (r0, 2)): U[2]}).subs({sp.Derivative(Uf(r0), r0): U[1]})
        e = e.replace(Uf, lambda arg: U[0])
        for n, fn in F.items():
            for k in (3, 2, 1):
                e = e.subs(sp.Derivative(fn, (x, k)), sp.Symbol(f"{n}_{k}"))
            e = e.subs(fn, sp.Symbol(f"{n}_0"))
        return e

    ders = [
        A,
        C,
        sp.Derivative(A, x),
        sp.Derivative(C, x),
        sp.Derivative(A, (x, 2)),
        sp.Derivative(C, (x, 2)),
    ]
    out = {}
    for nm, e in (("Ctt", Ctt), ("Att", Att), ("H", lin["tt"]), ("sigma_tx", lin["tx"])):
        c = sp.collect(
            sp.expand(e),
            ders + [sp.Derivative(A, t), sp.Derivative(C, t), sp.Derivative(C, t, x)],
            evaluate=False,
        )
        out[nm] = {str(k): str(to_sym(v)) for k, v in c.items()}
    return out


def main() -> int:
    import json

    P2, U, V = sv_background()
    lin, Ctt, Att, fns = linearise()
    eb = ellis_bronnikov_check(Ctt, Att, fns)
    rec = {
        "label": "derived; radial perturbations, GR + phantom + magnetic NED, delta Phi = 0 gauge",
        "sv_background": {"Phi_prime_squared": str(P2), "U_of_r": str(U), "V_of_x": str(V)},
        "ellis_bronnikov_residual": str(eb),
        "coefficients": coefficients(lin, Ctt, Att, fns),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s5_radial_equations.json").write_text(json.dumps(rec, indent=1) + "\n")
    print("Phi'^2 =", P2, " U(r) =", U, " V(x) =", V, " EB residual =", eb)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

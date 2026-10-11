"""Does the Hořava / khronometric "instantaneous" term give a controllable remote response?

Goal round 9, P3 (docs/goal_round9/README.md). Offline, symbolic; no data. Reuses the linearised
Einstein-aether (EA) machinery of scripts/goal_ea_response.py (round 8).

1. Controls [reproduced from Blas & Sibiryakov, arXiv:1110.2195v2, Sec. 2 and App. A]:
   - the decoupled khronon with an external source S^mu u_mu (their Eq. 10) has the propagator
     p_i p_j / (p^2 (alpha p0^2 - B p^2)), B = beta + lambda (their Eq. 13). Its instantaneous
     piece p_i p_j / (alpha p0^2 p^2) is transverse and traceless in position space (their
     Eq. 14) and gives no long-range term for a gradient source S_i = d_i sigma;
   - in EA the transverse aether propagator (delta_ij - p_i p_j/p^2)/((alpha + a0) p0^2 - a0 p^2)
     carries the opposite instantaneous piece, so the sum has none (their App. A).
2. New [derived]: the khronometric theory linearised about flat space with the aether at rest,
   with **ordinary matter coupled to the metric only** (no S^mu term), as EA with the transverse
   covariant aether removed in the action (u_perp = v_perp + h_0perp = 0, Jacobson
   arXiv:1001.4823v3). Its gauge-invariant tidal tensor E_ij = R_0i0j is compared with EA's.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import sympy as sp

_here = Path(__file__).resolve().parent
if "goal_ea_response" in sys.modules:
    ea = sys.modules["goal_ea_response"]
else:
    _spec = importlib.util.spec_from_file_location(
        "goal_ea_response", _here / "goal_ea_response.py"
    )
    ea = importlib.util.module_from_spec(_spec)
    sys.modules["goal_ea_response"] = ea
    _spec.loader.exec_module(ea)


def bs_controls() -> dict:
    """Momentum-space checks of Blas-Sibiryakov Eqs. (13)-(14) and App. A (3D, p = (px, py, pz))."""
    p0, a, b, a0 = sp.symbols("p0 alpha B alpha0", positive=True)
    px, py, pz, sig = sp.symbols("px py pz sigma")
    pv = [px, py, pz]
    p2 = px**2 + py**2 + pz**2
    # 1/(p^2 (a p0^2 - b p^2)) = [1/(a p0^2)] [1/p^2 + b/(a p0^2 - b p^2)]
    split = sp.simplify(
        1 / (p2 * (a * p0**2 - b * p2)) - (1 / (a * p0**2)) * (1 / p2 + b / (a * p0**2 - b * p2))
    )
    inst = sp.Matrix(3, 3, lambda i, j: pv[i] * pv[j] / (a * p0**2 * p2))
    # position space: p_i p_j / p^2 -> -d_i d_j (1/(4 pi r)); traceless and divergence-free
    # away from the source because nabla^2 (1/r) = 0 there
    xs = sp.symbols("x y z", real=True)
    r = sp.sqrt(sum(x**2 for x in xs))
    hess = sp.Matrix(3, 3, lambda i, j: sp.diff(1 / r, xs[i], xs[j]))
    trace = sp.simplify(hess.trace())
    div = [sp.simplify(sum(sp.diff(hess[i, j], xs[i]) for i in range(3))) for j in range(3)]
    # gradient source S_i = p_i sigma: the instantaneous exchange S.inst.S has no 1/p^2
    s = sp.Matrix([px * sig, py * sig, pz * sig])
    grad_exchange = sp.simplify((s.T * inst * s)[0])
    # EA: longitudinal (alpha -> alpha + a0) plus transverse instantaneous pieces cancel
    ap = a + a0
    long_inst = sp.Matrix(3, 3, lambda i, j: pv[i] * pv[j] / (ap * p0**2 * p2))
    trans = sp.Matrix(
        3, 3, lambda i, j: (sp.KroneckerDelta(i, j) - pv[i] * pv[j] / p2) / (ap * p0**2 - a0 * p2)
    )
    # instantaneous part of the transverse propagator: its p^2 -> 0-type piece, from
    # 1/(ap p0^2 - a0 p^2) = 1/(ap p0^2) + a0 p^2 / (ap p0^2 (ap p0^2 - a0 p^2))
    trans_inst = sp.Matrix(3, 3, lambda i, j: -pv[i] * pv[j] / (p2 * ap * p0**2))
    resid = sp.simplify(
        trans
        - trans_inst
        - sp.Matrix(
            3,
            3,
            lambda i, j: (
                sp.KroneckerDelta(i, j) / (ap * p0**2)
                + a0
                * (p2 * sp.KroneckerDelta(i, j) - pv[i] * pv[j])
                / (ap * p0**2 * (ap * p0**2 - a0 * p2))
            ),
        )
    )
    return {
        "partial_fraction_residual": str(split),
        "instantaneous_trace_position_space": str(trace),
        "instantaneous_divergence_position_space": [str(d) for d in div],
        "gradient_source_exchange": str(grad_exchange),
        "transverse_decomposition_residual": str(sp.simplify(resid)),
        "ea_instantaneous_sum": str(sp.simplify(long_inst + trans_inst)),
    }


def khronometric_vs_ea() -> dict:
    """Linear tidal response, matter coupled to the metric only. Khronometric = EA with the
    covariant transverse aether u_perp = v_perp + h_0perp set to zero in the action, so the
    h_0perp equation becomes E[h_0perp] - E[v_perp] (chain rule). Gauge h_0x = h_0y = h_0z = 0
    and unitary khronon gauge v_z = 0."""
    eqs, f, t = ea.field_equations()
    names = ea.NAMES
    e = dict(zip(names, eqs, strict=True))
    w, k = ea.W, ea.K
    fix = {f[names.index(n)]: 0 for n in ("vx", "vy", "vz", "h0z", "h0x", "h0y")}
    keep = [e[n] for n in names if n not in ("vx", "vy", "h0x", "h0y")]
    keep += [e["h0x"] - e["vx"], e["h0y"] - e["vy"]]
    unk = [x for x in f if x not in fix]
    t00, t0x, t0y, t0z, txx, txy, txz, tyy, tyz, tzz = t
    rho, jx, jy, sxx, syy, sxy = sp.symbols("rho jx jy Sxx Syy Sxy")
    src = {
        t00: rho,
        t0x: jx,
        t0y: jy,
        t0z: w * rho / k,
        txz: w * jx / k,
        tyz: w * jy / k,
        tzz: w**2 * rho / k**2,
        txx: sxx,
        tyy: syy,
        txy: sxy,
    }
    sols = sp.solve([q.subs(fix).subs(src) for q in keep], unk, dict=True)
    sol = sols[0]
    val = {x: sol.get(x, 0) for x in f}

    def h(a, b):
        name = {v: n for n, v in ea.INDEX.items()}[tuple(sorted((a, b)))]
        return val[f[names.index(name)]]

    kv = [-w, 0, 0, k]

    def riem(a, b, c, d):
        def dd(p, q):
            return -kv[p] * kv[q]

        return sp.Rational(1, 2) * (
            dd(b, c) * h(a, d) + dd(a, d) * h(b, c) - dd(a, c) * h(b, d) - dd(b, d) * h(a, c)
        )

    pairs = {"xx": (1, 1), "yy": (2, 2), "zz": (3, 3), "xy": (1, 2), "xz": (1, 3), "yz": (2, 3)}
    ek = {n: sp.factor(sp.simplify(riem(0, i, 0, j))) for n, (i, j) in pairs.items()}
    eea = ea.tidal_response()
    out = {"n_solutions": len(sols), "difference_c13_0": {}, "difference_general": {}}
    for n in pairs:
        d = sp.simplify(ek[n] - eea[n])
        out["difference_general"][n] = str(sp.factor(d))
        out["difference_c13_0"][n] = str(sp.factor(sp.simplify(d.subs(ea.C3, -ea.C1))))
    out["khronometric_Exz"] = str(ek["xz"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    out = {
        "provenance": "derived symbolically (no data): Blas-Sibiryakov controls and the linearised "
        "khronometric tidal response with matter coupled to the metric only",
        "bs_controls": bs_controls(),
        "khronometric_vs_ea": khronometric_vs_ea(),
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

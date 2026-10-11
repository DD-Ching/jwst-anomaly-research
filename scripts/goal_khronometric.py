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
from functools import lru_cache
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


@lru_cache(maxsize=1)
def khronometric_tidal() -> dict:
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
    if len(sols) != 1:
        raise RuntimeError(f"expected a unique linear solution, got {len(sols)}")
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
    return {n: sp.factor(sp.simplify(riem(0, i, 0, j))) for n, (i, j) in pairs.items()}


def khronometric_vs_ea() -> dict:
    """Helicity-frame (k = z) difference E_khronometric - E_EA, for every c_i and for c13 = 0."""
    ek = khronometric_tidal()
    eea = ea.tidal_response()
    pairs = ("xx", "yy", "zz", "xy", "xz", "yz")
    out = {
        "n_solutions": 1,  # khronometric_tidal raises otherwise
        "difference_c13_0": {},
        "difference_general": {},
    }
    for n in pairs:
        d = sp.simplify(ek[n] - eea[n])
        out["difference_general"][n] = str(sp.factor(d))
        out["difference_c13_0"][n] = str(sp.factor(sp.simplify(d.subs(ea.C3, -ea.C1))))
    out["khronometric_Exz"] = str(ek["xz"])
    return out


# ----------------------------------------------------------------------------- c13 != 0 (round 9b)
# Owner review of #164: GW170817 bounds |c13| <~ 1e-15 but does not set c13 = 0. Does the linear
# khronometric response to ordinary matter contain an instantaneous part when c13 != 0?

_SRC = sp.symbols("rho jx jy Sxx Syy Sxy")
KX, KY, KZ = sp.symbols("kx ky kz", real=True)
KVEC = (KX, KY, KZ)
K2 = KX**2 + KY**2 + KZ**2
_SIJ = sp.symbols("S11 S12 S13 S22 S23 S33", real=True)


def _smat():
    s11, s12, s13, s22, s23, s33 = _SIJ
    return sp.Matrix([[s11, s12, s13], [s12, s22, s23], [s13, s23, s33]])


def helicity_coefficients(e: dict) -> tuple[dict, dict]:
    """Rotation-covariant kernel from the helicity-frame (k = z) response. With khat = k/|k|,
    P = 1 - khat khat, V = P j, S^T = P_mn S_mn and S^TT = P S P - P S^T / 2:
    E_ij = a khat_i khat_j rho + b khat_i khat_j S^T + c P_ij rho + d P_ij S^T
           + v (khat_i V_j + khat_j V_i) + t S^TT_ij.
    Returns the six coefficients (functions of omega, k) and the frame residuals (0 when the
    ansatz reproduces every frame component)."""
    rho, jx, jy, sxx, syy, sxy = _SRC

    def coef(expr, sym):
        return sp.factor(sp.cancel(sp.diff(expr, sym)))

    a, b = coef(e["zz"], rho), coef(e["zz"], sxx)
    c = sp.factor(coef(e["xx"] + e["yy"], rho) / 2)
    d = sp.factor(coef(e["xx"] + e["yy"], sxx) / 2)
    t, v = coef(e["xy"], sxy), coef(e["xz"], jx)
    sig = sxx + syy
    model = {
        "xx": c * rho + d * sig + t * (sxx - syy) / 2,
        "yy": c * rho + d * sig - t * (sxx - syy) / 2,
        "zz": a * rho + b * sig,
        "xy": t * sxy,
        "xz": v * jx,
        "yz": v * jy,
    }
    res = {n: str(sp.cancel(e[n] - model[n])) for n in model}
    return {"a": a, "b": b, "c": c, "d": d, "v": v, "t": t}, res


def _sources():
    """A general conserved source through its stresses: j_i = k_m S_im / omega,
    rho = k_m k_n S_mn / omega^2 (any conserved T^{mu nu} at omega != 0)."""
    w = ea.W
    s = _smat()
    kv = sp.Matrix(KVEC)
    return s, kv, s * kv / w, (kv.T * s * kv)[0] / w**2


def _coefficients_at(co: dict, ps: dict | None) -> dict:
    """Coefficients at a parameter point; None is the GR limit c_i = eps r_i, eps -> 0 (EA's
    scalar block is 0/0 at c_i = 0)."""
    if ps is None:
        eps = sp.Symbol("eps")
        r = {ea.C1: eps / 10, ea.C2: eps / 7, ea.C3: eps / 20, ea.C4: eps / 5}
        return {n: sp.cancel(x.subs(r)).subs(eps, 0) for n, x in co.items()}
    return {n: sp.cancel(x.subs(ps)) for n, x in co.items()}


def kernel_3d(co: dict, ps: dict | None) -> sp.Matrix:
    """E_ij for a general wave vector (kx, ky, kz) and a general conserved source (`_sources`),
    at a parameter point. Written with k_i and k^2 only (v / k and a, b, c, d, t are even in k)."""
    k = ea.K
    s, kv, j, rho = _sources()
    pm = sp.eye(3) - kv * kv.T / K2
    vv = pm * j
    st = (pm * s).trace()
    stt = pm * s * pm - pm * st / 2
    c = _coefficients_at(co, ps)
    c["v"] = sp.cancel(c["v"] / k)
    c = {n: x.subs(k, sp.sqrt(K2)) for n, x in c.items()}
    return (
        (c["a"] * rho + c["b"] * st) * kv * kv.T / K2
        + (c["c"] * rho + c["d"] * st) * pm
        + c["v"] * (kv * vv.T + vv * kv.T)
        + c["t"] * stt
    )


def nonlocal_denominator_factors(m: sp.Matrix) -> list[str]:
    """Denominator factors free of omega. A bare k^2 = kx^2 + ky^2 + kz^2 is a spatially non-local
    term with no retardation (instantaneous). Factors omega^n (time integrals of a local source)
    and omega^2 - c^2 k^2 (retarded propagators) are causal."""
    found = set()
    for i in range(3):
        for j in range(i, 3):
            den = sp.denom(sp.cancel(sp.together(m[i, j])))
            for f, _ in sp.factor_list(den)[1]:
                if not f.has(ea.W) and f.free_symbols & set(KVEC):
                    found.add(str(f))
    return sorted(found)


# Exact rational test points; the kernel is rational in the c_i, so generic values suffice.
PARAM_SETS = {
    "GR limit": None,
    "c13=0": {
        ea.C1: sp.Rational(1, 10),
        ea.C2: sp.Rational(1, 7),
        ea.C3: -sp.Rational(1, 10),
        ea.C4: sp.Rational(1, 5),
    },
    "c13=3/20": {
        ea.C1: sp.Rational(1, 10),
        ea.C2: sp.Rational(1, 7),
        ea.C3: sp.Rational(1, 20),
        ea.C4: sp.Rational(1, 5),
    },
    "c13=1/1000": {
        ea.C1: sp.Rational(1, 10),
        ea.C2: sp.Rational(1, 7),
        ea.C3: -sp.Rational(1, 10) + sp.Rational(1, 1000),
        ea.C4: sp.Rational(1, 5),
    },
}


def c13_coupling():
    """A = 4 pi G c13^2 / (c14 (1 - c13)^2) and EA's c_V^2, in the c_i."""
    c13, c14 = ea.C1 + ea.C3, ea.C1 + ea.C4
    return 4 * sp.pi * ea.GS * c13**2 / (c14 * (1 - c13) ** 2), ea.literature_speeds()["V"]


def instantaneous_part() -> sp.Matrix:
    """I_ij = -2 A k_i k_j rho / k^2 (general direction, general conserved source)."""
    amp, _ = c13_coupling()
    _, kv, _, rho = _sources()
    return -2 * amp * kv * kv.T * rho / K2


def c13_check() -> dict:
    """1. Helicity ansatz for EA and khronometric; only v differs, by A omega k / (omega^2 -
    c_V^2 k^2): Delta E_ij = A omega (k_i P_jm + k_j P_im) j_m / (omega^2 - c_V^2 k^2) (owner).
    2. Khronometric depends on c1 only through c13 and c14 (Jacobson 1001.4823v3: c1 is redundant
    for a hypersurface-orthogonal aether), a check of the reduction.
    3. Split of Delta E with k.j = omega rho and omega^2/(k^2(omega^2 - v^2 k^2)) = 1/k^2 +
    v^2/(omega^2 - v^2 k^2): retarded(c_V) + I_ij, I_ij = -2 A k_i k_j rho / k^2.
    4. Bare-k^2 factors of the full 3D kernel (GR limit, EA, khronometric) and of khronometric - I.
    5. Pole bookkeeping of khronometric - I for a quadrupole source (S = -omega^2 q) at one k:
    which fronts follow I, and the static limit (equal to EA's)."""
    w, k = ea.W, ea.K
    co_ea, res_ea = helicity_coefficients(ea.tidal_response())
    co_kh, res_kh = helicity_coefficients(khronometric_tidal())
    amp, cv2 = c13_coupling()
    diff = {n: sp.factor(sp.cancel(co_kh[n] - co_ea[n])) for n in co_ea}
    out = {
        "frame_residuals_EA": res_ea,
        "frame_residuals_khronometric": res_kh,
        "khronometric_coefficients": {n: str(x) for n, x in co_kh.items()},
        "coefficient_difference_kh_minus_EA": {n: str(x) for n, x in diff.items()},
        "vector_difference_minus_owner_form": str(
            sp.cancel(diff["v"] - amp * w * k / (w**2 - cv2 * k**2))
        ),
    }
    c13s, c14s = sp.symbols("c13 c14", real=True)
    re = {ea.C3: c13s - ea.C1, ea.C4: c14s - ea.C1}
    out["khronometric_d_dc1_at_fixed_c13_c14"] = {
        n: str(sp.cancel(sp.diff(sp.cancel(x.subs(re)), ea.C1))) for n, x in co_kh.items()
    }
    s, kv, j, rho = _sources()
    pm = sp.eye(3) - kv * kv.T / K2
    owner = amp * w * (kv * (pm * j).T + (pm * j) * kv.T) / (w**2 - cv2 * K2)
    ret = (amp * w * (kv * j.T + j * kv.T) - 2 * amp * cv2 * kv * kv.T * rho) / (w**2 - cv2 * K2)
    inst = instantaneous_part()
    out["split_residual"] = [
        str(sp.cancel(sp.together(owner[i, jj] - ret[i, jj] - inst[i, jj])))
        for i in range(3)
        for jj in range(i, 3)
    ]
    bare = {}
    for pname, ps in PARAM_SETS.items():
        bare[f"EA {pname}"] = nonlocal_denominator_factors(kernel_3d(co_ea, ps))
        if ps is not None:
            mk = kernel_3d(co_kh, ps)
            bare[f"khronometric {pname}"] = nonlocal_denominator_factors(mk)
            bare[f"khronometric - I {pname}"] = nonlocal_denominator_factors(mk - inst.subs(ps))
    out["bare_k2_factors"] = bare
    # 5. one wave vector (|k| = 3) and one quadrupole q; S_mn = -omega^2 q_mn
    num = {KX: 1, KY: 2, KZ: 2}
    sq = {x: -(w**2) * v for x, v in zip(_SIJ, (3, -1, 2, 5, 1, -4), strict=True)}
    poles = {}
    for pname in ("c13=3/20", "c13=1/1000"):
        ps = PARAM_SETS[pname]
        mk, me = kernel_3d(co_kh, ps), kernel_3d(co_ea, ps)
        for idx in ((0, 1), (2, 2)):
            tot = sp.cancel(mk[idx].subs(num).subs(sq))
            ii = sp.cancel(inst[idx].subs(ps).subs(num).subs(sq))
            rest = sp.apart(sp.cancel(tot - ii), w)
            poles[f"{pname} E{idx[0] + 1}{idx[1] + 1}"] = {
                "I": str(ii),
                "khronometric_minus_I_terms": [str(x) for x in sp.Add.make_args(rest)],
                "static_khronometric_minus_EA": str(
                    sp.cancel(tot.subs(w, 0) - sp.cancel(me[idx].subs(num).subs(sq)).subs(w, 0))
                ),
            }
    out["pole_bookkeeping"] = poles
    return out


def d4_inv_r(qmat, x):
    """sum_mn q_mn d_i d_j d_m d_n (1/r) at x (3 x 3 numpy array)."""
    import numpy as np

    xs = sp.symbols("x y z", real=True)
    r = sp.sqrt(sum(v**2 for v in xs))
    out = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            expr = sum(
                qmat[m][n] * sp.diff(1 / r, xs[i], xs[j], xs[m], xs[n])
                for m in range(3)
                for n in range(3)
                if qmat[m][n]
            )
            out[i, j] = float(expr.subs(dict(zip(xs, x, strict=True))))
    return out


def position_space_check(n: int = 128, box: float = 16.0) -> dict:
    """Smooth, local, conserved, non-spherical source with its stresses:
    T^{ij} = q''_ij(t) f(x), T^{0i} = -d_j(q'_ij f), rho = d_i d_j(q_ij f), with
    f = (3465 / (512 pi R^3)) (1 - r^2/R^2)^4 for r < R (int f = 1, C^3) and q switched on smoothly
    from 0 (same prior state, decaying boundary conditions). The retarded parts vanish at the
    receiver before the fastest front, so the receiver sees only I_ij = -2 A d_i d_j del^-2 rho.
    Checks numerically (FFT, A = 1) that outside the source I_ij = (A / 2 pi) q_mn
    d_i d_j d_m d_n (1/r) (shell theorem: del^-2 f = -1/(4 pi r) outside R)."""
    import numpy as np

    big_r = 1.0
    h = box / n
    ax = (np.arange(n) - n // 2) * h
    x, y, z = np.meshgrid(ax, ax, ax, indexing="ij")
    rr = np.sqrt(x**2 + y**2 + z**2)
    f = np.where(rr < big_r, (1 - (rr / big_r) ** 2) ** 4, 0.0) * 3465 / (512 * np.pi * big_r**3)
    q = np.array([[0.3, 1.0, 0.0], [1.0, -0.5, 0.2], [0.0, 0.2, 0.2]])  # non-spherical
    kax = 2 * np.pi * np.fft.fftfreq(n, d=h)
    kk = np.meshgrid(kax, kax, kax, indexing="ij")
    k2 = kk[0] ** 2 + kk[1] ** 2 + kk[2] ** 2
    k2[0, 0, 0] = 1.0
    ft = np.fft.fftn(np.fft.ifftshift(f))
    rho_k = -sum(q[m, nn] * kk[m] * kk[nn] for m in range(3) for nn in range(3)) * ft
    out = {"n": n, "box_over_R": box, "points": []}
    for pt in ((4.0, 0.0, 0.0), (2.4, 2.4, 2.4), (0.0, 3.2, -2.0)):
        idx = tuple(int(round(c / h)) + n // 2 for c in pt)
        num = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                ik = -2 * kk[i] * kk[j] * rho_k / k2
                ik[0, 0, 0] = 0.0
                num[i, j] = np.fft.fftshift(np.fft.ifftn(ik).real)[idx]
        xyz = [ax[i] for i in idx]
        ana = d4_inv_r(q.tolist(), xyz) / (2 * np.pi)
        out["points"].append(
            {
                "x_over_R": [round(float(v), 3) for v in xyz],
                "max_rel_error": float(np.max(np.abs(num - ana)) / np.max(np.abs(ana))),
                "E_xx_numeric": float(num[0, 0]),
                "E_xx_analytic": float(ana[0, 0]),
            }
        )
    return out


def link_budget_c13(p=ea.P) -> dict:
    """Signal size of I_ij (model-internal linear response; ASSUMPTION-laden, generous).
    Outside the sender, I_ij = kappa' G I_mn(t) d_i d_j d_m d_n (1/r), kappa' = c13^2 /
    (c14 (1 - c13)^2), I_mn = int rho x_m x_n (energy and momentum conservation fix the monopole
    and the dipole's time dependence, so the quadrupole leads). Sender: round 8's mass and stroke
    moved non-spherically, Delta I = 2 M R stroke in I_xx - I_yy (ASSUMPTION); receiver on the x
    axis, angular factor 36 (d_x^4 - d_x^2 d_y^2 of 1/r = 36 / r^5). It lasts until the fastest
    front arrives, (r - R) / max(c_S, c_T). SNR = |I_xx| sqrt(tau) / ASD for a receive time tau
    inside that window. Returns kappa' and c14 (at |c13| = 1e-15) needed for the round-8 SNR."""
    import math

    ang = abs(d4_inv_r([[1, 0, 0], [0, -1, 0], [0, 0, 0]], (1.0, 0.0, 0.0))[0, 0])
    di = 2 * p.source_mass_kg * p.source_radius_m * p.source_stroke_m
    rows = []
    for d_ly in p.distances_ly:
        r = d_ly * ea.LY_M
        e_unit = ea.G_SI * di * ang / r**5  # I_xx / kappa'
        for tau in p.latencies_s:
            k_req = p.snr * p.detector_asd_s2_per_rthz / (e_unit * math.sqrt(tau))
            rows.append(
                {
                    "distance_ly": d_ly,
                    "receive_time_s": tau,
                    "window_needed_c_max_at_most": r / (ea.C_SI * tau),
                    "I_xx_per_kappa_s2": e_unit,
                    "kappa_required": k_req,
                    "c14_required_at_c13_1e-15": 1e-30 / k_req,
                }
            )
    return {"angular_factor": ang, "delta_I_kg_m2": di, "rows": rows}


# ------------------------------------------------------------------ PPN x receive window (round 9c)
# Owner review 3: can tuning (alpha -> 0, lambda + beta -> 0) give the round-8 receiver the kappa'
# it needs before the fastest front arrives, given the PPN bounds? Khronometric formulas only.
# alpha = c14, beta = c13, lambda = c2; kappa' = beta^2 / (alpha (1 - beta)^2).

# Weak-field PPN bounds as quoted by Yagi et al. 1311.7144v4 Sec. II (citing Will's review):
ALPHA1_MAX = 1e-4
ALPHA2_MAX = 1e-7
BETA_MAX = 1e-15  # GW170817, |c13| (round 8)
# Solar-system background for the perturbativity criterion (ASSUMPTIONs): Earth's orbital speed is
# a floor on its speed relative to any frame over a year; the CMB dipole speed 369.8 km/s is the
# usual guess for the aether frame; Phi = G M_sun / (1 AU c^2).
V_FLOOR = 2.98e4 / ea.C_SI
V_CMB = 3.698e5 / ea.C_SI
PHI_1AU = 9.87e-9


def khronometric_ppn():
    """Yagi et al. 1311.7144v4 Eqs. (48), (49), (118): alpha_1, alpha_2, c_s^2 (khronometric;
    exact in the couplings, weak-field 1PN)."""
    a, b, lam = sp.symbols("alpha beta lambda", real=True)
    a1 = 4 * (a - 2 * b) / (b - 1)
    a2 = (
        (a - 2 * b)
        * (-(b**2) + b * (a - 3) + a + lam * (-1 - 3 * b + 2 * a))
        / ((b - 1) * (lam + b) * (a - 2))
    )
    cs2 = (a - 2) * (b + lam) / (a * (b - 1) * (2 + b + 3 * lam))
    return (a, b, lam), a1, a2, cs2


def ppn_identity() -> dict:
    """General symbolic check (not parameter points): kappa' = (c_s^2 / 2)(alpha_2 - alpha_1 / 8)
    + (beta - alpha / 4) / (1 - beta)^2, and the factorised pieces used in the write-up."""
    (a, b, lam), a1, a2, cs2 = khronometric_ppn()
    kap = b**2 / (a * (1 - b) ** 2)
    rhs = cs2 / 2 * (a2 - a1 / 8) + (b - a / 4) / (1 - b) ** 2
    return {
        "identity_residual": str(sp.cancel(kap - rhs)),
        "alpha2_minus_alpha1_over_8": str(sp.factor(sp.cancel(a2 - a1 / 8))),
        "kappa_minus_second_term": str(sp.factor(sp.cancel(kap - (b - a / 4) / (1 - b) ** 2))),
        "singular_or_degenerate": [
            "alpha = 0 (kappa', c_s^2 diverge)",
            "lambda + beta = 0 (c_s^2 = 0; alpha_2 diverges unless alpha = 2 beta)",
            "alpha = 2 (c_s^2 = 0; alpha_2 diverges)",
            "beta = 1",
            "2 + beta + 3 lambda = 0 (c_s^2 diverges)",
            "c_s^2 < 0 (gradient instability; no bound derived)",
            "alpha = 2 beta (alpha_1 = alpha_2 = 0, kappa' = beta / (2 (1 - beta)^2))",
        ],
    }


def ppn_window_check(p=ea.P) -> dict:
    """For each round-8 distance and receive time tau:
    - the window X = D / (c tau) >= max(c_S, c_T) >= c_S;
    - the needed kappa' (link_budget_c13);
    - the PPN-allowed kappa' <= (X^2 / 2)(alpha2_max + alpha1_max / 8) + beta_max / (1 - beta_max)^2
      (alpha > 0, c_s^2 >= 0, either sign of beta; from ppn_identity);
    - the alpha_2 the needed kappa' implies through the same identity, and the perturbativity
      measure alpha_2 v^2 Phi at 1 AU (the PPN metric's preferred-frame term)."""
    import math

    coef = (ALPHA2_MAX + ALPHA1_MAX / 8) / 2
    second = BETA_MAX / (1 - BETA_MAX) ** 2
    lb = link_budget_c13(p)
    rows = []
    for r0 in lb["rows"]:
        x = r0["window_needed_c_max_at_most"]
        k_req = r0["kappa_required"]
        k_ppn = coef * x**2 + second
        a2_need = 2 * (k_req - second) / x**2 - ALPHA1_MAX / 8
        rows.append(
            {
                "distance_ly": r0["distance_ly"],
                "receive_time_s": r0["receive_time_s"],
                "X_max_cS": x,
                "kappa_required": k_req,
                "kappa_allowed_PPN": k_ppn,
                "log10_shortfall": math.log10(k_req / k_ppn),
                "alpha2_implied_by_required": a2_need,
                "alpha2_v2_phi_floor": a2_need * V_FLOOR**2 * PHI_1AU,
                "alpha2_v2_phi_cmb": a2_need * V_CMB**2 * PHI_1AU,
                # perturbativity alone (alpha_2 v^2 Phi < 1 at 1 AU, v = floor); no observed alpha_2
                "kappa_cap_perturbative_only": x**2 / 2 / (V_FLOOR**2 * PHI_1AU) + second,
            }
        )
    # receive time below which the PPN cap would reach the requirement (formal; flat-ASD assumption)
    thr = {}
    for d_ly in p.distances_ly:
        r = d_ly * ea.LY_M
        e_unit = ea.G_SI * lb["delta_I_kg_m2"] * lb["angular_factor"] / r**5
        # coef (r / c tau)^2 = snr asd / (e_unit sqrt(tau))
        #   ->  tau^(3/2) = coef r^2 e_unit / (c^2 snr asd)
        t32 = coef * r**2 * e_unit / (ea.C_SI**2 * p.snr * p.detector_asd_s2_per_rthz)
        thr[str(d_ly)] = float(t32 ** (2 / 3))
    return {
        "coef_cS2": coef,
        "beta_term": second,
        "rows": rows,
        "tau_threshold_s": thr,
        "assumptions": {
            "alpha1_max": ALPHA1_MAX,
            "alpha2_max": ALPHA2_MAX,
            "beta_max": BETA_MAX,
            "v_floor": V_FLOOR,
            "v_cmb": V_CMB,
            "phi_1au": PHI_1AU,
        },
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    out = {
        "provenance": "derived symbolically (no data): Blas-Sibiryakov controls and the linearised "
        "khronometric tidal response with matter coupled to the metric only",
        "bs_controls": bs_controls(),
        "khronometric_vs_ea": khronometric_vs_ea(),
        "c13_nonzero": c13_check(),
        "position_space_check": position_space_check(),
        "link_budget_c13": link_budget_c13(),
        "ppn_identity": ppn_identity(),
        "ppn_window_check": ppn_window_check(),
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

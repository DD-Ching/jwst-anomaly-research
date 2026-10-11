"""Einstein-aether (EA) spin-0 mode as a 1-bit channel: retarded tidal response to a matter source.

Goal round 8 (docs/goal_round8/README.md). Offline, symbolic; no data.

Linearised EA about flat space with the aether at rest (u = d_t), signature (-+++):

    S = (1/16 pi G) int sqrt(-g) [R - c1 (Du)^2 - c2 (div u)^2 - c3 Du.Du^T + c4 a^2] + S_matter,

with matter coupled to the metric only. Plane wave exp(i(k z - omega t)), aether constraint
imposed at linear order, gauge v^i = 0 = h_0z, a general conserved source T^{mu nu}. The receiver
reads the gauge-invariant tidal tensor E_ij = R_{0i0j} (geodesic deviation of a freely falling
detector at rest in the aether frame).

Positive controls [reproduced]:
- poles: c_T^2 = 1/(1 - c13), c_V^2 = (2c1 - c1^2 + c3^2)/(2 c14 (1 - c13)),
  c_S^2 = c123 (2 - c14)/(c14 (1 - c13)(2 + c13 + 3 c2)) (Jacobson & Mattingly gr-qc/0402005);
- static E_zz = 4 pi G_N rho with G_N = G/(1 - c14/2) (Carroll & Lim; Foster);
- GR tensor far field E_xy = -G omega^4 I_xy / r for c13 = 0.

New [derived] (linear order, c13 = 0, leading order in source size / wavelength, c = 1):
- a spherically pulsating source (Q = int rho r^2) has an exterior field made only of the c_S pole,
  with no instantaneous part:
  E_ij = -kappa G [delta_ij Q''''(t_S)/r + d_i d_j (Q''(t_S)/r)], t_S = t - r/c_S,
  kappa = c14 (c14 + 3 c2) / (3 c2 (2 - c14)^2)
  (SI: -kappa (G/c^2) [delta_ij Q''''/(c^2 r) + d_i d_j (Q''/r)]);
- any slow source (closed form found by the independent review, checked here per component):
  E^S_ij = -[G c14 / (2 (2 - c14))](delta_ij d_t^2 + d_i d_j)(d_t^2 I_kk + d_k d_l I_kl) acting on
  [ . ]_{t - r/c_S} / r; it reduces to the spherical result for I = delta Q / 3;
- far field -(c14/4) G omega^4 I_kk delta_ij / r + O(c14^2), equal to Yagi et al. Eq. (106);
- spin-1 with c13 != 0: kappa_V = c13^2 / (2 c14 c_V^4 (1 - c13)^2) <= 1 / ((1 - c13) c_V^2);
- a link budget: c_S = D / latency fixes c14 <= 2 c2 / ((2 + 3 c2) c_S^2), so kappa ~ 1/c_S^2.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache

import sympy as sp

G_SI = 6.67430e-11
C_SI = 2.99792458e8
LY_M = 9.4607304725808e15


@dataclass(frozen=True)
class Params:
    """Link-budget choices (ASSUMPTIONs, deliberately generous to the channel)."""

    c2_max: float = 0.095  # upper end of the post-GW170817 allowed c2 range (1802.04303)
    # 1802.04303 Eq. (3.23): region (i) 0 < c14 <= 2e-7 with c14 <~ c2 <~ 0.095 (c_S can be >> 1);
    # region (ii) 2e-6 <~ c14 <~ 2.5e-5 needs c2 - c14 <~ 2e-7, i.e. c_S ~ 1
    c14_max_region_i: float = 2e-7
    distances_ly: tuple = (20.0, 100.0)
    latencies_s: tuple = (1.0, 86400.0, 3.15576e7)
    # sender: a 1e12 kg shell (a large mountain) at R = 1 km pulsating radially by +-10 m at 1 Hz
    source_mass_kg: float = 1e12
    source_radius_m: float = 1e3
    source_stroke_m: float = 10.0
    source_freq_hz: float = 1.0
    # receiver: tidal-noise ASD, ~LIGO's best band expressed as omega^2 h (generous at 1 Hz)
    detector_asd_s2_per_rthz: float = 1e-20
    snr: float = 5.0


P = Params()

W, K = sp.symbols("omega k", positive=True)
C1, C2, C3, C4, GS = sp.symbols("c1 c2 c3 c4 G", real=True)
NAMES = ["h00", "h0x", "h0y", "h0z", "hxx", "hxy", "hxz", "hyy", "hyz", "hzz", "vx", "vy", "vz"]
INDEX = {
    "h00": (0, 0),
    "h0x": (0, 1),
    "h0y": (0, 2),
    "h0z": (0, 3),
    "hxx": (1, 1),
    "hxy": (1, 2),
    "hxz": (1, 3),
    "hyy": (2, 2),
    "hyz": (2, 3),
    "hzz": (3, 3),
}


def _hmat(s):
    h00, h0x, h0y, h0z, hxx, hxy, hxz, hyy, hyz, hzz = s[:10]
    return sp.Matrix(
        [[h00, h0x, h0y, h0z], [h0x, hxx, hxy, hxz], [h0y, hxy, hyy, hyz], [h0z, hxz, hyz, hzz]]
    )


@lru_cache(maxsize=1)
def field_equations():
    """Euler-Lagrange equations of the quadratic action in complex-amplitude form (a Hermitian
    density conj(F) M F); the matter term carries the calibrated factor 1/4 (checked against
    Phi = -G_N M / r and the GR quadrupole field)."""
    eta = sp.diag(-1, 1, 1, 1)
    f = sp.symbols(NAMES)
    fc = sp.symbols([n + "c" for n in NAMES])
    kv = [-W, 0, 0, K]  # d_a -> i k_a on the field, -i k_a on its conjugate
    h, hc = _hmat(f), _hmat(fc)
    hu, hcu = eta * h * eta, eta * hc * eta
    tr = lambda m: sum(eta[a, a] * m[a, a] for a in range(4))  # noqa: E731
    kk = sum(eta[a, a] * kv[a] ** 2 for a in range(4))
    div = [sum(sp.I * kv[a] * hu[a, b] for a in range(4)) for b in range(4)]
    divc = [sum(-sp.I * kv[a] * hcu[a, b] for a in range(4)) for b in range(4)]
    # Fierz-Pauli density for sqrt(-g) R at second order
    leh = -sp.Rational(1, 4) * kk * sum(hc[a, b] * hu[a, b] for a in range(4) for b in range(4))
    leh += sp.Rational(1, 2) * sum(divc[b] * eta[b, b] * div[b] for b in range(4))
    leh += -sp.Rational(1, 4) * (
        sum(divc[b] * sp.I * kv[b] * tr(h) for b in range(4))
        + sum(div[b] * (-sp.I) * kv[b] * tr(hc) for b in range(4))
    )
    leh += sp.Rational(1, 4) * kk * tr(hc) * tr(h)

    def dmat(s, sgn):
        hm = _hmat(s)
        v = s[10:]
        du = [hm[0, 0] / 2, v[0], v[1], v[2]]  # u^0 = 1 + h00/2 from u.u = -1
        kvv = kv if sgn > 0 else [-x for x in kv]

        def gam(m, a, b):
            return (
                sp.Rational(1, 2)
                * eta[m, m]
                * (sp.I * kvv[a] * hm[m, b] + sp.I * kvv[b] * hm[m, a] - sp.I * kvv[m] * hm[a, b])
            )

        return [[sp.I * kvv[a] * du[m] + gam(m, a, 0) for m in range(4)] for a in range(4)]

    d, dc = dmat(f, 1), dmat(fc, -1)
    lu = C1 * sum(dc[a][m] * eta[a, a] * eta[m, m] * d[a][m] for a in range(4) for m in range(4))
    lu += C2 * sum(dc[a][a] for a in range(4)) * sum(d[a][a] for a in range(4))
    lu += C3 * sum(dc[a][m] * d[m][a] for a in range(4) for m in range(4))
    lu += -C4 * sum(dc[0][m] * eta[m, m] * d[0][m] for m in range(4))
    lag = sp.expand((leh - lu) / (16 * sp.pi * GS))
    t = sp.symbols("T00 T0x T0y T0z Txx Txy Txz Tyy Tyz Tzz")
    tm = _hmat(list(t) + [0, 0, 0])
    eqs = []
    for n, xc in zip(NAMES, fc, strict=True):
        e = sp.diff(lag, xc)
        if n in INDEX:
            a, b = INDEX[n]
            e += sp.Rational(1, 4) * (1 if a == b else 2) * tm[a, b]
        eqs.append(sp.expand(e))
    return eqs, f, t


@lru_cache(maxsize=1)
def tidal_response():
    """E_ij(omega, k) for a general conserved source: rho, j_x, j_y, transverse stresses S_xx,
    S_yy, S_xy (the longitudinal components follow from conservation)."""
    eqs, f, t = field_equations()
    fix = {f[NAMES.index(n)]: 0 for n in ("vx", "vy", "vz", "h0z")}
    unk = [x for x in f if x not in fix]
    t00, t0x, t0y, t0z, txx, txy, txz, tyy, tyz, tzz = t
    rho, jx, jy, sxx, syy, sxy = sp.symbols("rho jx jy Sxx Syy Sxy")
    src = {
        t00: rho,
        t0x: jx,
        t0y: jy,
        t0z: W * rho / K,
        txz: W * jx / K,
        tyz: W * jy / K,
        tzz: W**2 * rho / K**2,
        txx: sxx,
        tyy: syy,
        txy: sxy,
    }
    sol = sp.solve([e.subs(fix).subs(src) for e in eqs], unk, dict=True)[0]
    val = {x: sol.get(x, 0) for x in f}

    def h(a, b):
        name = {v: n for n, v in INDEX.items()}[tuple(sorted((a, b)))]
        return val[f[NAMES.index(name)]]

    kv = [-W, 0, 0, K]

    def riem(a, b, c, d):
        dd = lambda p, q: -kv[p] * kv[q]  # noqa: E731
        return sp.Rational(1, 2) * (
            dd(b, c) * h(a, d) + dd(a, d) * h(b, c) - dd(a, c) * h(b, d) - dd(b, d) * h(a, c)
        )

    pairs = {"xx": (1, 1), "yy": (2, 2), "zz": (3, 3), "xy": (1, 2), "xz": (1, 3), "yz": (2, 3)}
    return {n: sp.factor(sp.simplify(riem(0, i, 0, j))) for n, (i, j) in pairs.items()}


def literature_speeds():
    c13, c14, c123 = C1 + C3, C1 + C4, C1 + C2 + C3
    return {
        "T": 1 / (1 - c13),
        "V": (2 * C1 - C1**2 + C3**2) / (2 * c14 * (1 - c13)),
        "S": c123 * (2 - c14) / (c14 * (1 - c13) * (2 + c13 + 3 * C2)),
    }


def controls() -> dict:
    """Pole positions against the published speeds, Newtonian limit and the GR tensor field."""
    e = tidal_response()
    rho, jx, jy, sxx, syy, sxy = sp.symbols("rho jx jy Sxx Syy Sxy")
    x = sp.symbols("x")
    lit = literature_speeds()
    poles = {}
    for name, comp, src in (
        ("S", "zz", {jx: 0, jy: 0, sxx: 0, syy: 0, sxy: 0}),
        ("V", "xz", {rho: 0, jy: 0, sxx: 0, syy: 0, sxy: 0}),
        ("T", "xy", {rho: 0, jx: 0, jy: 0, sxx: 0, syy: 0}),
    ):
        den = sp.fraction(sp.factor(e[comp].subs(src)))[1]
        roots = [r for r in sp.solve(den.subs(W, sp.sqrt(x) * K), x) if r != 0]
        poles[name] = [str(sp.simplify(r - lit[name])) for r in roots]
    static = sp.simplify(e["zz"].subs({jx: 0, jy: 0, sxx: 0, syy: 0, sxy: 0}).subs(W, 0))
    newton = sp.simplify(static - 4 * sp.pi * GS * rho / (1 - (C1 + C4) / 2))
    # GR quadrupole: E_xy -> 2 G omega^2 S_xy / r with S_xy = -omega^2 I_xy / 2, for c13 = 0
    exy = sp.simplify(e["xy"].subs(C3, -C1))
    tensor = sp.simplify(exy * (K**2 - W**2) / (8 * sp.pi * GS * W**2 * sxy))
    return {
        "pole_minus_literature_speed2": poles,
        "static_Ezz_minus_4piGN_rho": str(newton),
        "tensor_kernel_over_GR": str(tensor),
    }


def _c13_zero():
    c14 = sp.symbols("c14", positive=True)
    return c14, {C3: -C1, C4: c14 - C1}


def spherical_source() -> dict:
    """Pulsating spherical source, slow-source moments rho = -k^2 Q/6, S_xx = S_yy = -omega^2 Q/6.
    Splits E_ij(k) = E_xx delta_ij + (E_zz - E_xx) k_i k_j / k^2 into partial fractions in omega^2:
    the exterior field is the c_S-pole part only if no 1/k^2 (instantaneous) term remains."""
    e = tidal_response()
    rho, jx, jy, sxx, syy, sxy = sp.symbols("rho jx jy Sxx Syy Sxy")
    c14, sub = _c13_zero()
    q, w2 = sp.symbols("Q W")
    src = {rho: -(K**2) * q / 6, sxx: -(W**2) * q / 6, syy: -(W**2) * q / 6, jx: 0, jy: 0, sxy: 0}
    ezz = sp.factor(sp.simplify(e["zz"].subs(src).subs(sub).subs(W, sp.sqrt(w2))))
    exx = sp.factor(sp.simplify(e["xx"].subs(src).subs(sub).subs(W, sp.sqrt(w2))))
    cs2 = sp.factor(sp.simplify(literature_speeds()["S"].subs(sub)))
    diff = sp.factor(sp.simplify((ezz - exx) / K**2))
    pole_exx = sp.factor(sp.simplify(sp.limit(exx * (w2 - cs2 * K**2), w2, cs2 * K**2)))
    pole_diff = sp.factor(sp.simplify(sp.limit(diff * (w2 - cs2 * K**2), w2, cs2 * K**2)))
    alpha = sp.factor(sp.simplify(pole_diff / (q * K**2)))
    # E_ij = -(alpha/cS^4)(delta_ij dt^4 + d_i d_j dt^2) Q_ret/(4 pi r); kappa = alpha/(4 pi G cS^4)
    kappa = sp.factor(sp.simplify(alpha / (4 * sp.pi * GS * cs2**2)))
    # remainders after removing the pole must be polynomial in omega^2 (local contact terms)
    rem_exx = sp.simplify(exx - pole_exx / (w2 - cs2 * K**2))
    rem_diff = sp.simplify(diff - pole_diff / (w2 - cs2 * K**2))
    return {
        "c_S^2": str(cs2),
        "pole_Exx_over_alpha_cS2_k4_Q": str(sp.simplify(pole_exx / (alpha * cs2 * K**4 * q))),
        "kappa": str(kappa),
        "kappa_leading_c14_to_0": str(sp.series(kappa, c14, 0, 2).removeO()),
        "remainder_Exx": str(sp.factor(rem_exx)),
        "remainder_diff": str(sp.factor(rem_diff)),
        # an instantaneous (action-at-a-distance) piece would leave k in a denominator
        "no_instantaneous_term": not sp.denom(sp.together(rem_exx)).has(K)
        and not sp.denom(sp.together(rem_diff)).has(K),
    }


def far_field_general() -> dict:
    """c_S-pole far-field amplitude E^S_ij r / (G omega^4) for a slow quadrupole source, n = z."""
    e = tidal_response()
    rho, jx, jy, sxx, syy, sxy = sp.symbols("rho jx jy Sxx Syy Sxy")
    c14, sub = _c13_zero()
    ixx, iyy, izz, w2 = sp.symbols("Ixx Iyy Izz W")
    cs2 = sp.factor(sp.simplify(literature_speeds()["S"].subs(sub)))
    out = {}
    for comp in ("xx", "zz"):
        ex = sp.factor(
            sp.simplify(e[comp].subs({jx: 0, jy: 0, sxy: 0}).subs(sub).subs(W, sp.sqrt(w2)))
        )
        a = sp.factor(sp.simplify(sp.limit(ex * (w2 - cs2 * K**2), w2, cs2 * K**2)))
        onshell = {rho: -(K**2) * izz / 2, sxx: -w2 * ixx / 2, syy: -w2 * iyy / 2}
        amp = sp.simplify((-(a / cs2) / (4 * sp.pi)).subs(onshell).subs(K, sp.sqrt(w2 / cs2)))
        amp = sp.factor(sp.simplify(amp / (GS * w2**2)))
        out[comp] = {
            "exact": str(amp),
            "leading_c14_to_0": str(sp.simplify(sp.series(amp, c14, 0, 2).removeO())),
        }
    return out


def yagi_comparison() -> dict:
    """Published control: Yagi, Blas, Barausse & Yunes (1311.7144v4) Eq. (106) for weakly
    gravitating sources (sensitivity terms zero), breathing polarisation h_b = F/2 (2304.06801),
    transverse tidal E_xx = omega^2 h_b / 2 up to an overall sign convention:
    |E_xx| r / (G omega^4) = c14 [(3/2)(Z - 1) Q_nn + Z I_kk / 2] / (2 - c14), with Z from their
    Eqs. (42), (43), (107). Compares the I_kk and Q_nn = I_nn - I_kk/3 coefficients with ours."""
    c14, sub = _c13_zero()
    a1 = -8 * (C3**2 + C1 * C4) / (2 * C1 - C1**2 + C3**2)
    a2 = a1 / 2 - (C1 + 2 * C3 - C4) * (2 * C1 + 3 * C2 + C3 + C4) / (
        (2 - (C1 + C4)) * (C1 + C2 + C3)
    )
    cp = C1 + C3
    z = sp.factor(sp.simplify((((a1 - 2 * a2) * (1 - cp)) / (3 * (2 * cp - (C1 + C4)))).subs(sub)))
    yagi_ikk = sp.simplify(c14 * z / 2 / (2 - c14))
    yagi_qnn = sp.simplify(c14 * sp.Rational(3, 2) * (z - 1) / (2 - c14))
    ixx, iyy, izz = sp.symbols("Ixx Iyy Izz")
    ours = sp.sympify(far_field_general()["xx"]["exact"], locals={"c14": c14, "c2": C2})
    ikk, qnn = sp.symbols("Ikk Qnn")
    ours = sp.expand(-ours.subs({ixx: 0, iyy: ikk - qnn - ikk / 3, izz: qnn + ikk / 3}))
    ours_ikk = sp.simplify(ours.coeff(ikk))
    ours_qnn = sp.simplify(ours.coeff(qnn))
    return {
        "Z_c13_0": str(z),
        "Ikk_ours_minus_yagi": str(sp.simplify(ours_ikk - yagi_ikk)),
        "Qnn_ours_over_yagi": str(sp.simplify(ours_qnn / yagi_qnn)),
    }


def general_source_kernel() -> dict:
    """Any slow source, c13 = 0 (closed form found in the independent review, checked here):
    E^S_ij = -[G c14 / (2 (2 - c14))] (delta_ij d_t^2 + d_i d_j)(d_t^2 I_kk + d_k d_l I_kl)
    evaluated on [ . ]_{t - r/c_S} / r. In Fourier space its residue at omega^2 = c_S^2 k^2 is
    4 pi C c_S^2 (W delta_ij + k_i k_j)(W I_kk + k_k k_l I_kl), C = G c14 / (2 (2 - c14)).
    Compares that with the c_S-pole residue of the full response in the helicity frame (k = z)."""
    e = tidal_response()
    rho, jx, jy, sxx, syy, sxy = sp.symbols("rho jx jy Sxx Syy Sxy")
    c14, sub = _c13_zero()
    w2 = sp.symbols("W", positive=True)
    ixx, iyy, izz, ixy, ixz, iyz = sp.symbols("Ixx Iyy Izz Ixy Ixz Iyz")
    cs2 = sp.factor(sp.simplify(literature_speeds()["S"].subs(sub)))
    src = {
        rho: -(K**2) * izz / 2,
        jx: -sp.sqrt(w2) * K * ixz / 2,
        jy: -sp.sqrt(w2) * K * iyz / 2,
        sxx: -w2 * ixx / 2,
        syy: -w2 * iyy / 2,
        sxy: -w2 * ixy / 2,
    }
    cc = GS * c14 / (2 * (2 - c14))
    ikk = ixx + iyy + izz
    expect = {"xx": w2, "zz": w2 + K**2, "xy": 0, "xz": 0}
    out = {}
    for comp, tens in expect.items():
        ex = sp.factor(sp.simplify(e[comp].subs(sub).subs(W, sp.sqrt(w2)).subs(src)))
        res = sp.simplify(sp.limit(ex * (w2 - cs2 * K**2), w2, cs2 * K**2))
        pred = (4 * sp.pi * cc * cs2 * tens * (w2 * ikk + K**2 * izz)).subs(w2, cs2 * K**2)
        out[comp] = str(sp.simplify(res - pred))
    return out


def vector_sector_c13() -> dict:
    """Spin-1 mode with c13 != 0 (linear theory): far-field E_xz r / (G omega^4 I_xz) at the c_V
    pole, written with c_V, and a grid check of kappa_V <= 1 / ((1 - c13) c_V^2). The bound
    follows from 2 c14 c_V^2 (1 - c13) = 2 c1 (1 - c13) + c13^2 with c1 >= 0."""
    e = tidal_response()
    rho, jx, jy, sxx, syy, sxy = sp.symbols("rho jx jy Sxx Syy Sxy")
    c13, c14, cv = sp.symbols("c13 c14 c_V", positive=True)
    ixz = sp.symbols("Ixz")
    ex = e["xz"].subs({rho: 0, jy: 0, sxx: 0, syy: 0, sxy: 0}).subs(jx, -W * K * ixz / 2)
    ex = ex.subs({C3: c13 - C1, C4: c14 - C1})
    cv2 = sp.simplify(literature_speeds()["V"].subs({C3: c13 - C1, C4: c14 - C1}))
    res = sp.factor(sp.simplify(sp.limit(sp.factor(ex) * (W**2 - cv2 * K**2), W, sp.sqrt(cv2) * K)))
    amp = sp.simplify((-(res / cv2) / (4 * sp.pi)).subs(K, W / sp.sqrt(cv2)) / (GS * W**4 * ixz))
    c1sol = sp.solve(sp.Eq(cv2, cv**2), C1)[0]
    amp_cv = sp.factor(sp.simplify(amp.subs(C1, c1sol)))
    worst = 0.0
    for c13v in (1e-15, 1e-6, 1e-2):
        for c14v in (1e-25, 1e-12, 2e-7):
            for c1v in (c14v, 1e-3, 0.1):
                cvv2 = float(cv2.subs({c13: c13v, c14: c14v, C1: c1v}))
                a = abs(float(amp.subs({c13: c13v, c14: c14v, C1: c1v})))
                worst = max(worst, a * cvv2 * (1 - c13v))
    return {
        "far_field_Exz_r_over_G_w4_Ixz": str(amp_cv),
        "max_kappaV_times_cV2_(1-c13)_on_grid": worst,
    }


def link_budget(p: Params = P) -> list:
    """Largest spin-0 tidal signal at the receiver for a given latency, and the integration time
    for SNR p.snr. c_S = D/latency; c14 = 2 c2 (1 - c14/2)/((2 + 3 c2) c_S^2), capped at region (i);
    kappa = c14 (c14 + 3 c2) / (3 c2 (2 - c14)^2). Both near- and far-zone terms are kept."""
    omega = 2 * math.pi * p.source_freq_hz
    dq = 2 * p.source_mass_kg * p.source_radius_m * p.source_stroke_m  # Q amplitude, kg m^2
    rows = []
    for d_ly in p.distances_ly:
        r = d_ly * LY_M
        for lat in p.latencies_s:
            c_s = r / (C_SI * lat)  # in units of c
            c14 = min(2 * p.c2_max / ((2 + 3 * p.c2_max) * c_s**2), p.c14_max_region_i)
            # c2 that gives this c_S at this c14 (exact inversion of the c_S formula)
            c2 = 2 * c_s**2 * c14 / (2 - c14 - 3 * c_s**2 * c14)
            kappa = c14 * (c14 + 3 * c2) / (3 * c2 * (2 - c14) ** 2)
            far = G_SI * omega**4 * dq / (C_SI**4 * r)
            near = 2 * G_SI * omega**2 * dq / (C_SI**2 * r**3)  # |d_i d_j (Q''/r)| ~ 2 Q''/r^3
            e_sig = kappa * (far + near)
            t_int = (p.snr * p.detector_asd_s2_per_rthz / e_sig) ** 2 if e_sig > 0 else math.inf
            rows.append(
                {
                    "D_ly": d_ly,
                    "latency_s": lat,
                    "c_S_required_over_c": c_s,
                    "c14": c14,
                    "c2": c2,
                    "kappa": kappa,
                    "tidal_signal_s-2": e_sig,
                    "integration_time_s_for_snr": t_int,
                    "gr_quadrupole_scale_s-2": far,
                }
            )
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    out = {
        "provenance": "derived symbolically from the linearised Einstein-aether action (no data); "
        "the link budget uses the stated ASSUMPTION source and detector",
        "params": asdict(P),
        "controls": controls(),
        "spherical_source": spherical_source(),
        "far_field_general": far_field_general(),
        "general_source_kernel_residual": general_source_kernel(),
        "vector_sector_c13": vector_sector_c13(),
        "yagi_comparison": yagi_comparison(),
        "link_budget": link_budget(),
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

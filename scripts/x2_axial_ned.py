"""X2 S4: axial gravitational perturbations of SV with explicit matter (Yang & Lan 2603.20594).

Matter: magnetic nonlinear electrodynamics + phantom scalar (Bronnikov & Walia, arXiv:2112.13198).
In the axial (odd-parity) sector the scalar decouples and the gravitational variable H2 couples to
the EM variable H1 (Yang & Lan Eq. 56; G = c = M = 1, mu^2 = (l-1)(l+2), R = sqrt(r^2 + a^2),
f = 1 - 2/R):
  V22 = V_AF = f mu^2/R^2 + R (1/R)_{,r*r*}     (also the anisotropic-fluid single channel)
  V11 = (R-2)/(4 R^7) [2a^4(2mu^2+5) + a^2(20R + (8mu^2+17) r^2) + 6 r^2 R + (4mu^2+7) r^4]
  V12 = sqrt(3/(2 pi)) mu a (R-2)/R^(9/2),  V21 = 4 pi V12.
Since V12 V21 > 0, H1 -> sqrt(4 pi) H1 makes the matrix symmetric with off-diagonal
sqrt(4 pi) V12. That form, whose flux sum_i |d_t psi_i|^2 is conserved, is evolved here
(frequencies are unchanged).

psi_tt = psi_xx - V psi (x = r*) is integrated with second-order leapfrog on a uniform grid
covering both sides of the throat; the grid ends are far enough that their reflections arrive
after t_end.

  python scripts/x2_axial_ned.py   -> results/x2/s4_axial_ned.json   (several minutes)

Design frozen before computing: docs/state/units/X2_s4_matter.md. Derived quantities of a
conditional model: axial sector only, non-rotating, free evolution (no merger, no source
excitation eta).
"""

from __future__ import annotations

import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import x2_echo_delay as xe  # noqa: E402
import x2_scatter as xs  # noqa: E402

OUT = ROOT / "results" / "x2"
ELL = 2
MU2 = (ELL - 1) * (ELL + 2)
T_SUN = 4.925490947e-6
M_DET = 68.8  # GW150914, as in x2_s2_gw150914.py
RATIO_UL = 0.21
PMW = (0.2, 1.2)
EDGE = 0.015
# Yang & Lan supplementary tables (l = 2, M = 1), a = 2.20
PUBLISHED = {
    "a": 2.20,
    "ned_gr_led": complex(0.246140, -0.000111),
    "af": complex(0.308238, -0.000266),
}


# ----------------------------------------------------------------------------------- potentials
def v_af(r, a, mu2=MU2):
    R = np.sqrt(r * r + a * a)
    return (
        (R - 2) * (a * a * ((mu2 - 1) * R + 2) + r * r * ((mu2 + 2) * R - 6)) / (r * r + a * a) ** 3
    )


def v11(r, a, mu2=MU2):
    R = np.sqrt(r * r + a * a)
    br = 2 * a**4 * (2 * mu2 + 5) + a * a * (20 * R + (8 * mu2 + 17) * r * r)
    br = br + 6 * r * r * R + (4 * mu2 + 7) * r**4
    return (R - 2) / (4 * R**7) * br


def v12(r, a, mu2=MU2):
    R = np.sqrt(r * r + a * a)
    return math.sqrt(3 / (2 * math.pi)) * math.sqrt(mu2) * a * (R - 2) / R**4.5


def potentials(model, r, a):
    """(V_GW, V_EM, V_offdiag) on r for model 'ned', 'ned_off' (coupling removed), 'af', or an X2
    test-field proxy (1, 0, 'rw'; single channel)."""
    if model in ("ned", "ned_off"):
        off = math.sqrt(4 * math.pi) * v12(r, a) if model == "ned" else np.zeros_like(r)
        return v_af(r, a), v11(r, a), off
    if model == "af":
        return v_af(r, a), None, None
    return xs.potential_r(r, a, ELL, model), None, None


def gw_peak_x(model, a, h=0.01):
    """r* of the GW-channel barrier peak (r > 0)."""
    g = xs.Grid(a, ELL, 1, xmax=2000.0, h=h)
    v = potentials(model, g.r, a)[0]
    i = int(np.argmax(v))
    if i >= v.size - 2:
        raise ValueError("barrier peak not inside the grid")
    return float(g.x[i])


# ------------------------------------------------------------------------------------- evolution
def _grid_potentials(model, a, x):
    g = xs.Grid(a, ELL, 1, xmax=2500.0, h=0.01)
    ax = np.abs(x)
    inside = ax <= g.x[-1]
    # r ~ r* beyond 2500 M (log term < 0.6 %; V ~ 1e-6 there)
    r = np.where(inside, np.interp(ax, g.x, g.r), ax)
    return potentials(model, r, a)


def evolve(model, a, x1, x2, dx, t_end, init, observers, cfl=0.5):
    """Leapfrog on x in [-x1, x2]. init(x, dt) -> list of (psi(t=0), psi(t=dt)) per channel.
    Returns times and, per observer x, an array (n_channels, n_steps)."""
    x = np.arange(-x1, x2 + dx / 2, dx)
    vg, ve, vo = _grid_potentials(model, a, x)
    two = ve is not None
    dt = cfl * dx
    c2 = (dt / dx) ** 2
    (g_prev, g_cur), *em = init(x, dt)
    cg = 2.0 - 2.0 * c2 - dt * dt * vg
    if two:
        (e_prev, e_cur) = em[0] if em else (np.zeros_like(x), np.zeros_like(x))
        ce = 2.0 - 2.0 * c2 - dt * dt * ve
        co = dt * dt * vo
    idx = [int(np.argmin(np.abs(x - xo))) for xo in observers]
    nsteps = int(t_end / dt)
    rec = np.empty((len(idx), 2 if two else 1, nsteps))
    g_nxt = np.zeros_like(x)
    e_nxt = np.zeros_like(x) if two else None
    s = slice(1, -1)
    for n in range(nsteps):
        np.multiply(cg[s], g_cur[s], out=g_nxt[s])
        g_nxt[s] += c2 * (g_cur[2:] + g_cur[:-2])
        g_nxt[s] -= g_prev[s]
        if two:
            g_nxt[s] -= co[s] * e_cur[s]
            np.multiply(ce[s], e_cur[s], out=e_nxt[s])
            e_nxt[s] += c2 * (e_cur[2:] + e_cur[:-2])
            e_nxt[s] -= e_prev[s]
            e_nxt[s] -= co[s] * g_cur[s]
            e_prev, e_cur, e_nxt = e_cur, e_nxt, e_prev
        g_prev, g_cur, g_nxt = g_cur, g_nxt, g_prev
        for j, i in enumerate(idx):
            rec[j, 0, n] = g_cur[i]
            if two:
                rec[j, 1, n] = e_cur[i]
    t = (np.arange(nsteps) + 2) * dt
    return t, rec


# --------------------------------------------------------------------------------- mode fitting
def matrix_pencil(y, dt, order=12):
    """Complex frequencies w (psi ~ e^{-i w t}) and amplitudes of a sum of damped sinusoids."""
    n = y.size
    L = n // 3
    Y = np.array([y[i : i + L + 1] for i in range(n - L)])
    _, sv, vh = np.linalg.svd(Y, full_matrices=False)
    p = min(order, int(np.sum(sv > sv[0] * 1e-10)))
    v = vh[:p].conj().T
    v1, v2 = v[:-1], v[1:]
    z = np.linalg.eigvals(np.linalg.pinv(v1) @ v2)
    w = 1j * np.log(z) / dt
    Z = np.vander(z, n, increasing=True).T
    amp = np.linalg.lstsq(Z, y.astype(complex), rcond=None)[0]
    keep = w.real > 0
    return w[keep], np.abs(amp[keep]) * 2


def reproduce(dx, model="ned", t_fit=(350.0, 6000.0), sample=1.0):
    """Gaussian in the EM channel (H2 = 0) at a = 2.20, as in Yang & Lan Sec. 4; the coupling-off
    control puts the Gaussian in the GW channel instead."""
    a = PUBLISHED["a"]
    x_pk = gw_peak_x("af", a)
    xc, sig = x_pk + 10.0, 2.0
    t_end = t_fit[1] + 1.0
    X = t_end / 2 + 150.0

    def init(x, dt):
        gauss = np.exp(-((x - xc) ** 2) / sig**2)
        z = np.zeros_like(x)
        if model == "ned":
            return [(z, z.copy()), (gauss, gauss.copy())]
        return [(gauss, gauss.copy())]

    t, rec = evolve("ned" if model == "ned" else "af", a, X, X, dx, t_end, init, [50.0])
    step = max(1, int(round(sample / (t[1] - t[0]))))
    m = (t >= t_fit[0]) & (t <= t_fit[1])
    y = rec[0, 0, m][::step]
    w, amp = matrix_pencil(y, (t[1] - t[0]) * step)
    order = np.argsort(-amp)
    return [
        {"w_re": float(w[i].real), "w_im": float(w[i].imag), "amp": float(amp[i])} for i in order
    ]


def match(modes, target, rel_re=0.005, im_lo=0.67, im_hi=1.5):
    best = min(modes, key=lambda m: abs(m["w_re"] - target.real))
    ok = (
        abs(best["w_re"] / target.real - 1) <= rel_re
        and best["w_im"] < 0
        and im_lo <= best["w_im"] / target.imag <= im_hi
    )
    return best, bool(ok)


# ------------------------------------------------------------------------------ echo comparison
def echo_run(model, a, dx, w0=0.5307, tau=12.2):
    """Inward ringdown-shaped packet in the GW channel, 25 M inside that channel's peak (S3)."""
    x_pk = gw_peak_x(model, a)
    dte_4l = 4.0 * x_pk  # this model's own peak-to-peak tortoise distance
    x0 = x_pk - 25.0

    def p(u):
        return np.where(u > 0, np.exp(-u / tau) * np.sin(w0 * u), 0.0)

    def init(x, dt):
        u0 = x0 - x
        return [(p(u0), p(u0 - dt))]

    x_obs = 1500.0
    t_end = 3.0 * dte_4l + x_obs + 400.0
    t, rec = evolve(model, a, 6000.0, 7000.0, dx, t_end, init, [x_obs])
    dt = t[1] - t[0]
    uu = np.arange(0, 40 * tau, dt)
    pe = p(uu) ** 2
    ref_hrss = math.sqrt(np.trapezoid(pe, dx=dt))
    ref_c = float(np.trapezoid(uu * pe, uu) / np.trapezoid(pe, uu))
    t_ref = x_obs - x0 + ref_c  # free-flight arrival of the packet centroid
    gw = rec[0, 0]
    em = rec[0, 1] if rec.shape[1] > 1 else None
    pulses = []
    for k in (1, 2):
        tc = t_ref + k * dte_4l
        m = (t > tc - 0.45 * dte_4l) & (t < tc + 0.45 * dte_4l)
        e2 = gw[m] ** 2
        cen = float(np.trapezoid(t[m] * e2, t[m]) / np.trapezoid(e2, t[m]))
        d = {
            "hrss_ratio": math.sqrt(np.trapezoid(e2, t[m])) / ref_hrss,
            "arrival_after_free_flight_M": cen - t_ref,
            "rms_width_M": math.sqrt(
                np.trapezoid((t[m] - cen) ** 2 * e2, t[m]) / np.trapezoid(e2, t[m])
            ),
        }
        if em is not None:
            gdot, edot = np.gradient(gw[m], dt), np.gradient(em[m], dt)
            eg, ee = np.trapezoid(gdot**2, dx=dt), np.trapezoid(edot**2, dx=dt)
            d["em_energy_fraction"] = float(ee / (eg + ee))
        pulses.append(d)
    return {
        "model": str(model),
        "a_minus_2M": a - 2,
        "dx": dx,
        "dt_4L_M": dte_4l,
        "x_peak_M": x_pk,
        "pulses": pulses,
        "spacing_M": pulses[1]["arrival_after_free_flight_M"]
        - pulses[0]["arrival_after_free_flight_M"],
    }


def pmw_eta_max(run, m_det=M_DET):
    """eta_max from pulses inside the PMW (times after the direct signal, in seconds)."""
    tm = m_det * T_SUN
    t1 = run["pulses"][0]["arrival_after_free_flight_M"] * tm
    sp = run["spacing_M"] * tm
    times = [t1 + k * sp for k in range(6)]
    ks = [k for k, tk in enumerate(times) if PMW[0] + EDGE <= tk <= PMW[1] - EDGE]
    edge = [k for k, tk in enumerate(times) if abs(tk - PMW[0]) < EDGE or abs(tk - PMW[1]) < EDGE]
    p1, p2 = run["pulses"][0]["hrss_ratio"], run["pulses"][1]["hrss_ratio"]
    q = p2 / p1  # per-round-trip amplitude factor, extrapolated beyond pulse 2
    amps = [p1 * q**k for k in ks]
    F = math.sqrt(sum(x * x for x in amps)) if amps else 0.0
    return {
        "pulse_times_s": times[:4],
        "pulses_in_PMW": [k + 1 for k in ks],
        "pulses_on_edge": [k + 1 for k in edge],
        "F": F,
        "eta_max": RATIO_UL / F if F else None,
        "applicable_timing": bool(ks and ks[0] == 0 and len(ks) >= 2 and not edge),
    }


def _task(job):
    kind, *args = job
    return job, (reproduce(*args) if kind == "rep" else echo_run(*args))


def main() -> int:
    out = {
        "label": "derived; axial sector, NED magnetic + phantom scalar (Yang & Lan 2603.20594); "
        "conditional",
        "design": "docs/state/units/X2_s4_matter.md (frozen before computing)",
        "published_target": {k: str(v) for k, v in PUBLISHED.items()},
    }
    a = 2 + xe.a_for_delay(0.45 / (M_DET * T_SUN))  # S2's applicable GW150914 cell (dt = 0.45 s)
    jobs = [("rep", 0.1, "ned"), ("rep", 0.05, "ned"), ("rep", 0.1, "af")]
    jobs += [("echo", m, a, dx) for m in ("ned", 1, "rw") for dx in (0.2, 0.1)]
    # independent runs, one process each; the comparison is reported only if the reproduction passes
    with ProcessPoolExecutor(max_workers=min(len(jobs), os.cpu_count() or 1)) as ex:
        res = dict(ex.map(_task, jobs))
    rep = {}
    for dx in (0.1, 0.05):
        modes = res[("rep", dx, "ned")]
        best, ok = match(modes, PUBLISHED["ned_gr_led"])
        rep[f"ned_dx{dx}"] = {"modes": modes[:8], "closest": best, "pass": ok}
    ctl = res[("rep", 0.1, "af")]
    best, ok = match(ctl, PUBLISHED["af"])
    rep["control_coupling_off_af_dx0.1"] = {"modes": ctl[:8], "closest": best, "pass": ok}
    a0, a1 = rep["ned_dx0.1"]["closest"], rep["ned_dx0.05"]["closest"]
    rep["convergence_ok"] = bool(
        abs(a0["w_re"] / a1["w_re"] - 1) <= 1e-3 and abs(a0["w_im"] / a1["w_im"] - 1) <= 0.10
    )
    rep["pass"] = bool(rep["ned_dx0.05"]["pass"] and rep["convergence_ok"])
    out["reproduction"] = rep
    if rep["pass"]:
        cmp = {}
        for model in ("ned", 1, "rw"):
            r0, r1 = res[("echo", model, a, 0.2)], res[("echo", model, a, 0.1)]
            conv = all(
                abs(r0["pulses"][k]["hrss_ratio"] / r1["pulses"][k]["hrss_ratio"] - 1) <= 0.01
                and abs(
                    r0["pulses"][k]["arrival_after_free_flight_M"]
                    - r1["pulses"][k]["arrival_after_free_flight_M"]
                )
                <= 2.0
                for k in (0, 1)
            )
            cmp[str(model)] = {"dx0.2": r0, "dx0.1": r1, "converged": conv, "pmw": pmw_eta_max(r1)}
        out["comparison_gw150914_cell"] = cmp
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s4_axial_ned.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float)[:6000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

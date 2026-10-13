"""X2 S3: time-domain check of the frequency-domain echo amplitudes used in S2.

A ringdown-shaped wave packet starts just inside our barrier, moving towards the throat, on the full
symmetric SV double barrier (s = 1, l = 2). psi_tt = psi_xx - V(x) psi is integrated with a
second-order leapfrog scheme. An observer far outside records the echo pulses. Their h_rss
relative to the packet's own h_rss is compared with the frequency-domain prediction:
  pulse k: sqrt( int |P(w)|^2 |T|^2 |R|^(2(2k-1)) dw / int |P(w)|^2 dw ),
with P the packet spectrum and R, T from scripts/x2_scatter.py.

  python scripts/x2_s3_timedomain.py   -> results/x2/s3_timedomain.json  (about a minute)

Gate (docs/state/units/X2.md, written before the run): amplitudes within 10 % (pulses 1, 2) and
pulse spacing within 1 % of the geometric dt = 4L. The spacing half was added after review
(post hoc implementation of the pre-written criterion; it was not measured in the first run).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import x2_echo_delay as xe  # noqa: E402
import x2_scatter as xs  # noqa: E402

OUT = ROOT / "results" / "x2"


def run(a=None, s=1, dt_target_M=1181.0, w0=0.5307, tau=12.2, dx=0.2, cfl=0.5):
    if a is None:
        a = 2 + xe.a_for_delay(dt_target_M)
    g = xs.Grid(a, 2, s, xmax=2500.0, h=0.01)
    i_pk = int(np.argmax(g.V))
    x_pk = g.x[i_pk]
    X1, X2 = 6000.0, 7000.0
    x = np.arange(-X1, X2 + dx / 2, dx)
    ax = np.abs(x)
    V = np.where(ax <= g.x[-1], np.interp(ax, g.x, g.V), 6.0 / np.maximum(ax, 1.0) ** 2)
    dt = cfl * dx
    # inward-moving packet psi(x, 0) = p(x0 - x) with p(u) = exp(-u/tau) sin(w0 u) for u > 0,
    # i.e. a left-moving wave psi(x, t) = p(x0 - x - t); its time derivative is -p'(x0 - x - t)
    x0 = x_pk - 25.0

    def p(u):
        return np.where(u > 0, np.exp(-u / tau) * np.sin(w0 * u), 0.0)

    def dp(u):
        return np.where(u > 0, np.exp(-u / tau) * (w0 * np.cos(w0 * u) - np.sin(w0 * u) / tau), 0.0)

    u0 = x0 - x
    psi_prev = p(u0)
    psi = p(u0 - dt)  # psi at t = dt (exact free propagation; V is small at x0 - 0 over one step)
    x_obs = 1500.0
    i_obs = int(np.argmin(np.abs(x - x_obs)))
    t_end = 3.0 * xe.echo_delay(a) + x_obs + 400.0
    nsteps = int(t_end / dt)
    rec = np.empty(nsteps)
    c2 = (dt / dx) ** 2
    coef = 2.0 - 2.0 * c2 - dt * dt * V  # psi_next = coef psi + c2 (psi_l + psi_r) - psi_prev
    nxt = np.zeros_like(psi)
    for n in range(nsteps):
        np.multiply(coef[1:-1], psi[1:-1], out=nxt[1:-1])
        nxt[1:-1] += c2 * (psi[2:] + psi[:-2])
        nxt[1:-1] -= psi_prev[1:-1]
        psi_prev, psi, nxt = psi, nxt, psi_prev
        rec[n] = psi[i_obs]
    t = (np.arange(nsteps) + 2) * dt
    # pulse windows: arrival of pulse k at the observer ~ t_k = k*dt_echo + (x_obs - x0) offset
    dte = xe.echo_delay(a)
    t_direct = x_obs - x0  # straight-out travel time (reference)
    ref_hrss = math.sqrt(np.trapezoid(p(np.arange(0, 40 * tau, dt)) ** 2, dx=dt))
    pulses, t_arr, width = [], [], []
    for k in (1, 2):
        tc = t_direct + k * dte
        m = (t > tc - 0.45 * dte) & (t < tc + 0.45 * dte)
        e2 = rec[m] ** 2
        pulses.append(math.sqrt(np.trapezoid(e2, t[m])) / ref_hrss)
        tbar = float(np.trapezoid(t[m] * e2, t[m]) / np.trapezoid(e2, t[m]))  # energy centroid
        t_arr.append(tbar)
        width.append(
            float(math.sqrt(np.trapezoid((t[m] - tbar) ** 2 * e2, t[m]) / np.trapezoid(e2, t[m])))
        )
    # the reference packet's own energy centroid and rms width
    uu = np.arange(0, 40 * tau, dt)
    pe = p(uu) ** 2
    ref_c = float(np.trapezoid(uu * pe, uu) / np.trapezoid(pe, uu))
    ref_w = float(math.sqrt(np.trapezoid((uu - ref_c) ** 2 * pe, uu) / np.trapezoid(pe, uu)))
    # frequency-domain prediction with the packet spectrum
    u = np.arange(0, 60 * tau, 0.05)
    P = np.fft.rfft(p(u))
    wf = 2 * np.pi * np.fft.rfftfreq(u.size, d=0.05)
    keep = (wf > math.sqrt(max(g.v0, 1e-12)) * 1.001) & (wf < 3.0)
    R2, T2 = xs.single_barrier(g, wf[keep])
    P2 = np.abs(P[keep]) ** 2
    norm = np.sum(np.abs(P) ** 2)
    pred = [math.sqrt(np.sum(P2 * T2 * R2 ** (2 * k - 1)) / norm) for k in (1, 2)]
    return {
        "a_minus_2M": a - 2,
        "dt_echo_M": dte,
        "packet": {"w0": w0, "tau_M": tau},
        "grid": {"dx": dx, "dt": dt},
        "time_domain_pulse_hrss_ratio": pulses,
        "frequency_domain_prediction": pred,
        "relative_difference": [abs(pulses[i] / pred[i] - 1) for i in range(2)],
        "pulse_spacing_M": t_arr[1] - t_arr[0],
        "spacing_vs_4L_relative": (t_arr[1] - t_arr[0]) / dte - 1,
        "first_arrival_minus_direct_plus_4L_M": t_arr[0] - (t_direct + dte + ref_c),
        "pulse_rms_width_M": width,
        "packet_rms_width_M": ref_w,
    }


def main() -> int:
    out = {"label": "derived; s = 1, l = 2 test field; numerical consistency check"}
    out["run_dx0.2"] = run(dx=0.2)
    out["run_dx0.1"] = run(dx=0.1)
    out["pass_10pct"] = bool(max(out["run_dx0.1"]["relative_difference"]) <= 0.10)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s3_timedomain.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""X1 S4: model-free lower bound on the far-side chi^2 for the 2018-07-22 flare (checks S3).

Independent of the ray tracer, the orbit, eta, P, the inclination and light-travel delays. Every
far-side image has impact parameter b < b_c(a) (exact, scripts/wh_farside.py), so the spot's
total-light centroid u(t), a positive-weighted mean of image positions, lies in the disk
|u| <= r = b_c(a) theta_g around the true centre. With the same static light as S3 (flux
S = fs mean(F), fs <= min(0.5, 0.999 min F / mean F), offset |d| <= D):
  c(t) = c0 + [(F_t - S) u_t + S d] / F_t = c0 + q / F_t + v_t,
  q = S d, |q| <= S D;  v_t = (1 - S / F_t) u_t, |v_t| <= r (1 - S / F_t).
Every u_t is free in its disk, so the minimum chi^2 over (c0, S, q, v_t) is a lower bound on the
chi^2 of *any* far-side spot model under S3's static-light premises. The model is linear in
(c0, q, v_t) and every constraint is a second-order cone, so the problem is convex; SLSQP from
several starts finds the global minimum (start spread reported).

Gate (written in docs/state/units/X1.md before this run):
- consistency: chi2_LB(a) <= the S3 far-side chi^2 for every a and offset cap (a violation is an
  unexplained disagreement and S3 is then not reported as a finding);
- chi2_LB(a) - (best S3 ordinary chi^2) >= 9 (ASSUMPTION) for a <= 7 M means S3's disfavour holds
  without the tracer or the orbit model; below 9 it depends on them, and is stated so.

  python scripts/x1_s4.py   -> results/x1/s4_disk_bound.json

Every number is derived (our computation) under unestablished premises (X1.md).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import integrate, optimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh_farside as wf  # noqa: E402
import x1_s2  # noqa: E402
import x1_s3  # noqa: E402
import x1_tracer as xt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x1"
D_CAPS = (20.0, 40.0, 100.0, 200.0, 400.0)  # uas; 20 is S3's pre-written cap, the rest as in S3
N_STARTS = 6
DELTA_PASS = 9.0  # ASSUMPTION: Delta chi2 threshold for "holds model-free"


def disk_bound(d: dict, r: float, d_cap: float, fs_max: float = 0.5, seed: int = 0) -> dict:
    """Minimum chi^2 of c0 + q/F_t + v_t with |v_t| <= r (1 - S/F_t), |q| <= S d_cap."""
    y, f = d["xy"], d["flux"]
    w = 1.0 / d["err"] ** 2
    n = len(f)
    s_max = float(min(fs_max, 0.999 * f.min() / f.mean())) * f.mean()

    def unpack(x):
        return x[:2], x[2:4], x[4], x[5:].reshape(n, 2)

    def chi2(x):
        c0, q, _, v = unpack(x)
        res = y - c0 - q / f[:, None] - v
        return float((w * res**2).sum())

    def grad(x):
        c0, q, _, v = unpack(x)
        g_res = -2.0 * w * (y - c0 - q / f[:, None] - v)
        return np.concatenate([g_res.sum(0), (g_res / f[:, None]).sum(0), [0.0], g_res.ravel()])

    def cons(x):
        _, q, s, v = unpack(x)
        rho = r * (1.0 - s / f)
        return np.concatenate([rho**2 - (v**2).sum(1), rho, [(s * d_cap) ** 2 - q @ q]])

    def cons_jac(x):
        _, q, s, v = unpack(x)
        rho = r * (1.0 - s / f)
        j = np.zeros((2 * n + 1, x.size))
        idx = np.arange(n)
        j[idx, 4] = -2.0 * rho * r / f
        j[idx, 5 + 2 * idx] = -2.0 * v[:, 0]
        j[idx, 6 + 2 * idx] = -2.0 * v[:, 1]
        j[n + idx, 4] = -r / f
        j[2 * n, 2:4] = -2.0 * q
        j[2 * n, 4] = 2.0 * s * d_cap**2
        return j

    rng = np.random.default_rng(seed)
    c_w = (y * w).sum(0) / w.sum(0)
    best, chis = None, []
    for k in range(N_STARTS):
        s0 = 0.0 if k == 0 else rng.uniform(0.0, s_max)
        x0 = np.concatenate(
            [c_w + rng.normal(0, 20, 2) * (k > 0), np.zeros(2), [s0], np.zeros(2 * n)]
        )
        if k > 0:
            ang = rng.uniform(0, 2 * math.pi)
            x0[2:4] = 0.5 * s0 * d_cap * np.array([math.cos(ang), math.sin(ang)])
        bounds = [(None, None)] * 4 + [(0.0, s_max)] + [(None, None)] * (2 * n)
        rr = optimize.minimize(
            chi2,
            x0,
            jac=grad,
            method="SLSQP",
            bounds=bounds,
            constraints=[{"type": "ineq", "fun": cons, "jac": cons_jac}],
            options={"maxiter": 2000, "ftol": 1e-12},
        )
        viol = float(max(0.0, -cons(rr.x).min()))
        chis.append(rr.fun if viol < 1e-6 else math.inf)
        if viol < 1e-6 and (best is None or rr.fun < best.fun):
            best = rr
    if best is None:
        raise RuntimeError("no feasible SLSQP solution")
    c0, q, s, v = unpack(best.x)
    rho = r * (1.0 - s / f)
    on_edge = int((np.sqrt((v**2).sum(1)) > rho - 1e-3).sum())
    finite = [c for c in chis if math.isfinite(c)]
    return {
        "chi2_lb": float(best.fun),
        "start_spread": float(max(finite) - min(finite)),
        "n_starts_feasible": len(finite),
        "static_fraction": float(s / f.mean()),
        "static_offset_uas": float(math.hypot(*q) / s) if s > 1e-9 else 0.0,
        "centre_uas": c0.tolist(),
        "epochs_on_disk_edge": on_edge,
    }


def phi_far_hamilton(b: float, a: float, re: float, r0: float = 20.0) -> tuple[float, float]:
    """Far-side swept angle from Hamilton's equations (p_r evolves), separate from the tracer.

    H = (-1/f + f p_r^2 + b^2/A) / 2 with p_t = -1, p_phi = b, A = r^2 + a^2, f = 1 - 2/sqrt(A),
    integrated (DOP853) from coordinate r0 through the throat to r = -sqrt(re^2 - a^2); the range
    r0 -> infinity is added by its own quadrature in u = 1/r. Returns (phi, H at the end).
    """

    def rhs(_, y):
        r, _phi, pr = y
        aa = r * r + a * a
        f = 1.0 - 2.0 / math.sqrt(aa)
        fp = 2.0 * r / aa**1.5
        return [f * pr, b / aa, -0.5 * (fp / f**2 + fp * pr * pr - 2.0 * r * b * b / aa**2)]

    a0 = r0 * r0 + a * a
    f0 = 1.0 - 2.0 / math.sqrt(a0)
    r_t = -math.sqrt(re * re - a * a)

    def reached(_, y):
        return y[0] - r_t

    reached.terminal, reached.direction = True, -1.0
    pr0 = -math.sqrt((1.0 / f0 - b * b / a0) / f0)
    sol = integrate.solve_ivp(
        rhs, (0.0, 1e5), [r0, 0.0, pr0], method="DOP853", events=reached, rtol=1e-12, atol=1e-13
    )
    y = sol.y_events[0][0]
    aa = y[0] ** 2 + a * a
    f = 1.0 - 2.0 / math.sqrt(aa)

    def tail(u):
        r = 1.0 / u
        ar = r * r + a * a
        return b * r * r / (ar * math.sqrt(1.0 - b * b * (1.0 - 2.0 / math.sqrt(ar)) / ar))

    t = integrate.quad(tail, 0.0, 1.0 / r0, epsabs=1e-13, epsrel=1e-12)[0]
    return float(y[1] + t), float(0.5 * (-1.0 / f + f * y[2] ** 2 + b * b / aa))


def hamilton_check(re: float = 9.0) -> dict:
    """Tracer far-side swept angles vs phi_far_hamilton (S4 second route; tolerance 1e-8)."""
    rows = []
    for a in (2.1, 3.0, 5.0, 7.0):
        for x in (0.3, 0.9, 0.999):
            b = x * wf.b_crit(a)
            ph, h = phi_far_hamilton(b, a, re)
            ref = xt.phi_to_emitter(b, a, re, "far")[0]
            rows.append(
                {"a_M": a, "b_over_bc": x, "phi_hamilton": ph, "phi_tracer": ref, "H_end": h}
            )
    rel = max(abs(r["phi_hamilton"] / r["phi_tracer"] - 1.0) for r in rows)
    return {"emitter_R_M": re, "max_rel_diff": rel, "pass_1e-8": bool(rel < 1e-8), "rays": rows}


def main() -> int:
    d = x1_s3.load()
    s3 = json.loads((OUT / "s3_avg_fit.json").read_text())
    ord_best = min(r["chi2"] for r in s3["fits"] if r["model"] == "ordinary")
    far_s3 = {}
    for r in s3["fits"]:
        if r["model"] == "far":
            far_s3[r["a_M"]] = min(far_s3.get(r["a_M"], math.inf), r["chi2"])
    sens_s3 = {
        s["static_offset_max_uas"]: s["far_a2.1_chi2"] for s in s3["sensitivity_static_offset_P60"]
    }
    rows = []
    for a in x1_s2.A_FAR:
        r_uas = wf.b_crit(a) * x1_s2.THETA_G
        for cap in D_CAPS:
            b = disk_bound(d, r_uas, cap)
            s3_ref = far_s3[a] if cap == 20.0 else (sens_s3.get(cap) if a == 2.1 else None)
            rows.append(
                {
                    "a_M": a,
                    "b_c_M": wf.b_crit(a),
                    "disk_radius_uas": r_uas,
                    "static_offset_max_uas": cap,
                }
                | b
                | {
                    "s3_far_chi2": s3_ref,
                    "consistent_lb_le_s3": None
                    if s3_ref is None
                    else bool(b["chi2_lb"] <= s3_ref + 1e-6),
                    "delta_lb_minus_s3_ordinary_best": b["chi2_lb"] - ord_best,
                }
            )
    # disk radius at which the bound reaches the best ordinary chi^2 (cap 20 uas), for scale
    lo, hi = 0.0, 200.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if disk_bound(d, mid, 20.0)["chi2_lb"] > ord_best else (lo, mid)
    r_equal = 0.5 * (lo + hi)
    pre = [x for x in rows if x["static_offset_max_uas"] == 20.0]
    out = {
        "labels": "derived (our computation); data: processed numerical values (CDS"
        " J/A+A/618/L10 avg.dat); exploratory, no held-out flare (same flare as S3 and #172)",
        "method": "convex SOCP lower bound: spot centroid anywhere in the disk b < b_c(a) per"
        " epoch; S3 static-light nuisance; no tracer, orbit, P, eta, inclination or delay",
        "s3_ordinary_best_chi2": ord_best,
        "disk_radius_uas_where_lb_equals_ordinary_best_cap20": r_equal,
        "disk_radius_M_where_lb_equals_ordinary_best_cap20": r_equal / x1_s2.THETA_G,
        "theta_g_uas_per_M": x1_s2.THETA_G,
        "gate": {
            "all_consistent": all(x["consistent_lb_le_s3"] is not False for x in rows),
            "delta_threshold_assumption": DELTA_PASS,
            "holds_model_free_cap20": {
                str(x["a_M"]): bool(x["delta_lb_minus_s3_ordinary_best"] >= DELTA_PASS) for x in pre
            },
        },
        "rows": rows,
        "second_route_swept_angle": hamilton_check(),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s4_disk_bound.json").write_text(json.dumps(out, indent=1) + "\n")
    for x in rows:
        print(
            f"a={x['a_M']:.1f} r={x['disk_radius_uas']:.1f}uas"
            f" cap={x['static_offset_max_uas']:.0f} LB={x['chi2_lb']:.2f} S3={x['s3_far_chi2']}"
            f" dLB={x['delta_lb_minus_s3_ordinary_best']:.2f} fs={x['static_fraction']:.2f}"
            f" spread={x['start_spread']:.1e} edge={x['epochs_on_disk_edge']}"
        )
    print("r where LB = ordinary best:", r_equal, "uas")
    print(json.dumps(out["gate"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

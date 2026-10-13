"""X1 S4b: independent check of the tracer's per-image weights (point-source magnification).

The tracer weights each image by image_mag = b / (|dPhi/db| cos(al) sin(gamma)) (x1_tracer.py).
This script reaches the same quantity by another route, which shares only the metric:
- the swept angle Phi(b) and the photon direction at the emitter come from Hamilton's equations
  (p_r evolves; DOP853), not from the tracer's quadratures; cos(al) is read from p_r, not from
  the potential formula;
- the magnification is 1 / |det J|, J the 2-D Jacobian (central differences on the screen) of the
  map from the screen position to the ray's perpendicular offset from the emitter, built in 3-D
  from the hit point m = cos(Phi) o + sin(Phi) s_hat on the emitter sphere (areal radius R) and
  the ray direction there. This also tests the sin(gamma) and cos(al) factors.
image_mag is screen area per unit source solid angle on the emitter sphere; it is compared here per
unit cross-section, i.e. divided by R^2 (constant at fixed R; cancels in centroids and flux ratios).

Gate (written in docs/state/units/X1.md before this run): every far-side image at a = 2.1 and 5 M,
R = 9 M, i = 155 deg, 8 phases, within 0.5 % (ASSUMPTION, pre-written tolerance).

  python scripts/x1_s4b.py   -> results/x1/s4b_image_weights.json

Every number is derived (our computation) under unestablished premises (X1.md).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import integrate

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh_farside as wf  # noqa: E402
import x1_tracer as xt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x1"
TOL_REL = 5e-3  # ASSUMPTION: pre-written tolerance (X1.md "S4b gate")
R_EMIT, INCL, N_PHASE = 9.0, 155.0, 8
R0 = 100.0  # coordinate radius where the integration starts; beyond it, a quadrature tail
H_REL = (1e-4, 2e-4)  # screen step / b for the central differences (convergence pair)
FOLD_FRAC = 0.02  # near side: step at most this fraction of the distance to the fold
CASES = (  # (a, side, gate)
    (2.1, "far", True),
    (5.0, "far", True),
    (3.0, "far", False),
    (7.0, "far", False),
    (0.0, "near", False),
)


def hit_hamilton(b: float, a: float, re: float, side: str, branch: int = 0):
    """Swept angle, signed radial cosine and tangential sine of the backward ray at the emitter.

    H = (-1/f + f p_r^2 + b^2/A) / 2, p_t = -1, p_phi = b, A = r^2 + a^2, f = 1 - 2/sqrt(A); start
    at r = R0 inward. Far side: stop at r = -sqrt(re^2 - a^2). Near side: branch 0 stops at the
    first (inward) crossing of r = sqrt(re^2 - a^2), branch 1 at the outward one after the turn.
    cos(al) is the static-frame radial component, positive away from the throat. Returns None if
    the ray never reaches the emitter sphere on that branch.
    """

    def rhs(_, y):
        r, _phi, pr = y
        aa = r * r + a * a
        f = 1.0 - 2.0 / math.sqrt(aa)
        fp = 2.0 * r / aa**1.5
        return [f * pr, b / aa, -0.5 * (fp / f**2 + fp * pr * pr - 2.0 * r * b * b / aa**2)]

    rt = math.sqrt(re * re - a * a)
    target, direction = (-rt, -1.0) if side == "far" else (rt, -1.0 if branch == 0 else 1.0)

    def reached(_, y):
        return y[0] - target

    reached.terminal, reached.direction = True, direction

    def fallen(_, y):  # horizon (a < 2) or a far-side escape the event did not catch
        return math.sqrt(y[0] ** 2 + a * a) - 2.0 - 1e-6

    fallen.terminal = True
    a0 = R0 * R0 + a * a
    f0 = 1.0 - 2.0 / math.sqrt(a0)
    pr0 = -math.sqrt((1.0 / f0 - b * b / a0) / f0)
    sol = integrate.solve_ivp(
        rhs,
        (0.0, 1e5),
        [R0, 0.0, pr0],
        method="DOP853",
        events=[reached, fallen],
        rtol=1e-12,
        atol=1e-13,
    )
    if not len(sol.y_events[0]):
        return None
    r, phi, pr = sol.y_events[0][0]
    aa = r * r + a * a
    f = 1.0 - 2.0 / math.sqrt(aa)

    def tail(u):
        rr = 1.0 / u
        ar = rr * rr + a * a
        return b * rr * rr / (ar * math.sqrt(1.0 - b * b * (1.0 - 2.0 / math.sqrt(ar)) / ar))

    phi += integrate.quad(tail, 0.0, 1.0 / R0, epsabs=1e-13, epsrel=1e-12)[0]
    # dr/dlambda = f p_r; the areal radius grows away from the throat when r dr > 0
    cos_al = math.copysign(1.0, r) * f * pr
    sin_al = math.sqrt(f) * b / math.sqrt(aa)
    return float(phi), float(cos_al), float(sin_al)


def mag_jacobian(a, re, side, branch, s0, o, n, h_rel):
    """1 / |det J| of screen position -> perpendicular offset from the emitter (see module doc)."""
    b0 = float(np.linalg.norm(s0))
    sh0 = s0 / b0
    w = np.cross(o, sh0)  # second sky axis
    phi0, ca, sa = hit_hamilton(b0, a, re, side, branch)
    m0 = math.cos(phi0) * o + math.sin(phi0) * sh0
    t0 = -math.sin(phi0) * o + math.cos(phi0) * sh0
    k = ca * m0 + sa * t0
    p1 = np.cross(k, w)
    p1 /= np.linalg.norm(p1)
    p2 = np.cross(k, p1)

    def q(s):
        bb = float(np.linalg.norm(s))
        ph = hit_hamilton(bb, a, re, side, branch)[0]
        m = math.cos(ph) * o + math.sin(ph) * (s / bb)
        d = re * (m - n)
        return np.array([d @ p1, d @ p2])

    h = h_rel * b0
    if side == "near":  # stay on the branch: the near-side image curve folds at b_f
        h = min(h, FOLD_FRAC * (xt._fold_b(a, re) - b0))
    j = np.column_stack([(q(s0 + h * ax) - q(s0 - h * ax)) / (2 * h) for ax in (sh0, w)])
    return 1.0 / abs(np.linalg.det(j)), float(np.linalg.norm(m0 - n)), ca * ca + sa * sa - 1.0


def check_case(a: float, side: str) -> list[dict]:
    inc = math.radians(INCL)
    o = np.array([math.sin(inc), 0.0, math.cos(inc)])
    ex = np.array([math.cos(inc), 0.0, -math.sin(inc)])
    b_max = 3.0 * R_EMIT if side == "near" else wf.b_crit(a) * (1 - 1e-9)
    tables = xt._tables(a, R_EMIT, side, b_max)
    rows = []
    for ph in np.linspace(0.0, 2 * math.pi, N_PHASE, endpoint=False):
        n = np.array([math.cos(ph), math.sin(ph), 0.0])
        cg = float(n @ o)
        gamma = math.acos(cg)
        proj = n - cg * o
        e = proj / np.linalg.norm(proj) if np.linalg.norm(proj) > 1e-12 else ex
        images, _ = xt._images(a, R_EMIT, side, gamma, tables, 0)
        for bs, sgn, dphidb, br in images:
            # the tracer's weight is screen area per unit source solid angle on the emitter
            # sphere; per unit cross-section it is divided by R^2 (a constant at fixed R, so it
            # cancels in every centroid and flux ratio)
            mag_t = xt.image_mag(R_EMIT, bs, dphidb, gamma) / R_EMIT**2
            mags = []
            for hr in H_REL:
                mj, miss, null_err = mag_jacobian(a, R_EMIT, side, br, sgn * bs * e, o, n, hr)
                mags.append(mj)
            rows.append(
                {
                    "a_M": a,
                    "side": side,
                    "phase_rad": float(ph),
                    "image": "same" if sgn > 0 else "opposite",
                    "branch": br,
                    "b_M": bs,
                    "mag_tracer_over_R2": mag_t,
                    "mag_jacobian": mags[0],
                    "rel_diff": mags[0] / mag_t - 1.0,
                    "step_change_rel": mags[1] / mags[0] - 1.0,
                    "hit_miss": miss,
                    "null_err": null_err,
                }
            )
    return rows


def main() -> int:
    rows = []
    for a, side, _ in CASES:
        rows += check_case(a, side)
    gate = [r for r in rows if r["side"] == "far" and r["a_M"] in (2.1, 5.0)]
    summary = {}
    for a, side, g in CASES:
        sel = [r for r in rows if r["a_M"] == a and r["side"] == side]
        summary[f"{side}_a{a}"] = {
            "n_images": len(sel),
            "max_abs_rel_diff": max(abs(r["rel_diff"]) for r in sel),
            "max_abs_step_change": max(abs(r["step_change_rel"]) for r in sel),
            "max_hit_miss": max(r["hit_miss"] for r in sel),
            "gate_case": g,
        }
    out = {
        "stage": "X1 S4b image-weight check",
        "labels": "derived (our computation) under the X1 premises",
        "emitter_R_M": R_EMIT,
        "incl_deg": INCL,
        "n_phase": N_PHASE,
        "tolerance_rel": TOL_REL,
        "gate_pass": bool(max(abs(r["rel_diff"]) for r in gate) <= TOL_REL),
        "summary": summary,
        "images": rows,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s4b_image_weights.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"gate_pass": out["gate_pass"], "summary": summary}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

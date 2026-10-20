"""W3 negative-mass microlens inside a macro-lensed arc: gates and measurements.

Plan: results/w3_macro_shear/PLAN.md.

Everything is simulated / model_prediction. Writes results/w3_macro_shear/result.json and a figure.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

from jwst_anomaly import exotic_sim as es

OUT = Path(__file__).resolve().parents[1] / "results" / "w3_macro_shear"
HOLE_FRAC = 0.1  # ASSUMPTION (plan): "hole" = mu / mu_macro < 0.1


def as_image(t):
    ny, nx = t.meta["shape"]
    return (
        t["mu"].reshape(ny, nx),
        t["beta_x"].reshape(ny, nx)[0],
        t["beta_y"].reshape(ny, nx)[:, 0],
    )


def gate_g1():
    out = {}
    for sign in (1, -1):
        mu, bx, by = as_image(es.macro_magnification_map(1.0, 1.0, sign, 3.0, 3.0, 0.01))
        r = np.hypot(*np.meshgrid(bx, by))
        exact = es.total_magnification(r.ravel(), 1.0, sign).reshape(r.shape)
        sel = (
            (exact < 20) & (exact > 0.5) & (np.abs(r - 2.0) > 0.05)
        )  # avoid the caustic's own pixels
        bins = np.arange(0.0, 3.0, 0.1)
        idx = np.digitize(r[sel], bins)
        dev = []
        for k in np.unique(idx):
            m = idx == k
            dev.append(abs(mu[sel][m].mean() / exact[sel][m].mean() - 1))
        hole = (r < 1.95) & np.isfinite(exact)
        out[f"sign{sign:+d}"] = {
            "max_rel_dev_radial_profile": float(max(dev)),
            "mean_mu_inside_r1.95": float(mu[hole].mean()) if sign == -1 else None,
        }
    out["pass"] = all(
        v["max_rel_dev_radial_profile"] < 0.03 for v in out.values() if isinstance(v, dict)
    )
    return out


def caustic_points(lam1, lam2, sign, n=20000):
    """Critical curve from det J = 0 (polar roots on a radial grid), mapped to the source plane."""
    phi = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r = np.geomspace(0.05, 50, 20000)
    pts = []
    for p in phi:
        c, s = np.cos(p), np.sin(p)
        # J = A - sign / r^2 (I - 2 n n^T)
        a11 = lam1 - sign * (1 - 2 * c * c) / r**2
        a22 = lam2 - sign * (1 - 2 * s * s) / r**2
        a12 = sign * 2 * c * s / r**2
        det = a11 * a22 - a12 * a12
        k = np.nonzero(np.sign(det[:-1]) != np.sign(det[1:]))[0]
        for j in k:
            rr = r[j] - det[j] * (r[j + 1] - r[j]) / (det[j + 1] - det[j])
            tx, ty = rr * c, rr * s
            pts.append((lam1 * tx - sign * tx / rr**2, lam2 * ty - sign * ty / rr**2))
    return np.array(pts)


def gate_g2():
    lam1, lam2, sign = 0.7, 1.3, 1  # kappa 0, gamma 0.3
    t = es.macro_magnification_map(lam1, lam2, sign, 2.0, 2.0, 0.01)
    mu, bx, by = as_image(t)
    mac = 1 / abs(lam1 * lam2)
    pts = caustic_points(lam1, lam2, sign, 4000)
    ix = np.clip(((pts[:, 0] + 2.0) / 0.01).astype(int), 0, mu.shape[1] - 1)
    iy = np.clip(((pts[:, 1] + 2.0) / 0.01).astype(int), 0, mu.shape[0] - 1)
    inwin = (np.abs(pts[:, 0]) < 1.98) & (np.abs(pts[:, 1]) < 1.98)
    local_max = ndimage.maximum_filter(mu, size=5)
    frac_caustic_high = float(np.mean(local_max[iy[inwin], ix[inwin]] > 5 * mac))
    cmask = np.zeros_like(mu, bool)
    cmask[iy[inwin], ix[inwin]] = True
    near = ndimage.binary_dilation(cmask, iterations=2)
    high = mu > 20 * mac
    frac_high_near = float(near[high].mean()) if high.any() else float("nan")
    return {
        "lam1": lam1,
        "lam2": lam2,
        "n_caustic_points": int(inwin.sum()),
        "frac_caustic_points_with_local_mu_gt_5macro": frac_caustic_high,
        "frac_pixels_mu_gt_20macro_within_2px": frac_high_near,
        "pass": frac_caustic_high >= 0.9 and frac_high_near >= 0.9,
    }


def gate_g3():
    out = {}
    for sign in (1, -1):
        mu, _, _ = as_image(es.macro_magnification_map(0.7, 1.3, sign, 3.0, 3.0, 0.01))
        mac = 1 / (0.7 * 1.3)
        out[f"sign{sign:+d}"] = float(mu.mean() / mac - 1)
    out["pass"] = all(abs(v) < 0.02 for k, v in out.items() if k.startswith("sign"))
    out["note"] = "mis-specified as written (PLAN addendum 1); replaced by G3prime"
    return out


def gate_g3prime():
    """Disk integral of mu - 1 on the isolated lens vs pi (R sqrt(R^2 + 4 s) - R^2) (addendum 1)."""
    out = {}
    for sign in (1, -1):
        t = es.macro_magnification_map(1.0, 1.0, sign, 3.0, 3.0, 0.01)
        r = np.hypot(t["beta_x"], t["beta_y"])
        mu = np.asarray(t["mu"])
        for rad in (2.5, 3.0):
            got = float((mu[r < rad] - 1).sum() * 0.01**2)
            pred = float(np.pi * (rad * np.sqrt(rad * rad + 4 * sign) - rad * rad))
            out[f"sign{sign:+d}_R{rad}"] = {"got": got, "pred": pred, "rel": got / pred - 1}
    out["pass"] = all(abs(v["rel"]) < 0.02 for v in out.values())
    return out


# name: (lam1, lam2); lam2 = 1 means kappa = gamma (arc near the tangential critical curve)
CASES = {
    "isolated": (1.0, 1.0),
    "min_mu10": (0.1, 1.0),
    "min_mu100": (0.01, 1.0),
    "sad_mu10": (-0.1, 1.0),
    "sad_mu100": (-0.01, 1.0),
}


def measure(name, lam1, lam2, sign=-1):
    hx = max(3 * 2 * np.sqrt(abs(lam1)), 1.0)
    pix = 0.005 if abs(lam1) < 0.05 else 0.01
    t = es.macro_magnification_map(lam1, lam2, sign, hx, 3.0, pix)
    mu, bx, by = as_image(t)
    mac = 1 / abs(lam1 * lam2)
    rel = mu / mac
    hole = rel < HOLE_FRAC
    lab, _ = ndimage.label(hole)
    cy, cx = mu.shape[0] // 2, mu.shape[1] // 2
    central = lab == lab[cy, cx] if hole[cy, cx] else np.zeros_like(hole)
    pred = es.macro_hole_halfwidths(lam1, lam2)
    row = central[cy]
    col = central[:, cx]
    hw_x = float(row.sum() * pix / 2)
    hw_y = float(col.sum() * pix / 2)
    rho = 0.1
    k = int(round(rho / pix))
    yy, xx = np.mgrid[-k : k + 1, -k : k + 1]
    ker = (xx**2 + yy**2 <= k * k).astype(float)
    ker /= ker.sum()
    disk = ndimage.convolve(rel, ker, mode="nearest")
    tracks = {}
    for u0 in (0.0, 0.5, 1.0):
        iy = int(round((u0 + 3.0) / pix))
        ix = int(round((u0 + hx) / pix))
        tracks[f"along_x_u0_{u0}"] = {
            "min_point": float(rel[min(iy, mu.shape[0] - 1)].min()),
            "min_rho0.1": float(disk[min(iy, mu.shape[0] - 1)].min()),
        }
        if ix < mu.shape[1]:
            tracks[f"along_y_u0_{u0}"] = {
                "min_point": float(rel[:, ix].min()),
                "min_rho0.1": float(disk[:, ix].min()),
            }
    return {
        "lam1": lam1,
        "lam2": lam2,
        "mu_macro": float(1 / (lam1 * lam2)),
        "pixel": pix,
        "window_half_x": hx,
        "hole_area_thetaE2": float(central.sum() * pix * pix),
        "hole_area_all_lt0.1": float(hole.sum() * pix * pix),
        "hole_halfwidth_x": hw_x,
        "hole_halfwidth_y": hw_y,
        "pred_halfwidth_x": pred[0],
        "pred_halfwidth_y": pred[1],
        "tracks": tracks,
    }, (rel, bx, by)


def main():
    t0 = time.time()
    res = {"provenance": "simulated", "plan": "results/w3_macro_shear/PLAN.md"}
    res["G1"] = gate_g1()
    print("G1", res["G1"], round(time.time() - t0), file=sys.stderr)
    res["G2"] = gate_g2()
    print("G2", res["G2"], round(time.time() - t0), file=sys.stderr)
    res["G3"] = gate_g3()
    print("G3", res["G3"], round(time.time() - t0), file=sys.stderr)
    res["G3prime"] = gate_g3prime()
    print("G3prime", res["G3prime"], round(time.time() - t0), file=sys.stderr)
    if res["G1"]["pass"] and res["G2"]["pass"] and res["G3prime"]["pass"]:
        res["cases"] = {}
        fig, axes = plt.subplots(1, len(CASES), figsize=(4 * len(CASES), 4))
        for ax, (name, (l1, l2)) in zip(axes, CASES.items(), strict=True):
            res["cases"][name], (rel, bx, by) = measure(name, l1, l2)
            print(name, res["cases"][name], round(time.time() - t0), file=sys.stderr)
            ax.imshow(
                np.log10(np.clip(rel, 1e-3, 1e2)),
                origin="lower",
                extent=(bx[0], bx[-1], by[0], by[-1]),
                aspect="auto",
                cmap="magma",
                vmin=-3,
                vmax=2,
            )
            ax.set_title(f"{name}: λ1={l1}, λ2={l2}")
            ax.set_xlabel("β_x (θ_E), macro stretch axis")
        axes[0].set_ylabel("β_y (θ_E)")
        fig.suptitle("SIMULATED: log10 μ/μ_macro, negative-mass point lens in a macro model")
        fig.tight_layout()
        fig.savefig(OUT / "maps.png", dpi=70)
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "result.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()

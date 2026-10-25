"""RELICS Lenstool v1 amplitude over the high-magnification region (H_K re-test; D-034 bk)

Frozen plan: ``results/lens_consistency/relics_lenstool_mu/PLAN.md``. Same eight model pairs and
1" sky grid as ``relics_lenstool_batch.py``; the region S is now |mu_ref(z_s = 2)| >= 10, with mu
from the reference model's own deflection Jacobian. Gate T checks the Jacobian against each
reference's published kappa map first. Maps are streamed into memory and discarded.

    python scripts/relics_lenstool_mu.py     # -> summary.json
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import lens_consistency as lc
import numpy as np
import relics_lenstool_batch as rb

OUT = Path("results/lens_consistency/relics_lenstool_mu")
Z_S = 2.0
MU_MIN = 10.0  # PLAN step 2 (ASSUMPTION)
MU_SENS = (5.0, 20.0)  # sensitivity only
EDGE_ARCSEC = 3  # PLAN step 2
MIN_REGION_ARCSEC2 = 100.0  # PLAN step 5
MAX_REGION_FRACTION = 0.5  # gate V sub-region condition
GATE_T_TOL = 0.05  # median |kappa_derived - kappa_pub| (ASSUMPTION)
GATE_T_KAPPA_MIN = 0.3


def dls_ds(name: str, z_s: float = Z_S) -> float:
    from astropy.cosmology import FlatLambdaCDM

    m = lc.MODELS[name]
    h0, om = m["cosmology"]
    cosmo = FlatLambdaCDM(H0=h0, Om0=om)
    return float(
        cosmo.angular_diameter_distance_z1z2(m["z_lens"], z_s)
        / cosmo.angular_diameter_distance(z_s)
    )


def jacobian(ax: np.ndarray, ay: np.ndarray, wcs) -> tuple[np.ndarray, ...]:
    """d alpha / d theta along the pixel axes (arcsec / arcsec); NaN on the outer pixel ring."""
    scale = np.abs(wcs.proj_plane_pixel_scales()[0].to_value("arcsec"))
    dxx = np.gradient(ax, axis=1) / scale
    dxy = np.gradient(ax, axis=0) / scale
    dyx = np.gradient(ay, axis=1) / scale
    dyy = np.gradient(ay, axis=0) / scale
    for a in (dxx, dxy, dyx, dyy):
        a[[0, -1], :] = np.nan
        a[:, [0, -1]] = np.nan
    return dxx, dxy, dyx, dyy


def native_grid_world(shape, wcs):
    yy, xx = np.indices(shape)
    return wcs.pixel_to_world_values(xx.ravel(), yy.ravel())


def load_ref(name: str) -> dict:
    files = lc.MODELS[name]["files"]
    ax, wcs = rb.fetch(*files["alpha_x"])
    ay, _ = rb.fetch(*files["alpha_y"])
    kap, kwcs = rb.fetch(*files["kappa_map"])
    dxx, dxy, dyx, dyy = jacobian(ax, ay, wcs)
    kappa_d = 0.5 * (dxx + dyy)
    ra, dec = native_grid_world(ax.shape, wcs)
    kappa_p = rb.sample(kap, kwcs, ra, dec).reshape(ax.shape)
    sel = np.isfinite(kappa_d) & np.isfinite(kappa_p) & (kappa_p >= GATE_T_KAPPA_MIN)
    dev = float(np.median(np.abs(kappa_d[sel] - kappa_p[sel]))) if sel.any() else float("nan")
    s = dls_ds(name)
    det = (1 - s * dxx) * (1 - s * dyy) - (s * dxy) * (s * dyx)
    with np.errstate(divide="ignore"):
        mu = 1.0 / det
    return {
        "amp": np.hypot(ax, ay),
        "wcs": wcs,
        "mu": mu,
        "gate_t": {
            "median_abs_dkappa": dev,
            "n_pix": int(sel.sum()),
            "passed": bool(np.isfinite(dev) and dev <= GATE_T_TOL),
        },
        "dls_ds": s,
    }


def sample_nearest(arr: np.ndarray, wcs, ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    x, y = wcs.world_to_pixel_values(ra, dec)
    xi, yi = np.rint(x).astype(int), np.rint(y).astype(int)
    inside = (xi >= 0) & (yi >= 0) & (xi < arr.shape[1]) & (yi < arr.shape[0])
    out = np.full(ra.shape, np.nan)
    out[inside] = arr[yi[inside], xi[inside]]
    return out


def region_stats(a_lt, a_ref, mu, ok, mu_min: float) -> dict:
    sel = ok & (np.abs(mu) >= mu_min)
    n = int(sel.sum())
    if n == 0:
        return {"region_arcsec2": 0.0, "R": None}
    q = a_lt[sel] / a_ref[sel]
    return {
        "R": float(np.median(q)),
        "ratio_of_means": float(a_lt[sel].mean() / a_ref[sel].mean()),
        "q16_q84": [float(np.percentile(q, 16)), float(np.percentile(q, 84))],
        "region_arcsec2": float(n * rb.GRID_STEP**2),
        "region_fraction": float(n / ok.sum()),
    }


def measure(cluster: str) -> dict:
    from scipy.ndimage import binary_erosion

    _, lt_name, ref_name = rb.PAIRS[cluster]
    amp_lt, wcs_lt = rb.load_model(lt_name)
    ref = load_ref(ref_name)
    ra, dec = rb.sky_grid(amp_lt, wcs_lt)
    n = int(round(np.sqrt(ra.size)))
    a_lt = rb.sample(amp_lt, wcs_lt, ra, dec)
    a_ref = rb.sample(ref["amp"], ref["wcs"], ra, dec)
    mu = sample_nearest(ref["mu"], ref["wcs"], ra, dec)
    ok = np.isfinite(a_lt) & np.isfinite(a_ref) & np.isfinite(mu)
    ok_core = binary_erosion(
        ok.reshape(n, n), structure=np.ones((2 * EDGE_ARCSEC + 1,) * 2), border_value=0
    ).ravel()
    out = {
        "gate_t": ref["gate_t"],
        "dls_ds_z2": ref["dls_ds"],
        "overlap_core_arcsec2": float(ok_core.sum() * rb.GRID_STEP**2),
        "sensitivity": {str(m): region_stats(a_lt, a_ref, mu, ok_core, m)["R"] for m in MU_SENS},
    }
    out.update(region_stats(a_lt, a_ref, mu, ok_core, MU_MIN))
    if not ref["gate_t"]["passed"]:
        out["class"] = "unmeasured (gate T)"
    elif out["region_arcsec2"] < MIN_REGION_ARCSEC2:
        out["class"] = "unmeasured (|S| < 100 arcsec2)"
    else:
        out["class"] = rb.classify(out["R"])
    return out


def gate_v(rows: dict, factors: dict[str, float]) -> dict:
    out = {}
    for c, r in rows.items():
        if rb.PAIRS[c][0] != "validation":
            continue
        f = factors[c]
        expected = "offset" if f <= 0.8 else "near_1"
        measured = r["R"] is not None and not r["class"].startswith("unmeasured")
        out[c] = {
            "f_obs": f,
            "inv_R": 1.0 / r["R"] if measured else None,
            "expected_class": expected,
            "class_ok": r["class"] == expected,
            "tol_ok": measured and abs(1.0 / r["R"] - f) <= rb.GATE_V_TOL,
            "subregion_ok": measured and r["region_fraction"] <= MAX_REGION_FRACTION,
        }
    passed = all(v["class_ok"] and v["tol_ok"] and v["subregion_ok"] for v in out.values())
    return {"rows": out, "passed": passed}


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    factors = rb.known_factors()
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = dict(zip(rb.PAIRS, pool.map(measure, rb.PAIRS), strict=True))
    for c, r in rows.items():
        r.update(role=rb.PAIRS[c][0], lenstool=rb.PAIRS[c][1], reference=rb.PAIRS[c][2])
    gate = gate_v(rows, factors)
    summary = {
        "plan": str(OUT / "PLAN.md"),
        "rows": rows,
        "gate_v": gate,
        "h_k": rb.h_k(rows, gate),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for c, r in rows.items():
        frac = r.get("region_fraction")
        print(
            f"{c:9s} {r['role']:11s} T {r['gate_t']['median_abs_dkappa']:.3f} "
            f"R {r['R'] if r['R'] is None else round(r['R'], 3)} |S| {r['region_arcsec2']:.0f} "
            f"frac {frac if frac is None else round(frac, 3)} sens {r['sensitivity']} {r['class']}"
        )
    print("gate V:", gate["passed"], "| H_K:", summary["h_k"]["outcome"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""SMACS 0723: counter-images and time delays of AT 2024jfz in the Mahler+22 ICLv2 model.

AT 2024jfz is a point source seen only in the 2024-05-11 GO-4043 mosaics, 1.81" from catalogued
image 6.1 (``results/transients/smacs_arcdiff/RESULT.md``, h000). If its host is a lensed
background source, the model predicts counter-images that repeat the event at the model's time
delays. Plan: ``results/transients/smacs_at2024jfz/PLAN.md`` (frozen before the run).

* The TNS position is traced to the source plane at each z_s of a grid (ASSUMPTION) and re-imaged
  with ``lensmodel.find_images`` (analytic dPIE, best model; no MCMC spread).
* Fermat potential phi = |theta - beta|^2 / 2 - psi; psi differences are line integrals of the
  deflection (alpha = grad psi) along a straight path, checked against an L-shaped path.
* ``testable`` (the frozen rule) says whether a counter-image repeats the event inside a public
  JWST epoch bright enough to be seen.

Outputs (``--out``): ``result.json`` (model_prediction / derived).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import astropy.constants as const
import astropy.units as u
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

MODEL = "smacs0723-iclv2"
TNS_RADEC = (110.83416, -73.451752)  # AT 2024jfz, ../smacs_arcdiff/o4043_o12774/h000_tns.json
T0_MJD = 60441.82  # GO-4043 F115W epoch of the 2024 detection (vet/lightcurves.ecsv)
F0_NJY = 2579.0  # F090W forced flux in 2024 (vet/lightcurves.ecsv, 2578.64 rounded)
#: Public JWST NIRCam epochs covering the position (vet/lightcurves.ecsv, rounded).
EPOCHS_MJD = (59737.0, 60441.8, 61196.9, 61227.7)
Z_GRID = (0.6, 1.0, 1.425, 1.712766, 2.5, 4.0)  # ASSUMPTION; 1.712766 is system 6's model z
WINDOW_DAYS = (-2.0, 30.0)  # ASSUMPTION: the event stays near its 2024 flux for <= 30 observed d
MIN_FLUX_NJY = 10 * 17.0  # ASSUMPTION: 10x the largest blue-band forced error of any epoch
N_STEPS = 2000  # Simpson steps per path segment
PATH_TOL = 0.01  # straight vs L-shaped path: |difference| <= 1 % of |delta psi| (frozen)


def testable(
    delay_days: float,
    mu_ratio: float,
    t0: float = T0_MJD,
    epochs=EPOCHS_MJD,
    f0: float = F0_NJY,
    window=WINDOW_DAYS,
    min_flux: float = MIN_FLUX_NJY,
) -> list[float]:
    """Epochs (MJD) in which a counter-image repeating the event would be seen (frozen rule).

    An epoch qualifies if ``window[0] <= t_e - (t0 + delay) <= window[1]`` and the predicted flux
    ``f0 * |mu_ratio|`` is at least ``min_flux``. An empty list means "not testable".
    """
    if f0 * abs(mu_ratio) < min_flux:
        return []
    arrival = t0 + delay_days
    return [float(t) for t in epochs if window[0] <= t - arrival <= window[1]]


def psi_difference(deflect, p0, p1, n: int = N_STEPS, via=None) -> float:
    """psi(p1) - psi(p0) as the line integral of ``deflect(x, y) -> (ax, ay)`` (Simpson).

    ``via`` adds one corner point (an L-shaped path) for the path-independence check.
    """
    pts = [np.asarray(p0, float)] + ([np.asarray(via, float)] if via is not None else [])
    pts.append(np.asarray(p1, float))
    total = 0.0
    n += n % 2  # Simpson needs an even number of intervals
    w = np.ones(n + 1)
    w[1:-1:2], w[2:-1:2] = 4.0, 2.0
    for a, b in zip(pts[:-1], pts[1:], strict=True):
        t = np.linspace(0.0, 1.0, n + 1)
        x, y = a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])
        ax, ay = deflect(x, y)
        integrand = ax * (b[0] - a[0]) + ay * (b[1] - a[1])
        total += float(np.sum(w * integrand) / (3.0 * n))
    return total


def delay_factor_days(cosmo, z_l: float, z_s: float) -> float:
    """(1 + z_l) D_l D_s / (c D_ls) in days per arcsec^2."""
    dl = cosmo.angular_diameter_distance(z_l)
    ds = cosmo.angular_diameter_distance(z_s)
    dls = cosmo.angular_diameter_distance_z1z2(z_l, z_s)
    return float(
        ((1 + z_l) * dl * ds / (const.c * dls) * (u.arcsec.to(u.rad)) ** 2).to(u.day).value
    )


def solve(model, grid, z_s: float) -> dict:
    x0, y0 = (float(v[0]) for v in model.to_frame([TNS_RADEC[0]], [TNS_RADEC[1]]))
    f = float(model.dls_ds(z_s))

    def deflect(x, y):
        ax, ay = model.deflection_xy(x, y)
        return f * np.asarray(ax), f * np.asarray(ay)

    bx, by = (float(v[0]) for v in model.source_points([x0], [y0], z_s))
    imgs = lensmodel.find_images(model, grid, bx, by, z_s)
    out = {"z_s": z_s, "dls_ds": round(f, 5), "beta_xy": [round(bx, 4), round(by, 4)]}
    if len(imgs) == 0:
        return out | {"n_images": 0, "images": []}
    d0 = np.hypot(imgs["x"] - x0, imgs["y"] - y0)
    i0 = int(np.argmin(d0))
    tfac = delay_factor_days(model.cosmology, model.z_lens, z_s)
    p0 = (float(imgs["x"][i0]), float(imgs["y"][i0]))
    mu0 = float(imgs["magnification"][i0])

    def fermat(p, dpsi):
        return 0.5 * ((p[0] - bx) ** 2 + (p[1] - by) ** 2) - dpsi

    phi0 = fermat(p0, 0.0)
    rows, worst = [], 0.0
    for k in range(len(imgs)):
        p = (float(imgs["x"][k]), float(imgs["y"][k]))
        if k == i0:
            dpsi = 0.0
        else:
            dpsi = psi_difference(deflect, p0, p)
            alt = psi_difference(deflect, p0, p, via=(p[0], p0[1]))
            worst = max(worst, abs(alt - dpsi) / max(abs(dpsi), 1e-9))
        delay = (fermat(p, dpsi) - phi0) * tfac
        ratio = float(imgs["magnification"][k]) / mu0
        rows.append(
            {
                "ra": round(float(imgs["ra"][k]), 6),
                "dec": round(float(imgs["dec"][k]), 6),
                "sep_from_tns_arcsec": round(float(d0[k]), 3),
                "magnification": round(float(imgs["magnification"][k]), 3),
                "mu_ratio": round(ratio, 4),
                "delay_days": round(delay, 1),
                "reference_image": k == i0,
                "testable_epochs_mjd": [] if k == i0 else testable(delay, ratio),
                "next_window_mjd": None
                if k == i0 or F0_NJY * abs(ratio) < MIN_FLUX_NJY
                else [
                    round(T0_MJD + delay - WINDOW_DAYS[1], 1),
                    round(T0_MJD + delay - WINDOW_DAYS[0], 1),
                ],
            }
        )
    tested = any(r["testable_epochs_mjd"] for r in rows)
    return out | {
        "n_images": len(rows),
        "tns_to_reference_arcsec": round(float(d0[i0]), 3),
        "days_per_arcsec2": round(tfac, 3),
        "path_check_max_rel_diff": float(f"{worst:.3g}"),
        "path_check_pass": worst <= PATH_TOL,
        "outcome": "testable" if tested else "not testable",
        "images": rows,
    }


def posthoc(model, grid, out: Path) -> None:
    """Post hoc (not in the plan): a finer z_s grid against every MAST epoch in mast_epochs.json."""
    epochs = json.loads((out / "mast_epochs.json").read_text())["epochs"]
    days = sorted({round(float(e[3]), 1) for e in epochs})
    rows = []
    for z in np.round(np.arange(1.0, 2.51, 0.05), 3):
        r = solve(model, grid, float(z))
        for im in r["images"]:
            if im["reference_image"]:
                continue
            hits = testable(im["delay_days"], im["mu_ratio"], epochs=days)
            rows.append(
                {
                    "z_s": float(z),
                    "ra": im["ra"],
                    "dec": im["dec"],
                    "sep_from_tns_arcsec": im["sep_from_tns_arcsec"],
                    "mu_ratio": im["mu_ratio"],
                    "delay_days": im["delay_days"],
                    "testable_epochs_mjd": hits,
                }
            )
    doc = {
        "provenance": "model_prediction",
        "note": "post hoc (not in PLAN.md): z_s 1.00-2.50 step 0.05; epochs = unique days of "
        "mast_epochs.json (HST and JWST, any instrument); footprint coverage not checked",
        "epochs_mjd": days,
        "rows": rows,
    }
    (out / "posthoc.json").write_text(json.dumps(doc, indent=1) + "\n")
    for row in rows:
        if row["testable_epochs_mjd"]:
            print("posthoc", row)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path("results/transients/smacs_at2024jfz"))
    ap.add_argument("--half-width", type=float, default=60.0)
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--posthoc", action="store_true", help="finer z grid vs mast_epochs.json")
    args = ap.parse_args(argv)
    model, _files, _par = lc.load_model(MODEL)
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, args.half_width, args.step), args.half_width, args.step
    )
    if args.posthoc:
        posthoc(model, grid, args.out)
        return 0
    res = {
        "provenance": "model_prediction",
        "model": MODEL,
        "model_sha256": model.sha256,
        "tns_radec": list(TNS_RADEC),
        "t0_mjd": T0_MJD,
        "f0_njy": F0_NJY,
        "epochs_mjd": list(EPOCHS_MJD),
        "rule": {"window_days": list(WINDOW_DAYS), "min_flux_njy": MIN_FLUX_NJY},
        "runs": [solve(model, grid, z) for z in Z_GRID],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "result.json").write_text(json.dumps(res, indent=1) + "\n")
    for r in res["runs"]:
        print(r["z_s"], r.get("n_images"), r.get("outcome"), r.get("path_check_max_rel_diff"))
        for im in r["images"]:
            print("  ", im)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

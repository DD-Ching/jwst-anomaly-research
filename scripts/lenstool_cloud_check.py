"""Lenstool (conda-forge 8.15.5) vs the repository's dPIE port on the pinned CANUCS MACS0416 model.

Plan: ``results/lens_consistency/lenstool_cloud/PLAN.md`` (gate R, amended before the port run).
For every catalogued image, the port back-traces it and solves for all images on the 0.1" grid;
Lenstool's ``image.all`` (``runmode image 1`` on the same ``best.par``) is the reference.

    python scripts/lenstool_cloud_check.py --image-all <lenstool run dir>/image.all
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import lens_consistency as lc  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import linear_sum_assignment  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

MODEL = "macs0416-canucs"
GRID_STEP = 0.1
MATCH_ARCSEC = 0.1  # ASSUMPTION (plan): the port's solver grid step
FRAC_MIN = 0.95  # ASSUMPTION (plan)
FLAG = (64.048551, -24.064825)  # D115b extra image (#458)
FLAG_RADIUS = 0.7
FLAG_IMAGES = ("D115b.1", "D115b.2")
OUT = Path("results/lens_consistency/lenstool_cloud/result.json")


def read_image_all(path: str | Path) -> tuple[tuple[float, float], dict[str, np.ndarray]]:
    """Lenstool ``image.all``: reference (ra, dec) and, per input image id, its predicted (x, y)."""
    ref = None
    rows: dict[str, list[tuple[float, float]]] = {}
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "#REFERENCE":
            ref = (float(parts[2]), float(parts[3]))
            continue
        if parts[0].startswith("#"):
            continue
        rows.setdefault(parts[0], []).append((float(parts[1]), float(parts[2])))
    if ref is None:
        raise ValueError(f"{path}: no #REFERENCE line")
    return ref, {k: np.asarray(v, float) for k, v in rows.items()}


def match_pairs(ref_xy: np.ndarray, port_xy: np.ndarray) -> np.ndarray:
    """One-to-one (Hungarian) reference-to-port distances; NaN for unmatched reference images."""
    out = np.full(len(ref_xy), np.nan)
    if len(ref_xy) == 0 or len(port_xy) == 0:
        return out
    d = np.hypot(ref_xy[:, None, 0] - port_xy[None, :, 0], ref_xy[:, None, 1] - port_xy[None, :, 1])
    r, c = linear_sum_assignment(d)
    out[r] = d[r, c]
    return out


def gate_r(matched_dist: np.ndarray, flag_lenstool: dict, flag_port: dict) -> bool:
    """Frozen gate R: (a) >= FRAC_MIN of Lenstool images matched within MATCH_ARCSEC, and (b) every
    FLAG_IMAGES id has a flag image in both codes."""
    d = np.asarray(matched_dist, float)
    frac = float(np.mean(np.isfinite(d) & (d <= MATCH_ARCSEC))) if d.size else 0.0
    flags = all(flag_lenstool.get(i, False) and flag_port.get(i, False) for i in FLAG_IMAGES)
    return frac >= FRAC_MIN and flags


_STATE: dict = {}


def _solve(i: int) -> tuple[int, np.ndarray]:
    bt, model, grid = _STATE["bt"], _STATE["model"], _STATE["grid"]
    row = bt[i]
    pred = lensmodel.find_images(
        model, grid, float(row["beta_x"]), float(row["beta_y"]), float(row["z_used"])
    )
    return i, np.column_stack([np.asarray(pred["x"], float), np.asarray(pred["y"], float)])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--image-all", required=True)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args()

    ref, lt = read_image_all(args.image_all)
    files = lc.model_files(MODEL)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    if abs(model.ra0 - ref[0]) > 1e-7 or abs(model.dec0 - ref[1]) > 1e-7:
        raise SystemExit(f"frame mismatch: port {model.ra0, model.dec0} vs Lenstool {ref}")
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    half = lc.grid_half_width(model, images)
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, GRID_STEP), half, GRID_STEP
    )
    bt = lensmodel.backtrace_images(model, images, par["z_m_limit"])
    _STATE.update(bt=bt, model=model, grid=grid)
    with mp.get_context("fork").Pool(args.workers) as pool:
        port = dict(pool.map(_solve, range(len(bt))))

    fx, fy = model.to_frame(np.array([FLAG[0]]), np.array([FLAG[1]]))
    ids = [str(v) for v in bt["image_id"]]
    dists, per_image, missing_in_lenstool = [], {}, []
    n_port_unmatched = 0
    flag_lt, flag_port = {}, {}
    for i, iid in enumerate(ids):
        p = port[i]
        if iid not in lt:
            missing_in_lenstool.append(iid)
            continue
        ref_xy = lt[iid]
        d = match_pairs(ref_xy, p)
        dists.extend(d.tolist())
        n_port_unmatched += max(len(p) - int(np.isfinite(d).sum()), 0)
        per_image[iid] = {
            "n_lenstool": len(ref_xy),
            "n_port": len(p),
            "max_matched_arcsec": float(np.nanmax(d)) if np.isfinite(d).any() else None,
        }
        if iid in FLAG_IMAGES:
            flag_lt[iid] = bool(
                (np.hypot(ref_xy[:, 0] - fx[0], ref_xy[:, 1] - fy[0]) <= FLAG_RADIUS).any()
            )
            flag_port[iid] = bool(
                len(p) and (np.hypot(p[:, 0] - fx[0], p[:, 1] - fy[0]) <= FLAG_RADIUS).any()
            )
    d = np.asarray(dists, float)
    ok = np.isfinite(d)
    out = {
        "plan": "results/lens_consistency/lenstool_cloud/PLAN.md",
        "provenance": "model_prediction",
        "lenstool": "conda-forge lenstool 8.15.5 (linux-64), runmode image 1, pinned best.par",
        "grid_step_arcsec": GRID_STEP,
        "match_arcsec_assumption": MATCH_ARCSEC,
        "frac_min_assumption": FRAC_MIN,
        "n_catalogued_images": len(ids),
        "n_images_missing_in_image_all": len(missing_in_lenstool),
        "images_missing_in_image_all": missing_in_lenstool,
        "n_lenstool_predicted": int(d.size),
        "n_matched_within": int((ok & (d <= MATCH_ARCSEC)).sum()),
        "frac_matched_within": float(np.mean(ok & (d <= MATCH_ARCSEC))) if d.size else 0.0,
        "n_lenstool_unmatched": int((~ok).sum()),
        "n_port_unmatched": int(n_port_unmatched),
        "matched_dist_median_arcsec": float(np.median(d[ok])) if ok.any() else None,
        "matched_dist_p95_arcsec": float(np.percentile(d[ok], 95)) if ok.any() else None,
        "matched_dist_max_arcsec": float(d[ok].max()) if ok.any() else None,
        "flag_image_lenstool": flag_lt,
        "flag_image_port": flag_port,
        "gate_R": gate_r(d, flag_lt, flag_port),
        "per_image": per_image,
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "per_image"}, indent=1))


if __name__ == "__main__":
    main()

"""RXC J0600 581 under the RELICS GLAFIC v3 ``range/`` realisations (D-024 addendum 18).

Plan (frozen before the run): ``results/arcs/anti-shear-error-581-glafic/PLAN.md``. The
99 realisation maps (600 x 600, 0.4" px, pixel-unit deflections) are streamed into memory,
hashed and dropped; deflection in arcsec = pixel value x the map's own header pixel scale. The
route is validated on the best fit (V1 unit convention, V2 resolution, V3 centring) and by the
κ gate. Rules are addendum 13's functions, imported unchanged from ``arcs_shear_error.py``.
Predicted axes ``model_prediction``; classes and labels ``derived``.

    REQUESTS_CA_BUNDLE=<proxy + certifi bundle> python scripts/arcs_shear_error_581g.py --cache DIR
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import lensmodel

_DIR = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("arcs_shear_error", _DIR / "arcs_shear_error.py")
ase = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ase)
lc = ase.lc

OUT = ase.ROOT / "anti-shear-error-581-glafic"
BASE = ase.RELICS + "rxc0600m20/models/glafic/v3/"
BEST = BASE + "hlsp_relics_model_model_rxc0600m20_glafic_v3_{}.fits"
RANGE = BASE + "range/hlsp_relics_model_model_rxc0600m20_glafic-map{:03d}_v3_{}.fits"
ROW = 581
Z_LENS = ase.RXC["z_lens"]
DECIDING = ase.RXC["run"]  # the stored GLAFIC v3 run (best fit, 0.1" maps)

# frozen constants (ASSUMPTIONS, PLAN.md); deviations from addendum 13: 99 maps, gate map 098
N_RANGE = 99
RANGE_GATE_MAPS = (0, 49, 98)
V1_MAX_ARCSEC = 1e-3
V2_MAX_OFFSET_DEG = 2.0
V3_MAX_REL = 0.10


# ---- route functions (unit-tested offline) ----
def pixel_scale_arcsec(wcs: WCS) -> float:
    return float(abs(wcs.pixel_scale_matrix[1, 1]) * 3600.0)


def resample(data: np.ndarray, src: WCS, dst: WCS, shape: tuple[int, int]) -> np.ndarray:
    """Bilinear values of ``data`` (grid ``src``) at the pixel centres of grid ``dst``."""
    from scipy.ndimage import map_coordinates

    jj, ii = np.mgrid[0 : shape[0], 0 : shape[1]]
    ra, dec = dst.pixel_to_world_values(ii.ravel(), jj.ravel())
    x, y = src.world_to_pixel_values(ra, dec)
    out = map_coordinates(data, [y, x], order=1, mode="constant", cval=np.nan)
    return out.reshape(shape).astype(np.float32)


def v1_gate(max_abs: float) -> bool:
    return bool(np.isfinite(max_abs) and max_abs <= V1_MAX_ARCSEC)


def v2_gate(stored: list[str], new: list[str], d_offset) -> bool:
    d = np.asarray(d_offset, float)
    return bool(
        stored == new and len(d) > 0 and np.all(np.isfinite(d)) and d.max() <= V2_MAX_OFFSET_DEG
    )


def v3_stat(best_xy: np.ndarray, real_xy: np.ndarray) -> float:
    """Median over rows of |median-over-realisations α − best α| / |best α|.

    ``best_xy`` (n_rows, 2); ``real_xy`` (n_real, n_rows, 2)."""
    med = np.median(np.asarray(real_xy, float), axis=0)
    best = np.asarray(best_xy, float)
    return float(np.median(np.hypot(*(med - best).T) / np.hypot(*best.T)))


def v3_gate(stat: float) -> bool:
    return bool(np.isfinite(stat) and stat <= V3_MAX_REL)


# ---- I/O ----
def get(url: str) -> tuple[np.ndarray, WCS, fits.Header, str]:
    """Whole small FITS file into memory (retry on 429 / 5xx / connection errors)."""
    import time

    import requests

    for attempt in range(5):
        try:
            r = requests.get(url, timeout=180)
            if r.status_code == 200:
                body = r.content
                with fits.open(io.BytesIO(body)) as h:
                    hdr = h[0].header.copy()
                    data = np.asarray(h[0].data, np.float32)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", FITSFixedWarning)
                    wcs = WCS(hdr)
                return data, wcs, hdr, hashlib.sha256(body).hexdigest()
            if r.status_code != 429 and r.status_code < 500:
                r.raise_for_status()
        except (requests.ConnectionError, requests.Timeout):
            pass
        time.sleep(2.0 * 2**attempt)
    raise RuntimeError(f"{url}: failed after 5 attempts")


def model_from(ax, ay, wcs, centre, tag) -> lensmodel.MapLensModel:
    cosmo = FlatLambdaCDM(H0=ase.COSMO[0], Om0=ase.COSMO[1])
    return lensmodel.MapLensModel(ax, ay, wcs, Z_LENS, cosmo, source=tag, centre=centre)


def alpha_at(model, rows: Table) -> np.ndarray:
    x, y = model.to_frame(np.asarray(rows["ra"], float), np.asarray(rows["dec"], float))
    return np.stack(model.deflection_xy(x, y), axis=1)


def kappa_check(model, url: str, cache: Path) -> dict:
    data, wcs, _, sha = get(url)
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / url.rsplit("/", 1)[1]
    fits.PrimaryHDU(data, header=wcs.to_header()).writeto(path, overwrite=True)
    try:
        check = lc.map_check(model, path, "kappa", None)
    finally:
        path.unlink()
    check["map"] = path.name
    check["sha256"] = sha
    return check


def one_realisation(i: int, centre, row, ctrl, bg, ref_hdr, cache: Path) -> dict:
    ax, wcs, hx, sx = get(RANGE.format(i, "x-pixels-deflect"))
    ay, _, hy, sy = get(RANGE.format(i, "y-pixels-deflect"))
    for key in ("NAXIS1", "NAXIS2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD2_2"):
        for h in (hx, hy):
            if not np.isclose(float(h[key]), float(ref_hdr[key]), rtol=1e-9):
                raise ValueError(f"map{i:03d}: grid differs from map000 ({key})")
    s = pixel_scale_arcsec(wcs)
    model = model_from(ax * s, ay * s, wcs, centre, f"glafic-map{i:03d}_v3")
    rec = ase._rec(ase.reclass(model, row))
    rec.update(map=i, sha256_x=sx, sha256_y=sy)
    rec["control_retention"] = ase.retention(list(ase.reclass(model, ctrl)["orientation_class"]))
    rec["alpha_bg"] = alpha_at(model, bg).tolist()
    if i in RANGE_GATE_MAPS:
        check = kappa_check(model, RANGE.format(i, "kappa"), cache)
        rec["kappa_gate"] = {"kappa_check": check, "passes": ase.map_gate(check["median_rel_diff"])}
    return rec


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)

    row, ctrl = ase.split_rows(DECIDING, ROW)
    stored = Table.read(ase.ROOT / DECIDING)
    bg = stored[stored["strong_shear"] & (stored["population"] == "background")]

    # best fit: V1 unit convention; centre = best-fit reference pixel (as the stored run)
    bpx, bwcs, _, _ = get(BEST.format("x-pixels-deflect"))
    bpy, _, _, _ = get(BEST.format("y-pixels-deflect"))
    bax, _, _, sha_bx = get(BEST.format("x-arcsec-deflect"))
    bay, _, _, sha_by = get(BEST.format("y-arcsec-deflect"))
    sb = pixel_scale_arcsec(bwcs)
    v1_max = float(max(np.nanmax(np.abs(bpx * sb - bax)), np.nanmax(np.abs(bpy * sb - bay))))
    del bpx, bpy
    ref = (float(bwcs.wcs.crpix[0]) - 1.0, float(bwcs.wcs.crpix[1]) - 1.0)
    centre = tuple(float(v) for v in bwcs.pixel_to_world_values(*ref))
    best = model_from(bax, bay, bwcs, centre, "glafic v3 best fit")
    best_alpha = alpha_at(best, bg)

    # V2: best fit resampled onto the realisation grid
    _, rwcs, ref_hdr, _ = get(RANGE.format(0, "x-pixels-deflect"))
    shape = (int(ref_hdr["NAXIS2"]), int(ref_hdr["NAXIS1"]))
    coarse = model_from(
        resample(bax, bwcs, rwcs, shape),
        resample(bay, bwcs, rwcs, shape),
        rwcs,
        centre,
        "glafic v3 best fit resampled to 0.4 arcsec",
    )
    del bax, bay, best
    re = ase.reclass(coarse, bg)
    d_off = np.abs(np.asarray(re["offset_min"], float) - np.asarray(bg["offset_min"], float))
    v2 = v2_gate(list(bg["orientation_class"]), list(re["orientation_class"]), d_off)
    validation = {
        "best_pixel_scale_arcsec": sb,
        "range_pixel_scale_arcsec": pixel_scale_arcsec(rwcs),
        "V1_max_abs_arcsec": v1_max,
        "V1_passes": v1_gate(v1_max),
        "V2_n_rows": len(bg),
        "V2_classes_equal": list(bg["orientation_class"]) == list(re["orientation_class"]),
        "V2_max_offset_diff_deg": round(float(d_off.max()), 4),
        "V2_passes": v2,
        # post hoc (added after the first run; diagnostic, no rule change)
        "V2_rows_posthoc": [
            {
                "label": int(a),
                "stored_offset_min": round(float(b), 2),
                "resampled_offset_min": round(float(c), 2),
                "stored_g_min": round(float(d), 3),
                "resampled_g_min": round(float(e), 3),
            }
            for a, b, c, d, e in zip(
                bg["label"],
                bg["offset_min"],
                re["offset_min"],
                bg["g_min"],
                re["g_min"],
                strict=True,
            )
        ],
        "resampled_best_581": ase._rec(ase.reclass(coarse, row)),
        "sha256_best_arcsec": [sha_bx, sha_by],
    }
    print(json.dumps(validation, indent=1), flush=True)

    with ThreadPoolExecutor(args.workers) as pool:
        real = list(
            pool.map(
                lambda i: one_realisation(i, centre, row, ctrl, bg, ref_hdr, args.cache),
                range(N_RANGE),
            )
        )
    gates = {f"map{r['map']:03d}": r.pop("kappa_gate") for r in real if "kappa_gate" in r}
    rets = [r.pop("control_retention") for r in real]
    alphas = np.array([r.pop("alpha_bg") for r in real])
    v3 = v3_stat(best_alpha, alphas)
    validation.update(V3_median_rel=round(v3, 4), V3_passes=v3_gate(v3))
    r_med = float(np.median(rets))
    f = sum(r["class"] == "anti" for r in real) / len(real)
    valid = (
        validation["V1_passes"]
        and validation["V2_passes"]
        and validation["V3_passes"]
        and all(g["passes"] for g in gates.values())
        and len(gates) == len(RANGE_GATE_MAPS)
        and ase.control_gate(r_med)
    )
    label = ase.range_label(f, valid)
    labels = {"RXC J0600 581": label}

    def pct(k):
        return [round(float(v), 3) for v in np.percentile([r[k] for r in real], (16, 50, 84))]

    summary = {
        "plan": "results/arcs/anti-shear-error-581-glafic/PLAN.md",
        "params": {
            "kappa_max_rel": ase.KAPPA_MAX_REL,
            "control_min": ase.CONTROL_MIN,
            "robust_f": ase.ROBUST_F,
            "sensitive_f": ase.SENSITIVE_F,
            "range_gate_maps": list(RANGE_GATE_MAPS),
            "n_range": N_RANGE,
            "v1_max_arcsec": V1_MAX_ARCSEC,
            "v2_max_offset_deg": V2_MAX_OFFSET_DEG,
            "v3_max_rel": V3_MAX_REL,
            "cosmology": list(ase.COSMO),
            "z_lens": Z_LENS,
            "orientation_table": "lens_consistency defaults (aligned 30, anti 60, z margin 0.1)",
        },
        "provenance": {"axes": "model_prediction", "classes": "derived"},
        "validation": validation,
        "label": ROW,
        "deciding_run_row": ase._rec(row),
        "map_gates": gates,
        "control_n": len(ctrl),
        "control_retention_median": round(r_med, 3),
        "control_retention_range": [round(min(rets), 3), round(max(rets), 3)],
        "control_gate": ase.control_gate(r_med),
        "f_anti": round(f, 4),
        "f_aligned": round(sum(r["class"] == "aligned" for r in real) / len(real), 4),
        "label_rule3": label,
        "offset_min_pct": pct("offset_min"),
        "g_min_pct": pct("g_min"),
        "realisations": real,
        "labels": labels,
        "outcome": ase.outcome(labels),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in ("f_anti", "control_retention_median", "map_gates", "labels", "outcome")
            },
            default=str,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

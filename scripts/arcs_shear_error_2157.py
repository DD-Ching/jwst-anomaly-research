"""MACS0308 2157 under the RELICS Zitrin LTM-Gauss ``range/`` realisations (D-024 addendum 14).

Plan (frozen before the run): ``results/arcs/anti-shear-error-2157/PLAN.md``. Each 36 MB map is read
as one HTTP byte range covering the rows needed (no file is stored), sliced in columns, and turned
into a :class:`lensmodel.MapLensModel` with the header WCS sliced to the strip. The route is first
validated on the best fit (pixel × scale = arcsec maps; the stored Zitrin run's classes are
reproduced). Rules are addendum 13's functions, imported unchanged from ``arcs_shear_error.py``.
Predicted axes ``model_prediction``; classes and labels ``derived``.

    REQUESTS_CA_BUNDLE=<proxy bundle + certifi roots> python scripts/arcs_shear_error_2157.py --cache DIR
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import warnings
from collections.abc import Callable
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

OUT = ase.ROOT / "anti-shear-error-2157"
BASE = "https://archive.stsci.edu/hlsps/relics/macs0308p26/models/zitrin-ltm-gauss/"
BEST = BASE + "hlsp_relics_model_model_macs0308p26_zitrin-ltm-gauss_v1_{}.fits"
RANGE = BASE + "range/hlsp_relics_model_model_macs0308p26_zitrin-ltm-gauss-map{:03d}_v1_{}.fits"
ROW = 2157
Z_LENS = 0.356
DECIDING = "macs0308-relics/macs0308-relics-glafic/arcs.ecsv"
STORED_ZITRIN = "macs0308-relics/macs0308-relics-zitrin/arcs.ecsv"

# frozen thresholds (ASSUMPTIONS, PLAN.md)
MARGIN_PX = 16
V1_MAX_ARCSEC = 1e-3
V2_MAX_OFFSET_DEG = 0.5
HEADER_BYTES = 2880 * 4  # enough for these headers; checked for END below


# ---- strip reader (unit-tested offline) ----
def parse_header(head: bytes) -> tuple[fits.Header, int]:
    """FITS primary header and the byte offset of its data from the first header blocks."""
    text = head.decode("ascii", errors="replace")
    for k in range(0, len(text), 80):
        if text[k : k + 8].rstrip() == "END":
            offset = (k // 2880 + 1) * 2880
            return fits.Header.fromstring(text[: k + 80]), offset
    raise ValueError("END card not found in the bytes read")


def read_strip(get_range: Callable[[int, int], bytes], j0: int, j1: int) -> tuple[np.ndarray, WCS]:
    """Rows j0..j1-1 (0-based) of a 2-D BITPIX -32 primary image and the WCS sliced to them.

    ``get_range(start, end)`` returns bytes start..end inclusive (HTTP Range semantics)."""
    hdr, offset = parse_header(get_range(0, HEADER_BYTES - 1))
    if hdr["BITPIX"] != -32 or hdr["NAXIS"] != 2:
        raise ValueError("expected a 2-D float32 image")
    nx, ny = int(hdr["NAXIS1"]), int(hdr["NAXIS2"])
    j0, j1 = max(0, j0), min(ny, j1)
    row = 4 * nx
    raw = get_range(offset + j0 * row, offset + j1 * row - 1)
    if len(raw) != (j1 - j0) * row:
        raise ValueError(f"short read: {len(raw)} of {(j1 - j0) * row} bytes")
    data = np.frombuffer(raw, ">f4").reshape(j1 - j0, nx).astype(np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(hdr)
    return data, wcs.slice((slice(j0, j1), slice(None)))


def cut_columns(data: np.ndarray, wcs: WCS, i0: int, i1: int) -> tuple[np.ndarray, WCS]:
    i0, i1 = max(0, i0), min(data.shape[1], i1)
    return data[:, i0:i1], wcs.slice((slice(None), slice(i0, i1)))


def pixel_scale_arcsec(wcs: WCS) -> float:
    return float(abs(wcs.pixel_scale_matrix[1, 1]) * 3600.0)


def v1_gate(max_abs: float) -> bool:
    return bool(np.isfinite(max_abs) and max_abs <= V1_MAX_ARCSEC)


def v2_gate(stored: list[str], new: list[str], d_offset: np.ndarray) -> bool:
    d = np.asarray(d_offset, float)
    return bool(stored == new and len(d) > 0 and np.all(np.isfinite(d)) and d.max() <= V2_MAX_OFFSET_DEG)


# ---- I/O ----
class Http:
    def __init__(self, url: str):
        import requests

        self.url, self.session, self.digest = url, requests.Session(), hashlib.sha256()

    def __call__(self, start: int, end: int) -> bytes:
        import requests

        for attempt in range(4):
            try:
                r = self.session.get(self.url, headers={"Range": f"bytes={start}-{end}"}, timeout=180)
                r.raise_for_status()
                if r.status_code != 206:
                    raise RuntimeError(f"{self.url}: range not honoured ({r.status_code})")
                self.digest.update(r.content)
                return r.content
            except requests.RequestException:
                if attempt == 3:
                    raise
        raise RuntimeError("unreachable")


def bounds(tables: list[Table], wcs_full: WCS) -> tuple[int, int, int, int]:
    ra = np.concatenate([np.asarray(t["ra"], float) for t in tables])
    dec = np.concatenate([np.asarray(t["dec"], float) for t in tables])
    i, j = wcs_full.world_to_pixel_values(ra, dec)
    return (
        int(np.floor(j.min())) - MARGIN_PX,
        int(np.ceil(j.max())) + MARGIN_PX + 1,
        int(np.floor(i.min())) - MARGIN_PX,
        int(np.ceil(i.max())) + MARGIN_PX + 1,
    )


def strip(url: str, box) -> tuple[np.ndarray, WCS, str]:
    h = Http(url)
    data, wcs = read_strip(h, box[0], box[1])
    data, wcs = cut_columns(data, wcs, box[2], box[3])
    return data, wcs, h.digest.hexdigest()


def model_from(ax: np.ndarray, ay: np.ndarray, wcs: WCS, centre, tag: str) -> lensmodel.MapLensModel:
    return lensmodel.MapLensModel(
        ax, ay, wcs, Z_LENS, FlatLambdaCDM(H0=ase.COSMO[0], Om0=ase.COSMO[1]), source=tag, centre=centre
    )


def kappa_check(model, url: str, box, cache: Path) -> dict:
    data, wcs, _ = strip(url, box)
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / url.rsplit("/", 1)[1]
    fits.PrimaryHDU(data, header=wcs.to_header()).writeto(path, overwrite=True)
    try:
        check = lc.map_check(model, path, "kappa", None)
    finally:
        path.unlink()
    check["map"] = path.name
    return check


def one_realisation(i: int, box, centre, row: Table, ctrl: Table, cache: Path) -> dict:
    ax, wcs, sx = strip(RANGE.format(i, "x-pixels-deflect"), box)
    ay, _, sy = strip(RANGE.format(i, "y-pixels-deflect"), box)
    s = pixel_scale_arcsec(wcs)
    model = model_from(ax * s, ay * s, wcs, centre, f"zitrin-ltm-gauss-map{i:03d}")
    rec = ase._rec(ase.reclass(model, row))
    rec.update(map=i, sha256_x_strip=sx, sha256_y_strip=sy)
    rec["control_retention"] = ase.retention(list(ase.reclass(model, ctrl)["orientation_class"]))
    if i in ase.RANGE_GATE_MAPS:
        check = kappa_check(model, RANGE.format(i, "kappa"), box, cache)
        rec["kappa_gate"] = {"kappa_check": check, "passes": ase.map_gate(check["median_rel_diff"])}
    return rec


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)

    row, ctrl = ase.split_rows(DECIDING, ROW)
    stored = Table.read(ase.ROOT / STORED_ZITRIN)
    stored_bg = stored[stored["strong_shear"] & (stored["population"] == "background")]
    hdr, _ = parse_header(Http(RANGE.format(0, "x-pixels-deflect"))(0, HEADER_BYTES - 1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wfull = WCS(hdr)
    centre = tuple(float(v) for v in wfull.pixel_to_world_values(wfull.wcs.crpix[0] - 1, wfull.wcs.crpix[1] - 1))
    box = bounds([row, ctrl, stored_bg], wfull)

    # validation on the best fit
    bx, bwcs, _ = strip(BEST.format("x-pixels-deflect"), box)
    by, _, _ = strip(BEST.format("y-pixels-deflect"), box)
    axs, _, sha_ax = strip(BEST.format("x-arcsec-deflect"), box)
    ays, _, sha_ay = strip(BEST.format("y-arcsec-deflect"), box)
    s = pixel_scale_arcsec(bwcs)
    v1_max = float(max(np.nanmax(np.abs(bx * s - axs)), np.nanmax(np.abs(by * s - ays))))
    best = model_from(axs, ays, bwcs, centre, "zitrin-ltm-gauss best fit (strip)")
    re = ase.reclass(best, stored_bg)
    d_off = np.abs(np.asarray(re["offset_min"], float) - np.asarray(stored_bg["offset_min"], float))
    v2 = v2_gate(list(stored_bg["orientation_class"]), list(re["orientation_class"]), d_off)
    validation = {
        "strip_box_j0_j1_i0_i1": list(box),
        "pixel_scale_arcsec": s,
        "V1_max_abs_arcsec": v1_max,
        "V1_passes": v1_gate(v1_max),
        "V2_n_rows": len(stored_bg),
        "V2_classes_equal": list(stored_bg["orientation_class"]) == list(re["orientation_class"]),
        "V2_max_offset_diff_deg": round(float(d_off.max()), 4),
        "V2_passes": v2,
        "best_fit_2157": ase._rec(ase.reclass(best, row)),
        "sha256_best_arcsec_strips": [sha_ax, sha_ay],
    }
    print(json.dumps(validation, indent=1), flush=True)
    del bx, by, axs, ays, best

    with ThreadPoolExecutor(args.workers) as pool:
        real = list(pool.map(lambda i: one_realisation(i, box, centre, row, ctrl, args.cache), range(ase.N_RANGE)))
    gates = {f"map{r['map']:03d}": r.pop("kappa_gate") for r in real if "kappa_gate" in r}
    rets = [r.pop("control_retention") for r in real]
    r_med = float(np.median(rets))
    f = sum(r["class"] == "anti" for r in real) / len(real)
    valid = (
        validation["V1_passes"]
        and validation["V2_passes"]
        and all(g["passes"] for g in gates.values())
        and ase.control_gate(r_med)
    )
    label = ase.range_label(f, valid)
    labels = {"MACS0308 2157": label}

    def pct(k):
        return [round(float(v), 3) for v in np.percentile([r[k] for r in real], (16, 50, 84))]

    summary = {
        "plan": "results/arcs/anti-shear-error-2157/PLAN.md",
        "params": {
            "kappa_max_rel": ase.KAPPA_MAX_REL,
            "control_min": ase.CONTROL_MIN,
            "robust_f": ase.ROBUST_F,
            "sensitive_f": ase.SENSITIVE_F,
            "range_gate_maps": list(ase.RANGE_GATE_MAPS),
            "n_range": ase.N_RANGE,
            "margin_px": MARGIN_PX,
            "v1_max_arcsec": V1_MAX_ARCSEC,
            "v2_max_offset_deg": V2_MAX_OFFSET_DEG,
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
        "f_anti": f,
        "f_aligned": sum(r["class"] == "aligned" for r in real) / len(real),
        "label_rule3": label,
        "offset_min_pct": pct("offset_min"),
        "g_min_pct": pct("g_min"),
        "realisations": real,
        "labels": labels,
        "outcome": ase.outcome(labels),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("f_anti", "control_retention_median", "labels", "outcome")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())

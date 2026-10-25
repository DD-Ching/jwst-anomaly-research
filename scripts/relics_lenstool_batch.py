"""Does the RELICS Lenstool v1 amplitude offset follow the 2017-12 map batch? (H_K; D-034 bj)

Frozen plan: ``results/lens_consistency/relics_lenstool_batch/PLAN.md``. For eight clusters the
Lenstool v1 deflection field and a second model's are sampled on a common 1" sky grid; R is the
median |alpha_LT| / |alpha_ref| over the strong-lensing region |alpha_ref| >= 0.5 max. Gate V checks
R against the stored rule-7 factors of four clusters before the four test clusters are read. Maps
are streamed into memory, checked against their registered sha256 and discarded.

    python scripts/relics_lenstool_batch.py     # -> summary.json
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import lens_consistency as lc
import numpy as np

OUT = Path("results/lens_consistency/relics_lenstool_batch")
GRID_STEP = 1.0  # arcsec (PLAN step 1)
REGION_FRACTION = 0.5  # |alpha_ref| >= 0.5 max (PLAN step 2, ASSUMPTION)
OFFSET_MIN = 1.25  # PLAN step 4 classes (ASSUMPTION)
NEAR1 = (0.85, 1.15)
GATE_V_TOL = 0.12  # |1/R - f_obs| (PLAN gate V, ASSUMPTION)

#: cluster -> (role, Lenstool model, second model)
PAIRS: dict[str, tuple[str, str, str]] = {
    "macs0553": ("validation", "macs0553-relics-lenstool", "macs0553-relics-glafic"),
    "a3192": ("validation", "a3192-relics-lenstool", "a3192-relics-glafic"),
    "rxc2211": ("validation", "rxc2211-relics-lenstool", "rxc2211-relics-glafic"),
    "cl0152": ("validation", "cl0152-relics-lenstool", "cl0152-relics-glafic"),
    "rxc0600": ("test_batch", "rxc0600-relics-lenstool", "rxc0600-relics-glafic"),
    "macs0417": ("test_batch", "macs0417-relics-lenstool", "macs0417-canucs-maps"),
    "whl0137": ("test_other", "whl0137-relics-lenstool", "whl0137-relics-glafic"),
    "spt0615": ("test_other", "spt0615-relics-lenstool", "spt0615-relics-glafic"),
}


def known_factors() -> dict[str, float]:
    """Stored rule-7 medians of the validation clusters (derived; addendum bi summary)."""
    rows = json.loads(
        (Path("results/lens_consistency/relics_lenstool_zref/summary.json")).read_text()
    )
    return {c: float(r["f_obs"]) for c, r in rows["h_u"]["rows"].items()}


def fetch(url: str, sha256: str):
    """Whole FITS file into memory, hash-checked -> (2-D float array, celestial WCS)."""
    import requests
    from astropy.io import fits
    from astropy.wcs import WCS

    for attempt in range(4):
        try:
            r = requests.get(url, timeout=300)
            r.raise_for_status()
            break
        except requests.RequestException:
            if attempt == 3:
                raise
    raw = r.content
    if hashlib.sha256(raw).hexdigest() != sha256:
        raise RuntimeError(f"{url}: sha256 mismatch")
    if url.endswith(".gz"):
        raw = gzip.decompress(raw)
    with fits.open(io.BytesIO(raw)) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        return np.asarray(hdu.data, dtype=np.float64), WCS(hdu.header).celestial


def load_model(name: str):
    files = lc.MODELS[name]["files"]
    ax, wcs = fetch(*files["alpha_x"])
    ay, _ = fetch(*files["alpha_y"])
    if ax.shape != ay.shape:
        raise RuntimeError(f"{name}: alpha_x / alpha_y shapes differ")
    return np.hypot(ax, ay), wcs


def sample(amp: np.ndarray, wcs, ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """Bilinear |alpha| at sky positions; NaN outside the map."""
    from scipy.ndimage import map_coordinates

    x, y = wcs.world_to_pixel_values(ra, dec)
    inside = (x >= 0) & (y >= 0) & (x <= amp.shape[1] - 1) & (y <= amp.shape[0] - 1)
    out = np.full(ra.shape, np.nan)
    out[inside] = map_coordinates(amp, [y[inside], x[inside]], order=1)
    return out


def sky_grid(amp: np.ndarray, wcs, step: float = GRID_STEP):
    """1" tangent-plane grid covering the first map's footprint -> (ra, dec) arrays."""
    ny, nx = amp.shape
    ra0, dec0 = wcs.pixel_to_world_values((nx - 1) / 2, (ny - 1) / 2)
    cx, cy = wcs.pixel_to_world_values(
        np.array([0, nx - 1, 0, nx - 1]), np.array([0, 0, ny - 1, ny - 1])
    )
    dx = (cx - ra0) * np.cos(np.radians(dec0)) * 3600.0
    dy = (cy - dec0) * 3600.0
    half = float(max(np.abs(dx).max(), np.abs(dy).max()))
    t = np.arange(-half, half + step, step)
    gx, gy = np.meshgrid(t, t)
    dec = dec0 + gy / 3600.0
    ra = ra0 + gx / 3600.0 / np.cos(np.radians(dec0))
    return ra.ravel(), dec.ravel()


def ratio_stats(a_lt: np.ndarray, a_ref: np.ndarray, step: float = GRID_STEP) -> dict:
    """PLAN steps 2-3 on co-sampled magnitudes (NaN = outside a footprint)."""
    ok = np.isfinite(a_lt) & np.isfinite(a_ref)
    lt, ref = a_lt[ok], a_ref[ok]
    sel = ref >= REGION_FRACTION * ref.max()
    q = lt[sel] / ref[sel]
    return {
        "R": float(np.median(q)),
        "ratio_of_means": float(lt[sel].mean() / ref[sel].mean()),
        "q16_q84": [float(np.percentile(q, 16)), float(np.percentile(q, 84))],
        "region_arcsec2": float(sel.sum() * step**2),
        "overlap_arcsec2": float(ok.sum() * step**2),
        "region_fraction": float(sel.sum() / ok.sum()),
        "max_alpha_ref": float(ref.max()),
    }


def classify(r: float) -> str:
    if r >= OFFSET_MIN:
        return "offset"
    if NEAR1[0] <= r <= NEAR1[1]:
        return "near_1"
    return "intermediate"


def measure(cluster: str) -> dict:
    _, lt_name, ref_name = PAIRS[cluster]
    amp_lt, wcs_lt = load_model(lt_name)
    amp_ref, wcs_ref = load_model(ref_name)
    ra, dec = sky_grid(amp_lt, wcs_lt)
    stats = ratio_stats(sample(amp_lt, wcs_lt, ra, dec), sample(amp_ref, wcs_ref, ra, dec))
    stats["class"] = classify(stats["R"])
    return stats


def gate_v(rows: dict, factors: dict[str, float]) -> dict:
    out = {}
    for c, r in rows.items():
        if PAIRS[c][0] != "validation":
            continue
        f = factors[c]
        expected = "offset" if f <= 0.8 else "near_1"
        out[c] = {
            "f_obs": f,
            "inv_R": 1.0 / r["R"],
            "expected_class": expected,
            "class_ok": r["class"] == expected,
            "tol_ok": abs(1.0 / r["R"] - f) <= GATE_V_TOL,
        }
    return {"rows": out, "passed": all(v["class_ok"] and v["tol_ok"] for v in out.values())}


def h_k(rows: dict, gate: dict) -> dict:
    if not gate["passed"]:
        return {"outcome": "undetermined (gate V failed)"}
    batch = [rows[c]["class"] for c in rows if PAIRS[c][0] == "test_batch"]
    other = [rows[c]["class"] for c in rows if PAIRS[c][0] == "test_other"]
    if "near_1" in batch or "offset" in other:
        outcome = "rejected"
    elif all(b == "offset" for b in batch) and all(o == "near_1" for o in other):
        outcome = "supported"
    else:
        outcome = "undetermined (intermediate class)"
    return {"batch_classes": batch, "other_classes": other, "outcome": outcome}


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    factors = known_factors()
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = dict(zip(PAIRS, pool.map(measure, PAIRS), strict=True))
    for c, r in rows.items():
        r.update(role=PAIRS[c][0], lenstool=PAIRS[c][1], reference=PAIRS[c][2])
    gate = gate_v(rows, factors)
    summary = {"plan": str(OUT / "PLAN.md"), "rows": rows, "gate_v": gate, "h_k": h_k(rows, gate)}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for c, r in rows.items():
        print(
            f"{c:9s} {r['role']:11s} R {r['R']:.3f} "
            f"[{r['q16_q84'][0]:.2f}, {r['q16_q84'][1]:.2f}] {r['class']}"
        )
    print("gate V:", gate["passed"], "| H_K:", summary["h_k"]["outcome"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

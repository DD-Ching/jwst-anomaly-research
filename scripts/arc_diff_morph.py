"""Residual morphology of ``arc_diff.py`` hits: edges, streaks and PSF rosettes (D-027 addendum).

Outer MACS0416 tiles (2026-10-15) showed artifact classes that the bright-core flag misses (and
residuals at the mosaic edge):

- **streak**: a diffraction spike of a bright star, rotated between the epochs, leaves a long thin
  residual of one sign; hits are strung along it;
- **rosette**: a compact source (star, compact core) fainter than the bright-core level still leaves
  alternating-sign lobes around itself after PSF matching.

A point transient changes one PSF-sized spot. For each hit and band this measures, on the
significance map S of an ``arc_diff.py --maps`` file:

- ``extent_arcsec``: the length of the connected region of sign(S_hit) × S' ≥ ``--streak-level``
  that contains the hit, where S' is S smoothed by a Gaussian of ``--smooth-px`` (it joins the
  broken segments of a spike); the length is √(12 λ_max) of the region's pixel second moments
  (the full length of a uniform line);
- ``n_lobes``, ``both_signs``: connected regions of |S| ≥ ``--level`` (at least ``--min-lobe-px``
  pixels) within ``--lobe-arcsec`` of the hit, and whether both signs occur among them;
- ``edge_arcsec``: the distance to the nearest pixel without data (mosaic edge or a gap), or to the
  map border.

A hit is at the **edge** when ``edge_arcsec`` < ``--edge-arcsec`` in any band; **extended** when
``extent_arcsec`` > ``--max-extent-arcsec`` in the band of the largest |S_hit| (spike streaks, and
the smoothed lobes of bright rosettes); a **rosette** when ``n_lobes`` ≥ ``--min-lobes`` with both
signs in at least ``--min-rosette-bands`` bands; otherwise it passes. Every number is ``derived``;
the thresholds are ASSUMPTIONs. Validation: every Yan+2023 arc-transient position in the 20″
Warhol box passes in o211/o212; one (D21-S2, outside its discovery pair) is a rosette in o211/o213.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage


def _wcs(text: str) -> WCS:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(fits.Header.fromstring(text))


def component_extent(mask: np.ndarray, y: int, x: int) -> float:
    """√(12 λ_max) in pixels of the connected region of ``mask`` that contains (y, x); 0 if none."""
    if not mask[y, x]:
        return 0.0
    lab, _ = ndimage.label(mask, structure=np.ones((3, 3)))
    yy, xx = np.nonzero(lab == lab[y, x])
    if yy.size < 2:
        return 0.0
    cov = np.cov(np.vstack([xx, yy]).astype(float))
    return float(np.sqrt(12.0 * max(np.linalg.eigvalsh(cov)[-1], 0.0)))


def lobes(s: np.ndarray, y: float, x: float, r_px: float, level: float, min_px: int) -> list[int]:
    """Signs (+1/-1) of the regions of |S| ≥ level with ≥ min_px pixels touching the disk r_px."""
    yy, xx = np.mgrid[: s.shape[0], : s.shape[1]]
    disk = np.hypot(yy - y, xx - x) <= r_px
    out = []
    for sign in (1, -1):
        lab, n = ndimage.label(np.nan_to_num(sign * s) >= level, structure=np.ones((3, 3)))
        if not n:
            continue
        sizes = np.bincount(lab.ravel())
        for k in np.unique(lab[disk & (lab > 0)]):
            if sizes[k] >= min_px:
                out.append(sign)
    return out


def morph_band(maps, band: str, ra: float, dec: float, a) -> dict | None:
    """One band's morphology at (ra, dec); None when off the map or the S window is all NaN."""
    w = _wcs(str(maps[f"{band}_wcs"]))
    x0, y0 = (float(v) for v in w.all_world2pix(ra, dec, 0))
    scale = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600
    s = maps[f"{band}_s"]
    yi, xi = int(round(y0)), int(round(x0))
    r = max(1, int(round(a.peak_arcsec / scale)))
    if not (r <= yi < s.shape[0] - r and r <= xi < s.shape[1] - r):
        return None
    win = s[yi - r : yi + r + 1, xi - r : xi + r + 1]
    if not np.isfinite(win).any():
        return None
    k = int(np.nanargmax(np.abs(win)))
    py, px = yi - r + k // win.shape[1], xi - r + k % win.shape[1]
    s_hit = float(s[py, px])
    sgn = 1.0 if s_hit >= 0 else -1.0
    good = np.isfinite(s)
    sm = ndimage.gaussian_filter(np.where(good, s, 0.0), a.smooth_px)
    norm = ndimage.gaussian_filter(good.astype(float), a.smooth_px)
    sm = np.where(norm > 0.5, sm / np.maximum(norm, 1e-9), 0.0)
    # Smoothing lowers the noise of S by ≈ 1 / (2 √π σ) per pixel; rescale so the level is in σ.
    sm *= 2.0 * np.sqrt(np.pi) * a.smooth_px if a.smooth_px > 0 else 1.0
    ext = component_extent(sgn * sm >= a.streak_level, py, px) * scale
    pad = np.pad(good, 1, constant_values=False)
    edge = float(ndimage.distance_transform_edt(pad)[py + 1, px + 1]) * scale
    lb = lobes(s, py, px, a.lobe_arcsec / scale, a.level, a.min_lobe_px)
    return {
        "s_hit": s_hit,
        "extent_arcsec": ext,
        "n_lobes": len(lb),
        "both_signs": bool(1 in lb and -1 in lb),
        "edge_arcsec": edge,
    }


def classify(res: dict[str, dict], a) -> tuple[str, str | None]:
    """(verdict, top band) from the band of the largest |S_hit|."""
    if not res:
        return "not_covered", None
    top = max(res, key=lambda b: abs(res[b]["s_hit"]))
    r = res[top]
    if min(v["edge_arcsec"] for v in res.values()) < a.edge_arcsec:
        return "edge", top
    if r["extent_arcsec"] > a.max_extent_arcsec:
        return "extended", top
    n_ros = sum(v["n_lobes"] >= a.min_lobes and v["both_signs"] for v in res.values())
    if n_ros >= a.min_rosette_bands:
        return "rosette", top
    return "pass", top


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", type=Path, required=True, help="arc_diff.py --maps file")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--hits", type=Path, help="arc_diff.py hits.ecsv")
    src.add_argument("--known", type=Path, help="CSV with name, ra, dec (validation targets)")
    ap.add_argument("--frame-offset", type=float, nargs=2, default=[0.0, 0.0])
    ap.add_argument("--all-hits", action="store_true", help="include near_bright_core hits")
    ap.add_argument("--peak-arcsec", type=float, default=0.06)
    ap.add_argument("--level", type=float, default=3.0)
    ap.add_argument("--smooth-px", type=float, default=2.0)
    ap.add_argument("--streak-level", type=float, default=3.0)
    ap.add_argument("--max-extent-arcsec", type=float, default=3.0)
    ap.add_argument("--edge-arcsec", type=float, default=0.5)
    ap.add_argument("--lobe-arcsec", type=float, default=0.8)
    ap.add_argument("--min-lobe-px", type=int, default=4)
    ap.add_argument("--min-lobes", type=int, default=4)
    ap.add_argument("--min-rosette-bands", type=int, default=2)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.hits:
        t = Table.read(a.hits)
        if not a.all_hits and "near_bright_core" in t.colnames:
            t = t[~t["near_bright_core"]]
        targets = [(str(r["hit_id"]), float(r["ra"]), float(r["dec"])) for r in t]
    else:
        from arc_diff import load_known

        targets = [
            (k["name"], k["ra"], k["dec"]) for k in load_known(a.known, tuple(a.frame_offset))
        ]
    rows = []
    with np.load(a.maps) as m:
        bands = sorted({k.rsplit("_", 1)[0] for k in m.files if k.endswith("_wcs")})
        for name, ra, dec in targets:
            res = {b: morph_band(m, b, ra, dec, a) for b in bands}
            res = {b: r for b, r in res.items() if r is not None}
            verdict, top = classify(res, a)
            row = {"id": name, "ra": ra, "dec": dec, "verdict": verdict, "top_band": top}
            rows.append(row | {"bands": res})
    counts = {v: sum(r["verdict"] == v for r in rows) for v in sorted({r["verdict"] for r in rows})}
    out = {
        "provenance": "derived",
        "maps": a.maps.name,
        "thresholds_assumption": {
            k: getattr(a, k)
            for k in (
                "peak_arcsec",
                "level",
                "smooth_px",
                "streak_level",
                "max_extent_arcsec",
                "edge_arcsec",
                "lobe_arcsec",
                "min_lobe_px",
                "min_lobes",
                "min_rosette_bands",
            )
        },
        "counts": counts,
        "targets": rows,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(counts))


if __name__ == "__main__":
    main()

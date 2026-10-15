"""Compact-knot photometric residuals of an ``arc_diff.py --maps`` box (MACS0416 h044/h030 vetting).

For every local peak of the epoch-1 image (footprint ``--footprint-arcsec``) whose epoch-1
aperture flux is at least ``--min-njy``, the script measures the difference flux (``d`` = epoch 2
minus epoch 1 after matching) and the epoch-1 flux in the ``arc_diff_vet.py`` apertures and reports
the fractional residual ``fd / fe``. Targets (``--target name ra dec``) are measured the same way at
their positions and ranked against the knots of comparable flux (``--min-njy`` up to ``--max-njy``,
excluding knots within ``--exclude-arcsec`` of the target) in robust σ
(1.4826 × MAD around the knot median). An epoch-wide compact-core deficit (PSF, focus or kernel
mismatch) shows as a non-zero knot median; a transient shows as a target far outside the knot
distribution. Every number is ``derived``; the thresholds are ASSUMPTIONs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from arc_diff_vet import AP_RADIUS, NJY_PER_MJY, SR_PER_ARCSEC2, _wcs, disk_sum
from scipy import ndimage


def knot_peaks(e1: np.ndarray, size: int, pad: int) -> list[tuple[int, int]]:
    """Local maxima of ``e1`` (footprint ``size`` px) at least ``pad`` px from the edge."""
    finite = np.where(np.isfinite(e1), e1, -np.inf)
    peaks = (finite == ndimage.maximum_filter(finite, size=size)) & np.isfinite(e1)
    peaks[:pad], peaks[-pad:], peaks[:, :pad], peaks[:, -pad:] = False, False, False, False
    return [(int(y), int(x)) for y, x in zip(*np.nonzero(peaks), strict=True)]


def robust_z(value: float, ref: np.ndarray) -> tuple[float, float, float]:
    """(median, robust σ, z of ``value``) of the finite reference values; NaN when < 5 remain."""
    ref = np.asarray(ref, dtype=float)
    ref = ref[np.isfinite(ref)]
    if len(ref) < 5:
        return float("nan"), float("nan"), float("nan")
    med = float(np.median(ref))
    sig = 1.4826 * float(np.median(np.abs(ref - med)))
    return med, sig, (value - med) / sig if sig > 0 else float("nan")


def measure_band(maps, band: str, targets: list[tuple[str, float, float]], a) -> dict:
    w = _wcs(str(maps[f"{band}_wcs"]))
    scale = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600
    d, e1 = maps[f"{band}_d"], maps[f"{band}_e1img"]
    r_ap = AP_RADIUS[band] / scale
    to_njy = scale**2 * SR_PER_ARCSEC2 * NJY_PER_MJY
    fd, fe = disk_sum(d, r_ap) * to_njy, disk_sum(e1, r_ap) * to_njy
    # disk_sum counts NaN as 0, so an aperture touching a masked pixel is biased: flag it
    holed = disk_sum((~(np.isfinite(d) & np.isfinite(e1))).astype(float), r_ap) > 0
    pad = int(np.ceil(r_ap)) + 3
    size = 2 * int(round(a.footprint_arcsec / scale)) + 1
    knots, kpix = [], []
    for y, x in knot_peaks(e1, size, pad):
        if fe[y, x] >= a.min_njy and not holed[y, x]:
            kpix.append((y, x))
            ra, dec = (float(v) for v in w.all_pix2world(x, y, 0))
            knots.append(
                {
                    "ra": round(ra, 6),
                    "dec": round(dec, 6),
                    "e1_njy": round(float(fe[y, x]), 1),
                    "diff_njy": round(float(fd[y, x]), 2),
                }
            )
    e = np.array([k["e1_njy"] for k in knots])
    frac = np.array([k["diff_njy"] for k in knots]) / np.where(e > 0, e, np.nan)
    ok = np.isfinite(frac)
    comparable = frac[(e <= a.max_njy) & ok]
    bright = frac[(e > a.max_njy) & np.isfinite(frac)]
    out = {
        "n_knots": len(knots),
        "median_frac_comparable": float(np.median(comparable)) if len(comparable) else None,
        "n_comparable": int(len(comparable)),
        "median_frac_bright": float(np.median(bright)) if len(bright) else None,
        "n_bright": int(len(bright)),
        "targets": {},
        "knots": knots,
    }
    for name, ra, dec in targets:
        x0, y0 = (float(v) for v in w.all_world2pix(ra, dec, 0))
        yi, xi = int(round(y0)), int(round(x0))
        if not (pad <= yi < d.shape[0] - pad and pad <= xi < d.shape[1] - pad):
            out["targets"][name] = {"skipped": "outside the map interior"}
            continue
        f, f1 = float(fd[yi, xi]), float(fe[yi, xi])
        tf = f / f1 if f1 > 0 else float("nan")
        # The target's own knot (and close neighbours) are not part of its reference sample.
        far = np.array([np.hypot(y - y0, x - x0) * scale > a.exclude_arcsec for y, x in kpix])
        ref = frac[(e <= a.max_njy) & ok & far] if len(kpix) else comparable
        med, sig, z = robust_z(tf, ref)
        out["targets"][name] = {
            "diff_njy": round(f, 2),
            "e1_njy": round(f1, 1),
            "frac": round(tf, 4),
            "knot_median": round(med, 4),
            "knot_robust_sigma": round(sig, 4),
            "z_vs_knots": round(z, 2),
            "n_reference": int(len(ref)),
            "rank_from_top": int(np.sum(ref >= tf)) + 1 if np.isfinite(tf) else None,
            "aperture_has_nan": bool(holed[yi, xi]),
        }
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("maps", type=Path, help="arc_diff.py --maps .npz")
    ap.add_argument("--pair", required=True, help="label, e.g. o211_c004")
    ap.add_argument("--bands", nargs="+", default=["F150W", "F200W", "F277W", "F356W"])
    ap.add_argument("--target", nargs=3, action="append", default=[], metavar=("NAME", "RA", "DEC"))
    ap.add_argument("--footprint-arcsec", type=float, default=0.25)
    ap.add_argument("--min-njy", type=float, default=100.0)
    ap.add_argument("--max-njy", type=float, default=500.0, help="upper end of 'comparable' knots")
    ap.add_argument(
        "--exclude-arcsec",
        type=float,
        default=0.2,
        help="knots this close to a target are not its reference",
    )
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    targets = [(n, float(r), float(d)) for n, r, d in a.target]
    with np.load(a.maps) as maps:
        bands = {b: measure_band(maps, b, targets, a) for b in a.bands if f"{b}_d" in maps}
    summary = {
        "provenance": "derived",
        "pair": a.pair,
        "thresholds_assumption": {
            "footprint_arcsec": a.footprint_arcsec,
            "min_njy": a.min_njy,
            "max_njy_comparable": a.max_njy,
            "exclude_arcsec": a.exclude_arcsec,
            "aperture_radius_arcsec": AP_RADIUS,
        },
        "bands": bands,
    }
    a.out.write_text(json.dumps(summary, indent=1) + "\n")
    for b, r in bands.items():
        print(
            f"{b}: knots {r['n_knots']}, median frac comparable {r['median_frac_comparable']} "
            f"(n={r['n_comparable']}), bright {r['median_frac_bright']} (n={r['n_bright']})"
        )
        for n, t in r["targets"].items():
            if "skipped" in t:
                print(f"  {n}: {t['skipped']}")
                continue
            print(
                f"  {n}: frac {t['frac']:+.4f} z {t['z_vs_knots']:+.2f} rank {t['rank_from_top']}"
            )


if __name__ == "__main__":
    main()

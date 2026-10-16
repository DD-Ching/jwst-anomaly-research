"""Image-centroid check of the WCS that the forced-photometry tools use (``arc_diff.read_box``).

``reg_offsets.py`` compares the pipeline catalogues, whose ``sky_centroid`` comes from each
mosaic's full gwcs. ``hit_lightcurve.py`` instead reads a box with the FITS header WCS
(``read_box``). This script reads one box per mosaic, measures an iterated flux-weighted
centroid (radius ``R_ARCSEC``) at each control's catalogue position, converts it to the sky
with the box WCS, and reports the offset from the catalogue position of the same mosaic. A
large offset would mean that the tool, not the mosaic, is misregistered. All numbers are
``derived``; the radius and iteration count are ASSUMPTIONs.

    python scripts/reg_imgcheck.py --controls results/.../controls_1176_vs_2561_f356w.ecsv \\
        --mosaic 1176=s3://.../jw01176-o221_t010_nircam_clear-f356w_i2d.fits \\
        --mosaic 2561=s3://.../jw02561-o001_t003_nircam_clear-f356w_i2d.fits \\
        --center 3.5525 -30.3770 --half-arcsec 30 --out results/.../imgcheck_f356w.ecsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arc_diff import read_box  # noqa: E402

R_ARCSEC = 0.15  # ASSUMPTION: centroid radius (≈ PSF core in the LW bands)
N_ITER = 5  # ASSUMPTION: centroid iterations


def centroid(sci: np.ndarray, x: float, y: float, r_px: float) -> tuple[float, float]:
    """Iterated flux-weighted centroid in a circle of radius ``r_px`` (background: local median)."""
    for _ in range(N_ITER):
        x0, y0 = int(np.floor(x - 3 * r_px)), int(np.floor(y - 3 * r_px))
        x1, y1 = int(np.ceil(x + 3 * r_px)) + 1, int(np.ceil(y + 3 * r_px)) + 1
        if x0 < 0 or y0 < 0 or x1 > sci.shape[1] or y1 > sci.shape[0]:
            return float("nan"), float("nan")
        sub = sci[y0:y1, x0:x1]
        yy, xx = np.mgrid[y0:y1, x0:x1]
        d2 = (xx - x) ** 2 + (yy - y) ** 2
        bg = np.nanmedian(sub[d2 > (2 * r_px) ** 2])
        w = np.where(d2 <= r_px**2, sub - bg, 0.0)
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        if w.sum() <= 0:
            return float("nan"), float("nan")
        x, y = float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())
    return x, y


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--controls", type=Path, required=True)
    ap.add_argument("--mosaic", action="append", required=True, help="label=uri; first is 'cat'")
    ap.add_argument("--center", nargs=2, type=float, required=True)
    ap.add_argument("--half-arcsec", type=float, default=30.0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    t = Table.read(args.controls)
    ra0, dec0 = args.center
    cosd = np.cos(np.radians(dec0))
    near = (np.abs((t["ra"] - ra0) * cosd * 3600) < args.half_arcsec - 2) & (
        np.abs((t["dec"] - dec0) * 3600) < args.half_arcsec - 2
    )
    t = t[near]
    out = Table()
    out["ref_id"], out["ra"], out["dec"] = t["ref_id"], t["ra"], t["dec"]
    for k, spec in enumerate(args.mosaic):
        label, uri = spec.split("=", 1)
        box = read_box(uri, ra0, dec0, args.half_arcsec)
        # catalogue position in this mosaic: the reference position, plus (dx, dy) for the first
        cra = t["ra"] + (t["dx"] / 3600.0 / cosd if k == 0 else 0.0)
        cdec = t["dec"] + (t["dy"] / 3600.0 if k == 0 else 0.0)
        r_px = R_ARCSEC / box["scale"]
        ddx, ddy = [], []
        for a, d in zip(cra, cdec, strict=True):
            x, y = (float(v) for v in box["wcs"].all_world2pix(a, d, 0))
            cx, cy = centroid(box["sci"], x, y, r_px)
            if not np.isfinite(cx):
                ddx.append(np.nan)
                ddy.append(np.nan)
                continue
            sa, sd = (float(v) for v in box["wcs"].all_pix2world(cx, cy, 0))
            ddx.append((sa - a) * cosd * 3600)
            ddy.append((sd - d) * 3600)
        out[f"dx_{label}"], out[f"dy_{label}"] = ddx, ddy
        print(
            label,
            "n",
            int(np.isfinite(ddx).sum()),
            "median image−catalogue (arcsec)",
            round(float(np.nanmedian(ddx)), 4),
            round(float(np.nanmedian(ddy)), 4),
            "max",
            round(float(np.nanmax(np.hypot(ddx, ddy))), 3),
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.write(args.out, overwrite=True)


if __name__ == "__main__":
    main()

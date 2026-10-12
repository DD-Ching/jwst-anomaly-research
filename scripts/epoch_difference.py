"""Two-epoch difference of one source from the cutouts that ``transient_forced.py`` writes.

Epoch 2 is resampled onto epoch 1's pixel grid through the two WCSs (cubic spline, NaN outside).
Each epoch's local background, the median in a sky annulus around the target, is subtracted.
There is no PSF matching, so use the same filter and instrument in both epochs. The script reports:

* the flux-weighted centroid of epoch 1 and of the difference (positive pixels within
  ``--centroid-radius`` of the target; ASSUMPTION: 0.3″), and their separation;
* cumulative fluxes of epoch 1 and of the difference within radii around the epoch-1 centroid,
  and their ratio (how concentrated the change is).

Every number is ``derived``. ``--png`` writes epoch 1, epoch 2 and the difference side by side.

    python scripts/epoch_difference.py --ra 24.333856 --dec -8.426836 \\
        --epoch1 <out>/cutouts/F150W_e1/<stem>/c0000.fits \\
        --epoch2 <out>/cutouts/F150W_e2/<stem>/c0000.fits --png n0153_f150w.png
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage

RADII_ARCSEC = (0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
SKY_ANNULUS_ARCSEC = (0.6, 0.95)  # ASSUMPTION: outside the source, inside a 2″ cutout


def read_cutout(path: str | Path) -> tuple[np.ndarray, WCS]:
    """``(data, wcs)`` of the first image extension (SCI) of a cutout FITS file."""
    with fits.open(path) as hdul:
        hdu = next(h for h in hdul if h.data is not None)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            return np.asarray(hdu.data, dtype=float), WCS(hdu.header)


def _centroid(
    img: np.ndarray, mask: np.ndarray, xx: np.ndarray, yy: np.ndarray
) -> tuple[float, float]:
    w = np.nan_to_num(np.clip(img, 0.0, None)) * mask
    total = w.sum()
    if total <= 0:
        return float("nan"), float("nan")
    return float((w * xx).sum() / total), float((w * yy).sum() / total)


def difference(
    e1: np.ndarray,
    w1: WCS,
    e2: np.ndarray,
    w2: WCS,
    ra: float,
    dec: float,
    centroid_radius_arcsec: float = 0.3,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Measure the epoch-2 minus epoch-1 difference at (ra, dec); returns (results, images)."""
    ny, nx = e1.shape
    yy, xx = np.mgrid[:ny, :nx].astype(float)
    sky = w1.pixel_to_world(xx, yy)
    x2, y2 = w2.world_to_pixel(sky)
    e2r = ndimage.map_coordinates(e2, [y2, x2], order=3, mode="constant", cval=np.nan)
    scale = float(w1.proj_plane_pixel_scales()[0].to_value("arcsec"))

    x0, y0 = (float(v) for v in w1.world_to_pixel_values(ra, dec))
    r_target = np.hypot(xx - x0, yy - y0) * scale
    annulus = (r_target > SKY_ANNULUS_ARCSEC[0]) & (r_target < SKY_ANNULUS_ARCSEC[1])
    a = e1 - np.nanmedian(e1[annulus])
    b = e2r - np.nanmedian(e2r[annulus])
    d = b - a

    near = r_target < centroid_radius_arcsec
    c1 = _centroid(a, near, xx, yy)
    cd = _centroid(d, near, xx, yy)
    rr = np.hypot(xx - c1[0], yy - c1[1]) * scale
    cum1 = [float(np.nansum(a[rr < r])) for r in RADII_ARCSEC]
    cumd = [float(np.nansum(d[rr < r])) for r in RADII_ARCSEC]
    res = {
        "pixel_scale_arcsec": scale,
        "centroid1_offset_from_target_arcsec": float(np.hypot(c1[0] - x0, c1[1] - y0) * scale),
        "diff_centroid_offset_from_centroid1_arcsec": float(
            np.hypot(cd[0] - c1[0], cd[1] - c1[1]) * scale
        ),
        "radii_arcsec": list(RADII_ARCSEC),
        "cumflux_epoch1": cum1,
        "cumflux_diff": cumd,
        "diff_over_epoch1": [x / y if y else float("nan") for x, y in zip(cumd, cum1, strict=True)],
        "provenance": "derived",
        "assumptions": {
            "centroid_radius_arcsec": centroid_radius_arcsec,
            "sky_annulus_arcsec": list(SKY_ANNULUS_ARCSEC),
            "psf_matching": "none",
        },
    }
    return res, {"epoch1": a, "epoch2": b, "diff": d, "centroid1": np.array(c1)}


def _png(images: dict[str, np.ndarray], path: Path, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    v = float(np.nanpercentile(images["epoch1"], 99.5))
    fig, ax = plt.subplots(1, 3, figsize=(9, 3.2))
    for k, (key, lo, hi) in enumerate(
        (("epoch1", -0.1, 1.0), ("epoch2", -0.1, 1.0), ("diff", -0.3, 0.6))
    ):
        ax[k].imshow(images[key], origin="lower", cmap="gray", vmin=lo * v, vmax=hi * v)
        ax[k].plot(*images["centroid1"], "r+")
        ax[k].set_title(key, fontsize=9)
        ax[k].set_xticks([])
        ax[k].set_yticks([])
    fig.suptitle(title, fontsize=9)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--epoch1", type=Path, required=True, help="cutout FITS, epoch 1")
    ap.add_argument("--epoch2", type=Path, required=True, help="cutout FITS, epoch 2")
    ap.add_argument("--centroid-radius", type=float, default=0.3, help="arcsec (ASSUMPTION)")
    ap.add_argument("--png", type=Path)
    args = ap.parse_args(argv)
    e1, w1 = read_cutout(args.epoch1)
    e2, w2 = read_cutout(args.epoch2)
    res, images = difference(e1, w1, e2, w2, args.ra, args.dec, args.centroid_radius)
    res["inputs"] = {"epoch1": str(args.epoch1), "epoch2": str(args.epoch2)}
    json.dump(res, sys.stdout, indent=1)
    print()
    if args.png:
        _png(images, args.png, f"({args.ra:.6f}, {args.dec:.6f}); red + = epoch-1 centroid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

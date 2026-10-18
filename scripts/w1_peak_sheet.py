"""Contact sheet of W1 shear-screen peaks: a cutout of the F200W ``_i2d.fits`` (S3 section read,
not stored) around each summary's ``s_max_radec``, with the aperture and the used shape sticks.

    python scripts/w1_peak_sheet.py --summaries runs/z0.1/goodsn-o003/summary.json ... \
        --out results/w1_shear_goodsn/sheet_smax_flags.png

Pixels ``observed``; S and the peak ``derived``. Display only: no number is decided here.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import astropy.units as u  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.nddata import Cutout2D  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from astropy.wcs.utils import proj_plane_pixel_scales  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import inject_shear as ish  # noqa: E402
import w1_streak_mask as wsm  # noqa: E402

from jwst_anomaly import paths  # noqa: E402

HALF_ARCSEC = 20.0  # cutout half-width (ASSUMPTION, display only)


def cutout(uri: str, ra: float, dec: float, half: float = HALF_ARCSEC) -> Cutout2D:
    """Section of the SCI extension around (ra, dec), read remotely (trimmed at the edges)."""
    with fits.open(uri, use_fsspec=True, fsspec_kwargs={"anon": True}) as hdul:
        hdu = hdul["SCI"]
        pos = SkyCoord(ra, dec, unit="deg")
        return Cutout2D(hdu.section, pos, 2 * half * u.arcsec, wcs=WCS(hdu.header))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--summaries", nargs="+", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    n = len(args.summaries)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5.4), squeeze=False)
    for ax, path in zip(axes[0], args.summaries, strict=True):
        s = json.loads(path.read_text())
        ra, dec = s["screen"]["s_max_radec"]
        rel = ish.field_spec(s["field"])["catalog"]  # the stored path is machine-specific
        c = cutout(wsm.i2d_uri(rel), ra, dec)
        img = np.nan_to_num(np.asarray(c.data, float))
        lo, hi = np.percentile(img, [5, 99.5])
        ax.imshow(img, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        cat = Table.read(paths.data_root() / rel)
        sky = cat["sky_centroid"]
        x, y = c.wcs.world_to_pixel(sky)
        inside = (x > 0) & (y > 0) & (x < img.shape[1]) & (y < img.shape[0])
        ax.plot(x[inside], y[inside], "c+", ms=4, alpha=0.6)
        px, py = c.wcs.world_to_pixel(SkyCoord(ra, dec, unit="deg"))
        scale = float(np.mean(proj_plane_pixel_scales(c.wcs))) * 3600.0  # cutouts trim at edges
        r_ap = float(s["assumptions"]["screen"]["aperture_arcsec"])
        ax.add_patch(plt.Circle((px, py), r_ap / scale, fill=False, color="r"))
        ax.plot(px, py, "rx")
        sc = s["screen"]
        ax.set_title(
            f"{s['field']}: S {sc['s_max']:.2f}, p {sc['p_random_max']:.3f}, "
            f'max|S_x| {sc["b_max"]:.2f}\n{r_ap:g}" aperture (red), catalogue sources (cyan)',
            fontsize=9,
        )
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=90)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

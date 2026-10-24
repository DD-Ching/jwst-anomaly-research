"""Contact sheet of the Abell S1063 extra-image targets (``s1063_extra``) on GO 3293 F200W.

Reads ``results/lens_consistency/s1063_extra/predicted_abells1063.ecsv`` and cuts 6-arcsec stamps
from the pipeline i2d SCI extension by S3 range reads; writes ``sheet.png`` (display only, no
number is derived from it).

    python scripts/s1063_extra_sheet.py
"""

from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import s1063_extra as sx  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from astropy.wcs.utils import proj_plane_pixel_scales  # noqa: E402

HALF_ARCSEC = 3.0


def main() -> int:
    import s3fs

    t = Table.read(sx.OUT / f"predicted_{sx.FIELD}.ecsv")
    fs = s3fs.S3FileSystem(anon=True, default_block_size=1 << 20)
    key = sx.S3_3293.format(sx.OBS, sx.TILE, "i2d.fits")[5:]
    n = len(t)
    fig, axes = plt.subplots(1, n, figsize=(2.2 * n, 2.6), squeeze=False)
    warnings.simplefilter("ignore")
    with fs.open(key) as fh, fits.open(fh, lazy_load_hdus=True) as hl:
        h = hl["SCI"].header
        w = WCS(h)
        pix = float(proj_plane_pixel_scales(w.celestial).mean() * 3600.0)
        half = int(round(HALF_ARCSEC / pix))
        for ax, r in zip(axes[0], t, strict=True):
            x, y = (int(round(float(v))) for v in w.world_to_pixel_values(r["ra"], r["dec"]))
            d = np.asarray(hl["SCI"].section[y - half : y + half, x - half : x + half], float)
            lo, hi = np.nanpercentile(d, [5, 99.5])
            ax.imshow(d, origin="lower", cmap="gray_r", vmin=lo, vmax=hi)
            ax.add_patch(plt.Circle((half, half), r["search_arcsec"] / pix, fill=False, color="r"))
            ax.set_title(f"{r['system']} {r['image_class']}\n{r['class_local']}", fontsize=7)
            ax.set_xticks([])
            ax.set_yticks([])
    title = "Abell S1063 WSLAP+ (WLS) predicted targets, GO 3293 F200W, 6-arcsec stamps"
    fig.suptitle(title, fontsize=8)
    fig.tight_layout()
    fig.savefig(sx.OUT / "sheet.png", dpi=110)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

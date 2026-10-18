"""Contact sheet of the SPT0615 `arcs` strong-shear background rows of one model run, anti first.

4" F200W / F444W cutouts by S3 byte range (program 4212 o001 i2d); blue the observed axis
(``observed``), red the z = 2 predicted axis (``model_prediction``).

    python scripts/arcs_spt0615_sheet.py spt0615-relics-glafic \\
        results/arcs/spt0615-relics/background_glafic.png
"""

import sys

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS

plt.switch_backend("Agg")
model, png = sys.argv[1], sys.argv[2]
t = Table.read(f"results/arcs/spt0615-relics/{model}/arcs.ecsv")
b = t[t["strong_shear"] & (t["population"] == "background")]
order = np.argsort(b["orientation_class"] != "anti")
b = b[order]
url = "s3://stpubdata/jwst/public/jw04212/L3/t/o001/jw04212-o001_t001_nircam_clear-{}_i2d.fits"
bands = ["f200w", "f444w"]
fig, axes = plt.subplots(len(bands), len(b), figsize=(1.9 * len(b), 4.2), squeeze=False)
for j, band in enumerate(bands):
    with fits.open(url.format(band), use_fsspec=True, fsspec_kwargs={"anon": True}) as h:
        w = WCS(h["SCI"].header)
        for i, r in enumerate(b):
            c = SkyCoord(r["ra"], r["dec"], unit="deg")
            co = Cutout2D(h["SCI"].section, c, 4 * u.arcsec, wcs=w)
            d = np.nan_to_num(co.data)
            lo, hi = np.percentile(d, [5, 99.5])
            ax = axes[j, i]
            ax.imshow(np.arcsinh((d - lo) / max(hi - lo, 1e-9) * 10), origin="lower", cmap="gray_r")
            # PA east of north -> pixel axes via WCS offsets
            x0, y0 = co.wcs.world_to_pixel(c)
            for pa, col in ((r["pa_obs"], "tab:blue"), (r["pa_pred_z2"], "tab:red")):
                p2 = c.directional_offset_by(pa * u.deg, 1.2 * u.arcsec)
                x1, y1 = co.wcs.world_to_pixel(p2)
                ax.plot([2 * x0 - x1, x1], [2 * y0 - y1, y1], color=col, lw=0.8)
            ax.set_xticks([])
            ax.set_yticks([])
            if j == 0:
                ax.set_title(
                    f"{r['label']} {r['orientation_class']}\n"
                    f"m{r['mag']:.1f} g{r['g_min']:.2f} z{r['z_phot']:.1f}",
                    fontsize=6,
                )
        axes[j, 0].set_ylabel(band.upper(), fontsize=7)
fig.suptitle(
    f"SPT0615 arcs, {model}: strong-shear background rows "
    '(4" F200W/F444W; blue observed axis, red predicted z=2)',
    fontsize=7,
)
fig.tight_layout()
fig.savefig(png, dpi=74)  # keeps the PNG under 1 MB

"""Contact sheet and photo-z match diagnostic of the RXC J0018 `arcs` run (post hoc; no outcome
change).

Sheet: the strong-shear rows with a population (background, then member_or_foreground), 4" F200W
/ F444W cutouts by S3 byte range (program 6882 o001 i2d); blue the observed axis (``observed``), red
the z = 2 predicted axis (``model_prediction``). ``posthoc.json``: per population, the distance from
each strong-shear row to its nearest DJA photo-z row (the ``arcs`` match radius is 0.3").

    python scripts/arcs_rxc0018_sheet.py <rj0018p1626-grizli-v7.5-fix.eazypy.zout.fits>
"""

import json
import math
import sys
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS

plt.switch_backend("Agg")
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "rxc0018-relics"
MODEL = "rxc0018-relics-lenstool"
t = Table.read(OUT / MODEL / "arcs.ecsv")
s = t[t["strong_shear"]]
z = Table.read(sys.argv[1])
_, d, _ = SkyCoord(s["ra"], s["dec"], unit="deg").match_to_catalog_sky(
    SkyCoord(z["ra"], z["dec"], unit="deg")
)
post = {
    "photoz_rows": len(z),
    "strong_shear_rows": len(s),
    "nearest_photoz_arcsec": {
        pop: sorted(round(float(v), 2) for v in d.arcsec[s["population"] == pop])
        for pop in ("background", "member_or_foreground", "unknown")
    },
    "rows_without_photoz_within_0.3arcsec": int((d.arcsec > 0.3).sum()),
}
# Could the frozen rule say "replicates" here? A passing control (>= 70 % aligned) leaves at most
# n - ceil(0.7 n) anti rows; the best case is P(>= 1 anti) under the shape-noise expectation.
summ = json.loads((OUT / "summary.json").read_text())["runs"][MODEL]
n = summ["strong_shear_background"]["n"]
mu = summ["anti_expected_shape_noise_background"]["expected"]
post["reachability"] = {
    "n_background": n,
    "max_anti_with_passing_control": n - math.ceil(0.7 * n),
    "expected_anti": mu,
    "best_case_p_ge_1": round(1 - math.exp(-mu), 4),
    "replicates_reachable": (n - math.ceil(0.7 * n)) >= 1 and 1 - math.exp(-mu) < 0.05,
}
with open(OUT / "posthoc.json", "w") as f:
    f.write(json.dumps(post, indent=1) + "\n")
print(json.dumps(post))

b = s[s["population"] != "unknown"]
b = b[np.argsort(b["population"] != "background", kind="stable")]
url = "s3://stpubdata/jwst/public/jw06882/L3/t/o001/jw06882-o001_t001_nircam_clear-{}_i2d.fits"
bands = ["f200w", "f444w"]
fig, axes = plt.subplots(len(bands), len(b), figsize=(1.9 * len(b), 4.2), squeeze=False)
for j, band in enumerate(bands):
    with fits.open(url.format(band), use_fsspec=True, fsspec_kwargs={"anon": True}) as h:
        w = WCS(h["SCI"].header)
        for i, r in enumerate(b):
            c = SkyCoord(r["ra"], r["dec"], unit="deg")
            co = Cutout2D(h["SCI"].section, c, 4 * u.arcsec, wcs=w)
            img = np.nan_to_num(co.data)
            lo, hi = np.percentile(img, [5, 99.5])
            ax = axes[j, i]
            ax.imshow(
                np.arcsinh((img - lo) / max(hi - lo, 1e-9) * 10), origin="lower", cmap="gray_r"
            )
            x0, y0 = co.wcs.world_to_pixel(c)
            for pa, col in ((r["pa_obs"], "tab:blue"), (r["pa_pred_z2"], "tab:red")):
                x1, y1 = co.wcs.world_to_pixel(c.directional_offset_by(pa * u.deg, 1.2 * u.arcsec))
                ax.plot([2 * x0 - x1, x1], [2 * y0 - y1, y1], color=col, lw=0.8)
            ax.set_xticks([])
            ax.set_yticks([])
            if j == 0:
                ax.set_title(
                    f"{r['label']} {r['population'][:6]} {r['orientation_class']}\n"
                    f"m{r['mag']:.1f} g{r['g_min']:.2f} z{r['z_phot']:.2f}",
                    fontsize=6,
                )
        axes[j, 0].set_ylabel(band.upper(), fontsize=7)
fig.suptitle(
    f"RXC J0018 arcs, {MODEL}: strong-shear rows with a population "
    '(4" F200W/F444W; blue observed axis, red predicted z=2)',
    fontsize=7,
)
fig.tight_layout()
fig.savefig(OUT / "sheet_strong_shear.png", dpi=74)  # keeps the PNG under 1 MB

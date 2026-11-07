"""Contact sheet for SMACS f200w_3670 and lensed system 2 (results/smacs_jwst3/f200w_3670).

    python scripts/smacs_f3670_sheet.py <cutout dir>

Rows: f200w_3670 (12" and 4"), images 2.1, 2.2, 2.3 (Mahler+2023 ICLv2 arcs.dat, 4").
Columns: F150W, F200W, F444W. S3 byte-range cutouts of the jwst 3.0.0 program 2736 o001
mosaics (not stored). Pixels ``observed``.
"""

import glob
import sys
from pathlib import Path

import matplotlib
import numpy as np
from astropy.io import fits
from astropy.table import Table

from jwst_anomaly.cutouts import make_cutouts
from jwst_anomaly.pipeline import l3_image_uri

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CLOUD = {
    "s3_bucket": "stpubdata",
    "l3_key_pattern": "jwst/public/jw{proposal:05d}/L3/t/o{observation:03d}/{obs_id}_{suffix}",
}
PRE = "jw02736-o001_t001_nircam_clear-"
# (label, RA, Dec, size in arcsec); positions: pipeline catalog (3670) and ICLv2 arcs.dat (2.x)
ROWS = [
    ("f200w_3670", 110.836848, -73.458635, 12.0),
    ("f200w_3670", 110.836848, -73.458635, 4.0),
    ("2.3", 110.8364983, -73.4588136, 4.0),
    ("2.1", 110.8387288, -73.4510508, 4.0),
    ("2.2", 110.8407771, -73.4552122, 4.0),
]
BANDS = ["f150w", "f200w", "f444w"]


def main(out: str) -> None:
    fig, ax = plt.subplots(len(ROWS), len(BANDS), figsize=(3 * len(BANDS), 3 * len(ROWS)))
    for k, (lab, ra, dec, size) in enumerate(ROWS):
        tg = Table({"source_uid": [f"r{k}"], "ra": [ra], "dec": [dec]})
        for c, b in enumerate(BANDS):
            d_out = Path(f"{out}/{k}{b}")
            make_cutouts(l3_image_uri(CLOUD, PRE + b), tg, size_arcsec=size, out_dir=d_out)
            d = fits.getdata(glob.glob(f"{d_out}/**/*.fits", recursive=True)[0], "SCI")
            n = d.shape[0]
            v = np.nanpercentile(d, [5, 99.7])
            img = np.arcsinh((d - v[0]) / (v[1] - v[0]) * 10)
            ax[k, c].imshow(img, origin="lower", cmap="gray")
            ax[k, c].plot(n / 2, n / 2, "r+", ms=10)
            ax[k, c].set_title(f'{lab} {b} {size:g}"', fontsize=9)
            ax[k, c].axis("off")
    plt.tight_layout()
    plt.savefig(f"{out}/sheet.png", dpi=90)


if __name__ == "__main__":
    main(sys.argv[1])

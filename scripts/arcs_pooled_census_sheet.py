"""Contact sheet of the 6 pooled anti rows for the clump-pair census (looked at by eye).

Re-reads 8″ F200W S3 cutouts of the anti rows (temporary) and ``census.ecsv``:

    python scripts/arcs_pooled_census_sheet.py

Output: ``results/arcs/pooled-census/census_anti.png`` (pixels observed; titles derived).
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from scipy.ndimage import median_filter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arcs_pooled_census as pc  # noqa: E402

from jwst_anomaly import cutouts  # noqa: E402


def main() -> int:
    rows = pc.load_rows()
    cen = Table.read(pc.OUT / "census.ecsv")
    anti = rows[rows["orientation_class"] == "anti"]
    fig, ax = plt.subplots(1, len(anti), figsize=(3.6 * len(anti), 4.2), squeeze=False)
    with tempfile.TemporaryDirectory() as tmp:
        for axis, r in zip(ax[0], anti, strict=True):
            f, lab = str(r["field"]), str(r["label"])
            tg = Table({"source_uid": [lab], "ra": [r["ra"]], "dec": [r["dec"]]})
            c = cutouts.make_cutouts(pc.image_url(f), tg, size_arcsec=8.0, out_dir=Path(tmp) / lab)
            with fits.open(c["path"][0]) as h:
                ext = next(x for x in h if x.data is not None)
                img, w = ext.data.astype(float), WCS(ext.header)
            x, y = (int(v) for v in w.world_to_pixel_values(r["ra"], r["dec"]))
            k = 32  # half-size in px (2″ at 0.031″/px; ~4″ at 0.063″/px)
            s = img[max(y - k, 0) : y + k + 1, max(x - k, 0) : x + k + 1]
            s = s - median_filter(np.nan_to_num(s), 33)
            lo, hi = np.nanpercentile(s, [5, 99.7])
            axis.imshow(s, origin="lower", cmap="gray", vmin=lo, vmax=hi)
            m = cen[(cen["field"] == f) & (cen["label"] == lab)][0]
            axis.set_title(
                f"{f} {lab} g_min {m['g_min']:.2f}\n"
                f"{m['deblend_class']} gate {bool(m['gate_passed'])}"
                f' sep {m["sep_arcsec"]:.2f}" dBIC {m["dbic"]:.0f}',
                fontsize=8,
            )
            axis.set_xticks([])
            axis.set_yticks([])
    fig.tight_layout()
    fig.savefig(pc.OUT / "census_anti.png", dpi=80)
    return 0


if __name__ == "__main__":
    sys.exit(main())

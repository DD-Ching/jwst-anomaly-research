"""Contact sheet of the MACS1423 screen-anti rows for the clump-pair census (looked at by eye).

    python scripts/arcs_macs1423_census_sheet.py --cut-dir <census --cut-dir>/f200w/<image stem>

Output: ``results/arcs/macs1423-canucs/vetting/census_anti.png`` (pixels observed; titles derived).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from scipy.ndimage import median_filter  # noqa: E402

from jwst_anomaly import paths  # noqa: E402


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cut-dir", type=Path, required=True)
    p.add_argument(
        "--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs" / "macs1423-canucs"
    )
    a = p.parse_args(argv)
    j = json.loads((a.arcs_dir / "vetting" / "census.json").read_text())
    t = Table.read(a.arcs_dir / "macs1423-canucs-maps" / "bgsub_shapes.ecsv")
    labs = j["anti_sets"]["screen_anti"]
    fig, ax = plt.subplots(1, len(labs), figsize=(4 * len(labs), 4.6), squeeze=False)
    for axis, lab in zip(ax[0], labs, strict=True):
        with fits.open(a.cut_dir / f"{lab}.fits") as h:
            ext = next(x for x in h if x.data is not None)
            img, w = ext.data.astype(float), WCS(ext.header)
        r = t[[str(v) == lab for v in t["label"]]][0]
        x, y = (int(v) for v in w.world_to_pixel_values(r["ra"], r["dec"]))
        k = 32  # half-size in px (2″ at 0.031″/px)
        s = img[y - k : y + k + 1, x - k : x + k + 1]
        s = s - median_filter(np.nan_to_num(s), 33)  # ~1″ median, as the census
        lo, hi = np.nanpercentile(s, [5, 99.7])
        axis.imshow(s, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        rr = j["rows"][lab]
        axis.set_title(
            f'{lab} {rr["class"]} sep {rr["sep_arcsec"]:.2f}" dBIC {rr["dbic"]:.0f}\n'
            f"PA_sep {rr['pa_sep']:.0f}  pa_obs {r['pa_obs']:.0f}  pred_z2 {r['pa_pred_z2']:.0f}",
            fontsize=9,
        )
    fig.tight_layout()
    fig.savefig(a.arcs_dir / "vetting" / "census_anti.png", dpi=80)
    return 0


if __name__ == "__main__":
    sys.exit(main())

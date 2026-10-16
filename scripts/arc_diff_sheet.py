"""Contact sheet of ``arc_diff.py`` hits: epoch-1 image and significance map S per band.

Usage: ``arc_diff_sheet.py MAPS.npz MORPH.json OUT.png [RA,DEC ...]``. MAPS is an ``arc_diff.py
--maps`` file, MORPH an ``arc_diff_morph.py --out`` file (its targets and verdicts label the rows);
extra RA,DEC pairs are added as control rows (e.g. a steady reference source). Each stamp is 2″
wide; S is shown at ±10. For looking at hits before any conclusion, not a measurement.
"""

from __future__ import annotations

import json
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.wcs import WCS  # noqa: E402


def main(argv: list[str]) -> None:
    maps, morph, out = argv[:3]
    extra = [tuple(map(float, x.split(","))) for x in argv[3:]]
    z = np.load(maps, allow_pickle=True)
    bands = sorted({k.split("_")[0] for k in z.keys()})
    targets = json.load(open(morph))["targets"]
    rows = [(t["id"], t["ra"], t["dec"], t["verdict"]) for t in targets]
    rows += [("ctrl", r, d, "control") for r, d in extra]
    ncol = 2 * len(bands)
    fig, ax = plt.subplots(len(rows), ncol, figsize=(2.2 * ncol, 2.2 * len(rows)), squeeze=False)
    for i, (hid, ra, dec, verdict) in enumerate(rows):
        for j, b in enumerate(bands):
            w = WCS(fits.Header.fromstring(str(z[b + "_wcs"])))
            x, y = w.world_to_pixel_values(ra, dec)
            h = int(round(1.0 / (abs(w.pixel_scale_matrix[1, 1]) * 3600)))
            xi, yi = int(round(float(x))), int(round(float(y)))
            sl = (slice(max(yi - h, 0), yi + h + 1), slice(max(xi - h, 0), xi + h + 1))
            im, s = z[b + "_e1img"][sl], z[b + "_s"][sl]
            lo, hi = np.nanpercentile(im, [5, 99.5]) if im.size else (0.0, 1.0)
            ax[i, 2 * j].imshow(im, origin="lower", cmap="gray", vmin=lo, vmax=hi)
            ax[i, 2 * j + 1].imshow(s, origin="lower", cmap="RdBu_r", vmin=-10, vmax=10)
            ax[i, 2 * j].set_title(f"{hid} {verdict} {b} e1", fontsize=7)
            ax[i, 2 * j + 1].set_title(f"S (±10) {b}", fontsize=7)
            for a in ax[i, 2 * j : 2 * j + 2]:
                a.plot(xi - sl[1].start, yi - sl[0].start, "+", color="lime", ms=8)
                a.set_xticks([])
                a.set_yticks([])
    fig.tight_layout()
    fig.savefig(out, dpi=80)


if __name__ == "__main__":
    main(sys.argv[1:])

"""Contact sheet of the gate V2 disagreements (grid edge_pair != stamp_edge). Usage (repo root):
sheet_v2.py STAMP_ROOT   (STAMP_ROOT/<field>/stamps/e<k>/ from gate_v2.py); writes sheet_v2.png"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.table import Table, vstack  # noqa: E402

root, here = Path(sys.argv[1]), Path(__file__).parent
t = vstack(
    [Table.read(here / f"gate_v2/pairs_{f}.ecsv") for f in ("macs0416", "abell2744", "nexus")]
)
d = t[t["edge_pair"] != t["stamp_edge"]]
n = len(d)
if n == 0:
    sys.exit("no disagreements; no sheet")
cols = 8
fig, axs = plt.subplots(
    2 * ((n + cols - 1) // cols), cols, figsize=(cols * 1.5, 3.2 * ((n + cols - 1) // cols))
)
for a in axs.flat:
    a.axis("off")
for i, r in enumerate(d):
    hits = list((root / str(r["field"]) / "stamps" / f"e{r['epoch']}").rglob(f"{r['uid']}.fits"))
    row, col = 2 * (i // cols), i % cols
    if hits:
        with fits.open(hits[0]) as h:
            sci = np.nan_to_num(h["SCI"].data)
            wht = np.nan_to_num(h["WHT"].data)
        axs[row, col].imshow(
            sci,
            origin="lower",
            cmap="gray",
            vmin=np.percentile(sci, 5),
            vmax=np.percentile(sci, 99.5),
        )
        axs[row + 1, col].imshow(wht, origin="lower", cmap="viridis", vmin=0)
    tag = "grid only" if r["edge_pair"] else "stamp only"
    axs[row, col].set_title(
        f"{r['uid']} e{r['epoch']}\n{tag} g{r['grid_wht_rel']:.2f} s{r['stamp_wht_rel']:.2f}",
        fontsize=6,
    )
fig.suptitle('gate V2 disagreements: SCI (top) and WHT (bottom) 2" F200W stamps', fontsize=8)
fig.savefig(here / "sheet_v2.png", dpi=80, bbox_inches="tight")
print(n, "disagreements")

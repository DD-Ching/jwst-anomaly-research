"""Look-at sheet: each MACS0416 / Abell 2744 F200W weight grid with the injection sources marked
(red: edge pair in that epoch; white: covered, not edge). Usage (repo root): sheet.py OUT PNG
(OUT from run.py; needs the screen rebuilt, same inputs)."""

import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.table import Table  # noqa: E402

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

from jwst_anomaly import pipeline  # noqa: E402

_spec = importlib.util.spec_from_file_location("edge", Path(__file__).with_name("edge.py"))
ed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ed)

out, png = Path(sys.argv[1]), Path(sys.argv[2])
gaia_root = out.parent
p = ds.Params()
panels = []
for name in ("macs0416", "abell2744"):
    field = ds.load_field(ds.DEFAULT_CONFIG, name)
    cats, lc, scale, err, res, ordn = ds._screen(
        field, p, Table.read(gaia_root / name / "gaia.ecsv")
    )
    rel = np.load(out / name / "wht_rel.npy")
    inj = ds.monitored(lc, err, ordn, p) & ~ds.flagged(res)
    edge = ed.edge_pair(np.isfinite(lc["flux"][:, :, 0]) & np.isfinite(rel), rel)
    for k, e in enumerate(field["epochs"]):
        if "F200W" not in e["bands"]:
            continue
        grid, wcs = ed.wht_grid(pipeline.l3_image_uri(field["cloud"], ds.obs_id(e, "F200W")), 32)
        x, y = wcs.all_world2pix(lc["ra"], lc["dec"], 0)
        cov = inj & np.isfinite(lc["flux"][:, k, 0])
        panels.append((f"{name} e{k}", grid, x / 32, y / 32, cov & ~edge[:, k], cov & edge[:, k]))
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
for ax, (title, g, x, y, ok, bad) in zip(axes.flat, panels, strict=False):
    med = np.nanmedian(g[g > 0])
    ax.imshow(np.nan_to_num(g) / med, origin="lower", cmap="gray", vmin=0, vmax=2)
    ax.scatter(x[ok], y[ok], s=14, facecolors="none", edgecolors="cyan")
    ax.scatter(x[bad], y[bad], s=20, facecolors="none", edgecolors="red")
    ax.set_title(f"{title}: {int(bad.sum())} edge / {int(ok.sum() + bad.sum())}")
    ax.set_xticks([])
    ax.set_yticks([])
for ax in axes.flat[len(panels) :]:
    ax.axis("off")
fig.tight_layout()
fig.savefig(png, dpi=80)

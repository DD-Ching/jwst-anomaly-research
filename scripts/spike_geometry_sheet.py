"""16" F200W/second-band contact sheet for the spike-geometry check (results/spike_geometry).

    python scripts/spike_geometry_sheet.py <scratch dir holding out_*/runs> <cutout dir>

Rows: the known streak CEERS f200w_2842 and the one gate-passing flagged source per field from
``results/spike_geometry/summary.json``. S3 byte-range cutouts (not stored). Pixels ``observed``.
"""

import glob
import json
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

S = sys.argv[1]
out = sys.argv[2]
summ = json.load(open("results/spike_geometry/summary.json"))


CLOUD = {
    "s3_bucket": "stpubdata",
    "l3_key_pattern": "jwst/public/jw{proposal:05d}/L3/t/o{observation:03d}/{obs_id}_{suffix}",
}


def uri(obs):
    return l3_image_uri(CLOUD, obs)


def row(sample, uid):
    t = Table.read(glob.glob(f"{S}/out_*/runs/*/{sample}/sources.ecsv")[0])
    i = list(t["source_uid"]).index(uid)
    return t[[i]]


jobs = [
    (
        "ceers_t021",
        "jw01345-o001_t021_nircam_f200w_2842",
        "jw01345-o001_t021_nircam_clear-",
        ["f200w", "f356w"],
    ),
    (
        "ceers_t021",
        [r["source_uid"] for r in summ["ceers"]["flagged"] if r["gate_pass"]][0],
        "jw01345-o001_t021_nircam_clear-",
        ["f356w", "f200w"],
    ),
    (
        "smacs0723_nircam",
        [r["source_uid"] for r in summ["smacs"]["flagged"] if r["gate_pass"]][0],
        "jw02736-o001_t001_nircam_clear-",
        ["f444w", "f200w"],
    ),
]
fig, ax = plt.subplots(len(jobs), 2, figsize=(7, 3.5 * len(jobs)))
for k, (sample, uid, pre, bands) in enumerate(jobs):
    tg = row(sample, uid)["source_uid", "ra", "dec"]
    for c, b in enumerate(bands):
        ct = make_cutouts(uri(pre + b), tg, size_arcsec=16.0, out_dir=Path(f"{out}/{k}{b}"))
        f = glob.glob(f"{out}/{k}{b}/**/*.fits", recursive=True)[0]
        d = fits.getdata(f, "SCI")
        n = d.shape[0]
        v = np.nanpercentile(d, [5, 99.5])
        ax[k, c].imshow(np.arcsinh((d - v[0]) / (v[1] - v[0]) * 10), origin="lower", cmap="gray")
        ax[k, c].plot(n // 2, n // 2, "r+", ms=14)
        ax[k, c].set_title(f'{uid.split("_nircam_")[1]} {b} 16"', fontsize=9)
        ax[k, c].axis("off")
plt.tight_layout()
plt.savefig(f"{out}/sheet.png", dpi=80)

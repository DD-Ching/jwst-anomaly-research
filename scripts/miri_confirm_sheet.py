"""Contact sheet for ``miri_confirm.py``: single-band-only MIRI sources in F200W, F770W, F1800W.

Usage: ``python scripts/miri_confirm_sheet.py <run_dir> <out_dir>`` after ``miri_confirm.py`` wrote
``results/miri_confirm/summary.json``. Picks (seeded) 8 matched and 4 unmatched F770W, 4 unmatched
F1800W and 3 matched F1500W sources; cutouts are S3 byte-range reads from ``stpubdata``.
Visual check only (observed images); no number in RESULT.md depends on it.
"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

from jwst_anomaly import cutouts, viz

R, SP = Path(sys.argv[1]), Path(sys.argv[2])
SP.mkdir(parents=True, exist_ok=True)
s = Table.read(R / "smacs0723_miri/sources.ecsv")
g = json.load(open("results/miri_confirm/summary.json"))["groups"]["single_band_only"]
matched = set(g["matched_uids"])
q = Table.read(R / "smacs0723_miri/quality.ecsv")
only = np.asarray(np.ma.filled(q["quality_reason"].astype(str), "")) == "single_band"
anchor = np.asarray(s["anchor_band"]).astype(str)
uid = np.asarray(s["source_uid"]).astype(str)
rng = np.random.default_rng(0)
pick = []
for band, want_match, k in [
    ("F770W", True, 8),
    ("F770W", False, 4),
    ("F1800W", False, 4),
    ("F1500W", True, 3),
]:
    m = only & (anchor == band) & (np.isin(uid, list(matched)) == want_match)
    idx = np.flatnonzero(m)
    pick += list(rng.choice(idx, min(k, len(idx)), replace=False))
t = s[pick]
tg = Table({"source_uid": t["source_uid"], "ra": t["ra"], "dec": t["dec"]})
tg["rank"] = np.arange(1, len(tg) + 1)
print(
    [
        (r, u, anchor[i], u in matched)
        for r, (i, u) in enumerate(zip(pick, uid[pick], strict=True), 1)
    ]
)
imgs = [
    "jwst/public/jw02736/L3/t/o002/jw02736-o002_t001_miri_f770w_i2d.fits",
    "jwst/public/jw02736/L3/t/o002/jw02736-o002_t001_miri_f1800w_i2d.fits",
    "jwst/public/jw02736/L3/t/o001/jw02736-o001_t001_nircam_clear-f200w_i2d.fits",
]
cs = []
for key in imgs:
    c = cutouts.make_cutouts(
        "s3://stpubdata/" + key,
        tg,
        size_arcsec=3.0,
        out_dir=SP / "cut",
        storage_options={"anon": True},
    )
    c["rank"] = tg["rank"]
    cs.append(c)
allc = vstack(cs, metadata_conflicts="silent")
viz.contact_sheet(
    allc,
    SP / "contact_sheet_rgb.png",
    title=(
        "SMACS MIRI single-band-only: 1-8 F770W matched, 9-12 F770W unmatched, "
        "13-16 F1800W unmatched, 17-19 F1500W matched"
    ),
)

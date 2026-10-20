"""Copy of a field's lightcurves.ecsv in which only its D-039 `on_sky` flag targets are unvetoed flags
(PLAN.md "Screen" / "Measurement"). Usage (repo root): make_targets.py FIELD OUT
(rewrites OUT/lightcurves.ecsv; writes OUT/targets_match.json with the weight_edge-removed targets)
"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

field, out = sys.argv[1], Path(sys.argv[2])
N = {"macs0416": 12, "abell2744": 8}[field]
DET_EPOCH = {"macs0416": 2, "abell2744": 1}[field]  # the only F200W catalogue epoch of the targets
v = Table.read(f"results/dimming/persistence_run/{field}/verdicts.ecsv")
on_sky = sorted(
    str(u)
    for u, k, s in zip(v["uid"], v["kind"], v["verdict"], strict=True)
    if k == "flag" and s == "on_sky"
)
assert len(on_sky) == N, on_sky
t = Table.read(out / "lightcurves.ecsv")
uid = np.asarray(t["uid"], str)
is_t = np.isin(uid, on_sky)
assert is_t.sum() == N
det = np.stack(
    [np.asarray(t[c], bool) for c in sorted(t.colnames) if c.startswith("F200W_e") and c.endswith("_det")],
    axis=1,
)
other = np.delete(det, DET_EPOCH, axis=1)
masked = ~np.isfinite(np.asarray(t[f"F200W_e{DET_EPOCH}_flux"], float))  # weight-edge term
ok = ~other.any(axis=1) & (det[:, DET_EPOCH] | masked)
assert bool(np.all(ok[is_t])), "uid does not match the persistence-run target"
flagged = np.asarray(t["flagged"], bool)
single = np.asarray(t["single_epoch"], bool)
veto = np.asarray(t["catalogue_veto"], bool)
# still a flag vetoed by single_epoch alone -> measured; otherwise explained by the weight-edge term
measured = is_t & flagged & ~(veto & ~single)
removed = sorted(set(uid[is_t]) - set(uid[measured]))
t["catalogue_veto"] = np.where(measured, False, veto | flagged)
t.write(out / "lightcurves.ecsv", overwrite=True)
info = {"field": field, "n_targets": N, "measured": int(measured.sum()), "weight_edge": removed}
(out / "targets_match.json").write_text(json.dumps(info, indent=1) + "\n")
print(json.dumps(info))

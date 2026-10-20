"""Copy of a field's lightcurves.ecsv in which only its D-039 `on_sky` flag targets are unvetoed
flags (PLAN.md "Screen" / "Measurement"). Usage (repo root):
make_targets.py FIELD OUT NOTERM_OUT
OUT holds the production screen (coverage term on), NOTERM_OUT the same screen with
`--no-weight-coverage`. Rewrites OUT/lightcurves.ecsv; writes OUT/targets_match.json.
"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

field, out, noterm = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
N = {"macs0416": 12, "abell2744": 8}[field]
DET_EPOCH = {"macs0416": 2, "abell2744": 1}[field]  # the only F200W catalogue epoch of the targets
OTHER_VETOES = ("near_star", "bright_neighbour", "edge_proxy", "blended", "sharp_artifact")

v = Table.read(f"results/dimming/persistence_run/{field}/verdicts.ecsv")
on_sky = sorted(
    str(u)
    for u, k, s in zip(v["uid"], v["kind"], v["verdict"], strict=True)
    if k == "flag" and s == "on_sky"
)
assert len(on_sky) == N, on_sky


def eligible(t: Table) -> np.ndarray:
    """Flagged, vetoed by `single_epoch` alone (the persistence-run selection)."""
    other = np.zeros(len(t), bool)
    for c in OTHER_VETOES:
        other |= np.asarray(t[c], bool)
    return np.asarray(t["flagged"], bool) & np.asarray(t["single_epoch"], bool) & ~other


def f200w_det(t: Table) -> np.ndarray:
    n = 1 + max(int(c[7:].split("_")[0]) for c in t.colnames if c.startswith("F200W_e"))
    return np.stack([np.asarray(t[f"F200W_e{k}_det"], bool) for k in range(n)], axis=1)


# uid check on the screen the persistence run used (no coverage term): every target is eligible
# and F200W-catalogued in its detection epoch only
t0 = Table.read(noterm / "lightcurves.ecsv")
is0 = np.isin(np.asarray(t0["uid"], str), on_sky)
assert is0.sum() == N
only = np.zeros(f200w_det(t0).shape[1], bool)
only[DET_EPOCH] = True
assert bool(np.all(f200w_det(t0)[is0] == only)), "uid does not match the persistence-run target"
assert bool(np.all(eligible(t0)[is0]))

t = Table.read(out / "lightcurves.ecsv")
uid = np.asarray(t["uid"], str)
assert np.array_equal(uid, np.asarray(t0["uid"], str))
is_t = np.isin(uid, on_sky)
measured = is_t & eligible(t)
removed = sorted(set(uid[is_t]) - set(uid[measured]))
# the two screens differ only by the coverage term, so a target eligible without it and not with it
# is removed by the term (checked: its pairs are masked, NaN flux, in at least one epoch)
for u in removed:
    k = int(np.flatnonzero(uid == u)[0])
    ne = f200w_det(t).shape[1]
    assert any(
        np.isfinite(t0[f"F200W_e{e}_flux"][k]) and not np.isfinite(t[f"F200W_e{e}_flux"][k])
        for e in range(ne)
    ), u
t["catalogue_veto"] = np.where(measured, False, np.asarray(t["catalogue_veto"]) | t["flagged"])
t.write(out / "lightcurves.ecsv", overwrite=True)
info = {"field": field, "n_targets": N, "measured": int(measured.sum()), "weight_edge": removed}
(out / "targets_match.json").write_text(json.dumps(info, indent=1) + "\n")
print(json.dumps(info))

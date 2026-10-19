"""Copy of the NEXUS lightcurves.ecsv in which only the 12 D-039 `on_sky` targets are unvetoed flags
(PLAN.md "Measurement"). Usage (repo root): make_targets.py SCREEN_OUT  (rewrites lightcurves.ecsv there)
"""

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

out = Path(sys.argv[1])
v = Table.read("results/dimming/persistence_run/nexus/verdicts.ecsv")
on_sky = {
    str(u)
    for u, k, s in zip(v["uid"], v["kind"], v["verdict"], strict=True)
    if k == "flag" and s == "on_sky"
}
assert len(on_sky) == 12, on_sky
t = Table.read(out / "lightcurves.ecsv")
is_t = np.isin(np.asarray(t["uid"], str), sorted(on_sky))
assert is_t.sum() == 12 and bool(np.all(t["flagged"][is_t]))
# the targets were vetoed by single_epoch alone (persistence_run selection)
assert not bool(
    np.any(np.asarray(t["catalogue_veto"])[is_t] & ~np.asarray(t["single_epoch"])[is_t])
)
t["catalogue_veto"] = np.where(
    is_t, False, np.asarray(t["catalogue_veto"]) | np.asarray(t["flagged"])
)
t.write(out / "lightcurves.ecsv", overwrite=True)
print("targets", int((t["flagged"] & ~t["catalogue_veto"]).sum()))

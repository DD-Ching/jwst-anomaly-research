"""Real flags the swap clears must be the ones already tested: the compact single-epoch vanish flags
of the production screen (coverage term on) against the persistence-run targets (`kind == flag`).
Usage (repo root): check_targets.py GAIA_ROOT CACHE_ROOT FIELD  -> runs/<field>/targets_check.json
"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

gaia_root, cache_root, name = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
p = ds.Params()
field = ds.load_field(ds.DEFAULT_CONFIG, name)
gaia = Table.read(gaia_root / name / "gaia.ecsv")
cats, lc, scale, err, res, ordn = ds._screen(field, p, gaia, cache_root / name / "wht_rel.npz")
static_veto = ds.catalogue_veto({**ordn, "single_epoch": np.zeros_like(ordn["single_epoch"])})
t = res["vanish"] & ordn["single_epoch"] & ~static_veto & ordn["point_like"]
now = {f"{name}-d{i:05d}" for i in np.flatnonzero(t)}
v = Table.read(f"results/dimming/persistence_run/{name}/verdicts.ecsv")
v = v[v["kind"] == "flag"]
verdict = {str(u): str(s) for u, s in zip(v["uid"], v["verdict"], strict=True)}
out = {
    "field": name,
    "n_targets_term_on": len(now),
    "n_tested_targets": len(verdict),
    "not_tested": sorted(now - set(verdict)),
    "verdicts_of_term_on_targets": {u: verdict.get(u, "not_tested") for u in sorted(now)},
}
path = Path("results/dimming/swap_inject/runs") / name / "targets_check.json"
path.write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out))

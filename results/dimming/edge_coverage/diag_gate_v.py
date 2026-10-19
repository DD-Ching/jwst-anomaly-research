"""Post hoc diagnostic of the failed gate V: per NEXUS edge target and epoch, proxy coverage,
catalogue detection, neighbour count and wht_rel. Usage (repo root): diag_gate_v.py OUT
(writes results/dimming/edge_coverage/gate_v_diagnostic.json)."""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

out = Path(sys.argv[1])
p = ds.Params()
field = ds.load_field(ds.DEFAULT_CONFIG, "nexus")
cats, lc, scale, err, res, ordn = ds._screen(field, p, Table.read(out.parent / "nexus/gaia.ecsv"))
rel = np.load(out / "nexus" / "wht_rel.npy")
tg = Table.read("results/dimming/nexus_o014/targets.ecsv")
tg = tg[np.asarray(tg["cut_edge"]) | np.asarray(tg["cut_low_weight"])]
rows = {}
for u in np.asarray(tg["uid"], str):
    i = int(u.split("-d")[1])
    rows[u] = {
        "covered": [bool(x) for x in np.isfinite(lc["flux"][i, :, 0])],
        "detected": [bool(x) for x in lc["det"][i, :, 0]],
        "neighbours": [int(x) for x in lc["neighbours"][i, :, 0]],
        "wht_rel": [None if not np.isfinite(x) else round(float(x), 3) for x in rel[i]],
    }
Path("results/dimming/edge_coverage/gate_v_diagnostic.json").write_text(
    json.dumps({"post_hoc": True, "targets": rows}, indent=1) + "\n"
)
for u, r in rows.items():
    print(
        u,
        "cov",
        "".join("1" if c else "0" for c in r["covered"]),
        "det",
        "".join("1" if c else "0" for c in r["detected"]),
        "rel",
        r["wht_rel"],
    )

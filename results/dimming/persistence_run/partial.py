"""Per-(source, epoch) verdicts from the stored per-exposure tables of run.py, for fields whose run
did not finish (RESULT.md "Deviations"). Usage (repo root): partial.py FIELD GAIA_ROOT"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402
import persistence_check as pc  # noqa: E402

name, groot = sys.argv[1], Path(sys.argv[2])
d = Path("results/dimming/persistence_run") / name
p = ds.Params()
field = ds.load_field(ds.DEFAULT_CONFIG, name)
cats, lc, scale, err, res, ordn = ds._screen(field, p, Table.read(groot / name / "gaia.ecsv"))
static_veto = ds.catalogue_veto({**ordn, "single_epoch": np.zeros_like(ordn["single_epoch"])})
target = ds.flagged(res) & ordn["single_epoch"] & ~static_veto & ordn["point_like"] & res["vanish"]
rows = []
for f in sorted(d.glob("exposures_e*.ecsv")):
    k = int(f.stem.rsplit("_e", 1)[1])
    e = Table.read(f)
    for uid in sorted(set(e["uid"])):
        s = e[e["uid"] == uid]
        i = int(uid.rsplit("-d", 1)[1])
        snr, sus, npr = (np.asarray(s[c]) for c in ("snr", "suspect", "n_prior"))
        rows.append(
            {
                "uid": uid,
                "epoch": k,
                "kind": "flag" if target[i] else "monitored",
                "verdict": pc.classify(snr, sus, npr),
                "n_exposures": len(s),
                "mosaic_snr": float(lc["flux"][i, k, 0] / err[i, k, 0]),
            }
        )
t = Table(rows=rows)
t.meta.update(provenance="derived", source=f"results/dimming/persistence_run/partial.py {name}")
t.write(d / "verdicts.ecsv", overwrite=True)
(d / "counts.json").write_text(
    json.dumps(
        {
            "field": name,
            "n_targets": int(target.sum()),
            "epochs_done": sorted({int(r["epoch"]) for r in rows}),
            "complete": False,
        },
        indent=1,
    )
    + "\n"
)
print(
    name,
    dict(
        zip(
            *np.unique([f"{r['kind']}:{r['verdict']}" for r in rows], return_counts=True),
            strict=True,
        )
    ),
)

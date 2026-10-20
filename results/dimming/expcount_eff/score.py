"""Score the frozen exposure-count efficiency model (PLAN.md) on the NEXUS monitored controls.
Usage (repo root): score.py  (reads results/dimming/persistence_run/nexus/verdicts.ecsv)"""

import json
import math
from pathlib import Path

import numpy as np
from astropy.table import Table

TRAIN = {0, 1, 3, 4}
HELD = {5, 6, 7}
SWAP = {"le4": 0.412, "gt4": 0.961}  # PLAN.md, frozen


def wilson(k: int, n: int, z: float = 1.96) -> list[float] | None:
    if not n:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 3), round(c + h, 3)]


t = Table.read("results/dimming/persistence_run/nexus/verdicts.ecsv")
m = t[(t["kind"] == "monitored") & (t["mosaic_snr"] >= 10)]
out = {}
for split, ep in (("train", TRAIN), ("held_out", HELD)):
    s = m[np.isin(m["epoch"], sorted(ep))]
    out[split] = {
        "epochs": sorted({int(e) for e in s["epoch"]}),
        "persistence": int((s["verdict"] == "persistence").sum()),
    }
    for cls, sel in (("le4", s["n_exposures"] <= 4), ("gt4", s["n_exposures"] > 4)):
        q = s[sel]
        k, n = int((q["verdict"] == "on_sky").sum()), len(q)
        row = {"on_sky": k, "n": n, "wilson95": wilson(k, n)}
        row["verdicts"] = {
            str(v): int(c)
            for v, c in zip(*np.unique(q["verdict"], return_counts=True), strict=True)
        }
        if split == "held_out":
            if n < 10:
                row["rule1"] = "not tested (< 10 pairs)"
            else:
                row["rule1"] = "holds" if row["wilson95"][1] >= SWAP[cls] else "fails"
        out[split][cls] = row
    out[split]["persistence_wilson95"] = wilson(out[split]["persistence"], len(s))
Path("results/dimming/expcount_eff/score.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out, indent=1))

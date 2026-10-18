"""E1 + EP post hoc diagnostic (not in the frozen plan): EP-BAT / EP-GBM long-lag cells with joint
detections removed first, as D-074 addendum 6 "Revisit if" asks for satellite x satellite channels.

A joint EP event has a same-direction BAT or GBM pair at |lag| < 1 h (the control pairs); those EP
events and their BAT / GBM partners are removed from both catalogues, then the same 18-cell counts
and jit null are recomputed. Writes ``results/e1_ep/diag_no_joint.json``.

  python scripts/e1_ep_diag.py [--jit 10000] [--cpu 4]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import e1_ep as ep  # noqa: E402
import e1_events as e1  # noqa: E402
import e1_swift as sw  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jit", type=int, default=10000)
    ap.add_argument("--cpu", type=int, default=1)
    a = ap.parse_args(argv)
    p = ep.P
    ev = sw.load(sw.fetch(False))
    ev = {"GBM": ev["GBM"], "BAT": ev["BAT"], "EP": ep.load_ep()}
    drop = {k: set() for k in ev}
    for other in ("BAT", "GBM"):
        t = e1.close_pairs(ev["EP"], ev[other], p, False, 3600.0)
        t = t[np.abs(t["lag_s"]) < 3600.0]
        drop["EP"] |= set(t["a"])
        drop[other] |= set(t["b"])
    kept = {k: v[~np.isin(v["name"], sorted(drop[k]))] for k, v in ev.items()}
    samples = {k: en.Sample.from_table(v) for k, v in kept.items()}
    obs = ep.count_all(samples, p)
    jit = ep.run_null("jit", a.jit, samples, a.cpu, 6_000_000)
    mask = ep.test_mask(p)
    pa = en.analytic_p(obs, jit)
    labels = en.lag_labels(p)
    cells = {}
    for c, (x, y) in enumerate(ep.CHANNELS):
        for k in range(len(labels)):
            for j, cl in enumerate(e1.CLASSES[:2]):
                if mask[c, k, j]:
                    cells[f"{x}-{y} {labels[k]} {cl}"] = {
                        "obs": int(obs[c, k, j]),
                        "jit_mean": round(float(jit[:, c, k, j].mean()), 2),
                        "jit_sd": round(float(jit[:, c, k, j].std()), 2),
                        "analytic_p": float(pa[c, k, j]),
                    }
    pmin = min(v["analytic_p"] for v in cells.values())
    n_cells = int(mask.sum())
    out = {
        "provenance": "derived",
        "post_hoc": True,
        "removed": {k: len(v) for k, v in drop.items()},
        "removed_ep": sorted(drop["EP"]),
        "events_kept": {k: len(v) for k, v in kept.items()},
        "n_jit": a.jit,
        "cells": cells,
        "min_cell_p": pmin,
        "tested_cells": n_cells,
        "min_cell_bonferroni": min(1.0, pmin * n_cells),
    }
    with open(ep.OUT / "diag_no_joint.json", "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: out[k] for k in ("removed", "min_cell_p", "min_cell_bonferroni")}))
    print(json.dumps({k: v for k, v in cells.items() if "GBM" in k or "BAT" in k}, indent=0))
    return 0


if __name__ == "__main__":
    sys.exit(main())

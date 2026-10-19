"""Apply the frozen PLAN.md rules to run.py outputs. Usage (repo root): summarise.py OUT_ROOT"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

BINS = [(0, 10), (10, 20), (20, 50), (50, np.inf)]  # mosaic F200W S/N (PLAN.md)
MIN_ON_SKY, MAX_PERSIST, MIN_PAIRS = 0.8, 0.02, 10  # ASSUMPTION (PLAN.md rule 1)


def efficiency(mon: Table) -> list[dict]:
    snr = np.asarray(mon["mosaic_snr"], float)
    v = np.asarray(mon["verdict"]).astype(str)
    out = []
    for lo, hi in BINS:
        m = (snr >= lo) & (snr < hi)
        n = int(m.sum())
        out.append(
            {
                "snr_lo": lo,
                "snr_hi": None if np.isinf(hi) else hi,
                "n_pairs": n,
                **{
                    f"n_{k}": int((v[m] == k).sum())
                    for k in ("on_sky", "persistence", "inconclusive", "undetected")
                },
                "on_sky_fraction": float((v[m] == "on_sky").mean()) if n else None,
            }
        )
    return out


def usable(eff: list[dict], persist_rate: float) -> bool:
    """Rule 1: on_sky >= 0.8 in every S/N >= 20 bin with >= 10 pairs, and persistence <= 2 %."""
    hi = [b for b in eff if b["snr_lo"] >= 20 and b["n_pairs"] >= MIN_PAIRS]
    return (
        bool(hi)
        and all(b["on_sky_fraction"] >= MIN_ON_SKY for b in hi)
        and persist_rate <= MAX_PERSIST
    )


if __name__ == "__main__":
    root = Path(sys.argv[1])
    t = vstack(
        [Table.read(f) for f in sorted(root.glob("*/verdicts.ecsv"))], metadata_conflicts="silent"
    )
    mon = t[t["kind"] == "monitored"]
    tg = t[t["kind"] == "flag"]
    eff = efficiency(mon)
    pr = (
        float(np.mean(np.asarray(mon["verdict"]).astype(str) == "persistence"))
        if len(mon)
        else float("nan")
    )
    res = {
        "provenance": "derived",
        "monitored_pairs": len(mon),
        "monitored_persistence_rate": pr,
        "efficiency_by_mosaic_snr": eff,
        "rule1_usable": usable(eff, pr),
        "targets": [
            {
                k: (r[k].item() if hasattr(r[k], "item") else str(r[k]))
                for k in (
                    "uid",
                    "epoch",
                    "verdict",
                    "n_exposures",
                    "n_detected",
                    "n_suspect",
                    "n_clean",
                    "mosaic_snr",
                )
            }
            for r in tg
        ],
        "targets_to_vet": [str(r["uid"]) for r in tg if r["verdict"] == "on_sky"],
        "n_failed_reads": int(  # n_failed is per (field, epoch), repeated on each of its rows
            sum(
                {
                    (str(r["uid"]).rsplit("-d", 1)[0], int(r["epoch"])): int(r["n_failed"])
                    for r in t
                }.values()
            )
        ),
    }
    (root / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "targets"}, indent=1))

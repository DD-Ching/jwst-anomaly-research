"""W3 short-event recovery study (owner directive V2.1, priority 1): stage 1, a frozen sample.

Runs the production injection chain (``w3_moa.run_inject``, CHAIN_VERSION unchanged) for t_E = 3 d
only, on one field, with a seed that no limit uses, and writes a tracked table with a frozen
``split`` column (dev / validation, by crc32 of the row index and seed) fixed before any analysis.
Later stages tune classifiers on ``dev`` only and report recovery on ``validation`` once.

    python scripts/w3_short_te.py --field gb12 --per-cell 1000 --seed 3003 --procs 4
    python scripts/w3_short_te.py --field gb12 --seed 3003 --diagnose   # dev split only
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import w3_moa as wm  # noqa: E402
from astropy.table import Table  # noqa: E402

OUT = ROOT / "results" / "w3_moa" / "short_te"


SPIKE_SIGMA = 3.0  # ASSUMPTION: spike / umbra epochs: injected flux beyond 3 median errors
JOB_KEYS = ("event_id", "tE", "rho", "u0", "t0", "Is")  # the inputs that define one injection


def split_of(row, seed: int) -> str:
    """Frozen assignment (ASSUMPTION: 50/50) from a SHA-256 of the injection's own inputs and the
    seed. It does not depend on the table's row order (pool completion order) or on any outcome."""
    key = f"{seed}|" + "|".join(repr(float(row[k])) if k != "event_id" else str(row[k])
                                for k in JOB_KEYS)  # fmt: skip
    return "validation" if hashlib.sha256(key.encode()).digest()[0] & 1 else "dev"


def sampling_diagnostics(tab: Table, field: str) -> Table:
    """Per flagged injection (``derived``): how well the injected signal is sampled. Spike nights
    have an epoch where the injected flux rises above 3 median errors; umbra nights one where it
    falls below −3 median errors (the jackknife's own feature threshold)."""
    import numpy as np

    wm.set_field(field)
    pre = wm.read_prescreen()
    flagged = tab[np.asarray(tab["flagged"], bool)]
    arrays = wm.load_arrays(list(dict.fromkeys(flagged["event_id"])), aux=False, pre=pre)
    rows = []
    for r in flagged:
        t, f, sf = arrays[r["event_id"]][:3]
        fs = float(wm.moa.mag_to_counts(float(r["Is"]), wm.moa.parse_event_id(r["event_id"])[1]))
        prm = {k: float(r[k]) for k in ("tE", "rho", "u0", "t0")}
        sig = wm.injected_signal(t, "W3", fs, prm)
        thr = SPIKE_SIGMA * float(np.median(sf))
        nt = wm.nights(t)
        up, down = sig > thr, sig < -thr
        tc = prm["t0"]
        rows.append(
            {
                "event_id": r["event_id"],
                "split": r["split"],
                "rho": prm["rho"],
                "u0": prm["u0"],
                "Is": float(r["Is"]),
                "failed_test": str(r["failed_test"] or "")
                if not np.ma.is_masked(r["failed_test"])
                else "",
                "recovered": bool(r["recovered"]),
                "spike_nights": int(np.unique(nt[up]).size),
                "spike_before": bool((up & (t < tc)).any()),
                "spike_after": bool((up & (t > tc)).any()),
                "umbra_nights": int(np.unique(nt[down]).size),
                "feature_nights": int(np.unique(nt[up | down]).size),
            }
        )
    out = Table(rows=rows)
    out.meta.update(
        provenance="derived",
        source=f"scripts/w3_short_te.py --diagnose: injected-signal sampling of flagged {field} "
        "t_E = 3 d injections (simulated signals on real carriers)",
    )
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--field", default="gb12")
    ap.add_argument("--per-cell", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=3003)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--diagnose", action="store_true", help="sampling diagnostics of the sample")
    ap.add_argument("--resplit", action="store_true", help="rewrite the split column (job key)")
    a = ap.parse_args(argv)
    if a.resplit:
        path = OUT / f"injections_{a.field}_te3_seed{a.seed}.ecsv.gz"
        tab = Table.read(path, format="ascii.ecsv")
        tab["split"] = [split_of(r, a.seed) for r in tab]
        wm.w3.write_ecsv_gz(tab, path)
        print(
            f"re-split {path}: {sum(tab['split'] == 'dev')} dev, "
            f"{sum(tab['split'] != 'dev')} validation"
        )
        return 0
    if a.diagnose:
        tab = Table.read(
            OUT / f"injections_{a.field}_te3_seed{a.seed}.ecsv.gz", format="ascii.ecsv"
        )
        d = sampling_diagnostics(tab, a.field)
        out = OUT / f"diagnostics_{a.field}_te3_seed{a.seed}.ecsv"
        d.write(out, overwrite=True)
        print(f"wrote {out}: {len(d)} flagged injections")
        return 0
    wm.set_field(a.field)
    prod = wm.out_dir() / f"injections_{a.field}.ecsv"
    keep = prod.with_suffix(".production.ecsv")
    had_prod = prod.exists()
    if keep.exists():
        raise SystemExit(f"{keep} exists (an interrupted run?): restore it to {prod} first")
    if had_prod:
        shutil.copy2(prod, keep)  # run_inject writes to the production path
    te_grid = wm.INJ_TE
    wm.INJ_TE = (3.0,)
    try:
        path = wm.run_inject(a.procs, a.per_cell, 0, a.seed, False, "lf")
        tab = Table.read(path)
    finally:
        wm.INJ_TE = te_grid
        if had_prod:
            shutil.move(keep, prod)
        else:
            prod.unlink(missing_ok=True)  # never leave the study table as production injections
    tab["split"] = [split_of(r, a.seed) for r in tab]
    tab.meta["study"] = "W3 short-event recovery, stage 1 (t_E = 3 d; frozen dev/validation split)"
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"injections_{a.field}_te3_seed{a.seed}.ecsv.gz"
    wm.w3.write_ecsv_gz(tab, out)
    print(f"wrote {out}: {len(tab)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

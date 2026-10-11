"""W3 short-event recovery study (owner directive V2.1, priority 1): stage 1, a frozen sample.

Runs the production injection chain (``w3_moa.run_inject``, CHAIN_VERSION unchanged) for t_E = 3 d
only, on one field, with a seed that no limit uses, and writes a tracked table with a frozen
``split`` column (dev / validation, by crc32 of the row index and seed) fixed before any analysis.
Later stages tune classifiers on ``dev`` only and report recovery on ``validation`` once.

    python scripts/w3_short_te.py --field gb12 --per-cell 1000 --seed 3003 --procs 4
"""

from __future__ import annotations

import argparse
import shutil
import sys
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import w3_moa as wm  # noqa: E402
from astropy.table import Table  # noqa: E402

OUT = ROOT / "results" / "w3_moa" / "short_te"


def split_of(i: int, seed: int) -> str:
    """Frozen assignment (ASSUMPTION: 50/50): independent of every outcome column."""
    return "validation" if zlib.crc32(f"{seed}:{i}".encode()) % 2 else "dev"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--field", default="gb12")
    ap.add_argument("--per-cell", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=3003)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args(argv)
    wm.set_field(a.field)
    prod = wm.out_dir() / f"injections_{a.field}.ecsv"
    keep = prod.with_suffix(".production.ecsv")
    if prod.exists():
        shutil.copy2(prod, keep)  # run_inject writes to the production path
    wm.INJ_TE = (3.0,)
    try:
        path = wm.run_inject(a.procs, a.per_cell, 0, a.seed, False, "lf")
        tab = Table.read(path)
    finally:
        if keep.exists():
            shutil.move(keep, prod)
    tab["split"] = [split_of(i, a.seed) for i in range(len(tab))]
    tab.meta["study"] = "W3 short-event recovery, stage 1 (t_E = 3 d; frozen dev/validation split)"
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"injections_{a.field}_te3_seed{a.seed}.ecsv.gz"
    wm.w3.write_ecsv_gz(tab, out)
    print(f"wrote {out}: {len(tab)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

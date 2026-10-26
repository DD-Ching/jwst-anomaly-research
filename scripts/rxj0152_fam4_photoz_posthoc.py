"""Post hoc (not in the plan): rxj0152_fam4_photoz with the controls' median offset removed.

Shifts the RELICS positions by minus the median (RELICS - SLICE) control offset stored in
``summary.json``, re-matches, and writes ``posthoc.json``. Deviation from the frozen plan
(step 2: no offset correction); it does not change the frozen outcome.

    python scripts/rxj0152_fam4_photoz_posthoc.py --bpz <RELICS cl0152m13 catalogue>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rxj0152_fam4_photoz as f  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bpz", type=Path, required=True)
    ap.add_argument(
        "--out", type=Path, default=f.OUT, help="directory of the frozen run's summary.json"
    )
    args = ap.parse_args(argv)
    summary = json.loads((args.out / "summary.json").read_text())
    dra, ddec = summary["control_median_offset_arcsec"]
    cat = f.read_bpz(args.bpz)
    cat["Dec"] = np.asarray(cat["Dec"], float) - ddec / 3600
    cat["RA"] = np.asarray(cat["RA"], float) - dra / 3600 / np.cos(np.radians(cat["Dec"]))
    sl = Table.read(f.SLICE)
    controls = [
        {"image": str(r["image_id"]), "z": float(r["z"]), "match": f.match(r["ra"], r["dec"], cat)}
        for r in sl
        if str(r["system"]) in f.CONTROL_FAMILIES
    ]
    segs = []
    for s in summary["segments"]:
        m = f.match(s["ra"], s["dec"], cat)
        segs.append({"segment": s["segment"], "match": m, "class": f.classify(m)})
    gp = f.gate_p(controls)
    out = {
        "post_hoc": True,
        "offset_removed_arcsec": [dra, ddec],
        "gate_P": gp,
        "controls": controls,
        "segments": segs,
        "outcome_if_frozen_rules": f.outcome(gp["pass"], [s["class"] for s in segs]),
    }
    (args.out / "posthoc.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"gate_P": gp, "classes": [s["class"] for s in segs]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Detailed review of the first ``--capacity`` entries of a batch review queue (PLAN §5).

For each reviewed entry: a contact sheet (``arc_diff_sheet.py`` on the box's screen maps, that hit
only), forced photometry in every covering public mosaic (``hit_lightcurve.py``), and the class rule,
applied by code to the two search epochs (the box's epoch-1 programme and VENUS 6882):

- ``edge``: the hit lies within 1″ of the box edge (box half-width 12″);
- per band, z = (f2 - f1) / sqrt(e1² + e2²) from the 0.15″ forced photometry;
- ``no significant change``: no search band with |z| >= 3;
- ``real change``: both bands with the same sign and at least one with |z| >= 3;
- ``undetermined``: anything else (e.g. one band only, or opposite signs).

The by-eye sheet verdict comes from ``--sheet-verdicts`` (JSON {hit key: text}); a verdict starting
with "veto" makes the entry ``artifact``. Entries beyond capacity are written to ``pending.json``
unchanged; nothing is dropped. Every number is ``derived``.

    python scripts/batch_review.py --batch DIR --scratch SCRATCH --capacity 10 [--sheet-verdicts V.json]
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def key(e: dict) -> str:
    return f"{e['box_id']}:{e['hit_id']}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--batch", type=Path, required=True)
    ap.add_argument("--boxes", type=Path, required=True)
    ap.add_argument("--scratch", type=Path, required=True)
    ap.add_argument("--capacity", type=int, default=10)
    ap.add_argument("--sheet-verdicts", type=Path)
    ap.add_argument("--skip-photometry", action="store_true", help="reuse lightcurves.ecsv")
    a = ap.parse_args(argv)
    from astropy.table import Table

    queue = json.loads((a.batch / "review_queue.json").read_text())["entries"]
    boxes = {b["box_id"]: b for b in json.loads(a.boxes.read_text())["chosen"]}
    sel, pend = queue[: a.capacity], queue[a.capacity :]
    out = a.batch / "review"
    out.mkdir(parents=True, exist_ok=True)
    for e in sel:
        d = a.batch / e["box_id"] / "screen"
        morph = json.loads((d / "morph_all.json").read_text())
        one = {**morph, "targets": [t for t in morph["targets"] if t["id"] == e["hit_id"]]}
        tmp = out / f"morph_{e['box_id']}_{e['hit_id']}.json"
        tmp.write_text(json.dumps(one) + "\n")
        subprocess.run(
            [sys.executable, str(HERE / "arc_diff_sheet.py"),
             str(a.scratch / e["box_id"] / "screen" / "maps.npz"), str(tmp),
             str(out / f"sheet_{e['box_id']}_{e['hit_id']}.png")],
            check=False, capture_output=True,
        )  # fmt: skip
        tmp.unlink()
    if sel and not a.skip_photometry:
        cmd = [sys.executable, str(HERE / "hit_lightcurve.py")]
        for e in sel:
            cmd += ["--pos", key(e).replace(":", "__"), str(e["ra"]), str(e["dec"])]
        subprocess.run([*cmd, "--threads", "8", "--out", str(out)], check=False)
    lc = Table.read(out / "lightcurves.ecsv") if (out / "lightcurves.ecsv").exists() else None
    verdicts = json.loads(a.sheet_verdicts.read_text()) if a.sheet_verdicts else {}
    rows = []
    for e in sel:
        b = boxes[e["box_id"]]
        c = math.cos(math.radians(b["dec"]))
        dx, dy = (e["ra"] - b["ra"]) * c * 3600, (e["dec"] - b["dec"]) * 3600
        edge = max(abs(dx), abs(dy)) > 11.0
        z, phot = {}, {}
        if lc is not None:
            m = lc[lc["name"] == key(e).replace(":", "__")]
            for band in ("F150W", "F444W"):
                o1 = b["e1"][band]["obs_id"]
                o2 = b["e2"][band]["obs_id"]
                r1 = m[m["obs_id"] == o1]
                r2 = m[m["obs_id"] == o2]
                if len(r1) and len(r2):
                    f1, e1_ = float(r1["flux_njy"][0]), float(r1["err_njy"][0])
                    f2, e2_ = float(r2["flux_njy"][0]), float(r2["err_njy"][0])
                    z[band] = round((f2 - f1) / math.hypot(e1_, e2_), 2)
                    phot[band] = [round(f1, 1), round(e1_, 1), round(f2, 1), round(e2_, 1)]
        zs = [v for v in z.values() if np.isfinite(v)]
        if edge:
            cls = "artifact (edge)"
        elif len(zs) < 2:
            cls = "undetermined (photometry missing)"
        elif max(abs(v) for v in zs) < 3:
            cls = "no significant change"
        elif all(v > 0 for v in zs) or all(v < 0 for v in zs):
            cls = "real change"
        else:
            cls = "undetermined"
        sv = verdicts.get(key(e))
        if sv and sv.lower().startswith("veto"):
            cls = f"artifact (sheet: {sv})"
        rows.append(
            {**e, "edge": edge, "z": z, "phot_njy_f1_e1_f2_e2": phot, "sheet": sv, "class": cls}
        )
    (out / "review.json").write_text(
        json.dumps({"label": "derived; PLAN §5 classes", "reviewed": rows}, indent=1) + "\n"
    )
    (a.batch / "pending.json").write_text(
        json.dumps(
            {"label": "review-list entries beyond capacity: pending, not cleared", "entries": pend},
            indent=1,
        )
        + "\n"
    )
    for r in rows:
        print(r["rank"], key(r), r["max_abs_sigma"], r["z"], r["class"])
    print("pending", len(pend))
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""No-injection control at the positions of morphology-vetoed injections (2026-10-16).

Implements ``results/transients/veto_control/PLAN.md`` (written before the measurement). An
``arc_diff.py --inject-grid`` star is added to epoch 2 with one sign (positive), so an opposite-sign
(negative) lobe that vetoes it is either already in the no-injection difference map at that
position (H1, neighbour residuals) or made by the injection itself (H2). ``arc_diff_morph.py
--known`` on the no-injection maps measures the lobes; this script builds its target list and
compares vetoed with passed injected hits. Every number is ``derived``; the thresholds of the
decision rule are ASSUMPTIONs fixed in the plan.

    python scripts/arc_veto_control.py targets --runs <dir> ... --out targets.csv
    python scripts/arc_veto_control.py compare --runs <dir> ... --control morph_noinj.json \\
        --maps base.npz --field smacs --out compare.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

LEVEL = 3.0  # |S| of a lobe (arc_diff_morph Params.level)
H1_MIN_VETOED = 0.70
H1_MIN_EXCESS = 0.30
H2_MAX_VETOED = 0.30
VETOES = ("rosette", "dipole")


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "1", "yes"}


def injected_hits(run: Path) -> list[dict]:
    """Injected stars of one run that are hits and not near a bright core, with their verdict."""
    targets = json.loads((run / "morph.json").read_text())["targets"]
    verdict = {t["id"]: t["verdict"] for t in targets}
    s_hit = {
        t["id"]: max((abs(b["s_hit"]) for b in t["bands"].values()), default=0.0) for t in targets
    }
    rows = []
    with open(run / "injected.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["hit"] and not _truthy(r["near_bright_core"]):
                if r["hit"] not in verdict:
                    raise ValueError(f"{run.name}/{r['hit']} has no morphology verdict")
                rows.append(
                    {
                        "name": f"{run.name}:{r['name']}",
                        "ra": float(r["ra"]),
                        "dec": float(r["dec"]),
                        "mag": r["mag"],
                        "verdict": verdict[r["hit"]],
                        "s_hit": s_hit[r["hit"]],
                    }
                )
    return rows


def has_negative_lobe(bands: dict) -> bool:
    """A negative lobe (|S| >= LEVEL, morph lobe definition) near the position in any band."""
    return any(p <= -LEVEL for b in bands.values() for p in b.get("lobe_peaks", []))


def centre_negative(bands: dict) -> bool:
    """Post hoc: the no-injection S at the position itself is <= -LEVEL in some band (then the
    counted 'lobe' may be the region under the star rather than one around it)."""
    return any(b.get("s_hit", 0.0) <= -LEVEL for b in bands.values())


def negative_ratio(bands: dict, s_hit: float) -> float:
    """Post hoc: strongest pre-existing negative lobe over the injected hit's top |S|.

    Region peaks are global (a region touching the disk may peak outside it) and the two values
    may come from different bands, so this is a coarse indicator only."""
    neg = [-p for b in bands.values() for p in b.get("lobe_peaks", []) if p <= -LEVEL]
    return round(max(neg) / s_hit, 3) if neg and s_hit > 0 else 0.0


def decide(vetoed: list[bool], passed: list[bool]) -> dict:
    """Plan decision rule on per-position 'pre-existing negative lobe' flags."""
    fv = sum(vetoed) / len(vetoed) if vetoed else None
    fp = sum(passed) / len(passed) if passed else None
    if fv is None:
        outcome = "undetermined"
    elif fv >= H1_MIN_VETOED and fp is not None and fv - fp >= H1_MIN_EXCESS:
        outcome = "H1"
    elif fv <= H2_MAX_VETOED:
        outcome = "H2"
    else:
        outcome = "mixed"
    return {
        "n_vetoed": len(vetoed),
        "n_passed": len(passed),
        "frac_vetoed": None if fv is None else round(fv, 3),
        "frac_passed": None if fp is None else round(fp, 3),
        "outcome": outcome,
    }


def compare(runs: list[Path], control: Path, maps: Path | None, field: str | None) -> dict:
    ctrl = {t["id"]: t for t in json.loads(control.read_text())["targets"]}
    mask = None
    if maps is not None:
        from arc_chain import _pixel, load_mask

        mask, w, *_ = load_mask(maps, field)
    elif field is not None:
        raise ValueError("--field needs --maps (the mask strata would be empty)")
    rows = []
    for run in runs:
        for r in injected_hits(run):
            if r["name"] not in ctrl:
                raise ValueError(f"{r['name']} missing from {control.name}; re-run --known")
            c = ctrl[r["name"]]
            neg = has_negative_lobe(c["bands"])
            row = r | {
                "control_verdict": c["verdict"],
                "noinj_negative_lobe": neg,
                "noinj_centre_negative_post_hoc": centre_negative(c["bands"]),
                "noinj_neg_ratio_post_hoc": negative_ratio(c["bands"], r["s_hit"]),
            }
            if mask is not None:
                p = _pixel(mask.shape, w, r["ra"], r["dec"])
                row["in_mask"] = None if p is None else bool(mask[p[1], p[0]])
            row["noinj_lobe_peaks"] = {b: v.get("lobe_peaks", []) for b, v in c["bands"].items()}
            rows.append(row)
    out: dict = {"provenance": "derived", "control": control.name, "strata": {}}
    for mag in sorted({r["mag"] for r in rows}) + ["all"]:
        for where in ("all", "out_mask", "in_mask"):
            sel = [
                r
                for r in rows
                if (mag == "all" or r["mag"] == mag)
                and (where == "all" or r.get("in_mask") is (where == "in_mask"))
            ]
            v = [r["noinj_negative_lobe"] for r in sel if r["verdict"] in VETOES]
            p = [r["noinj_negative_lobe"] for r in sel if r["verdict"] == "pass"]
            vet = [r for r in sel if r["verdict"] in VETOES]
            pas = [r for r in sel if r["verdict"] == "pass"]
            split = {
                f"post_hoc_{k}_centre_negative": sum(r["noinj_centre_negative_post_hoc"] for r in g)
                for k, g in (("vetoed", vet), ("passed", pas))
            } | {
                f"post_hoc_{k}_lobe_without_negative_centre": sum(
                    r["noinj_negative_lobe"] and not r["noinj_centre_negative_post_hoc"] for r in g
                )
                for k, g in (("vetoed", vet), ("passed", pas))
            }
            out["strata"][f"{mag}/{where}"] = (
                decide(v, p)
                | split
                | {
                    "post_hoc_neg_ratio_vetoed": sorted(
                        r["noinj_neg_ratio_post_hoc"] for r in sel if r["verdict"] in VETOES
                    ),
                    "post_hoc_neg_ratio_passed": sorted(
                        r["noinj_neg_ratio_post_hoc"] for r in sel if r["verdict"] == "pass"
                    ),
                }
            )
    out["rows"] = rows
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("targets")
    t.add_argument("--runs", type=Path, nargs="+", required=True)
    t.add_argument("--out", type=Path, required=True)
    c = sub.add_parser("compare")
    c.add_argument("--runs", type=Path, nargs="+", required=True)
    c.add_argument("--control", type=Path, required=True)
    c.add_argument("--maps", type=Path, required=True)
    c.add_argument("--field", required=True)
    c.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "targets":
        rows = [r for run in a.runs for r in injected_hits(run)]
        a.out.parent.mkdir(parents=True, exist_ok=True)
        with open(a.out, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=["name", "ra", "dec", "mag", "verdict", "s_hit"])
            wr.writeheader()
            wr.writerows(rows)
        print(len(rows))
        return 0
    out = compare(a.runs, a.control, a.maps, a.field)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["strata"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

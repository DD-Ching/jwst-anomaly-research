"""Fading / vanishing calibration of ``arc_diff.py`` with pre-matching injections (PLAN.md).

``make``: the fixed injection sample for one seed. Positions come from ``arc_diff.grid_job`` on the
**epoch-1** box only (pixels above the 90th percentile, ``--sep`` apart); classes cycle
(unchanged, appear, partial, vanish) and magnitudes cycle independently over ``--mags``:

- unchanged: m → m;  appear: none → m;  partial: m → m + 0.75 (flux × 0.5);  vanish: m → none.

Writes ``name,ra,dec,mag1,mag2,cls,mag`` for ``arc_diff.py --inject-pre``.

``eval``: per class × magnitude, the stages of PLAN.md §3 from one or more runs (each an
``injected_pre.csv`` + ``hits.ecsv`` + ``arc_diff_morph.py --all-hits`` JSON):
hit (≥ 2 bands), not near a bright core, morph ``pass``, correct sign, review list
(hit ∧ unflagged ∧ pass) and review with the correct sign; for ``unchanged`` the false-hit count.
Binomial 68 % intervals are Wilson intervals. Every number is ``derived``.

    python scripts/fade_calib.py make --uri <e1 F150W i2d> --ra 15.73 --dec -49.253 \
        --seed 20 --out s20.csv
    python scripts/fade_calib.py eval --run RUN_DIR MORPH.json [--run …] --out eval.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

CLASSES = ("unchanged", "appear", "partial", "vanish")
FADE_CLASSES = ("partial", "vanish")
PARTIAL_DMAG = 0.75  # ASSUMPTION: flux × 0.5 (2.5 log10 2 = 0.753)


def class_mags(cls: str, m: float) -> tuple[float | None, float | None]:
    """(epoch-1, epoch-2) AB magnitudes of one injected source; None = no source."""
    return {
        "unchanged": (m, m),
        "appear": (None, m),
        "partial": (m, m + PARTIAL_DMAG),
        "vanish": (m, None),
    }[cls]


def sample_rows(positions, mags) -> list[dict]:
    rows = []
    for i, (ra, dec) in enumerate(positions):
        cls, m = CLASSES[i % len(CLASSES)], mags[i % len(mags)]
        m1, m2 = class_mags(cls, m)
        rows.append(
            {
                "name": f"pre{i:03d}",
                "ra": f"{ra:.8f}",
                "dec": f"{dec:.8f}",
                "mag1": "" if m1 is None else f"{m1:.3f}",
                "mag2": "" if m2 is None else f"{m2:.3f}",
                "cls": cls,
                "mag": f"{m:.2f}",
            }
        )
    return rows


def sign_ok(cls: str, hit_sign: str, s_by_band: dict[str, float | None]) -> bool:
    """Correct sign: '-' for fades ('mixed' if the largest |S| is negative), '+' for appear."""
    if cls in FADE_CLASSES:
        if hit_sign == "-":
            return True
        if hit_sign == "mixed":
            v = [x for x in s_by_band.values() if x is not None]
            return bool(v) and min(v, key=lambda t: -abs(t)) < 0
        return False
    if cls == "appear":
        if hit_sign == "+":
            return True
        if hit_sign == "mixed":
            v = [x for x in s_by_band.values() if x is not None]
            return bool(v) and min(v, key=lambda t: -abs(t)) > 0
        return False
    return True


def wilson(k: int, n: int, z: float = 1.0) -> list[float] | None:
    if n == 0:
        return None
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [round(max(c - h, 0.0), 3), round(min(c + h, 1.0), 3)]


def stage_rows(run_dir: Path, morph_path: Path) -> list[dict]:
    verdict = {t["id"]: t["verdict"] for t in json.loads(morph_path.read_text())["targets"]}
    out = []
    with open(run_dir / "injected_pre.csv") as f:
        for r in csv.DictReader(f):
            hit = r["hit"] or None
            bands = {
                k: (float(r[k]) if r.get(k) not in (None, "", "None") else None)
                for k in r
                if k.isupper() and k.startswith("F")
            }
            unflag = hit is not None and r["near_bright_core"] != "True"
            mp = hit is not None and verdict.get(hit) == "pass"
            so = hit is not None and sign_ok(r["cls"], r["hit_sign"], bands)
            out.append(
                {
                    "cls": r["cls"],
                    "mag": r["mag"],
                    "hit": hit is not None,
                    "unflagged": unflag,
                    "morph_pass": mp,
                    "sign_ok": so,
                    "review": unflag and mp,
                    "review_sign": unflag and mp and so,
                }
            )
    return out


def summarize(rows: list[dict]) -> dict:
    keys = ("hit", "unflagged", "morph_pass", "sign_ok", "review", "review_sign")
    out: dict = {}
    for cls in CLASSES:
        for mag in sorted({r["mag"] for r in rows}):
            sel = [r for r in rows if r["cls"] == cls and r["mag"] == mag]
            if not sel:
                continue
            n = len(sel)
            cell = {"n": n}
            for k in keys:
                c = sum(r[k] for r in sel)
                cell[k] = c
                cell[k + "_frac"] = round(c / n, 3)
                cell[k + "_68"] = wilson(c, n)
            out.setdefault(cls, {})[mag] = cell
    unch = [r for r in rows if r["cls"] == "unchanged"]
    fh = sum(r["hit"] for r in unch)
    out["unchanged_false_hits"] = {
        "n": len(unch),
        "false_hits": fh,
        "frac": round(fh / len(unch), 3) if unch else None,
        "frac_68": wilson(fh, len(unch)),
    }
    return out


def pass_rule(summary: dict, mag: str = "26.00") -> dict:
    """PLAN.md §4: vanish at 26.0 AB review-with-sign ≥ 0.80 and unchanged false hits ≤ 0.10."""
    v = summary.get("vanish", {}).get(mag)
    u = summary["unchanged_false_hits"]
    ok_v = bool(v and v["review_sign_frac"] >= 0.80)
    ok_u = bool(u["frac"] is not None and u["frac"] <= 0.10)
    return {
        "vanish_review_sign_frac": v["review_sign_frac"] if v else None,
        "unchanged_false_hit_frac": u["frac"],
        "vanish_ok": ok_v,
        "unchanged_ok": ok_u,
        "passed": ok_v and ok_u,
    }


def rank_boxes(ra, dec, half, step, exclude, min_images):
    """Candidate box centres on a ``step``-arcsec grid over the image list, ranked by the number of
    catalogued images inside (descending), dropping boxes that overlap ``exclude`` boxes
    (list of (ra, dec)) or hold fewer than ``min_images``. Coverage only; no image data used."""
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    d0 = float(np.median(dec))
    cosd = math.cos(math.radians(d0))
    x, y = (ra - ra.mean()) * cosd * 3600, (dec - d0) * 3600
    gx = np.arange(x.min(), x.max() + step, step)
    gy = np.arange(y.min(), y.max() + step, step)
    ex = [((r - ra.mean()) * cosd * 3600, (dd - d0) * 3600) for r, dd in exclude]
    out = []
    for cy in gy:
        for cx in gx:
            if any(abs(cx - a) < 2 * half and abs(cy - b) < 2 * half for a, b in ex):
                continue
            n = int(((abs(x - cx) <= half) & (abs(y - cy) <= half)).sum())
            if n >= min_images:
                out.append((n, float(ra.mean() + cx / 3600 / cosd), float(d0 + cy / 3600)))
    out.sort(key=lambda t: -t[0])
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make")
    m.add_argument("--uri", required=True)
    m.add_argument("--ra", type=float, required=True)
    m.add_argument("--dec", type=float, required=True)
    m.add_argument("--half-arcsec", type=float, default=12.0)
    m.add_argument("--n", type=int, default=80)
    m.add_argument("--sep", type=float, default=1.5)
    m.add_argument("--margin", type=float, default=1.1)
    m.add_argument("--mags", type=float, nargs="+", default=[26.0, 27.0, 27.5])
    m.add_argument("--seed", type=int, required=True)
    m.add_argument("--out", type=Path, required=True)
    e = sub.add_parser("eval")
    e.add_argument("--run", nargs=2, action="append", required=True, metavar=("DIR", "MORPH"))
    e.add_argument("--out", type=Path, required=True)
    bx = sub.add_parser("boxes", help="choose extra boxes from coverage only (PLAN §6)")
    bx.add_argument("--field", default="elgordo")
    bx.add_argument("--uri", nargs="+", required=True, help="every i2d that must cover the box")
    bx.add_argument("--exclude", type=float, nargs=2, action="append", default=[])
    bx.add_argument("--half-arcsec", type=float, default=12.0)
    bx.add_argument("--step", type=float, default=2.0)
    bx.add_argument("--min-images", type=int, default=5)
    bx.add_argument("--min-finite", type=float, default=0.95)
    bx.add_argument("--n", type=int, default=2)
    bx.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    if a.cmd == "boxes":
        from arc_chain import image_positions
        from arc_diff import read_box

        ra, dec = image_positions(a.field)
        chosen, checked = [], []
        for n, r, d in rank_boxes(ra, dec, a.half_arcsec, a.step, a.exclude, a.min_images):
            if any(
                abs((r - c["ra"]) * math.cos(math.radians(d)) * 3600) < 2 * a.half_arcsec
                and abs((d - c["dec"]) * 3600) < 2 * a.half_arcsec
                for c in chosen
            ):
                continue
            fin = []
            for u in a.uri:
                try:
                    b = read_box(u, r, d, a.half_arcsec)
                    full = (2 * math.ceil(a.half_arcsec / b["scale"]) + 1) ** 2
                    fin.append(round(float(np.isfinite(b["sci"]).sum()) / full, 4))
                except (ValueError, OSError):  # outside the mosaic, or a failed read
                    fin.append(0.0)
            row = {"ra": round(r, 6), "dec": round(d, 6), "n_images": n, "finite": fin}
            checked.append(row)
            if min(fin) >= a.min_finite:
                chosen.append(row)
            if len(chosen) == a.n:
                break
        out = {
            "label": "derived; chosen from image lists and coverage only, before any difference",
            "rule": vars(a) | {"out": str(a.out)},
            "chosen": chosen,
            "checked": checked,
        }
        a.out.write_text(json.dumps(out, indent=1) + "\n")
        print(json.dumps(chosen))
        return 0
    if a.cmd == "make":
        from arc_diff import grid_job

        pos = grid_job(a.uri, a.ra, a.dec, a.half_arcsec, a.n, a.sep, a.margin, a.seed)
        rows = sample_rows(pos, a.mags)
        with open(a.out, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)
        print(f"{len(rows)} rows -> {a.out}")
        return 0
    rows = [r for d, mp in a.run for r in stage_rows(Path(d), Path(mp))]
    summary = summarize(rows)
    out = {
        "label": "derived; injected sources are simulations",
        "runs": [d for d, _ in a.run],
        "summary": summary,
        "pass_rule": pass_rule(summary),
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["pass_rule"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

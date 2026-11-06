"""Match the SMACS 0723 candidate records (jwst 2.0.1 runs) to a jwst 3.0.0 reference-sample run.

Plan: ``results/smacs_jwst3/PLAN.md``. For each ``docs/candidates/jw02736-*.md`` record, read
its recorded position and rank, find the nearest merged source of the new run (``sources.ecsv``)
and its rank (``scores.ecsv``), and write ``results/smacs_jwst3/candidates.csv`` plus per-sample
counts in ``summary.json``.

Usage: python scripts/smacs_jwst3_compare.py <run_dir> [--radius 0.15]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

REPO = Path(__file__).resolve().parents[1]
POS = re.compile(r"RA ([0-9]+\.[0-9]+), Dec [−-]([0-9]+\.[0-9]+)")
FEATURES = ("ref_mag", "ref_sharpness", "ref_roundness", "color_f150w_f200w", "color_f200w_f277w")
RANK = re.compile(r"rank (\d+), ensemble ([0-9]+\.[0-9]+)")


def read_records(folder: Path) -> list[dict]:
    """Recorded position, rank and ensemble score of each SMACS record (first occurrence)."""
    out = []
    for path in sorted(folder.glob("jw02736-*.md")):
        text = path.read_text(encoding="utf-8")
        pos, rank = POS.search(text), RANK.search(text)
        if pos is None:
            continue
        out.append(
            {
                "record": path.stem,
                "ra": float(pos.group(1)),
                "dec": -float(pos.group(2)),
                "old_rank": int(rank.group(1)) if rank else None,
                "old_score": float(rank.group(2)) if rank else None,
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--sample", default="smacs0723_nircam")
    ap.add_argument("--radius", type=float, default=0.15, help="match radius, arcsec (ASSUMPTION)")
    ap.add_argument("--out", type=Path, default=REPO / "results/smacs_jwst3/candidates.csv")
    args = ap.parse_args()

    sample_dir = args.run_dir / args.sample
    sources = Table.read(sample_dir / "sources.ecsv")
    scores = Table.read(sample_dir / "scores.ecsv")
    star_scores = args.run_dir / f"{args.sample}-stars" / "scores.ecsv"
    stars = Table.read(star_scores) if star_scores.exists() else None
    feats = Table.read(sample_dir / "features.ecsv")
    quality = Table.read(sample_dir / "quality.ecsv")
    feat_of = {str(r["source_uid"]): r for r in feats}
    qual_of = {str(r["source_uid"]): r for r in quality}
    rank_of = {str(r["source_uid"]): (int(r["rank"]), float(r["score"])) for r in scores}
    star_rank = {str(r["source_uid"]): int(r["rank"]) for r in stars} if stars is not None else {}

    cat = SkyCoord(sources["ra"], sources["dec"], unit="deg")
    rows = []
    for rec in read_records(REPO / "docs/candidates"):
        c = SkyCoord(rec["ra"], rec["dec"], unit="deg")
        sep = c.separation(cat).to(u.arcsec).value
        i = int(sep.argmin())
        uid = str(sources["source_uid"][i])
        matched = sep[i] <= args.radius
        new_rank, new_score = rank_of.get(uid, (None, None)) if matched else (None, None)
        rows.append(
            {
                **rec,
                "new_uid": uid if matched else "",
                "sep_arcsec": round(float(sep[i]), 3),
                "n_within_1arcsec": int((sep <= 1.0).sum()),
                "new_rank": new_rank,
                "new_score": None if new_score is None else round(new_score, 4),
                "new_rank_star": star_rank.get(uid) if matched else None,
                "n_ranked": len(scores),
                "quality_ok": bool(qual_of[uid]["quality_ok"]) if matched else None,
                "quality_reason": str(qual_of[uid]["quality_reason"]) if matched else "",
                **{
                    k: (None if not matched else round(float(feat_of[uid][k]), 3)) for k in FEATURES
                },
            }
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {"run_id": args.run_dir.name, "match_radius_arcsec": args.radius, "samples": {}}
    for sdir in sorted(p for p in args.run_dir.iterdir() if (p / "sources.ecsv").exists()):
        q = Table.read(sdir / "quality.ecsv")
        f = Table.read(sdir / "features.ecsv")
        ok = np.asarray(q["quality_ok"], dtype=bool)
        entry = {"n_merged": len(q), "n_quality_ok": int(ok.sum())}
        if (sdir / "scores.ecsv").exists():
            entry["n_ranked"] = len(Table.read(sdir / "scores.ecsv"))
        if "ref_sharpness" in f.colnames:
            sharp = np.ma.filled(np.ma.asarray(f["ref_sharpness"], dtype=float), np.nan)
            entry["ref_sharpness_nan_frac_gated"] = round(float(np.isnan(sharp[ok]).mean()), 3)
        summary["samples"][sdir.name] = entry
    out_json = args.out.with_name("summary.json")
    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    for r in rows:
        print(r)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

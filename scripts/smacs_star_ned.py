"""NED objects near the quality-passing star stratum of a run (TASKS "Then" item 4).

Usage: python scripts/smacs_star_ned.py <run_dir> [--sample smacs0723_nircam] [--radius 0.5]
       [--out results/smacs_jwst3/star_ned/matches.csv]

Reads ``<sample>/populations.ecsv`` (star basis per source) and ``<sample>-stars/scores.ecsv``
(the ranked, quality-passing stars), queries NED (``crossmatch.query_matches``, services
``ned``) within ``--radius`` and writes one row per star-NED pair, plus stars without a match
(empty NED columns). NED content is observed; ``star_basis`` and the separation are derived.
Galaxy-like NED types are an ASSUMPTION listed in ``GALAXY_TYPES``.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import crossmatch

REPO = Path(__file__).resolve().parents[1]
# NED object types read as "a galaxy-like object" (ASSUMPTION; NED type list:
# https://ned.ipac.caltech.edu/help/ui/nearposn-list_objecttypes). NED "*" is star or point source.
GALAXY_TYPES = frozenset(
    {"G", "G_Lens", "GPair", "GTrpl", "GGroup", "QSO", "Q_Lens", "RadioS", "IrS"}
)


def _s(v: object) -> str:
    return "" if np.ma.is_masked(v) else str(v)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--sample", default="smacs0723_nircam")
    ap.add_argument("--radius", type=float, default=0.5, help="NED radius, arcsec (ASSUMPTION)")
    ap.add_argument("--out", type=Path, default=REPO / "results/smacs_jwst3/star_ned/matches.csv")
    args = ap.parse_args()

    pops = Table.read(args.run_dir / args.sample / "populations.ecsv")
    basis = {str(u): str(b) for u, b in zip(pops["source_uid"], pops["star_basis"], strict=True)}
    stars = Table.read(args.run_dir / f"{args.sample}-stars" / "scores.ecsv")
    srcs = Table.read(args.run_dir / args.sample / "sources.ecsv")
    pos = {
        str(u): (float(r), float(d))
        for u, r, d in zip(srcs["source_uid"], srcs["ra"], srcs["dec"], strict=True)
    }
    uids = [str(u) for u in stars["source_uid"]]
    targets = Table(
        {"source_uid": uids, "ra": [pos[u][0] for u in uids], "dec": [pos[u][1] for u in uids]}
    )
    targets.meta.update(provenance="derived", source=f"ranked star stratum of {args.run_dir.name}")
    m = crossmatch.query_matches(targets, radius_arcsec=args.radius, services=("ned",))
    if m.meta.get("services_failed"):
        raise SystemExit(f"NED failed: {m.meta['services_failed']}")
    rank = {str(u): int(r) for u, r in zip(stars["source_uid"], stars["rank"], strict=True)}

    rows = []
    matched = set()
    for r in m:
        uid = str(r["source_uid"])
        matched.add(uid)
        t = _s(r["match_type"])
        rows.append(
            {
                "source_uid": uid,
                "star_rank": rank[uid],
                "star_basis": basis[uid],
                "ra": f"{pos[uid][0]:.6f}",
                "dec": f"{pos[uid][1]:.6f}",
                "ned_name": _s(r["match_id"]),
                "ned_type": t,
                "galaxy_like": t in GALAXY_TYPES,
                "sep_arcsec": f"{float(r['sep_arcsec']):.3f}",
                "ned_redshift": _s(r["redshift"]),
            }
        )
    for uid in uids:
        if uid not in matched:
            rows.append(
                {
                    "source_uid": uid,
                    "star_rank": rank[uid],
                    "star_basis": basis[uid],
                    "ra": f"{pos[uid][0]:.6f}",
                    "dec": f"{pos[uid][1]:.6f}",
                    "ned_name": "",
                    "ned_type": "",
                    "galaxy_like": False,
                    "sep_arcsec": "",
                    "ned_redshift": "",
                }
            )
    rows.sort(key=lambda d: (d["star_rank"], float(d["sep_arcsec"] or "inf")))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    summary: dict[str, dict[str, int]] = {}
    for b in sorted(set(basis[u] for u in uids)):
        us = [u for u in uids if basis[u] == b]
        gal = {d["source_uid"] for d in rows if d["galaxy_like"] and d["star_basis"] == b}
        gal03 = {
            d["source_uid"]
            for d in rows
            if d["galaxy_like"] and d["star_basis"] == b and float(d["sep_arcsec"]) <= 0.3
        }
        summary[b] = {
            "n_stars": len(us),
            "n_ned_galaxy_radius": len(gal),
            "n_ned_galaxy_0p3": len(gal03),
        }
    out = {
        "run_id": args.run_dir.name,
        "radius_arcsec": args.radius,
        "ned_backend": m.meta.get("backends_used"),
        "query_utc": m.meta.get("query_utc"),
        "by_basis": summary,
    }
    args.out.with_name("summary.json").write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()

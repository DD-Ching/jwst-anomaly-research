"""Abell 2744 system 700.1 (UNCOVER 78): public NIRSpec MSA coverage of the arc ends.

Frozen plan: ``results/lens_consistency/a2744_s700_nirspec/PLAN.md``. MAST CAOM footprints
(NIRSpec/MSA, 1" around each knot) and DJA v4.6 rows within 0.5" of each knot; ``decide`` gives
the outcome.

    python scripts/a2744_s700_nirspec.py [--dja PATH]
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import a2744_uncover as au  # noqa: E402
import macs0647_specz as m6  # noqa: E402

from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

OUT = Path("results/lens_consistency/a2744_s700_nirspec")
# Frozen rules (PLAN.md; ASSUMPTIONs)
RELEASED_AFTER = "2026-10-20"
FOOTPRINT_RADIUS = "1 arcsec"
Z_KNOWN = 2.3151  # UNCOVER v2.0 table, system 78, Price+24 NIRSpec
KNOWN_PROGRAMME = "2561"
END_A = ("78.11", "78.21", "78.31")
END_B = ("78.12", "78.22", "78.32")


def decide(gate_k: bool, end_a: bool, end_b: bool, new_public: list[str]) -> str:
    if not gate_k:
        return "not validated"
    if end_a and end_b:
        return "both ends measured"
    if new_public:
        return "new input, not in DJA"
    return "no new input"


def knots() -> Table:
    t = Table.read(fetch_catalog(*au.UNCOVER_V2["images"]))
    ids = np.array([str(i) for i in t["ID"]])
    k = t[np.char.startswith(ids, "78.")]
    out = Table({"image_id": [str(i) for i in k["ID"]], "ra": k["RA"], "dec": k["Dec"]})
    return out


def footprints(k: Table) -> dict[str, list[dict]]:
    from astroquery.mast import Observations

    cut = Time(RELEASED_AFTER).mjd
    out: dict[str, list[dict]] = {}
    for r in k:
        o = Observations.query_criteria(
            coordinates=f"{r['ra']} {r['dec']}",
            radius=FOOTPRINT_RADIUS,
            obs_collection="JWST",
            instrument_name="NIRSPEC/MSA",
        )
        progs: dict[str, dict] = {}
        for x in o:
            p = str(x["proposal_id"])
            rel = float(x["t_obs_release"]) if np.isfinite(x["t_obs_release"]) else None
            d = progs.setdefault(
                p, {"programme": p, "n_rows": 0, "planned": True, "release_min": None}
            )
            d["n_rows"] += 1
            if int(x["calib_level"]) >= 0:
                d["planned"] = False
            if rel is not None and (d["release_min"] is None or rel < d["release_min"]):
                d["release_min"] = rel
        for d in progs.values():
            m = d.pop("release_min")
            d["release_min"] = Time(m, format="mjd").iso[:10] if m is not None else None
            d["new_public"] = bool((not d["planned"]) and m is not None and m > cut)
        out[str(r["image_id"])] = sorted(progs.values(), key=lambda d: d["programme"])
    return out


def dja_path(arg: str | None) -> tuple[Path, Path | None]:
    if arg:
        return Path(arg), None
    tmp = Path(tempfile.mkdtemp(prefix="dja46_"))
    gz = tmp / m6.TABLE["name"]
    with urllib.request.urlopen(m6.TABLE["url"], timeout=600) as r, open(gz, "wb") as fh:
        shutil.copyfileobj(r, fh)
    csv = tmp / "table.csv"
    with gzip.open(gz, "rb") as fi, open(csv, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    gz.unlink()
    return csv, tmp


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dja", help="local uncompressed DJA v4.6 csv (default: stream and delete)")
    args = ap.parse_args(argv)
    k = knots()
    fp = footprints(k)
    csv, tmp = dja_path(args.dja)
    try:
        rows, n = m6.match_rows(csv, k)
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)
    spec = {i: [r for r in v if m6._is_spec(r)] for i, v in rows.items()}
    known_z = any(
        abs(float(r["z_best"]) - Z_KNOWN) <= m6.Z_AGREE * (1 + Z_KNOWN)
        for v in spec.values()
        for r in v
        if r.get("z_best") not in (None, "")
    )
    known_prog = any(
        d["programme"] == KNOWN_PROGRAMME and not d["planned"] for v in fp.values() for d in v
    )
    gate_k = known_prog and known_z
    end_a = any(spec.get(i) for i in END_A)
    end_b = any(spec.get(i) for i in END_B)
    new_public = sorted({d["programme"] for v in fp.values() for d in v if d["new_public"]})
    planned = sorted({d["programme"] for v in fp.values() for d in v if d["planned"]})
    outcome = decide(gate_k, end_a, end_b, new_public)
    summary = {
        "plan": "PLAN.md",
        "provenance": "observed (archive metadata: MAST CAOM footprints, DJA v4.6 catalogue rows)",
        "params": {
            "released_after": RELEASED_AFTER,
            "footprint_radius": FOOTPRINT_RADIUS,
            "dja_match_arcsec": m6.MATCH_RADIUS_ARCSEC,
            "spec_grade": m6.SPEC_GRADE,
            "z_agree": m6.Z_AGREE,
            "z_known": Z_KNOWN,
        },
        "dja_table": m6.TABLE,
        "dja_rows_read": n,
        "knots": [
            {"image_id": str(r["image_id"]), "ra": float(r["ra"]), "dec": float(r["dec"])}
            for r in k
        ],
        "footprints": fp,
        "dja_rows": rows,
        "gate_k": {
            "known_programme_returned": bool(known_prog),
            "known_z_recovered": bool(known_z),
            "passed": bool(gate_k),
        },
        "end_a_measured": bool(end_a),
        "end_b_measured": bool(end_b),
        "new_public_programmes": new_public,
        "planned_programmes": planned,
        "outcome": outcome,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(
        json.dumps(
            {
                k_: summary[k_]
                for k_ in (
                    "gate_k",
                    "end_a_measured",
                    "end_b_measured",
                    "new_public_programmes",
                    "planned_programmes",
                    "outcome",
                )
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

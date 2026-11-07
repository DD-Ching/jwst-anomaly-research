"""Which open SMACS 0723 records fall inside the program 2736 NIRISS and MIRI footprints.

Answers TASKS "Later: NIRISS imaging of program 2736" and whether a second instrument at the
2022 epoch could test the 2736 NIRCam position of ``f200w_384`` (D-017 addendum 3). Reuses the
geometry of ``scripts/grism_4043_coverage.py``; record positions come from
``results/grism_4043/coverage.json`` (candidates.csv plus the 2736 F200W catalog positions of
``1816`` and ``384``). Footprints are the MAST ``s_region`` polygons of the level-3 products,
queried live with ``astroquery.mast`` (one query; the obs_ids and polygon counts are stored).

Optionally (``--noirot``) each record is matched to the Noirot+2023 NIRISS WFSS redshift catalogue
(``redshift_catalog.txt`` from https://niriss.github.io/smacs0723; sha256 checked against
``NOIROT_SHA256``). ``Q_FLAG`` is 1 = multiple features, 2 = single feature (paper Table B1).

Being inside the footprint is necessary, not sufficient, for a useful measurement (ASSUMPTION).

    python scripts/smacs_niriss_miri_coverage.py --noirot <dir>/redshift_catalog.txt \\
        --out results/smacs_niriss_miri_coverage/summary.json
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grism_4043_coverage import coverage, parse_polygons  # noqa: E402

INSTRUMENTS = ("NIRISS/IMAGE", "NIRISS/WFSS", "MIRI/IMAGE")
RECORDS = Path("results/grism_4043/coverage.json")
NOIROT_SHA256 = "c1218230df02a6c1d396b0a3981d07672703441ca9dc76b935c200613264dfc3"
MATCH_RADIUS_ARCSEC = 0.5  # ASSUMPTION: same-object match radius for the redshift catalogue


def footprints(table) -> dict[str, list[np.ndarray]]:
    """Level-3 rows of ``table`` → {obs_id: polygons}."""
    fp: dict[str, list[np.ndarray]] = {}
    for row in table:
        if int(row["calib_level"]) == 3:
            fp.setdefault(str(row["obs_id"]), []).extend(parse_polygons(str(row["s_region"])))
    return dict(sorted(fp.items()))


def covered_by(cov: dict) -> list[str]:
    return sorted(k for k, v in cov.items() if v["inside"])


def noirot_matches(path: Path, recs: list[dict]) -> list[dict]:
    """Nearest Noirot+2023 entry per record and whether it lies within ``MATCH_RADIUS_ARCSEC``."""
    import hashlib

    from astropy.coordinates import SkyCoord
    from astropy.io import ascii

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != NOIROT_SHA256:
        raise ValueError(f"{path}: sha256 {digest} != pinned {NOIROT_SHA256}")
    t = ascii.read(path, format="basic", comment="#")
    cat = SkyCoord(t["RA"], t["DEC"], unit="deg")
    out = []
    for r in recs:
        sep = SkyCoord(r["ra"], r["dec"], unit="deg").separation(cat).arcsec
        i = int(np.argmin(sep))

        def val(col: str, i: int = i):
            v = t[col][i]
            return None if np.ma.is_masked(v) else (v.item() if hasattr(v, "item") else v)

        out.append(
            {
                "nearest_id": int(t["ID"][i]),
                "sep_arcsec": round(float(sep[i]), 3),
                "matched": bool(sep[i] <= MATCH_RADIUS_ARCSEC),
                "z_grism": val("Z_GRISM"),
                "zerr_grism": val("ZERR_GRISM"),
                "q_flag": val("Q_FLAG"),
                "z_spec": val("Z_SPEC"),
                "ref_spec": val("REF_SPEC"),
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--noirot", type=Path, help="Noirot+2023 redshift_catalog.txt")
    args = ap.parse_args()

    from astroquery.mast import Observations

    t = Observations.query_criteria(proposal_id="2736", instrument_name=list(INSTRUMENTS))
    fp = footprints(t)
    recs = json.loads(RECORDS.read_text())["records"]
    rows = []
    for r in recs:
        cov = coverage(r["ra"], r["dec"], fp)
        rows.append(
            {
                "name": r["name"],
                "ra": r["ra"],
                "dec": r["dec"],
                "source": r["source"],
                "covered_by": covered_by(cov),
                "coverage": cov,
            }
        )
    if args.noirot:
        for r, m in zip(rows, noirot_matches(args.noirot, rows), strict=True):
            r["noirot2023"] = m
    summary = {
        "query": {
            "service": "MAST Observations.query_criteria (astroquery)",
            "proposal_id": "2736",
            "instrument_name": list(INSTRUMENTS),
            "utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        },
        "label": "observed footprints; derived inside/outside",
        "n_polygons": {k: len(v) for k, v in fp.items()},
        "records_source": str(RECORDS),
        "n_records": len(rows),
        "n_covered_any": sum(bool(r["covered_by"]) for r in rows),
        "noirot2023": None
        if not args.noirot
        else {
            "url": "https://niriss.github.io/smacs0723/redshift_catalog.txt",
            "sha256": NOIROT_SHA256,
            "match_radius_arcsec": MATCH_RADIUS_ARCSEC,
            "q_flag_meaning": {"1": "multiple features", "2": "single feature"},
        },
        "records": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=1) + "\n")
    for r in rows:
        print(r["name"].split("_", 3)[-1], len(r["covered_by"]), r.get("noirot2023"))


if __name__ == "__main__":
    main()

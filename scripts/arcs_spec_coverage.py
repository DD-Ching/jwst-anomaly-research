"""`arcs`: public spectral coverage of the open anti rows (D-024 addenda 13-14 handoff).

Plan (frozen before any query): ``results/arcs/spec-coverage/PLAN.md``. Metadata only: ESO ObsCore
(TAP, footprint ``INTERSECTS``) and MAST CAOM (3" cone, footprint containment tested here). Gate Q
(controls) runs first; targets get a class only when it passes. Archive rows are ``observed``;
containment, classes and the outcome are ``derived``; "footprint contains the point = coverage" and
the instrument groups are ASSUMPTIONs.

    python scripts/arcs_spec_coverage.py
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "spec-coverage"
ESO_TAP = "https://archive.eso.org/tap_obs"
MAST_RADIUS = 3.0  # arcsec cone before the footprint test

TARGETS = (
    ("RXC J0600", 581, 90.0304080782502, -20.147417436311873, True),
    ("MACS0308", 2157, 47.23285635900545, 26.755093856038467, True),
    (
        "WHL0137",
        2194,
        24.341443283127845,
        -8.454895554127202,
        False,
    ),  # explained; reported, no label
)
CONTROL_POS = ("Abell 2744 core", 3.5862, -30.4002)
CONTROL_NEG = ("empty sky", 90.5, -35.0)

# Instrument groups (ASSUMPTION). Optical/NIR IFUs and slitless modes give a spectrum of every
# covered source; MSA / fixed slits only of the targets they were placed on.
ESO_IFU = {
    "MUSE",
    "KMOS",
    "SINFONI",
    "ERIS",
    "VIMOS",
}  # VIMOS counted only for its IFU mode (below)
SLITLESS_RE = re.compile(r"WFSS|GRISM|G102|G141|G800L|G280|NIS/WFSS|NRC_WFSS|GR150", re.I)
MSA_RE = re.compile(r"MSA|NRS_FIXEDSLIT|FIXED ?SLIT|NRS_IFU|IFU", re.I)


def parse_polygons(s_region: str) -> list[np.ndarray]:
    """Return the POLYGON vertex arrays (deg) of an STC-S string; CIRCLE becomes a 64-gon."""
    polys: list[np.ndarray] = []
    if not s_region:
        return polys
    tokens = s_region.replace(",", " ").split()
    i = 0
    while i < len(tokens):
        tok = tokens[i].upper()
        if tok in ("POLYGON", "CIRCLE"):
            j = i + 1
            # skip frame / refpos words
            while j < len(tokens) and not _isnum(tokens[j]):
                j += 1
            nums = []
            while j < len(tokens) and _isnum(tokens[j]):
                nums.append(float(tokens[j]))
                j += 1
            if tok == "POLYGON" and len(nums) >= 6:
                polys.append(np.array(nums[: len(nums) // 2 * 2]).reshape(-1, 2))
            elif tok == "CIRCLE" and len(nums) >= 3:
                ra0, dec0, r = nums[:3]
                t = np.linspace(0, 2 * np.pi, 64, endpoint=False)
                cosd = np.cos(np.radians(dec0))
                polys.append(np.c_[ra0 + r * np.cos(t) / cosd, dec0 + r * np.sin(t)])
            i = j
        else:
            i += 1
    return polys


def _isnum(s: str) -> bool:
    try:
        float(s)
    except ValueError:
        return False
    return True


def contains(poly: np.ndarray, ra: float, dec: float) -> bool:
    """Point-in-polygon on the tangent plane at (ra, dec) (even-odd rule; small footprints)."""
    dra = (poly[:, 0] - ra + 180.0) % 360.0 - 180.0
    x = dra * np.cos(np.radians(dec))
    y = poly[:, 1] - dec
    inside = False
    n = len(x)
    for k in range(n):
        x1, y1, x2, y2 = x[k], y[k], x[(k + 1) % n], y[(k + 1) % n]
        if (y1 > 0) != (y2 > 0):
            xc = x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
            if xc > 0:
                inside = not inside
    return inside


def region_contains(s_region: str, ra: float, dec: float) -> bool:
    return any(contains(p, ra, dec) for p in parse_polygons(s_region))


def mode_group(archive: str, instrument: str, mode: str, dataproduct: str) -> str:
    """'ifu_slitless', 'slit' or 'other' (ASSUMPTION; see module constants)."""
    ins = (instrument or "").upper()
    mode = mode or ""
    if archive == "ESO":
        if ins == "VIMOS":
            return "ifu_slitless" if (dataproduct or "").lower() == "cube" else "slit"
        if ins in ESO_IFU:
            return "ifu_slitless"
        if (dataproduct or "").lower() == "spectrum":
            return "slit"
        return "other"  # ALMA and other radio / imaging cubes: not an optical/NIR spectrum
    if SLITLESS_RE.search(f"{ins} {mode}"):
        return "ifu_slitless"
    if "NRS_IFU" in mode.upper() or ins.upper().startswith("NIRSPEC/IFU"):
        return "ifu_slitless"
    if MSA_RE.search(f"{ins} {mode}") or ins.upper().startswith(("NIRSPEC", "STIS", "COS")):
        return "slit"
    return "other"


def classify(groups: list[str]) -> str:
    """Frozen per-target class from the groups of the public, covering observations."""
    if "ifu_slitless" in groups:
        return "spectrum_available"
    if "slit" in groups:
        return "msa_footprint_only"
    return "none"


def eso_rows(ra: float, dec: float, now_iso: str) -> list[dict]:
    import pyvo

    q = (
        "SELECT obs_id, dp_id, instrument_name, dataproduct_type, obs_release_date, target_name, "
        "s_region, t_exptime, em_min, em_max, proposal_id FROM ivoa.ObsCore WHERE "
        "dataproduct_type IN ('cube','spectrum') AND "
        f"INTERSECTS(s_region, POINT('ICRS',{ra},{dec}))=1"
    )
    t = pyvo.dal.TAPService(ESO_TAP).search(q).to_table()
    out = []
    for r in t:
        rel = str(r["obs_release_date"])
        out.append(
            {
                "archive": "ESO",
                "obs_id": str(r["obs_id"]),
                "dp_id": str(r["dp_id"]),
                "instrument": str(r["instrument_name"]),
                "mode": "",
                "dataproduct": str(r["dataproduct_type"]),
                "release": rel,
                "public": rel[:10] <= now_iso[:10],
                "target": str(r["target_name"]),
                "proposal": str(r["proposal_id"]),
                "exptime_s": float(r["t_exptime"]) if r["t_exptime"] is not np.ma.masked else None,
                "em_min_m": float(r["em_min"]) if r["em_min"] is not np.ma.masked else None,
                "em_max_m": float(r["em_max"]) if r["em_max"] is not np.ma.masked else None,
                "contains": True,  # INTERSECTS on the server
            }
        )
    return out


def mast_rows(ra: float, dec: float, now_mjd: float) -> list[dict]:
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Observations

    c = SkyCoord(ra * u.deg, dec * u.deg)
    t = Observations.query_criteria(
        coordinates=c, radius=MAST_RADIUS * u.arcsec, dataproduct_type=["spectrum", "cube"]
    )
    out = []
    for r in t:
        rel = float(r["t_obs_release"]) if r["t_obs_release"] is not np.ma.masked else np.inf
        out.append(
            {
                "archive": "MAST",
                "obs_id": str(r["obs_id"]),
                "instrument": str(r["instrument_name"]),
                "mode": str(r["filters"]),
                "collection": str(r["obs_collection"]),
                "dataproduct": str(r["dataproduct_type"]),
                "release_mjd": rel,
                "public": rel <= now_mjd,
                "target": str(r["target_name"]),
                "proposal": str(r["proposal_id"]),
                "exptime_s": float(r["t_exptime"]) if r["t_exptime"] is not np.ma.masked else None,
                "contains": region_contains(str(r["s_region"]), ra, dec),
            }
        )
    return out


def summarize(rows: list[dict]) -> dict:
    cov = [r for r in rows if r["contains"] and r["public"]]
    for r in rows:
        r["group"] = mode_group(r["archive"], r["instrument"], r.get("mode", ""), r["dataproduct"])
    groups = [r["group"] for r in cov]
    by = {}
    for r in cov:
        key = f"{r['archive']}:{r['instrument']}:{r['group']}"
        by[key] = by.get(key, 0) + 1
    return {
        "n_rows": len(rows),
        "n_public_covering": len(cov),
        "by_instrument": by,
        "class": classify(groups),
    }


def gate_q(pos: dict, neg: dict) -> dict:
    bi = pos["by_instrument"]
    muse = any(k.startswith("ESO:MUSE:") for k in bi)
    niriss = any(k.startswith("MAST:NIRISS") and k.endswith("ifu_slitless") for k in bi)
    return {
        "pos_muse": muse,
        "pos_niriss_wfss": niriss,
        "neg_none": neg["class"] == "none",
        "passed": bool(muse and niriss and neg["class"] == "none"),
    }


def main(argv: list[str] | None = None) -> int:
    from astropy.time import Time

    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    now = Time.now()
    now_iso, now_mjd = now.isot, now.mjd

    def query(ra, dec):
        rows = eso_rows(ra, dec, now_iso) + mast_rows(ra, dec, now_mjd)
        return rows, summarize(rows)

    _, pos = query(*CONTROL_POS[1:])
    _, neg = query(*CONTROL_NEG[1:])
    gate = gate_q(pos, neg)
    targets = []
    for field, label, ra, dec, labelled in TARGETS:
        rows, s = query(ra, dec)
        cls = (
            s["class"] if (gate["passed"] and labelled) else ("undetermined" if labelled else None)
        )
        targets.append(
            {
                "field": field,
                "label": label,
                "ra": ra,
                "dec": dec,
                "labelled": labelled,
                "summary": s,
                "row_class": cls,
                "covering": [r for r in rows if r["contains"] and r["public"]],
            }
        )
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {
        "plan": "results/arcs/spec-coverage/PLAN.md",
        "queried_utc": now_iso,
        "provenance": {
            "observed": "ESO ObsCore and MAST CAOM metadata rows",
            "derived": "footprint containment, groups, classes, outcome",
            "assumption": "footprint containment = coverage; instrument groups; MAST 3 arcsec cone",
        },
        "controls": {
            "positive": {"name": CONTROL_POS[0], **pos},
            "negative": {"name": CONTROL_NEG[0], **neg},
        },
        "gate_q": gate,
        "targets": targets,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    print(
        json.dumps(
            {
                "gate_q": gate,
                **{
                    f"{t['field']} {t['label']}": [t["row_class"], t["summary"]["by_instrument"]]
                    for t in targets
                },
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

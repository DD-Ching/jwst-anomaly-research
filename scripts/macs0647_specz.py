"""MACS0647: DJA NIRSpec spectroscopic redshifts for the Meena+2023 families, and the CLASH
LTM-Gauss maps at those redshifts.

Frozen plan: ``results/lens_consistency/macs0647_specz/PLAN.md``. Rules S1-S4 are the functions
below. The DJA table (``TABLE["url"]``) is streamed, not stored (cloud disk rule):

    curl -sS <TABLE url> | gunzip -c | python scripts/macs0647_specz.py /dev/stdin

"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import m0417_images as mi
import numpy as np
from astropy.table import Table

OUT = Path("results/lens_consistency/macs0647_specz")
IMAGES = Path("results/lens_consistency/macs0647_inputs/meena23_images.ecsv")
MODEL = "macs0647-clash-ltm"  # the usable model of D-034 addendum w
TABLE = {
    "name": "dja_msaexp_emission_lines_v4.6.csv.gz",
    "url": "https://s3.amazonaws.com/msaexp-nirspec/extractions/dja_msaexp_emission_lines_v4.6.csv.gz",
    "last_modified": "2026-09-15T13:01:52Z",
    "bytes": 171878708,
    "sha256": "ee64bb460f417d3fcf96fc53818c06ca0e27d784833063ae8eb48feb18c6ddf3",
    "n_rows": 99478,
    "accessed": "2026-10-24",
}
# Frozen rules (PLAN.md; ASSUMPTIONs)
MATCH_RADIUS_ARCSEC = 0.5
SPEC_GRADE = 3
Z_AGREE = 0.01  # in units of (1 + z)
S3_MIN_IMAGES = 4
S3_MIN_MATCHED = 0.8
MATCH_ARCSEC = 3.0


def _agree(z: list[float]) -> bool:
    """All values within Z_AGREE (1 + z) of each other (full spread, not distance to the median)."""
    z = np.asarray(z, float)
    return bool((np.max(z) - np.min(z)) / (1 + np.median(z)) <= Z_AGREE)


def image_z(rows: list[dict]) -> tuple[str, float | None]:
    """S1: the image's spec-z from its matched DJA rows (grade 3 only)."""
    z = [float(r["z_best"]) for r in rows if str(r.get("grade", "")).strip() == str(SPEC_GRADE)]
    if not z:
        return "none", None
    if not _agree(z):
        return "ambiguous", None
    return "spec_z", float(np.median(z))


def family_z(z_by_image: dict[str, float]) -> tuple[str, float | None]:
    """S2: the family's spec-z from its images' S1 values."""
    z = list(z_by_image.values())
    if not z:
        return "none", None
    if not _agree(z):
        return "inconsistent", None
    return "spec_z", float(np.median(z))


def gate_s3(d: np.ndarray) -> dict:
    """S3: control families' certain images at their spec-z; NaN = no prediction (unmatched)."""
    d = np.asarray(d, float)
    n = len(d)
    frac = float(np.mean(np.isfinite(d) & (d <= MATCH_ARCSEC))) if n else 0.0
    return {
        "n_images": n,
        "frac_within": frac,
        "passes": bool(n >= S3_MIN_IMAGES and frac >= S3_MIN_MATCHED),
    }


def readings(d: np.ndarray, s3_passes: bool) -> list[str]:
    """S4: per held-out certain image with a family spec-z."""
    if not s3_passes:
        return ["not_decided"] * len(d)
    d = np.asarray(d, float)
    return ["no_flag" if (np.isfinite(x) and x <= MATCH_ARCSEC) else "flag" for x in d]


def sep_arcsec(ra1, dec1, ra2, dec2) -> np.ndarray:
    ra1, dec1, ra2, dec2 = (np.radians(np.asarray(v, float)) for v in (ra1, dec1, ra2, dec2))
    c = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600


KEEP = (
    "file",
    "srcid",
    "ra",
    "dec",
    "grating",
    "root",
    "sn50",
    "zline",
    "z_best",
    "z_phot",
    "grade",
)


def match_rows(table_csv: Path, images: Table) -> tuple[dict[str, list[dict]], int]:
    """DJA rows within MATCH_RADIUS_ARCSEC of each image (all Meena+2023 rows), and the number
    of rows read."""
    ira, idec = np.asarray(images["ra"], float), np.asarray(images["dec"], float)
    ids = np.asarray(images["image_id"], str)
    out: dict[str, list[dict]] = {}
    n = 0
    with open(table_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            n += 1
            try:
                ra, dec = float(row["ra"]), float(row["dec"])
            except (TypeError, ValueError):
                continue
            s = sep_arcsec(ra, dec, ira, idec)
            for k in np.where(s <= MATCH_RADIUS_ARCSEC)[0]:
                r = {c: row.get(c) for c in KEEP}
                r["sep_arcsec"] = round(float(s[k]), 3)
                out.setdefault(str(ids[k]), []).append(r)
    return out, n


def shared_sources(matches: dict[str, list[dict]], spec_images: set[str]) -> dict[str, list[str]]:
    """DJA ``srcid``s whose rows set the S1 z of more than one image (not independent spectra)."""
    by_src: dict[str, set[str]] = {}
    for i, rows in matches.items():
        if i in spec_images:
            for r in rows:
                if str(r.get("grade", "")).strip() == str(SPEC_GRADE):
                    by_src.setdefault(str(r["srcid"]), set()).add(i)
    return {s: sorted(v) for s, v in by_src.items() if len(v) > 1}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("table_csv", type=Path)
    ap.add_argument("--grid-step", type=float, default=0.25)
    args = ap.parse_args(argv)
    allimg = Table.read(IMAGES, format="ascii.ecsv")
    allimg["image_id"] = np.asarray(allimg["image_id"], str)
    allimg["system"] = np.asarray(allimg["system"], str)
    matches, n_rows = match_rows(args.table_csv, allimg)
    certain = allimg[~np.asarray(allimg["uncertain"], bool)]
    s1 = {i: image_z(matches.get(i, [])) for i in np.asarray(allimg["image_id"], str)}
    fam = np.asarray(certain["system"], str)
    s2 = {}
    for f in dict.fromkeys(fam):
        zi = {
            i: s1[i][1]
            for i in np.asarray(certain["image_id"], str)[fam == f]
            if s1[i][0] == "spec_z"
        }
        s2[f] = (*family_z(zi), sorted(zi))
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "table": TABLE,
        "rows_read": n_rows,
        "rules": {
            "MATCH_RADIUS_ARCSEC": MATCH_RADIUS_ARCSEC,
            "SPEC_GRADE": SPEC_GRADE,
            "Z_AGREE": Z_AGREE,
            "S3_MIN_IMAGES": S3_MIN_IMAGES,
            "S3_MIN_MATCHED": S3_MIN_MATCHED,
            "MATCH_ARCSEC": MATCH_ARCSEC,
        },
        "matches": matches,
        "S1": {i: v[0] for i, v in s1.items() if v[0] != "none" or i in matches},
        "S1_z": {i: v[1] for i, v in s1.items() if v[1] is not None},
        "S2": {f: {"status": v[0], "z": v[1], "images": v[2]} for f, v in s2.items()},
        "shared_sources": shared_sources(matches, {i for i, v in s1.items() if v[0] == "spec_z"}),
    }
    zfam = {f: v[1] for f, v in s2.items() if v[0] == "spec_z"}
    if not zfam:
        summary["outcome"] = "no input"
    else:
        t = certain[np.isin(fam, list(zfam))]
        t["z"] = [zfam[f] for f in np.asarray(t["system"], str)]
        t["a"] = np.ones(len(t))
        solver = mi.Solver(MODEL, t, args.grid_step)
        ip = solver.image_plane(t)
        order = {i: k for k, i in enumerate(np.asarray(ip["image_id"], str))}
        ip = ip[[order[i] for i in np.asarray(t["image_id"], str)]]
        d = np.asarray(ip["dtheta_arcsec"], float)
        held = np.asarray(t["new_in_jwst"], bool)  # as macs0647_clash_heldout
        s3 = gate_s3(d[~held])
        r = readings(d[held], s3["passes"])
        summary["model"] = MODEL
        summary["model_sha256"] = solver.model.sha256
        summary["grid_step"] = args.grid_step
        summary["per_image"] = {
            str(i): {
                "system": str(s),
                "z_spec": float(z),
                "dtheta_arcsec": (None if not np.isfinite(x) else round(float(x), 3)),
                "n_predicted": int(n),
            }
            for i, s, z, x, n in zip(
                t["image_id"], t["system"], t["z"], d, ip["n_predicted"], strict=True
            )
        }
        summary["S3"] = s3
        # Post hoc (not in PLAN.md): each DJA source counted once; a family whose S2 z comes only
        # from sources an earlier family already used is left out of the control.
        used: set[str] = set()
        once = []
        for f, v in s2.items():
            if v[0] != "spec_z":
                continue
            src = {
                str(r["srcid"])
                for i in v[2]
                for r in matches.get(i, [])
                if str(r.get("grade", "")).strip() == str(SPEC_GRADE)
            }
            if src <= used:
                continue
            used |= src
            once.append(f)
        sel = ~held & np.isin(np.asarray(t["system"], str), once)
        summary["post_hoc_S3_sources_once"] = {
            **gate_s3(d[sel]),
            "families": sorted({str(x) for x in np.asarray(t["system"], str)[sel]}),
        }
        summary["S4"] = dict(zip(np.asarray(t["image_id"], str)[held], r, strict=True))
        flags = [i for i, v in summary["S4"].items() if v == "flag"]
        if not s3["passes"]:
            summary["outcome"] = "not decided"
        elif not held.any():
            summary["outcome"] = "no held-out input"
        else:
            summary["outcome"] = "flags" if flags else "no flag"
        summary["flags"] = flags
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(
        json.dumps(
            {k: summary.get(k) for k in ("S2", "S3", "S4", "outcome", "per_image")}, indent=1
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

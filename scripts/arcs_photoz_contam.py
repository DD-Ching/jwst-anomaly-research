"""`arcs`: photo-z contamination of the 3 single-component pooled anti rows (D-024 add. 9 handoff).

Plan (frozen before the run): ``results/arcs/photoz-contam/PLAN.md``. For each row: (1) a
spectroscopic redshift from NED (TAP ``NEDTAP.objdir``, 1.0") or DJA v4.6 (grade 3, 0.5"); (2) the
row's F606W - F814W offset from the cluster red sequence fitted to RELICS BPZ members; (3) the
RELICS BPZ zbmin as a second photo-z. ``label`` and ``outcome`` apply the frozen rules.
Catalogue values are ``observed``; fits, offsets, labels and the outcome are ``derived``;
thresholds are ASSUMPTIONs.

    python scripts/arcs_photoz_contam.py --relics-dir DIR --dja PATH
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arcs_spt0615 import read_bpz  # noqa: E402
from macs0647_specz import TABLE as DJA_TABLE  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "photoz-contam"
ARCS = Path(__file__).resolve().parents[1] / "results" / "arcs"
NED_TAP = "https://ned.ipac.caltech.edu/tap"

# Rows (field, label, primary-model run, RELICS cluster code, z_cl = lens_consistency z_lens).
ROWS = (
    ("WHL0137", 2194, "whl0137-relics/whl0137-relics-glafic", "whl0137m08", 0.566),
    ("RXC J0600", 581, "rxc0600-relics/rxc0600-relics-glafic", "rxc0600m20", 0.43),
    ("MACS0308", 2157, "macs0308-relics/macs0308-relics-glafic", "macs0308p26", 0.356),
)
# Frozen rules (PLAN.md; ASSUMPTIONs)
NED_RADIUS = 1.0  # arcsec
DJA_RADIUS = 0.5  # arcsec
RELICS_RADIUS = 1.0  # arcsec
DZ_MEMBER = 0.03  # z_spec <= z_cl + DZ_MEMBER (1 + z_cl) -> member_or_foreground
RS_DZ = 0.05  # fitting set |zb - z_cl| <= RS_DZ (1 + z_cl)
RS_MAG = 23.0  # F814W limit
STEL_MAX = 0.8
RS_NSIG = 2.0
BG_ZMIN = 0.1  # validation background set: zbmin > z_cl + BG_ZMIN
GATE_RECALL = 0.8
GATE_FBG = 0.3


def sep_arcsec(ra1, dec1, ra2, dec2) -> np.ndarray:
    d = np.deg2rad
    c = np.sin(d(dec1)) * np.sin(d(dec2)) + np.cos(d(dec1)) * np.cos(d(dec2)) * np.cos(
        d(np.asarray(ra1) - np.asarray(ra2))
    )
    return np.rad2deg(np.arccos(np.clip(c, -1, 1))) * 3600.0


def member_cut(z_cl: float) -> float:
    return z_cl + DZ_MEMBER * (1 + z_cl)


def spec_class(z_spec: float | None, z_cl: float) -> str | None:
    """Test 1: None without a spectroscopic redshift."""
    if z_spec is None or not np.isfinite(z_spec):
        return None
    return "member_or_foreground" if z_spec <= member_cut(z_cl) else "background"


def fit_red_sequence(mag: np.ndarray, col: np.ndarray, n_iter: int = 3, clip: float = 3.0) -> dict:
    """Straight line colour(mag), sigma-clipped; sigma is the clipped rms scatter."""
    keep = np.isfinite(mag) & np.isfinite(col)
    p = np.array([0.0, np.nanmedian(col[keep])])
    sig = float(np.nanstd(col[keep]))
    for _ in range(n_iter):
        p = np.polyfit(mag[keep], col[keep], 1)
        r = col - np.polyval(p, mag)
        sig = float(np.std(r[keep]))
        new = np.isfinite(mag) & np.isfinite(col) & (np.abs(r) <= clip * sig)
        if new.sum() < 5:
            break
        keep = new
    return {"slope": float(p[0]), "intercept": float(p[1]), "sigma": sig, "n_used": int(keep.sum())}


def rs_offset(rs: dict, mag: float, col: float) -> float:
    """Colour offset from the red sequence in units of sigma."""
    return float((col - (rs["slope"] * mag + rs["intercept"])) / rs["sigma"])


def gate(recall: float, f_bg: float) -> bool:
    return bool(recall >= GATE_RECALL and f_bg <= GATE_FBG)


def label(spec: str | None, rs: str | None, photoz_disagrees: bool | None) -> str:
    """Per-row label (PLAN.md). rs: red_sequence_like / off_sequence / None (undetermined).
    photoz_disagrees: None when there is no second photo-z."""
    if spec == "member_or_foreground":
        return "contaminated"
    if spec == "background":
        return "background_spec"
    if rs == "red_sequence_like" and photoz_disagrees is not False:
        return "contaminated"
    if rs == "off_sequence" and photoz_disagrees is not True:
        return "background_photometric"
    return "undetermined"


def outcome(labels: list[str]) -> str:
    n_c = sum(lab == "contaminated" for lab in labels)
    n_b = sum(lab.startswith("background") for lab in labels)
    if n_c >= 2:
        return (
            "open part of the pooled anti excess explained by photo-z contamination (exploratory)"
        )
    if n_c == 0 and n_b >= 2:
        return "photo-z contamination not supported; the rows stay open"
    return "mixed / undetermined"


def colour_mag(b: Table) -> tuple[np.ndarray, np.ndarray]:
    m6 = np.asarray(b["f606w_mag"], float)
    m8 = np.asarray(b["f814w_mag"], float)
    bad = (m6 <= 0) | (m6 >= 90) | (m8 <= 0) | (m8 >= 90)
    m6, m8 = np.where(bad, np.nan, m6), np.where(bad, np.nan, m8)
    return m8, m6 - m8


def red_sequence(b: Table, z_cl: float) -> dict:
    mag, col = colour_mag(b)
    zb, zmin = np.asarray(b["zb"], float), np.asarray(b["zbmin"], float)
    gal = np.asarray(b["stel"], float) < STEL_MAX
    fit_set = gal & (mag < RS_MAG) & (np.abs(zb - z_cl) <= RS_DZ * (1 + z_cl)) & np.isfinite(col)
    rs = fit_red_sequence(mag[fit_set], col[fit_set])
    off = (col - (rs["slope"] * mag + rs["intercept"])) / rs["sigma"]
    within = np.abs(off) <= RS_NSIG
    bg_set = gal & (mag < RS_MAG) & (zmin > z_cl + BG_ZMIN) & np.isfinite(col)
    rs.update(
        n_fit_set=int(fit_set.sum()),
        recall=float(within[fit_set].mean()) if fit_set.any() else float("nan"),
        n_bg_set=int(bg_set.sum()),
        f_bg=float(within[bg_set].mean()) if bg_set.any() else float("nan"),
    )
    rs["gate_passed"] = gate(rs["recall"], rs["f_bg"])
    return rs


def ned_redshifts(ra: float, dec: float) -> list[dict]:
    import pyvo

    q = (
        "SELECT prefname, ra, dec, z, zunc, zflag, prefphytype, z_bibcode FROM NEDTAP.objdir WHERE "
        f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {NED_RADIUS / 3600.0})) = 1"
    )
    t = pyvo.dal.TAPService(NED_TAP).run_sync(q).to_table()
    out = []
    for r in t:
        z = r["z"]
        out.append(
            {
                "prefname": str(r["prefname"]),
                "sep_arcsec": round(float(sep_arcsec(ra, dec, float(r["ra"]), float(r["dec"]))), 3),
                "z": None if np.ma.is_masked(z) else float(z),
                "zflag": "" if np.ma.is_masked(r["zflag"]) else str(r["zflag"]).strip(),
                "type": str(r["prefphytype"]).strip(),
                "z_bibcode": "" if np.ma.is_masked(r["z_bibcode"]) else str(r["z_bibcode"]),
            }
        )
    return out


def ned_spec_z(rows: list[dict]) -> float | None:
    """Nearest NED redshift whose flag is not photometric (NED zflag 'PHOT' marks photo-z)."""
    ok = [r for r in rows if r["z"] is not None and "PHOT" not in r["zflag"].upper()]
    return min(ok, key=lambda r: r["sep_arcsec"])["z"] if ok else None


def dja_rows(table_csv: Path, targets: list[tuple[float, float]]) -> list[list[dict]]:
    ra0 = np.array([t[0] for t in targets])
    de0 = np.array([t[1] for t in targets])
    out: list[list[dict]] = [[] for _ in targets]
    with open(table_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                ra, dec = float(row["ra"]), float(row["dec"])
            except (TypeError, ValueError):
                continue
            if np.min(np.abs(de0 - dec)) > 0.01:
                continue
            s = sep_arcsec(ra, dec, ra0, de0)
            for k in np.where(s <= DJA_RADIUS)[0]:
                out[k].append(
                    {
                        c: row.get(c)
                        for c in ("file", "srcid", "ra", "dec", "grating", "z_best", "grade")
                    }
                    | {"sep_arcsec": round(float(s[k]), 3)}
                )
    return out


def dja_spec_z(rows: list[dict]) -> float | None:
    ok = [r for r in rows if str(r.get("grade", "")).strip() == "3"]
    return float(min(ok, key=lambda r: r["sep_arcsec"])["z_best"]) if ok else None


def row_position(run: str, lab: int) -> dict:
    t = Table.read(ARCS / run / "arcs_anti.ecsv")
    r = t[np.asarray(t["label"]) == lab][0]
    return {
        k: float(r[k])
        for k in ("ra", "dec", "mag", "z_phot", "z160", "z840", "g_min", "ellipticity")
    }


def hst_offset(run: str) -> tuple[float, float]:
    """Stored BPZ-to-JWST shift (arcsec) where the field run has one, else (0, 0)."""
    p = ARCS / run.split("/")[0] / "photoz_offset.json"
    if p.exists():
        d = json.loads(p.read_text())
        if d.get("shift_applied"):
            return float(d["dra_arcsec"]), float(d["ddec_arcsec"])
    return 0.0, 0.0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--relics-dir", type=Path, required=True)
    ap.add_argument("--dja", type=Path, required=True, help="uncompressed DJA v4.6 csv")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    pos = [row_position(run, lab) for _, lab, run, _, _ in ROWS]
    dja = dja_rows(args.dja, [(p["ra"], p["dec"]) for p in pos])
    rows, rs_all = [], {}
    for (field, lab, run, code, z_cl), p, drows in zip(ROWS, pos, dja, strict=True):
        b = read_bpz(args.relics_dir / f"hlsp_relics_hst_acs-wfc3ir_{code}_multi_v1_cat.txt")
        rs = red_sequence(b, z_cl)
        rs_all[field] = rs
        dra, ddec = hst_offset(run)
        bra = np.asarray(b["RA"], float) + dra / 3600.0 / np.cos(np.deg2rad(p["dec"]))
        bde = np.asarray(b["Dec"], float) + ddec / 3600.0
        s = sep_arcsec(p["ra"], p["dec"], bra, bde)
        k = int(np.argmin(s))
        matched = bool(s[k] <= RELICS_RADIUS)
        mag, col = colour_mag(b)
        rel = {
            "id": int(b["id"][k]),
            "sep_arcsec": round(float(s[k]), 3),
            "matched": matched,
            "f814w_mag": float(mag[k]),
            "f606w_f814w": float(col[k]),
            "zb": float(b["zb"][k]),
            "zbmin": float(b["zbmin"][k]),
            "zbmax": float(b["zbmax"][k]),
            "odds": float(b["odds"][k]),
            "stel": float(b["stel"][k]),
        }
        ned = ned_redshifts(p["ra"], p["dec"])
        z_ned, z_dja = ned_spec_z(ned), dja_spec_z(drows)
        z_spec = z_dja if z_dja is not None else z_ned
        spec = spec_class(z_spec, z_cl)
        if matched and np.isfinite(col[k]) and np.isfinite(mag[k]):
            # post hoc (not in the plan): offset reported even where the gate failed; never labels
            rel["posthoc_rs_offset_sigma"] = round(rs_offset(rs, float(mag[k]), float(col[k])), 3)
        rs_cls = None
        if matched and rs["gate_passed"] and np.isfinite(col[k]) and np.isfinite(mag[k]):
            rel["rs_offset_sigma"] = rs_offset(rs, float(mag[k]), float(col[k]))
            rs_cls = (
                "red_sequence_like" if abs(rel["rs_offset_sigma"]) <= RS_NSIG else "off_sequence"
            )
        # second photo-z only where the stored photo-z was DJA eazy (MACS0308's is this BPZ)
        second = code != "macs0308p26" and matched
        disagrees = bool(rel["zbmin"] <= member_cut(z_cl)) if second else None
        rows.append(
            {
                "field": field,
                "label": lab,
                "z_cl": z_cl,
                "member_cut": round(member_cut(z_cl), 4),
                "stored": p,
                "relics": rel,
                "ned": ned,
                "dja": drows,
                "z_spec_ned": z_ned,
                "z_spec_dja": z_dja,
                "spec_class": spec,
                "rs_class": rs_cls,
                "photoz_disagrees": disagrees,
                "row_label": label(spec, rs_cls, disagrees),
            }
        )
    res = {
        "plan": "results/arcs/photoz-contam/PLAN.md",
        "provenance": {
            "observed": "RELICS HLSP catalogues (DOI 10.17909/T9SP45), NED objdir, DJA v4.6 rows, "
            "stored arcs_anti.ecsv rows",
            "derived": "red-sequence fits, offsets, labels, outcome",
            "assumption": "z_cl, radii, DZ_MEMBER, RS cuts, 2 sigma, gate thresholds",
        },
        "dja_table": DJA_TABLE,
        "red_sequence": rs_all,
        "rows": rows,
        "outcome": outcome([r["row_label"] for r in rows]),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    for r in rows:
        print(
            r["field"],
            r["label"],
            r["spec_class"],
            r["rs_class"],
            r["photoz_disagrees"],
            r["row_label"],
            r["relics"],
        )
    print(json.dumps(rs_all, indent=1))
    print(res["outcome"])
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())

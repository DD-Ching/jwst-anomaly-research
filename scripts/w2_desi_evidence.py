"""W2 DESI pairs evidence audit (results/w2_desi_pairs/evidence/PLAN.md): identity checks.

Metadata only; no spectrum is retrieved or read, and #472's statistic is not run.

  identity  the 13 rejected pairs with DESI spectra of both components (#472 pairs.json): their
            desi_dr1.zpix rows by TARGETID, Gaia DR3 sources within 5" and LS DR10 Tractor sources
            within 6" of the system position, and the match of each DESI component to them
  build     evidence.json: identity flags + pinned catalogue rows (D-064 inputs, Hennawi et al.
            2006) + the hand-curated literature.json

    python scripts/w2_desi_evidence.py identity
    python scripts/w2_desi_evidence.py build
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w2_desi_pairs as wd  # noqa: E402
import w12_niq as niq  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "w2_desi_pairs"
OUT = RES / "evidence"
CACHE = ROOT / "data" / "cache" / "w2_desi_pairs"
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
GAIA_R, LS_R = 5.0, 6.0  # search radii around the system position (arcsec)
MATCH_R = 0.5  # a DESI fibre position counts as one Gaia / Tractor source within this (arcsec)
CLOSE = 2.0  # #472 PLAN section 4 close-pair list (fibres are 1.5")


def tap_csv(url, query, retries=3):
    last = ""
    for k in range(retries):
        try:
            r = requests.post(
                url,
                data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query},
                timeout=300,
            )
            if r.status_code == 200 and "," in r.text.splitlines()[0]:
                lines = r.text.strip().splitlines()
                cols = [c.strip('"') for c in lines[0].split(",")]
                return [dict(zip(cols, ln.split(","), strict=True)) for ln in lines[1:]]
            last = f"HTTP {r.status_code}: {r.text[:300]}"
        except requests.RequestException as exc:
            last = type(exc).__name__
        time.sleep(5 * (k + 1))
    raise RuntimeError(f"TAP failed ({last}): {query[:120]}")


def targets():
    sy = json.loads((RES / "pairs.json").read_text())["systems"]
    return [s for s in sy if s["role"] == "rejected" and s.get("targetid_a")]


def nearest(ra, dec, rows, rakey="ra", deckey="dec"):
    best = None
    for r in rows:
        d = wd.sep_arcsec(ra, dec, float(r[rakey]), float(r[deckey]))
        if best is None or d < best[0]:
            best = (d, r)
    return best


def cmd_identity(args):
    tg = targets()
    ids = [s[k] for s in tg for k in ("targetid_a", "targetid_b")]
    zpix = tap_csv(
        wd.TAP,
        "SELECT targetid, survey, program, mean_fiber_ra, mean_fiber_dec, z, zwarn, spectype "
        f"FROM desi_dr1.zpix WHERE targetid IN ({','.join(ids)})",
    )
    circles = " OR ".join(
        "CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {s['ra']}, {s['dec']}, {GAIA_R / 3600}))=1"
        for s in tg
    )
    gaia = tap_csv(
        GAIA_TAP,
        "SELECT source_id, ra, dec, phot_g_mean_mag, parallax, parallax_over_error "
        f"FROM gaiadr3.gaia_source WHERE {circles}",
    )
    boxes = " OR ".join(  # #472's zpix box clauses (RA wrap included), on Tractor's columns
        c.replace("mean_fiber_ra", "ra").replace("mean_fiber_dec", "dec")
        for s in tg
        for c in wd._clauses((s["ra"], s["dec"], LS_R))
    )
    tractor = tap_csv(
        wd.TAP,
        "SELECT ra, dec, type, mag_g, mag_r, mag_z, ref_cat, ref_id FROM ls_dr10.tractor "
        f"WHERE brick_primary = 1 AND ({boxes})",
    )
    out = {
        "label": "observed metadata (desi_dr1.zpix and ls_dr10.tractor via NOIRLab Data Lab TAP; "
        "gaiadr3.gaia_source via ESA Gaia TAP); matches derived",
        "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "radii_arcsec": {"gaia": GAIA_R, "tractor": LS_R, "match": MATCH_R},
        "systems": {},
    }
    for s in tg:
        rec = {"ra": s["ra"], "dec": s["dec"], "sep_cat": s["sep_cat"], "sep_desi": s["sep_desi"]}
        g = [
            r
            for r in gaia
            if wd.sep_arcsec(s["ra"], s["dec"], float(r["ra"]), float(r["dec"])) <= GAIA_R
        ]
        t = [
            r
            for r in tractor
            if wd.sep_arcsec(s["ra"], s["dec"], float(r["ra"]), float(r["dec"])) <= LS_R
        ]
        rec["gaia"] = [
            {k: r[k] for k in ("source_id", "ra", "dec", "phot_g_mean_mag", "parallax_over_error")}
            for r in g
        ]
        rec["tractor"] = t
        for side in ("a", "b"):
            tid = s[f"targetid_{side}"]
            zr = [r for r in zpix if r["targetid"] == tid]
            ra, dec = float(zr[0]["mean_fiber_ra"]), float(zr[0]["mean_fiber_dec"])
            ng, nt = nearest(ra, dec, g), nearest(ra, dec, t)
            rec[side] = {
                "targetid": tid,
                "zpix": zr,
                "gaia_nearest": None
                if ng is None
                else {"source_id": ng[1]["source_id"], "d": round(ng[0], 2)},
                "tractor_nearest": None
                if nt is None
                else {
                    "d": round(nt[0], 2),
                    **{k: nt[1][k] for k in ("type", "mag_g", "mag_r", "mag_z")},
                },
            }
        out["systems"][s["name"]] = rec
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "identity.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        f"identity.json: {len(tg)} systems, {len(zpix)} zpix rows, {len(gaia)} Gaia, "
        f"{len(tractor)} Tractor"
    )


def flags(rec):
    """Identity flags (derived; descriptive, no classification)."""
    a, b = rec["a"], rec["b"]
    ga = a["gaia_nearest"] if a["gaia_nearest"] and a["gaia_nearest"]["d"] <= MATCH_R else None
    gb = b["gaia_nearest"] if b["gaia_nearest"] and b["gaia_nearest"]["d"] <= MATCH_R else None
    ta = (
        a["tractor_nearest"]
        if a["tractor_nearest"] and a["tractor_nearest"]["d"] <= MATCH_R
        else None
    )
    tb = (
        b["tractor_nearest"]
        if b["tractor_nearest"] and b["tractor_nearest"]["d"] <= MATCH_R
        else None
    )
    mags = [float(x["mag_g"]) for x in (ta, tb) if x and x["mag_g"] not in ("", "nan")]
    return {
        "gaia_a": ga["source_id"] if ga else None,
        "gaia_b": gb["source_id"] if gb else None,
        "both_match_distinct_gaia": bool(ga and gb and ga["source_id"] != gb["source_id"]),
        "tractor_match_a": ta,
        "tractor_match_b": tb,
        "n_tractor_6arcsec": len(rec["tractor"]),
        "dg_tractor": round(abs(mags[0] - mags[1]), 2) if len(mags) == 2 else None,
        "sep_desi_minus_cat": round(rec["sep_desi"] - rec["sep_cat"], 2),
        "close_le_2arcsec": rec["sep_desi"] <= CLOSE,
    }


def catalogue_rows(name, rec):
    """Pinned catalogue rows (the D-064 inputs and Hennawi et al. 2006) at this system."""
    rows = {}
    for f, key in (
        ("J_MNRAS_520_3305_table1.tsv", "Lemon2023"),
        ("J_AJ_135_496_table3.tsv", "SQLS-DR3 table3"),
        ("J_AJ_140_403_table2.tsv", "SQLS-DR5 table2"),
        ("J_AJ_143_119_table3.tsv", "SQLS-DR7 table3"),
    ):
        t = niq.read_tsv(CACHE / f)
        pos = (
            [(float(r["RAJ2000"]), float(r["DEJ2000"])) for r in t]
            if "RAJ2000" in t.colnames
            else [
                niq.sdss_name_radec(r["SDSS"]) for r in t
            ]  # SQLS rows: position from the SDSS name
        )
        hit = [
            r
            for r, (ra, dec) in zip(t, pos, strict=True)
            if wd.sep_arcsec(rec["ra"], rec["dec"], ra, dec) < 3
        ]
        if hit:
            rows[key] = [{c: str(r[c]) for c in t.colnames if str(r[c])} for r in hit]
    t = niq.read_tsv(CACHE / "J_AJ_131_1_binqso.tsv")
    pos = [
        (float(z["mean_fiber_ra"]), float(z["mean_fiber_dec"]))
        for z in (rec["a"]["zpix"][0], rec["b"]["zpix"][0])
    ]
    hit = []
    for r in t:
        for c in ("1", "2"):
            ra, dec = niq.sdss_name_radec(
                "J" + r[f"RA{c}"].replace(" ", "") + r[f"DE{c}"].replace(" ", "")
            )
            if any(wd.sep_arcsec(ra, dec, p[0], p[1]) < 1.0 for p in pos):
                hit.append(r)
                break
    if hit:
        rows["Hennawi2006"] = [{c: str(r[c]) for c in t.colnames if str(r[c])} for r in hit]
    return rows


def cmd_build(args):
    ident = json.loads((OUT / "identity.json").read_text())
    lit = json.loads((OUT / "literature.json").read_text())
    out = {
        "label": "derived (identity.json, pinned catalogue rows, literature.json); "
        "no spectrum read",
        "systems": {},
    }
    for name, rec in ident["systems"].items():
        out["systems"][name] = {
            "position": {
                "ra": rec["ra"],
                "dec": rec["dec"],
                "sep_cat": rec["sep_cat"],
                "sep_desi": rec["sep_desi"],
            },
            "desi": {
                s: {
                    k: rec[s]["zpix"][0][k]
                    for k in ("targetid", "mean_fiber_ra", "mean_fiber_dec", "z", "program")
                }
                for s in ("a", "b")
            },
            "identity": flags(rec),
            "catalogue": catalogue_rows(name, rec),
            "literature": lit["systems"].get(name),
        }
    missing = sorted(set(ident["systems"]) - set(lit["systems"]))
    if missing:
        raise SystemExit(f"literature.json lacks {missing}")
    (OUT / "evidence.json").write_text(json.dumps(out, indent=1) + "\n")
    print(f"evidence.json: {len(out['systems'])} systems")


ESO_TAP = "https://archive.eso.org/tap_obs/sync"
KOA_TAP = "https://koa.ipac.caltech.edu/TAP/sync"
KOA_TABLES = ("koa_nirc2", "koa_osiris", "koa_kcwi", "koa_lris", "koa_esi", "koa_mosfire")
ARCHIVE_R = 3.0  # the target position must fall inside the observation's footprint (arcsec)
KOA_R = 60.0  # KOA lists pointing centres only (arcsec)
KOA_EXTRA = {  # instrument set-up columns (imager vs spectrograph, filter, slit, exposure)
    "koa_osiris": ", object, imagmode, ifilter, itime",
    "koa_lris": ", object, slitname, grisname, graname, elaptime",
}
CADC_TAP = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync"


def _tap_any(url, query):
    """Rows of a sync TAP CSV query, or the error text (an archive that fails is 'not checked')."""
    try:
        r = requests.get(
            url,
            params={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query},
            timeout=180,
        )
    except requests.RequestException as exc:
        return f"error: {type(exc).__name__}"
    if r.status_code != 200:
        return f"error: HTTP {r.status_code}"
    import csv
    import io

    return list(csv.DictReader(io.StringIO(r.text)))


def cmd_archives(args):
    """Public observations at <= 3 chosen positions (metadata only): MAST (all missions), ESO,
    KOA and CADC."""
    from astroquery.mast import Observations

    ident = json.loads((OUT / "identity.json").read_text())["systems"]
    out = {
        "label": "observed archive metadata (MAST CAOM via astroquery, ESO ObsCore TAP, KOA TAP, "
        "CADC ObsCore TAP); "
        "no data product downloaded",
        "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "radii_arcsec": {"mast_eso_footprint": ARCHIVE_R, "koa_pointing": KOA_R},
        "systems": {},
    }
    keep_mast = (
        "obs_collection",
        "instrument_name",
        "filters",
        "dataproduct_type",
        "calib_level",
        "t_min",
        "t_exptime",
        "em_min",
        "em_max",
        "obs_id",
        "proposal_id",
        "dataRights",
        "target_name",
    )
    for name in args.names:
        rec = ident[name]
        ra, dec = rec["ra"], rec["dec"]
        res = {}
        try:
            t = Observations.query_region(f"{ra} {dec}", radius=f"{ARCHIVE_R} arcsec")
            res["mast"] = [{k: str(r[k]) for k in keep_mast if k in t.colnames} for r in t]
        except Exception as exc:  # noqa: BLE001 - recorded as 'not checked'
            res["mast"] = f"error: {type(exc).__name__}"
        res["eso"] = _tap_any(
            ESO_TAP,
            "SELECT dp_id, obs_collection, instrument_name, dataproduct_type, dataproduct_subtype, "
            "t_min, t_exptime, em_min, em_max, s_resolution, obs_release_date, target_name "
            "FROM ivoa.ObsCore "
            f"WHERE INTERSECTS(s_region, CIRCLE('ICRS', {ra}, {dec}, {ARCHIVE_R / 3600}))=1",
        )
        res["koa"] = {}
        cone = (
            "WHERE CONTAINS(POINT('J2000', ra, dec), "
            f"CIRCLE('J2000', {ra}, {dec}, {KOA_R / 3600}))=1"
        )
        for tb in KOA_TABLES:
            rows = _tap_any(
                KOA_TAP,
                "SELECT koaid, date_obs, ra, dec, targname, koaimtyp, progid, progpi, progtitl, "
                "propint"
                f"{KOA_EXTRA.get(tb, '')} FROM {tb} {cone}",
            )
            if isinstance(rows, str):  # a table without one of the extra columns
                rows = _tap_any(
                    KOA_TAP, f"SELECT koaid, date_obs, ra, dec, targname, progid FROM {tb} {cone}"
                )
            res["koa"][tb] = rows
        res["cadc"] = _tap_any(  # CFHT, Gemini, Subaru (via CADC), JCMT and others
            CADC_TAP,
            "SELECT obs_collection, instrument_name, dataproduct_type, t_min, t_exptime, em_min, "
            "em_max, "
            "s_resolution, obs_id, obs_release_date FROM ivoa.ObsCore "
            f"WHERE INTERSECTS(s_region, CIRCLE('ICRS', {ra}, {dec}, {ARCHIVE_R / 3600}))=1",
        )
        out["systems"][name] = res
        n_m = len(res["mast"]) if isinstance(res["mast"], list) else res["mast"]
        n_e = len(res["eso"]) if isinstance(res["eso"], list) else res["eso"]
        n_k = {k: (len(v) if isinstance(v, list) else v) for k, v in res["koa"].items()}
        n_c = len(res["cadc"]) if isinstance(res["cadc"], list) else res["cadc"]
        print(name, "MAST", n_m, "ESO", n_e, "KOA", n_k, "CADC", n_c)
    (OUT / "archives.json").write_text(json.dumps(out, indent=1) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("identity").set_defaults(func=cmd_identity)
    sub.add_parser("build").set_defaults(func=cmd_build)
    p = sub.add_parser("archives")
    p.add_argument("names", nargs="+")
    p.set_defaults(func=cmd_archives)
    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

"""W3 (negative-mass point lens, n = 1, eps < 0) feasibility trial on ZTF DR24 light curves.

Plan: ``results/w3_ztf/PLAN.md`` (fixed before any light curve was read). Data: IRSA's ZTF DR24 HATS
collection on S3 (``s3://ipac-irsa-ztf/ztf/enhanced/dr24/``, anonymous), read column-wise with
pyarrow; nothing is stored except the small tracked tables under ``results/w3_ztf/`` and a
session cache.

Stages (``python scripts/w3_ztf.py <stage>``):

- ``sample``: metadata-only star selection, groups and the development/holdout split;
- ``fetch``: the light-curve partition's rows for the sample and for the comparison stars (cache);
- ``coverage``: good epochs per (field, rcid, filter) and per-epoch limits from comparison stars;
- ``inject``: W3 and ordinary-control injections on development carriers, through pre-screen and
  vetting;
- ``search``: pre-screen + vetting on one half (``--half dev|holdout``).

Every number is ``derived`` (``simulated`` for injections); thresholds are ASSUMPTIONs in the plan.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import zlib
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

BUCKET = "ipac-irsa-ztf/ztf/enhanced/dr24"
LC_PART = f"{BUCKET}/lc/hats/ztf_dr24_lc-hats/dataset/Norder=5/Dir=0/Npix=3839/part0.snappy.parquet"
OBJ_PART = (
    f"{BUCKET}/objects/hats/ztf_dr24_objects-hats/dataset/Norder=4/Dir=0/Npix=959"
    "/part0.snappy.parquet"
)
LC_ORDER, LC_PIX = 5, 3839
OUT = Path(__file__).resolve().parents[1] / "results" / "w3_ztf"
CACHE = Path(os.environ.get("JWST_ANOMALY_DATA", "data")) / "cache" / "w3_ztf"
SAMPLE = OUT / "sample.json"
VDIR = OUT / "v3"
VERSION = 3
EDGE_MARGIN = (
    10.0  # arcsec inside an exposure's quadrant footprint for a non-detection to count (v4)
)


def configure(pix: int, vdir: str | None = None) -> None:
    """v4: another light-curve partition of the same order-4 pixel (independent validation). v5:
    ``vdir`` (e.g. "v5") writes the v4 method's outputs for the v1-v3 sample in its own folder."""
    global LC_PIX, LC_PART, OBJ_PART, CACHE, SAMPLE, VDIR, VERSION, N_MAX
    if pix == 3839 and not (vdir and int(vdir.lstrip("v")) >= 11):  # v11: own sample in 3839
        if vdir:
            VDIR, VERSION = OUT / vdir, int(vdir.lstrip("v"))
        return
    LC_PIX = pix
    LC_PART = LC_PART.replace("Npix=3839", f"Npix={pix}")
    CACHE = CACHE / str(pix)
    VDIR = OUT / (vdir or "v4")
    SAMPLE = VDIR / "sample.json"
    VERSION = int((vdir or "v4").lstrip("v"))
    if pix // 4 != 959:  # v6: another sky area; objects from the order-4 parent pixel
        OBJ_PART = OBJ_PART.replace("Npix=959", f"Npix={pix // 4}")
        N_MAX = N_MAX6


def oid_quad(oid: int) -> tuple[int, int, int]:
    """(field, rcid, filterid) encoded in a ZTF DR objectid: field, filter, ccd (2), qid (1), 8
    digits."""
    s = str(int(oid))
    field, filt, ccd, q = int(s[:-12]), int(s[-12]), int(s[-11:-9]), int(s[-9])
    return field, (ccd - 1) * 4 + q - 1, filt


def ambiguous(members) -> bool:
    """Two objectids of one group in the same quadrant and filter: two sources or a duplicate."""
    keys = [oid_quad(m) for m in members]
    return len(keys) != len(set(keys))


def in_quad(expo, ra: float, dec: float, margin: float = EDGE_MARGIN) -> np.ndarray:
    """Per exposure row: (ra, dec) inside the quadrant footprint (IRSA corners ra1..dec4, from each
    image's WCS) by at least ``margin`` arcsec. The corners are projected gnomonically about the
    star, so the great-circle edges are straight lines (v4 review fix). A non-finite corner
    counts as outside; a table without corner columns counts as inside (pre-v4)."""
    cols = getattr(expo, "colnames", None) or list(getattr(expo, "columns", []))
    if "ra1" not in cols:
        return np.ones(len(expo), bool)
    cra = np.radians(np.column_stack([np.asarray(expo[f"ra{i}"], float) for i in range(1, 5)]))
    cde = np.radians(np.column_stack([np.asarray(expo[f"dec{i}"], float) for i in range(1, 5)]))
    a0, d0 = np.radians(ra), np.radians(dec)
    dra = cra - a0
    cosc = np.sin(d0) * np.sin(cde) + np.cos(d0) * np.cos(cde) * np.cos(dra)
    x = np.cos(cde) * np.sin(dra) / cosc * 206264.806
    y = (np.cos(d0) * np.sin(cde) - np.sin(d0) * np.cos(cde) * np.cos(dra)) / cosc * 206264.806
    # order the corners by angle around their centre (counter-clockwise)
    ang = np.arctan2(y - y.mean(1, keepdims=True), x - x.mean(1, keepdims=True))
    o = np.argsort(ang, axis=1)
    x, y = np.take_along_axis(x, o, 1), np.take_along_axis(y, o, 1)
    x2, y2 = np.roll(x, -1, 1), np.roll(y, -1, 1)
    # signed distance of the origin (the star) to each edge, positive inside
    with np.errstate(invalid="ignore", divide="ignore"):
        d = (x * y2 - x2 * y) / np.hypot(x2 - x, y2 - y)
    return np.all(d >= margin, axis=1) & np.isfinite(d).all(axis=1) & (cosc > 0).all(axis=1)


# plan section 4
MAG_RANGE = (14.5, 17.5)
SHARP_MAX = 0.3
MIN_GOODOBS = 300
ISO_RADIUS = 5.0  # arcsec
ISO_DMAG = 2.5
GROUP_RADIUS = 1.5  # arcsec
N_MAX = 1000
N_MAX6 = 5000  # v6: up to 5,000 groups, split into calibration and search halves
SPLIT6 = "v6"  # v6 split salt: search if crc32("v6:" + key) is even, else calibration
SEED = 7
COMP_MAG = (14.0, 18.0)
# plan sections 3, 5, 6
GOOD_FRAC = 0.5
EPOCH_TOL = 0.0005  # d
DIM_MAG = 1.0
MIN_RUN = 3
MIN_NIGHTS = 3
MIN_BRACKET = 20
COMMON_FRAC = 0.2
INTER = (0.3, 1.0)
MAX_INTER = 1


def s3():
    import s3fs

    return s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-west-2"})


def holdout(key: int) -> bool:
    return zlib.crc32(str(int(key)).encode()) % 2 == 0


def in_partition(ra, dec) -> np.ndarray:
    import astropy.units as u
    from astropy_healpix import HEALPix

    hp = HEALPix(nside=2**LC_ORDER, order="nested")
    return hp.lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg) == LC_PIX


# ----------------------------------------------------------------------------------------- sample


def select(obj, old=None) -> tuple[list[dict], dict]:
    """Groups of the sample from an objects table (pandas), plan section 4; v4 also drops ambiguous
    groups and groups within 2 x GROUP_RADIUS of an ``old`` (ra, dec) position."""
    from scipy.spatial import cKDTree

    obj = obj[in_partition(obj["ra"].to_numpy(), obj["dec"].to_numpy())].reset_index(drop=True)
    dec0 = float(np.median(obj["dec"]))
    xy = np.column_stack(
        [obj["ra"].to_numpy() * np.cos(np.radians(dec0)) * 3600, obj["dec"].to_numpy() * 3600]
    )
    tree = cKDTree(xy)
    exclude = None
    if old is not None and len(old):
        old = np.asarray(old, float)
        exclude = cKDTree(
            np.column_stack([old[:, 0] * np.cos(np.radians(dec0)) * 3600, old[:, 1] * 3600])
        )
    ok = (
        (obj["fid"] == 2)
        & obj["refmag"].between(*MAG_RANGE)
        & (obj["refsharp"].abs() <= SHARP_MAX)
        & (obj["transient"] == 0)
        & (obj["ngoodobs"] >= MIN_GOODOBS)
    ).to_numpy()
    refmag = obj["refmag"].to_numpy()
    cand = np.flatnonzero(ok)
    iso = []
    for i in (
        cand
    ):  # isolated: no other star (beyond GROUP_RADIUS) within ISO_RADIUS brighter than m + ISO_DMAG
        nb = tree.query_ball_point(xy[i], ISO_RADIUS)
        if not any(
            np.hypot(*(xy[j] - xy[i])) > GROUP_RADIUS
            and np.isfinite(refmag[j])
            and refmag[j] < refmag[i] + ISO_DMAG
            for j in nb
        ):
            iso.append(i)
    iso = sorted(iso, key=lambda i: -int(obj["ngoodobs"][i]))
    used = set()
    groups = []
    n_amb = n_old = 0
    for i in iso:
        if i in used:
            continue
        mem = [j for j in tree.query_ball_point(xy[i], GROUP_RADIUS) if j not in used]
        used.update(mem)
        if VERSION >= 4:  # v4: no ambiguous groups, no physical source of the v1-v3 sample
            if ambiguous([int(obj["oid"][j]) for j in mem]):
                n_amb += 1
                continue
            if exclude is not None and len(exclude.query_ball_point(xy[i], 2 * GROUP_RADIUS)):
                n_old += 1
                continue
        groups.append(
            {
                "key": int(obj["oid"][i]),
                "ra": float(obj["ra"][i]),
                "dec": float(obj["dec"][i]),
                "refmag_r": float(refmag[i]),
                "members": [int(obj["oid"][j]) for j in mem],
                # v4 merge check: farthest member and nearest non-member object (arcsec)
                "member_sep_max": float(max(np.hypot(*(xy[j] - xy[i])) for j in mem)),
                "nearest_other": min(  # None: no other object within ISO_RADIUS
                    (
                        float(np.hypot(*(xy[j] - xy[i])))
                        for j in tree.query_ball_point(xy[i], ISO_RADIUS)
                        if j not in mem
                    ),
                    default=None,
                ),
            }
        )
    n_all = len(groups)
    rng = np.random.default_rng(SEED)
    pick = rng.permutation(n_all)[:N_MAX]
    groups = [groups[k] for k in sorted(pick)]
    for g in groups:
        g["holdout"] = holdout(g["key"]) if VERSION < 4 else False  # v4: all are validation
        if VERSION >= 6:  # v6: "holdout" is the search group, the rest the calibration group
            salt = SPLIT6 if VERSION < 11 else f"v{VERSION}"  # v11: its own split salt
            g["holdout"] = zlib.crc32(f"{salt}:{int(g['key'])}".encode()) % 2 == 0
    info = {
        "objects_in_area": int(len(obj)),
        "r_band_passing_cuts": int(ok.sum()),
        "isolated": len(iso),
        "groups_before_draw": int(n_all),
        "groups": len(groups),
        "holdout": int(sum(g["holdout"] for g in groups)),
        "ambiguous_groups_dropped": n_amb,
        "groups_dropped_near_old_sample": n_old,
    }
    return groups, info


def stage_sample(_a) -> None:
    import pyarrow.parquet as pq

    cols = ["oid", "ra", "dec", "field", "ccdid", "qid", "fid", "refmag", "refsharp"]
    cols += ["ngoodobs", "transient"]
    t0 = time.time()
    obj = pq.ParquetFile(s3().open(OBJ_PART)).read(columns=cols).to_pandas()
    old = None
    if VERSION >= 4:
        old = [(g["ra"], g["dec"]) for g in json.loads((OUT / "sample.json").read_text())["groups"]]
        if VERSION >= 11:  # v11: the v4 sample (Npix 3838) is used too
            old += [
                (g["ra"], g["dec"])
                for g in json.loads((OUT / "v4" / "sample.json").read_text())["groups"]
            ]
    groups, info = select(obj, old)
    if VERSION >= 6:
        info["calibration_groups"] = int(sum(not g["holdout"] for g in groups))
        info["search_groups"] = int(sum(g["holdout"] for g in groups))
    info["read_s"] = round(time.time() - t0, 1)
    info["lc_partition"] = LC_PART
    SAMPLE.parent.mkdir(parents=True, exist_ok=True)
    SAMPLE.write_text(
        json.dumps({"label": "derived; plan section 4", "info": info, "groups": groups}) + "\n"
    )
    excl = {m for g in groups for m in g["members"]} if VERSION >= 6 else set()
    n = save_comparisons(obj, excl)
    print(json.dumps(info), n)


def save_comparisons(obj, exclude=frozenset()) -> int:
    """Comparison stars for coverage: all objects in the area in the comparison magnitude range
    (v6: the control group, i.e. without any calibration or search group member)."""
    comp = obj[in_partition(obj["ra"].to_numpy(), obj["dec"].to_numpy())]
    comp = comp[comp["refmag"].between(*COMP_MAG) & (comp["transient"] == 0)]
    comp = comp[~comp["oid"].isin(list(exclude))]
    CACHE.mkdir(parents=True, exist_ok=True)
    np.save(CACHE / "comparison_oids.npy", comp["oid"].to_numpy())
    return len(comp)


def stage_comparisons(_a) -> None:
    """v5: only the comparison stars (5 object columns), for an existing sample."""
    import pyarrow.parquet as pq

    cols = ["oid", "ra", "dec", "refmag", "transient"]
    obj = pq.ParquetFile(s3().open(OBJ_PART)).read(columns=cols).to_pandas()
    excl = set()
    if VERSION >= 6:  # v7 diagnosis: rebuild the v6 control group (no sample member)
        excl = {m for g in json.loads(SAMPLE.read_text())["groups"] for m in g["members"]}
    print(save_comparisons(obj, excl))


# ------------------------------------------------------------------------------------------ fetch


def stage_fetch(_a) -> None:
    """Rows of the light-curve partition for sample members and comparison stars, one row group at a
    time (each ~250 MB in memory; the partition is never stored)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    groups = json.loads(SAMPLE.read_text())["groups"]
    members = np.array([m for g in groups for m in g["members"]], dtype=np.int64)
    comp = np.load(CACHE / "comparison_oids.npy")
    want = pa.array(np.unique(np.concatenate([members, comp])))
    f = pq.ParquetFile(s3().open(LC_PART))
    cols = ["objectid", "filterid", "fieldid", "rcid", "objra", "objdec", "lightcurve"]
    parts, t0 = [], time.time()
    for i in range(f.metadata.num_row_groups):
        t = f.read_row_group(i, columns=cols)
        parts.append(t.filter(pc.is_in(t["objectid"], value_set=want)))
        print(
            f"row group {i}: kept {parts[-1].num_rows} of {t.num_rows}, {time.time() - t0:.0f} s",
            flush=True,
        )
    tab = pa.concat_tables(parts)
    pq.write_table(tab, CACHE / "lc_rows.parquet")
    found = set(tab["objectid"].to_numpy())
    info = {
        "rows": tab.num_rows,
        "sample_members_found": int(sum(m in found for m in members)),
        "sample_members": int(len(members)),
        "read_s": round(time.time() - t0, 1),
    }
    (OUT / "fetch_info.json" if VERSION < 4 else VDIR / "fetch_info.json").write_text(
        json.dumps({"label": "derived", **info}) + "\n"
    )
    print(info)


# --------------------------------------------------------------------------------------- coverage


def load_rows():
    import pyarrow.parquet as pq

    t = pq.read_table(CACHE / "lc_rows.parquet")
    lc = t["lightcurve"].combine_chunks()
    rows = {}
    for i, o in enumerate(t["objectid"].to_numpy()):
        r = lc[i].as_py()
        rows[int(o)] = {
            "filt": int(t["filterid"][i].as_py()),
            "quad": (int(t["fieldid"][i].as_py()), int(t["rcid"][i].as_py())),
            "ra": float(t["objra"][i].as_py()),
            "dec": float(t["objdec"][i].as_py()),
            "t": np.asarray(r["hmjd"], float),
            "m": np.asarray(r["mag"], float),
            "e": np.asarray(r["magerr"], float),
            "f": np.asarray(r["catflags"], np.int64),
        }
    return rows


def PALOMAR():  # noqa: N802 (constant-like factory; no network lookup of the site registry)
    from astropy.coordinates import EarthLocation

    return EarthLocation.from_geodetic(-116.8597, 33.3563, 1712.0)


MATCH_TOL = 0.001  # d: light-curve hmjd to the exposure's heliocentric mid-time (v3)


def exposure_times(expo, ra: float, dec: float) -> np.ndarray:
    """Heliocentric MJD of each exposure's mid-time (obsjd = start) at (ra, dec). v8: half of the
    ``exptime`` column when present (private exposures last up to 300 s), else 30-s exposures."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time

    cols = getattr(expo, "colnames", None) or list(getattr(expo, "columns", []))
    half = np.asarray(expo["exptime"], float) / 2.0 if "exptime" in cols else 15.0
    t = Time(np.asarray(expo["obsjd"], float) + half / 86400.0, format="jd", scale="utc")
    ltt = t.light_travel_time(
        SkyCoord(ra * u.deg, dec * u.deg),
        kind="heliocentric",
        location=PALOMAR(),
    )
    return (t.utc + ltt).mjd


def match_exposures(t, et) -> np.ndarray:
    """Index into the (sorted) exposure times ``et`` of each time, -1 if none within MATCH_TOL."""
    t = np.asarray(t, float)
    if not len(et):
        return np.full(len(t), -1)
    if len(et) == 1:  # v8: the general branch below would return index -1 for a single exposure
        return np.where(np.abs(t - et[0]) <= MATCH_TOL, 0, -1)
    i = np.clip(np.searchsorted(et, t), 1, len(et) - 1)
    near = np.where(np.abs(t - et[i - 1]) <= np.abs(t - et[i]), i - 1, i)
    return np.where(np.abs(t - et[near]) <= MATCH_TOL, near, -1)


FILL: dict = {}  # v10: (field, rcid, filter) -> expids filled from psfcat (set by merge_psfcat)
PSF_R_MATCH, PSF_R_AMB = (
    1.5,
    3.0,
)  # arcsec: unique match / second source = ambiguous (v9 plan section 2)


def psfcat_match(sra, sdec, ra, dec) -> np.ndarray:
    """v10: per star (ra, dec) the index of the nearest psfcat source within PSF_R_MATCH arcsec,
    -1 if none (missing), -2 if a second source lies within PSF_R_AMB (ambiguous)."""
    sra, sdec = np.asarray(sra, float), np.asarray(sdec, float)
    out = np.full(len(ra), -1)
    for i, (a, d) in enumerate(zip(ra, dec, strict=True)):
        sep = np.hypot((sra - a) * np.cos(np.radians(d)), sdec - d) * 3600.0
        if not len(sep):
            continue
        k = np.argsort(sep)[:2]
        if sep[k[0]] > PSF_R_MATCH:
            continue
        out[i] = -2 if len(k) > 1 and sep[k[1]] <= PSF_R_AMB else k[0]
    return out


def merge_psfcat(rows, cats, expo) -> dict:
    """v10: add psfcat detections of fill exposures to ``rows`` (in place) and set FILL.

    ``cats`` maps (field, rcid, filter, expid) to (ra, dec, mag + MAGZP, sigmag, flags) arrays of
    one validated psfcat file. A star gets a point only where DR24 has no detection of it in that
    exposure (DR24 is never overwritten); flags != 0 or an ambiguous match enter as a flagged
    detection (catflags 1). A missing star gets nothing (state 4 in ``coverage``). Returns
    per-key counts."""
    fid = {"zg": 1, "zr": 2, "zi": 3}
    erow = {
        (
            int(r["field"]),
            (int(r["ccdid"]) - 1) * 4 + int(r["qid"]) - 1,
            fid[str(r["filtercode"])],
            int(r["expid"]),
        ): i
        for i, r in enumerate(expo)
    }
    by_q: dict = {}
    for o, r in rows.items():
        by_q.setdefault((*r["quad"], r["filt"]), []).append(o)
    log: dict = {}
    FILL.clear()
    for key, (sra, sdec, smag, serr, sflag) in cats.items():
        q, expid = key[:3], key[3]
        FILL.setdefault(q, set()).add(expid)
        oids = by_q.get(q, [])
        if not oids or key not in erow:
            continue
        e = expo[erow[key] : erow[key] + 1]
        # mid-time once per file at the stars' mean position, as coverage() does per quadrant
        ra0, de0 = np.mean([rows[o]["ra"] for o in oids]), np.mean([rows[o]["dec"] for o in oids])
        t = float(exposure_times(e, ra0, de0)[0])
        idx = psfcat_match(sra, sdec, [rows[o]["ra"] for o in oids], [rows[o]["dec"] for o in oids])
        c = log.setdefault(q, {"added_clean": 0, "added_flagged": 0, "missing": 0, "dr24_kept": 0})
        for o, j in zip(oids, idx, strict=True):
            r = rows[o]
            if len(r["t"]) and np.min(np.abs(r["t"] - t)) <= MATCH_TOL:
                c["dr24_kept"] += 1
                continue
            if j == -1:
                c["missing"] += 1
                continue
            m, err = (np.nan, np.nan) if j == -2 else (float(smag[j]), float(serr[j]))
            # ambiguous, flagged or a non-finite magnitude: a flagged detection (magnitude unused)
            bad = j == -2 or int(sflag[j]) != 0 or not np.isfinite(m)
            c["added_flagged" if bad else "added_clean"] += 1
            r["t"], r["m"] = np.append(r["t"], t), np.append(r["m"], m)
            r["e"], r["f"] = np.append(r["e"], err), np.append(r["f"], 1 if bad else 0)
    return log


def star_state(r, idx, n_exp):
    """Per exposure index: 1 detected with catflags 0 and not dim, 2 detected and dim (>= DIM_MAG
    fainter than the median), 3 flagged detection; 0 none. Duplicates in one exposure count once."""
    ok = r["f"] == 0
    med = float(np.median(r["m"][ok])) if ok.any() else np.nan
    st = np.zeros(n_exp, np.int8)
    mag = np.full(n_exp, np.nan)
    rank = {0: 0, 3: 1, 2: 2, 1: 3}  # clean beats dim beats flagged beats none
    for k, m, good in zip(idx, r["m"], ok, strict=True):
        if k < 0:
            continue
        v = 3 if not good else (2 if m >= med + DIM_MAG else 1)
        if rank[v] > rank[int(st[k])]:
            st[k], mag[k] = v, m
    return st, med, mag


def coverage(rows, comp: set[int], members: set[int], expo=None) -> dict:
    """v3: per (field, rcid, filter) the official exposures (IRSA metadata), each matched to the
    light-curve times; per exposure the number of comparison stars in play (detections before and
    after it), detected clean, and detected but dim, each star counted once (plan v3 A)."""
    from astropy.table import Table

    expo = expo if expo is not None else Table.read(CACHE / "exposures.ecsv")
    fid = {"zg": 1, "zr": 2, "zi": 3}
    by_q: dict = {}
    for o, r in rows.items():
        by_q.setdefault((*r["quad"], r["filt"]), []).append(o)
    cov = {}
    ccd = np.asarray(expo["ccdid"])
    qid = np.asarray(expo["qid"])
    fld = np.asarray(expo["field"])
    flt = np.array([fid[str(x)] for x in expo["filtercode"]])
    for q, oids in by_q.items():
        field, rcid, filt = q
        sel = (fld == field) & (ccd == rcid // 4 + 1) & (qid == rcid % 4 + 1) & (flt == filt)
        e = expo[sel]
        if not len(e):
            continue
        comps = [o for o in oids if o in comp]
        if not comps:
            continue
        ra = float(np.mean([rows[o]["ra"] for o in comps]))
        dec = float(np.mean([rows[o]["dec"] for o in comps]))
        et = exposure_times(e, ra, dec)
        order = np.argsort(et)
        e, et = e[order], et[order]
        n = len(et)
        inplay = np.zeros(n, int)
        clean = np.zeros(n, int)
        states, spans, covered, unmatched = {}, {}, {}, 0
        for o in oids:
            r = rows[o]
            idx = match_exposures(r["t"], et)
            unmatched += int((idx < 0).sum())
            st, med, mag = star_state(r, idx, n)
            lo, hi = (r["t"].min(), r["t"].max()) if len(r["t"]) else (np.inf, -np.inf)
            play = (et > lo) & (et < hi)
            inside = in_quad(e, r["ra"], r["dec"])  # v4: footprint of each exposure
            if q in FILL:  # v10: no psfcat row in a fill exposure is unknown, never a non-detection
                st[(st == 0) & inside & np.isin(np.asarray(e["expid"]), list(FILL[q]))] = 4
            states[o], spans[o], covered[o] = (st, med, mag), play, inside
            if o in comp:
                inplay += play & (inside | (st > 0))
                clean += play & (st == 1)
        info = np.asarray(e["infobits"])
        good = (info == 0) & (inplay > 0) & (clean >= GOOD_FRAC * np.maximum(inplay, 1))
        cov[q] = {
            "t": et,
            "expid": np.asarray(e["expid"]),
            "maglimit": np.asarray(e["maglimit"], float),
            "infobits": info,
            "good": good,
            "inplay": inplay,
            "clean": clean,
            "states": states,
            "spans": spans,
            "covered": covered,
            "comp": set(comps),
            "unmatched_detections": unmatched,
        }
    return cov


# --------------------------------------------------------------------------------------- series


def series(group, rows, cov, filt):
    """v3: time-ordered (t, dmag, is_nd, limit_dmag, quad, expid, base_mag, dimfrac) of a group in
    one filter over good exposures where the star is in play; dimfrac excludes the group's own
    members from the comparison counts."""
    out, flagged = [], set()
    mem = set(group["members"])
    for o in group["members"]:
        r = rows.get(o)
        if r is None or r["filt"] != filt:
            continue
        q = (*r["quad"], filt)
        c = cov.get(q)
        if c is None or o not in c["states"]:
            continue
        st, med, mag = c["states"][o]
        if (st == 1).sum() < 10:
            continue
        own_play = np.zeros(len(c["t"]), int)
        own_clean = np.zeros(len(c["t"]), int)
        cv = c.get("covered", {})
        ones = np.ones(len(c["t"]), bool)
        for m in mem & c["comp"]:
            sm = c["states"][m][0]
            own_play += c["spans"][m] & (cv.get(m, ones) | (sm > 0))
            own_clean += c["spans"][m] & (sm == 1)
        play_n = c["inplay"] - own_play
        clean_n = c["clean"] - own_clean
        # v4: a non-detection counts only where the star is inside the exposure's footprint
        mine = c["spans"][o] & c["good"] & (cv.get(o, ones) | (st > 0))
        for k in np.flatnonzero(mine):
            if st[k] == 4:
                continue  # v10: missing from a psfcat fill exposure: unknown for this star only
            if st[k] == 3:
                flagged.add(int(c["expid"][k]))  # v4: the exposure is unknown for the whole group
                continue  # flagged detection: unknown
            nd = st[k] == 0
            dm = 0.0 if nd else float(mag[k] - med)  # the magnitude star_state chose
            frac = (
                float((play_n[k] - clean_n[k]) / play_n[k]) if play_n[k] > 0 else float("nan")
            )  # unknown: fails
            out.append(
                (
                    float(c["t"][k]),
                    dm,
                    nd,
                    float(c["maglimit"][k] - med),
                    q,
                    int(c["expid"][k]),
                    med,
                    frac,
                )
            )
    out = [x for x in out if not (x[2] and x[5] in flagged)]
    out.sort(key=lambda x: x[0])
    return dedupe(out)


def dedupe(out):
    """v4: one entry per exposure (expid); a detection beats a non-detection."""
    best: dict = {}
    for x in out:
        k = x[5]
        if k not in best or (best[k][2] and not x[2]):
            best[k] = x
    return sorted(best.values(), key=lambda x: x[0])


# --------------------------------------------------------------------------------- pre-screen, vet


def dim_state(s):
    return np.array([x[2] or x[1] >= DIM_MAG for x in s], bool)


def runs(d):
    """(start, end) index pairs of consecutive True runs."""
    out, i = [], 0
    while i < len(d):
        if d[i]:
            j = i
            while j + 1 < len(d) and d[j + 1]:
                j += 1
            out.append((i, j))
            i = j + 1
        else:
            i += 1
    return out


def flagged_runs(s):
    """Plan section 5: >= MIN_RUN consecutive D epochs on >= MIN_NIGHTS nights, >= MIN_BRACKET
    baseline epochs on each side."""
    d = dim_state(s)
    t = np.array([x[0] for x in s])
    out = []
    for i, j in runs(d):
        nights = len(
            set(np.floor(t[i : j + 1] + 0.2).astype(int))
        )  # v2: night boundary ~local noon
        if j - i + 1 >= MIN_RUN and nights >= MIN_NIGHTS:
            if (~d[:i]).sum() >= MIN_BRACKET and (~d[j + 1 :]).sum() >= MIN_BRACKET:
                out.append((i, j))
    return out


def vet(s, run, other, cov, dimfrac_override=None) -> dict:
    """Plan section 6; returns each test's outcome."""
    i, j = run
    t = np.array([x[0] for x in s])
    dm = np.array([x[1] for x in s])
    nd = np.array([x[2] for x in s], bool)
    d = dim_state(s)
    res = {}
    # 1 common mode
    fr = []
    for x in s[i : j + 1]:
        k = x[5]
        f = (dimfrac_override or {}).get(k)
        fr.append(f if f is not None else (x[7] if len(x) > 7 else 0.0))
    res["common_mode_frac"] = float(np.mean(fr))
    res["common_mode_ok"] = res["common_mode_frac"] < COMMON_FRAC
    # 2 sharp edges
    span = max(t[j] - t[i], 1.0)
    win = ((t >= t[i] - span) & (t < t[i])) | ((t > t[j]) & (t <= t[j] + span))
    inter = win & ~nd & (dm >= INTER[0]) & (dm < INTER[1])
    res["n_intermediate"] = int(inter.sum())
    res["sharp_ok"] = res["n_intermediate"] <= MAX_INTER
    # 3 one-off: no other D-run of >= MIN_RUN epochs
    others = [r for r in runs(d) if r != (i, j) and r[1] - r[0] + 1 >= MIN_RUN]
    res["n_other_runs"] = len(others)
    res["oneoff_ok"] = not others
    # 4 achromatic coverage
    if other:
        to = np.array([x[0] for x in other])
        inside = (to >= t[i]) & (to <= t[j])
        dio = dim_state(other)[inside]
        res["other_filter_epochs"] = int(inside.sum())
        res["achromatic_ok"] = bool(dio.all()) if inside.any() else True
    else:
        res["other_filter_epochs"] = 0
        res["achromatic_ok"] = True
    res["pass"] = all(res[k] for k in ("common_mode_ok", "sharp_ok", "oneoff_ok", "achromatic_ok"))
    return res


def screen(s, other, cov, dimfrac_override=None):
    """All flagged runs and their vetting; the light curve passes if any run passes."""
    out = []
    for r in flagged_runs(s):
        v = vet(s, r, other, cov, dimfrac_override)
        out.append({"t_start": s[r[0]][0], "t_end": s[r[1]][0], "n": r[1] - r[0] + 1, **v})
    return out


def context(groups_half):
    """Rows and coverage. v6: only the control group and the members of ``groups_half`` are kept
    (the search group is never loaded during calibration); a usage record is appended to
    ``VDIR / "source_usage.jsonl"``."""
    rows = load_rows()
    comp = set(np.load(CACHE / "comparison_oids.npy").tolist())
    members = {m for g in json.loads(SAMPLE.read_text())["groups"] for m in g["members"]}
    if VERSION >= 6:
        if groups_half is None:  # every group (only after the validation passed; logged)
            groups_half = json.loads(SAMPLE.read_text())["groups"]
        use = {m for g in groups_half for m in g["members"]}
        rows = {o: r for o, r in rows.items() if o in comp or o in use}
        search = {
            m
            for g in json.loads(SAMPLE.read_text())["groups"]
            if g["holdout"]
            for m in g["members"]
        }
        rec = {
            "stage": sys.argv[1] if len(sys.argv) > 1 else "?",
            "groups": len(groups_half),
            "search_groups_in_use": int(sum(g["holdout"] for g in groups_half)),
            "rows_loaded": len(rows),
            "control_rows": int(sum(o in comp for o in rows)),
            "search_member_rows": int(sum(o in search for o in rows)),
        }
        VDIR.mkdir(parents=True, exist_ok=True)
        with open(VDIR / "source_usage.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
    cov = coverage(rows, comp, members)
    return rows, cov


# ------------------------------------------------------------------------------------- injection


def apply_dmag(s, ddm, cov=None, rng=None, nd_force=None):
    """New series: detections shifted by ddm (mag); epochs fainter than their limit become ND.
    With ``cov`` and ``rng``, noise is added for the fainter level from the epoch's magerr-mag
    relation (plan section 7; v2)."""
    out = []
    for x, dd, f in zip(
        s, ddm, nd_force if nd_force is not None else [False] * len(s), strict=True
    ):
        t, dm, nd, lim, q, k, base, *rest = x
        if nd:
            out.append(x)
            continue
        new = dm + dd
        if f or not np.isfinite(new) or new > lim:
            out.append((t, 0.0, True, lim, q, k, base, *rest))
            continue
        if rng is not None and dd != 0:
            # v3: magerr ~ 1.0857 / SNR with SNR = 5 at the exposure's maglimit (sky-limited)
            e_old = 0.2171 * 10 ** (0.4 * (dm - lim))
            e_new = 0.2171 * 10 ** (0.4 * (new - lim))
            new += rng.normal(0, np.sqrt(max(e_new**2 - e_old**2, 0.0)))
        out.append((t, new, False, lim, q, k, base, *rest))
    return out


def w3_dmag(t, t0, te, u0, rho, blend):
    from jwst_anomaly import exotic_sim as es

    f = np.asarray(
        es.inject_light_curve(t, 1.0, t0, te, u0, n=1.0, sign=-1, rho=rho, blend=blend)["flux"]
    )
    with np.errstate(divide="ignore"):
        return -2.5 * np.log10(np.clip(f, 0, None))


def control_model(kind, t, rng):
    """An ordinary-control light-curve model dmag(t) with parameters drawn once (v2: callable, so
    the other filter sees the same event)."""
    if kind == "eclipse":
        p, dur, ph0 = rng.uniform(1, 50), rng.uniform(0.1, 2.0), None
        ph0 = rng.uniform(0, p)
        return lambda x: np.where(((np.asarray(x) - ph0) % p) < dur, 1.5, 0.0)
    if kind == "fade":
        t1, td, tr = rng.uniform(t.min(), t.max()), rng.uniform(5, 20), rng.uniform(30, 200)

        def fade(x):
            x = np.asarray(x)
            return np.where(
                x < t1 - td,
                0.0,
                np.where(x < t1, 3.0 * (x - (t1 - td)) / td, 3.0 * np.exp(-(x - t1) / tr)),
            )

        return fade
    if kind == "pspl":
        t0, te, u0 = rng.uniform(t.min(), t.max()), rng.uniform(3, 100), rng.uniform(0, 1)

        def pspl(x):
            u = np.hypot(u0, (np.asarray(x) - t0) / te)
            return -2.5 * np.log10((u**2 + 2) / (u * np.sqrt(u**2 + 4)))

        return pspl
    raise ValueError(kind)


def stage_inject(a) -> None:
    """v2 (after review): independent seeds (W3 cell c: 1100 + c; controls 1200 + k), noise from the
    exposure's maglimit (v3), a hit must overlap the umbra t0 +- tE sqrt(4 - u0^2) and must not be a
    run already flagged in the uninjected carrier; controls use one event model in both filters.
    The common-mode control sets the comparison-star dim fraction to 0.5 at the injected epochs, so
    its rejection by the common-mode test holds by construction."""
    groups = [g for g in json.loads(SAMPLE.read_text())["groups"] if not g["holdout"]]
    rows, cov = context(groups)
    carriers = []
    for g in groups:
        s = series(g, rows, cov, 2)
        if len(s) >= 2 * MIN_BRACKET + MIN_RUN:
            o = series(g, rows, cov, 1)
            base = {(x["t_start"], x["t_end"]) for x in screen(s, o, cov)}
            carriers.append((g, s, o, base))
    res = {"version": VERSION, "carriers": len(carriers), "w3": [], "controls": {}}
    n = a.n
    cells = [
        (r, b, te)
        for r, b in ((0.01, 1.0), (0.1, 1.0), (0.01, 0.7), (0.1, 0.7))
        for te in (3.0, 10.0, 30.0, 100.0)
    ]
    for c, (rho, blend, te) in enumerate(cells):
        if a.cells == "calibrated" and (te == 3.0 or (rho == 0.1 and te == 100.0)):
            continue  # v6: only the 10 cells that reached 0.5 in v4 (seeds keep their cell index)
        rng = np.random.default_rng(1100 + c)
        rec = flag = 0
        used = set()
        for _ in range(n):
            g, s, o, base = carriers[rng.integers(len(carriers))]
            used.add(g["key"])
            t = np.array([x[0] for x in s])
            t0, u0 = rng.uniform(t.min(), t.max()), rng.uniform(0, 2)
            half = te * np.sqrt(4 - u0**2)
            si = apply_dmag(s, w3_dmag(t, t0, te, u0, rho, blend), cov, rng)
            to = np.array([x[0] for x in o]) if o else np.array([])
            oi = apply_dmag(o, w3_dmag(to, t0, te, u0, rho, blend), cov, rng) if o else o
            hit = [
                x
                for x in screen(si, oi, cov)
                if x["t_start"] <= t0 + half
                and x["t_end"] >= t0 - half
                and (x["t_start"], x["t_end"]) not in base
            ]
            flag += bool(hit)
            rec += any(x["pass"] for x in hit)
        res["w3"].append(
            {"rho": rho, "blend": blend, "tE": te, "n": n, "flagged": flag, "recovered": rec}
        )
        res["w3"][-1]["distinct_stars"] = len(used)
        print(res["w3"][-1], flush=True)
    for kk, kind in enumerate(("eclipse", "fade", "pspl", "common")):
        rng = np.random.default_rng(1200 + kk)
        flag = passed = 0
        used = set()
        for _ in range(n):
            g, s, o, base = carriers[rng.integers(len(carriers))]
            used.add(g["key"])
            t = np.array([x[0] for x in s])
            over = None
            if kind == "common":
                k0 = rng.integers(MIN_BRACKET, len(s) - MIN_BRACKET - 5)
                ddm = np.zeros(len(s))
                ddm[k0 : k0 + 5] = 3.0
                over = {s[k][5]: 0.5 for k in range(k0, k0 + 5)}
                si, oi = apply_dmag(s, ddm, cov, rng), o
            else:
                model = control_model(kind, t, rng)
                si = apply_dmag(s, model(t), cov, rng)
                oi = apply_dmag(o, model(np.array([x[0] for x in o])), cov, rng) if o else o
            sc = [x for x in screen(si, oi, cov, over) if (x["t_start"], x["t_end"]) not in base]
            flag += bool(sc)
            passed += any(x["pass"] for x in sc)
        res["controls"][kind] = {"n": n, "flagged": flag, "passed_vetting": passed}
        res["controls"][kind]["distinct_stars"] = len(used)
        print(kind, res["controls"][kind], flush=True)
    VDIR.mkdir(parents=True, exist_ok=True)
    name = "injections_dev_v3.json" if VERSION < 4 else f"injections_v{VERSION}.json"
    (VDIR / name).write_text(json.dumps({"label": "simulated", **res}, indent=1) + "\n")


def stage_search(a) -> None:
    want = a.half == "holdout"
    groups = [
        g
        for g in json.loads(SAMPLE.read_text())["groups"]
        if a.half == "all" or g["holdout"] == want
    ]
    rows, cov = context(groups)
    out, n_ser = [], 0
    for g in groups:
        for filt, other in ((2, 1), (1, 2)):
            s = series(g, rows, cov, filt)
            if len(s) < 2 * MIN_BRACKET + MIN_RUN:
                continue
            n_ser += 1
            for x in screen(s, series(g, rows, cov, other), cov):
                out.append({"key": g["key"], "ra": g["ra"], "dec": g["dec"], "filter": filt, **x})
    summ = {
        "version": VERSION,
        "half": a.half,
        "groups": len(groups),
        "series_searched": n_ser,
        "flagged_runs": len(out),
        "flagged_groups": len({x["key"] for x in out}),
        "passed_groups": sorted({x["key"] for x in out if x["pass"]}),
    }
    VDIR.mkdir(parents=True, exist_ok=True)
    (VDIR / f"search_{a.half}_v{VERSION}.json").write_text(
        json.dumps({"label": "derived", "summary": summ, "runs": out}, indent=1, default=float)
        + "\n"
    )
    print(json.dumps(summ))


def stage_exposures(_a) -> None:
    """Official exposure metadata (IRSA TAP ztf.ztf_current_meta_sci) for every quadrant of the
    fetched light curves (v3); cached, with a small tracked summary."""
    import io
    import urllib.parse
    import urllib.request

    import pyarrow.parquet as pq
    from astropy.table import Table, vstack

    t = pq.read_table(CACHE / "lc_rows.parquet", columns=["fieldid", "rcid"])
    quads = set(zip(t["fieldid"].to_numpy().tolist(), t["rcid"].to_numpy().tolist(), strict=True))
    parts, nbytes = [], 0
    for field in sorted({f for f, _ in quads}):
        rc = sorted(r for f, r in quads if f == field)
        cond = " OR ".join(f"(ccdid={r // 4 + 1} AND qid={r % 4 + 1})" for r in rc)
        sql = (
            "SELECT field,ccdid,qid,filtercode,obsjd,exptime,expid,infobits,maglimit,seeing,"
            "filefracday,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4 "
            f"FROM ztf.ztf_current_meta_sci WHERE field={field} AND ({cond})"
        )
        url = "https://irsa.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode(
            {"QUERY": sql, "FORMAT": "csv"}
        )
        for k in range(4):  # v6: bounded retry on resets (TAP, 2026-10-18)
            try:
                raw = urllib.request.urlopen(url, timeout=300).read()
                break
            except OSError:
                if k == 3:
                    raise
                time.sleep(10 * (k + 1))
        nbytes += len(raw)
        parts.append(Table.read(io.BytesIO(raw), format="csv"))
    expo = vstack(parts)
    expo.write(CACHE / "exposures.ecsv", overwrite=True)
    VDIR.mkdir(parents=True, exist_ok=True)
    summ = {"label": "observed (IRSA TAP)", "rows": len(expo), "quadrants": len(quads)}
    summ["tap_bytes"] = nbytes
    (VDIR / "exposures_summary.json").write_text(json.dumps(summ) + "\n")
    print(summ)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "stage", choices=["sample", "comparisons", "fetch", "exposures", "inject", "search"]
    )
    ap.add_argument("--half", choices=["dev", "holdout", "all"], default="dev")
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--cells", choices=["all", "calibrated"], default="all")
    ap.add_argument("--pix", type=int, default=3839, help="light-curve partition (v4: 3838)")
    ap.add_argument("--vdir", default=None, help="output folder for the v1-v3 sample (v5)")
    a = ap.parse_args(argv)
    if a.vdir and a.stage == "sample" and a.pix == 3839:
        ap.error("--vdir reuses the v1-v3 sample; the sample stage would overwrite it")
    if a.pix not in (3838, 3839) and not a.vdir:
        ap.error("another partition needs --vdir (v4 outputs would be overwritten)")
    configure(a.pix, a.vdir)
    fn = {
        "sample": stage_sample,
        "comparisons": stage_comparisons,
        "fetch": stage_fetch,
        "exposures": stage_exposures,
        "inject": stage_inject,
        "search": stage_search,
    }
    fn[a.stage](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

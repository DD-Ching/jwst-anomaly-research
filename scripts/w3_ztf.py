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

# plan section 4
MAG_RANGE = (14.5, 17.5)
SHARP_MAX = 0.3
MIN_GOODOBS = 300
ISO_RADIUS = 5.0  # arcsec
ISO_DMAG = 2.5
GROUP_RADIUS = 1.5  # arcsec
N_MAX = 1000
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


def select(obj) -> tuple[list[dict], dict]:
    """Groups of the sample from an objects table (pandas), plan section 4."""
    from scipy.spatial import cKDTree

    obj = obj[in_partition(obj["ra"].to_numpy(), obj["dec"].to_numpy())].reset_index(drop=True)
    dec0 = float(np.median(obj["dec"]))
    xy = np.column_stack(
        [obj["ra"].to_numpy() * np.cos(np.radians(dec0)) * 3600, obj["dec"].to_numpy() * 3600]
    )
    tree = cKDTree(xy)
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
    for i in iso:
        if i in used:
            continue
        mem = [j for j in tree.query_ball_point(xy[i], GROUP_RADIUS) if j not in used]
        used.update(mem)
        groups.append(
            {
                "key": int(obj["oid"][i]),
                "ra": float(obj["ra"][i]),
                "dec": float(obj["dec"][i]),
                "refmag_r": float(refmag[i]),
                "members": [int(obj["oid"][j]) for j in mem],
            }
        )
    n_all = len(groups)
    rng = np.random.default_rng(SEED)
    pick = rng.permutation(n_all)[:N_MAX]
    groups = [groups[k] for k in sorted(pick)]
    for g in groups:
        g["holdout"] = holdout(g["key"])
    info = {
        "objects_in_area": int(len(obj)),
        "r_band_passing_cuts": int(ok.sum()),
        "isolated": len(iso),
        "groups_before_draw": int(n_all),
        "groups": len(groups),
        "holdout": int(sum(g["holdout"] for g in groups)),
    }
    return groups, info


def stage_sample(_a) -> None:
    import pyarrow.parquet as pq

    cols = ["oid", "ra", "dec", "field", "ccdid", "qid", "fid", "refmag", "refsharp"]
    cols += ["ngoodobs", "transient"]
    t0 = time.time()
    obj = pq.ParquetFile(s3().open(OBJ_PART)).read(columns=cols).to_pandas()
    groups, info = select(obj)
    info["read_s"] = round(time.time() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "sample.json").write_text(
        json.dumps({"label": "derived; plan section 4", "info": info, "groups": groups}) + "\n"
    )
    # comparison stars for coverage (all objects in the area in the comparison magnitude range)
    comp = obj[in_partition(obj["ra"].to_numpy(), obj["dec"].to_numpy())]
    comp = comp[comp["refmag"].between(*COMP_MAG) & (comp["transient"] == 0)]
    CACHE.mkdir(parents=True, exist_ok=True)
    np.save(CACHE / "comparison_oids.npy", comp["oid"].to_numpy())
    print(json.dumps(info), len(comp))


# ------------------------------------------------------------------------------------------ fetch


def stage_fetch(_a) -> None:
    """Rows of the light-curve partition for sample members and comparison stars, one row group at a
    time (each ~250 MB in memory; the partition is never stored)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    groups = json.loads((OUT / "sample.json").read_text())["groups"]
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
    (OUT / "fetch_info.json").write_text(json.dumps({"label": "derived", **info}) + "\n")
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
            "t": np.asarray(r["hmjd"], float),
            "m": np.asarray(r["mag"], float),
            "e": np.asarray(r["magerr"], float),
            "f": np.asarray(r["catflags"], np.int64),
        }
    return rows


CLUSTER_GAP = 0.003  # d: hmjd of one exposure spreads ~1e-4 d across a quadrant (v2 fix)


def clusters(times) -> np.ndarray:
    """Median times of exposure clusters (sorted times split at gaps > CLUSTER_GAP)."""
    t = np.sort(np.asarray(times, float))
    if not len(t):
        return t
    cut = np.flatnonzero(np.diff(t) > CLUSTER_GAP) + 1
    return np.array([np.median(c) for c in np.split(t, cut)])


def assign(t, centres) -> np.ndarray:
    """Index of the exposure cluster of each time, -1 if none within CLUSTER_GAP."""
    t = np.asarray(t, float)
    if not len(centres):
        return np.full(len(t), -1)
    i = np.clip(np.searchsorted(centres, t), 1, len(centres) - 1)
    left = np.where(np.abs(t - centres[i - 1]) <= np.abs(t - centres[i]), i - 1, i)
    return np.where(np.abs(t - centres[left]) <= CLUSTER_GAP, left, -1)


def coverage(rows, comp: set[int], members: set[int]) -> dict:
    """Per (field, rcid, filter): exposure clusters, good epochs, per-epoch 5-sigma limit and
    magerr relation, and the fraction of comparison stars in the dim state per epoch (plan
    section 3; limit: deviation; v2: exposure clusters instead of fixed time bins, comparison
    stars include sample members as the plan says)."""
    by_q: dict = {}
    for o, r in rows.items():
        if o in comp:
            by_q.setdefault((*r["quad"], r["filt"]), []).append(r)
    cov = {}
    for q, rs in by_q.items():
        n = len(rs)
        cen = clusters(np.concatenate([r["t"] for r in rs]))
        det: dict = {}
        dim: dict = {}
        mags: dict = {}
        for r in rs:
            ok = r["f"] == 0
            med = np.median(r["m"][ok]) if ok.any() else np.nan
            for k, m, e, good in zip(assign(r["t"], cen), r["m"], r["e"], ok, strict=True):
                if good and k >= 0:
                    det[k] = det.get(k, 0) + 1
                    mags.setdefault(k, []).append((m, e))
                    if m >= med + DIM_MAG:
                        dim[k] = dim.get(k, 0) + 1
        good = {int(k) for k, c in det.items() if c >= GOOD_FRAC * n}
        lim, err = {}, {}
        for k in good:
            me = np.array(mags[k])
            b, a = (
                np.polyfit(me[:, 0], np.log10(me[:, 1]), 1)
                if len(me) >= 10 and np.ptp(me[:, 0]) > 0.5
                else (np.nan, np.nan)
            )
            ok_fit = np.isfinite(b) and b > 0
            lim[k] = float(np.clip((np.log10(0.2171) - a) / b, 18.0, 22.5)) if ok_fit else 20.5
            err[k] = (float(a), float(b)) if ok_fit else (np.log10(0.02) - 0.2 * 17, 0.2)
        cov[q] = {
            "n_comp": n,
            "centres": cen,
            "good": good,
            "limit": lim,
            "err": err,
            "dimfrac": {k: 1.0 - (det[k] - dim.get(k, 0)) / n for k in good},
        }
    return cov


# --------------------------------------------------------------------------------------- series


def series(group, rows, cov, filt):
    """Time-ordered (t, dmag, is_nd, limit_dmag, quad, epoch, base_mag) of a group in one filter
    over good epochs; dmag relative to each objectid's median (field zero points differ)."""
    out = []
    for o in group["members"]:
        r = rows.get(o)
        if r is None or r["filt"] != filt:
            continue
        q = (*r["quad"], filt)
        c = cov.get(q)
        if c is None:
            continue
        ok = r["f"] == 0
        if ok.sum() < 10:
            continue
        med = float(np.median(r["m"][ok]))
        seen = {}
        for k, m, good in zip(assign(r["t"], c["centres"]), r["m"], ok, strict=True):
            if k >= 0:
                seen[int(k)] = m - med if good else None  # None: flagged detection (unknown)
        for k in c["good"]:
            if k in seen and seen[k] is None:
                continue
            nd = k not in seen
            t = float(c["centres"][k])
            out.append((t, 0.0 if nd else seen[k], nd, c["limit"][k] - med, q, k, med))
    out.sort(key=lambda x: x[0])
    return out


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
        k, q = x[5], x[4]
        f = (dimfrac_override or {}).get(k)
        fr.append(f if f is not None else cov[q]["dimfrac"].get(k, 0.0))
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
    rows = load_rows()
    comp = set(np.load(CACHE / "comparison_oids.npy").tolist())
    members = {
        m for g in json.loads((OUT / "sample.json").read_text())["groups"] for m in g["members"]
    }
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
        t, dm, nd, lim, q, k, base = x
        if nd:
            out.append(x)
            continue
        new = dm + dd
        if f or not np.isfinite(new) or new > lim:
            out.append((t, 0.0, True, lim, q, k, base))
            continue
        if cov is not None and rng is not None and dd != 0:
            a, b = cov[q]["err"][k]
            e_old, e_new = 10 ** (a + b * (base + dm)), 10 ** (a + b * (base + new))
            new += rng.normal(0, np.sqrt(max(e_new**2 - e_old**2, 0.0)))
        out.append((t, new, False, lim, q, k, base))
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
    epoch's magerr relation, a hit must overlap the umbra t0 +- tE sqrt(4 - u0^2) and must not be a
    run already flagged in the uninjected carrier; controls use one event model in both filters.
    The common-mode control sets the comparison-star dim fraction to 0.5 at the injected epochs, so
    its rejection by the common-mode test holds by construction."""
    groups = [
        g for g in json.loads((OUT / "sample.json").read_text())["groups"] if not g["holdout"]
    ]
    rows, cov = context(groups)
    carriers = []
    for g in groups:
        s = series(g, rows, cov, 2)
        if len(s) >= 2 * MIN_BRACKET + MIN_RUN:
            o = series(g, rows, cov, 1)
            base = {(x["t_start"], x["t_end"]) for x in screen(s, o, cov)}
            carriers.append((g, s, o, base))
    res = {"version": 2, "carriers": len(carriers), "w3": [], "controls": {}}
    n = a.n
    cells = [
        (r, b, te)
        for r, b in ((0.01, 1.0), (0.1, 1.0), (0.01, 0.7), (0.1, 0.7))
        for te in (3.0, 10.0, 30.0, 100.0)
    ]
    for c, (rho, blend, te) in enumerate(cells):
        rng = np.random.default_rng(1100 + c)
        rec = flag = 0
        for _ in range(n):
            g, s, o, base = carriers[rng.integers(len(carriers))]
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
        print(res["w3"][-1], flush=True)
    for kk, kind in enumerate(("eclipse", "fade", "pspl", "common")):
        rng = np.random.default_rng(1200 + kk)
        flag = passed = 0
        for _ in range(n):
            g, s, o, base = carriers[rng.integers(len(carriers))]
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
        print(kind, res["controls"][kind], flush=True)
    (OUT / "injections_dev_v2.json").write_text(
        json.dumps({"label": "simulated", **res}, indent=1) + "\n"
    )


def stage_search(a) -> None:
    want = a.half == "holdout"
    groups = [
        g for g in json.loads((OUT / "sample.json").read_text())["groups"] if g["holdout"] == want
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
        "version": 2,
        "half": a.half,
        "groups": len(groups),
        "series_searched": n_ser,
        "flagged_runs": len(out),
        "flagged_groups": len({x["key"] for x in out}),
        "passed_groups": sorted({x["key"] for x in out if x["pass"]}),
    }
    (OUT / f"search_{a.half}_v2.json").write_text(
        json.dumps({"label": "derived", "summary": summ, "runs": out}, indent=1, default=float)
        + "\n"
    )
    print(json.dumps(summ))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=["sample", "fetch", "inject", "search"])
    ap.add_argument("--half", choices=["dev", "holdout"], default="dev")
    ap.add_argument("--n", type=int, default=200)
    a = ap.parse_args(argv)
    fn = {
        "sample": stage_sample,
        "fetch": stage_fetch,
        "inject": stage_inject,
        "search": stage_search,
    }
    fn[a.stage](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

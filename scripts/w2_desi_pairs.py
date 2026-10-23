"""W2: one quasar seen twice, or two quasars? DESI DR1 spectra of rejected lensed-quasar pairs
(``results/w2_desi_pairs/PLAN.md``, frozen before any spectrum was retrieved).

Stages (each writes into ``--out``, default ``results/w2_desi_pairs``):

    python scripts/w2_desi_pairs.py pairs     # sample (D-064 inputs + Hennawi binaries), DESI pairs
    python scripts/w2_desi_pairs.py fetch --roles control,binary   # SPARCL spectra, line windows
    python scripts/w2_desi_pairs.py validate  # control gate (PLAN section 4)
    python scripts/w2_desi_pairs.py fetch --roles rejected         # only after the gate passed
    python scripts/w2_desi_pairs.py apply     # targets, refused unless the gate passed

Spectra are ``observed`` (DESI DR1 via SPARCL, anonymous public access); the statistic and
classes are ``derived``. "Consistent with one source" is necessary for a dark lens, never
evidence of one.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w12_niq as niq  # noqa: E402

C_KMS = 299792.458
LINES = {"CIV": 1549.06, "CIII]": 1908.73, "MgII": 2798.75, "Hbeta": 4862.68, "Halpha": 6564.61}
WMIN, WMAX = 3650.0, 9800.0
CONT_V = (9000.0, 15000.0)  # |v| range of the continuum windows (km/s)
LINE_V = 6000.0  # |v| of the profile window (km/s)
MIN_PIX, MIN_SNR = 50, 10.0
T_D = 3.0  # PLAN section 4
GATE = {"min_controls": 6, "pos_frac": 0.80, "neg_frac": 0.70}
SEP_TOL, POS_TOL, EXTRA = 0.5, 1.0, 1.5  # arcsec
HEN_THETA_MAX, HEN_DV_MAX = 10.0, 2000.0
TAP = "https://datalab.noirlab.edu/tap/sync"
SPARCL = "https://astrosparcl.datalab.noirlab.edu/api"


# ---------------------------------------------------------------- statistic (PLAN section 3)
def line_profile(wave, flux, ivar, z, rest):
    """Profile p = F / C - 1 and its sigma on |v| <= LINE_V (continuum C); NaN where unusable."""
    v = C_KMS * (wave / (rest * (1 + z)) - 1)
    good = ivar > 0
    cw = good & (np.abs(v) >= CONT_V[0]) & (np.abs(v) <= CONT_V[1])
    if cw.sum() < 10 or not ((v[cw] < 0).any() and (v[cw] > 0).any()):
        return None
    coef = np.polyfit(v[cw], flux[cw], 1, w=np.sqrt(ivar[cw]))
    cont = np.polyval(coef, v)
    win = np.abs(v) <= LINE_V
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(good & (cont > 0), flux / cont - 1, np.nan)
        sp = np.where(good & (cont > 0), 1 / np.sqrt(ivar) / cont, np.nan)
    return v[win], p[win], sp[win]


def line_chi2(pa, sa, pb, sb):
    """Fit p_a ~ s p_b (s >= 0) by minimising chi2 with sigma^2 = sa^2 + s^2 sb^2; chi2_nu, s, N."""
    ok = np.isfinite(pa) & np.isfinite(pb) & np.isfinite(sa) & np.isfinite(sb)
    pa, sa, pb, sb = pa[ok], sa[ok], pb[ok], sb[ok]
    n = len(pa)
    if n < 2:
        return np.nan, np.nan, n
    grid = np.linspace(
        0, 10, 1001
    )  # s >= 0 (PLAN); an equivalent-width ratio above 10 is not a single source
    chi = np.array([np.sum((pa - s * pb) ** 2 / (sa**2 + s**2 * sb**2)) for s in grid])
    i = int(np.argmin(chi))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    fine = np.linspace(lo, hi, 201)
    chf = np.array([np.sum((pa - s * pb) ** 2 / (sa**2 + s**2 * sb**2)) for s in fine])
    j = int(np.argmin(chf))
    return float(chf[j] / (n - 1)), float(fine[j]), n


def line_snr(p, sp):
    ok = np.isfinite(p) & np.isfinite(sp)
    return float(np.sum(p[ok]) / np.sqrt(np.sum(sp[ok] ** 2))) if ok.any() else np.nan


def velocity_offset(pa, pb, dv_pix):
    """Cross-correlation peak of p_a and p_b over +-3000 km/s (reported only)."""
    ok = np.isfinite(pa) & np.isfinite(pb)
    a, b = np.where(ok, pa, 0.0), np.where(ok, pb, 0.0)
    kmax = int(3000 / dv_pix)
    lags = np.arange(-kmax, kmax + 1)
    cc = [
        np.sum(a[max(0, k) : len(a) + min(0, k)] * b[max(0, -k) : len(b) - max(0, k)]) for k in lags
    ]
    return float(lags[int(np.argmax(cc))] * dv_pix)


def pair_statistic(wave, fa, iva, fb, ivb, z):
    """PLAN section 3: per-line chi2_nu with a free scale; D = mean over counted lines."""
    out = {"lines": {}, "D": None, "dv_kms": None}
    ds, dvs = [], []
    for name, rest in LINES.items():
        lo, hi = rest * (1 + z) * (1 - CONT_V[1] / C_KMS), rest * (1 + z) * (1 + CONT_V[1] / C_KMS)
        if lo < WMIN or hi > WMAX:
            continue
        a, b = line_profile(wave, fa, iva, z, rest), line_profile(wave, fb, ivb, z, rest)
        if a is None or b is None:
            out["lines"][name] = {"used": False, "why": "no continuum window"}
            continue
        v, pa, sa = a
        _, pb, sb = b
        chi, s, n = line_chi2(pa, sa, pb, sb)
        snr_a, snr_b = line_snr(pa, sa), line_snr(pb, sb)
        used = n >= MIN_PIX and snr_a >= MIN_SNR and snr_b >= MIN_SNR and np.isfinite(chi)
        dv = velocity_offset(pa, pb, float(np.median(np.diff(v)))) if n >= MIN_PIX else None
        out["lines"][name] = {
            "used": bool(used),
            "chi2_nu": chi,
            "s": s,
            "n_pix": n,
            "snr_a": snr_a,
            "snr_b": snr_b,
            "dv_kms": dv,
        }
        if used:
            ds.append(chi)
            dvs.append(dv)
    if ds:
        out["D"] = float(np.mean(ds))
        out["dv_kms"] = float(np.median(dvs))
    return out


# ------------------------------------------------- sample and DESI pairs (PLAN section 2)
def tap(query, retries=3):
    last = ""
    for k in range(retries):
        try:
            r = requests.post(
                TAP,
                data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query},
                timeout=300,
            )
            if r.status_code == 200 and r.text.startswith("targetid"):
                lines = r.text.strip().splitlines()
                cols = lines[0].split(",")
                return [dict(zip(cols, ln.split(","), strict=True)) for ln in lines[1:]]
            last = f"HTTP {r.status_code}: {r.text[:300]}"
        except requests.RequestException as exc:
            last = type(exc).__name__
        time.sleep(5 * (k + 1))
    raise RuntimeError(f"TAP failed ({last}): {query[:120]}")


def _box(ra, dec, rad_arcsec):
    d = rad_arcsec / 3600.0
    dr = d / math.cos(math.radians(dec))
    return ra - dr, ra + dr, dec - d, dec + d


ZPIX = []  # rows prefetched by zpix_prefetch (batched box ORs; Data Lab takes no uploads)


def _clauses(bx):
    a, b, c, d = _box(*bx)
    spans = [(a, b)] if a >= 0 and b <= 360 else [(a % 360, 360), (0, b % 360)]  # RA wrap
    return [
        f"(mean_fiber_ra BETWEEN {lo} AND {hi} AND mean_fiber_dec BETWEEN {c} AND {d})"
        for lo, hi in spans
    ]


def zpix_prefetch(boxes, per_query=20):
    """One TAP query per batch of boxes (ra, dec, radius arcsec): ORed BETWEEN clauses."""
    for i in range(0, len(boxes), per_query):
        clauses = [cl for bx in boxes[i : i + per_query] for cl in _clauses(bx)]
        q = (
            "SELECT targetid, mean_fiber_ra, mean_fiber_dec, z, zwarn, spectype "
            "FROM desi_dr1.zpix WHERE " + " OR ".join(clauses)
        )
        ZPIX.extend(tap(q))


def zpix_near(ra, dec, rad_arcsec):
    """QSO rows with zwarn 0 inside the box, one per TARGETID, from the prefetched rows."""
    a, b, c, d = _box(ra, dec, rad_arcsec)
    half = 0.5 * (b - a)
    rows = [
        r
        for r in ZPIX
        if abs((float(r["mean_fiber_ra"]) - ra + 180) % 360 - 180) <= half
        and c <= float(r["mean_fiber_dec"]) <= d
        and r["spectype"] == "QSO"
        and int(r["zwarn"]) == 0
    ]
    uniq = {}
    for r in rows:  # one row per TARGETID (zpix can hold several survey / program coadds)
        uniq.setdefault(r["targetid"], r)
    return list(uniq.values())


def sep_arcsec(ra1, de1, ra2, de2):
    d1, d2 = math.radians(de1), math.radians(de2)
    c = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(
        math.radians(ra1 - ra2)
    )
    return math.degrees(math.acos(min(1.0, max(-1.0, c)))) * 3600


def best_pair(rows, sep_cat):
    best = None
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            s = sep_arcsec(
                float(rows[i]["mean_fiber_ra"]),
                float(rows[i]["mean_fiber_dec"]),
                float(rows[j]["mean_fiber_ra"]),
                float(rows[j]["mean_fiber_dec"]),
            )
            if abs(s - sep_cat) <= SEP_TOL and (best is None or abs(s - sep_cat) < best[0]):
                best = (abs(s - sep_cat), rows[i], rows[j], s)
    return best


def cmd_pairs(args):
    out = Path(args.out)
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    pins = niq.pinned()
    tables = {}
    for src in niq.INPUTS:
        f, _ = niq.fetch(src, cache, pins)
        tables[src] = niq.read_tsv(f)
    fb, _ = niq.fetch(niq.HENNAWI, cache, pins)
    binq = niq.read_tsv(fb)
    s = niq.build_sample(tables)
    systems = []
    for r in s:
        sep = float(r["sep_cat"])
        if np.isfinite(sep):
            systems.append(
                {
                    "role": "control" if r["group"] == "control" else "rejected",
                    "name": str(r["name"]),
                    "ra": float(r["ra"]),
                    "dec": float(r["dec"]),
                    "sep_cat": sep,
                    "comment": str(r["comment"]),
                }
            )
    ra1, de1 = niq._coords(binq["RA1"], binq["DE1"])
    ra2, de2 = niq._coords(binq["RA2"], binq["DE2"])
    th = np.array([niq._float(x) for x in binq["theta"]])
    dv = np.array([niq._float(x) for x in binq["dV"]])
    for i in range(len(binq)):
        if (
            all(np.isfinite([ra1[i], de1[i], ra2[i], de2[i], th[i], dv[i]]))
            and th[i] <= HEN_THETA_MAX
            and abs(dv[i]) < HEN_DV_MAX
        ):
            systems.append(
                {
                    "role": "binary",
                    "name": f"{binq['Name1'][i]}+{binq['Name2'][i]}".replace(" ", ""),
                    "ra": float(ra1[i]),
                    "dec": float(de1[i]),
                    "sep_cat": float(th[i]),
                    "comment": f"Hennawi2006 dV={dv[i]:.0f}",
                    "ra2": float(ra2[i]),
                    "dec2": float(de2[i]),
                }
            )
    boxes = [(sy["ra"], sy["dec"], sy["sep_cat"] + EXTRA) for sy in systems]
    boxes += [(sy["ra2"], sy["dec2"], POS_TOL + 0.5) for sy in systems if sy["role"] == "binary"]
    zpix_prefetch(boxes)
    rows = []
    for sy in systems:
        q = zpix_near(sy["ra"], sy["dec"], sy["sep_cat"] + EXTRA)
        if sy["role"] == "binary":
            a = [
                r
                for r in q
                if sep_arcsec(
                    sy["ra"], sy["dec"], float(r["mean_fiber_ra"]), float(r["mean_fiber_dec"])
                )
                <= POS_TOL
            ]
            q2 = zpix_near(sy["ra2"], sy["dec2"], POS_TOL + 0.5)
            b = [
                r
                for r in q2
                if sep_arcsec(
                    sy["ra2"], sy["dec2"], float(r["mean_fiber_ra"]), float(r["mean_fiber_dec"])
                )
                <= POS_TOL
            ]
            pair = (
                (
                    0,
                    a[0],
                    b[0],
                    sep_arcsec(
                        float(a[0]["mean_fiber_ra"]),
                        float(a[0]["mean_fiber_dec"]),
                        float(b[0]["mean_fiber_ra"]),
                        float(b[0]["mean_fiber_dec"]),
                    ),
                )
                if a and b
                else None
            )
        else:
            pair = best_pair(q, sy["sep_cat"])
        rec = {
            **{k: sy[k] for k in ("role", "name", "ra", "dec", "sep_cat", "comment")},
            "n_qso": len(q),
            "status": "pair" if pair else "unknown (no DESI pair)",
        }
        if pair:
            _, pa, pb, sdesi = pair
            rec.update(
                targetid_a=pa["targetid"],
                targetid_b=pb["targetid"],
                z_a=float(pa["z"]),
                z_b=float(pb["z"]),
                sep_desi=float(sdesi),
            )
        rows.append(rec)
    json.dump(
        {
            "label": "derived (D-064 pinned inputs, Hennawi et al. 2006, desi_dr1.zpix via "
            "NOIRLab Data Lab TAP)",
            "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "systems": rows,
        },
        open(out / "pairs.json", "w"),
        indent=1,
    )
    (out / "pairs.json").open("a").write("\n")
    for role in ("control", "binary", "rejected"):
        n = [r for r in rows if r["role"] == role]
        print(role, len(n), "with a DESI pair:", sum(r["status"] == "pair" for r in n))


# ---------------------------------------------------------------- spectra (SPARCL)
def sparcl_find(ra, dec, rad_arcsec):
    d = rad_arcsec / 3600.0
    dr = d / math.cos(math.radians(dec))
    body = {
        "outfields": ["sparcl_id", "specid", "ra", "dec", "redshift", "spectype", "data_release"],
        "search": [
            ["ra", ra - dr, ra + dr],
            ["dec", dec - d, dec + d],
            ["data_release", "DESI-DR1"],
        ],
    }
    return _post(f"{SPARCL}/find/?limit=100", body, 180)[1:]


def _post(url, body, timeout, retries=4):
    """SPARCL POST with back-off on 5xx / connection errors (a refusal is raised, never skipped)."""
    for k in range(retries):
        try:
            r = requests.post(url, json=body, timeout=timeout)
            if r.status_code < 500:
                r.raise_for_status()
                return r.json()
        except requests.ConnectionError:
            pass
        time.sleep(10 * 2**k)
    raise RuntimeError(f"SPARCL unavailable after {retries} tries: {url}")


def sparcl_get(ids):
    j = _post(
        f"{SPARCL}/spectras/?include=sparcl_id,flux,ivar,wavelength,specid,redshift,mask"
        "&format=json&dataset_list=DESI-DR1",
        ids,
        600,
    )
    if not j[0].get("status", {}).get("success"):
        raise RuntimeError(f"SPARCL refused: {j[0]}")
    return j[1:]


def windows_mask(wave, z):
    m = np.zeros(len(wave), bool)
    for rest in LINES.values():
        lo, hi = rest * (1 + z) * (1 - CONT_V[1] / C_KMS), rest * (1 + z) * (1 + CONT_V[1] / C_KMS)
        m |= (wave >= lo) & (wave <= hi)
    return m


PAIR_COLS = ("wave_a", "flux_a", "ivar_a", "wave_b", "flux_b", "ivar_b")


def save_pair(dst, arrays):
    """Line windows of both spectra as text (float32 values, round-trip exact via %.9g)."""
    import gzip

    n = len(arrays["wave_a"])
    with gzip.open(dst, "wt", compresslevel=9) as f:
        f.write(",".join(PAIR_COLS) + "\n")
        for i in range(n):
            f.write(",".join(f"{float(arrays[c][i]):.9g}" for c in PAIR_COLS) + "\n")


def load_pair(path):
    """Inverse of save_pair: float32 arrays, as stored."""
    t = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return {c: t[:, k].astype(np.float32) for k, c in enumerate(PAIR_COLS)}


def cmd_fetch(args):
    out = Path(args.out)
    roles = set(args.roles.split(","))
    if "rejected" in roles and not json.load(open(out / "gate.json"))["passed"]:
        raise SystemExit("gate not passed: targets are not retrieved (PLAN section 4)")
    sdir = out / "spectra"
    sdir.mkdir(exist_ok=True)
    P = json.load(open(out / "pairs.json"))
    flog = out / "fetch_log.json"
    log = json.load(open(flog))["entries"] if flog.exists() else []

    def save_log():
        json.dump(
            {
                "label": "SPARCL DESI-DR1 retrievals (observed spectra; line windows kept)",
                "entries": log,
            },
            open(flog, "w"),
            indent=1,
        )
        flog.open("a").write("\n")

    for sy in P["systems"]:
        if sy["role"] not in roles or sy["status"] != "pair":
            continue
        dst = sdir / f"{sy['role']}_{sy['name'].replace('/', '_')}.csv.gz"
        if dst.exists():
            continue
        recs = sparcl_find(sy["ra"], sy["dec"], sy["sep_cat"] + EXTRA + 1.0)
        want = {sy["targetid_a"]: "a", sy["targetid_b"]: "b"}
        hits = [r for r in recs if str(r["specid"]) in want]
        spec = sparcl_get([r["sparcl_id"] for r in hits]) if hits else []
        best = {}
        for sp in spec:  # one spectrum per TARGETID: the most good pixels (PLAN section 2)
            tid = str(sp["specid"])
            ngood = int(np.sum((np.asarray(sp["ivar"]) > 0) & (np.asarray(sp["mask"]) == 0)))
            if tid in want and (tid not in best or ngood > best[tid][0]):
                best[tid] = (ngood, sp, sp["sparcl_id"])
        if set(best) != set(want):
            log.append(
                {
                    "name": sy["name"],
                    "role": sy["role"],
                    "status": "unknown (SPARCL spectrum missing)",
                }
            )
            save_log()
            continue
        z = 0.5 * (sy["z_a"] + sy["z_b"])
        arrays, ids = {}, {}
        for tid, side in want.items():
            _, sp, sid = best[tid]
            w = np.asarray(sp["wavelength"], float)
            m = windows_mask(w, z)
            iv = np.where(np.asarray(sp["mask"]) == 0, np.asarray(sp["ivar"], float), 0.0)
            arrays[f"wave_{side}"] = w[m].astype(np.float32)
            arrays[f"flux_{side}"] = np.asarray(sp["flux"], float)[m].astype(np.float32)
            arrays[f"ivar_{side}"] = iv[m].astype(np.float32)
            ids[side] = sid
        save_pair(dst, arrays)
        log.append(
            {
                "name": sy["name"],
                "role": sy["role"],
                "status": "ok",
                "sparcl_ids": ids,
                "bytes": dst.stat().st_size,
                "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )
        save_log()
        time.sleep(1)
    print(
        sum(e["status"] == "ok" for e in log),
        "fetched;",
        sum(e["status"] != "ok" for e in log),
        "missing",
    )


def analyse(out, roles):
    P = json.load(open(out / "pairs.json"))
    res = []
    for sy in P["systems"]:
        if sy["role"] not in roles:
            continue
        rec = {
            k: sy.get(k)
            for k in ("role", "name", "sep_cat", "sep_desi", "z_a", "z_b", "comment", "status")
        }
        f = out / "spectra" / f"{sy['role']}_{sy['name'].replace('/', '_')}.csv.gz"
        if sy["status"] != "pair":
            res.append(rec)
            continue
        if not f.exists():
            rec["status"] = "unknown (SPARCL spectrum missing)"
            res.append(rec)
            continue
        d = load_pair(f)
        wa, wb = d["wave_a"].astype(float), d["wave_b"].astype(float)
        if len(wa) != len(wb) or np.max(np.abs(wa - wb)) > 1e-3:
            raise RuntimeError(f"{sy['name']}: wavelength grids differ")
        st = pair_statistic(
            wa,
            d["flux_a"].astype(float),
            d["ivar_a"].astype(float),
            d["flux_b"].astype(float),
            d["ivar_b"].astype(float),
            0.5 * (sy["z_a"] + sy["z_b"]),
        )
        rec.update(st)
        rec["status"] = (
            "unknown (no usable line)"
            if st["D"] is None
            else ("two different spectra" if st["D"] > T_D else "consistent with one source")
        )
        res.append(rec)
    return res


def wilson(k, n, zz=1.96):
    if n == 0:
        return None
    p = k / n
    den = 1 + zz * zz / n
    c = (p + zz * zz / (2 * n)) / den
    h = zz * math.sqrt(p * (1 - p) / n + zz * zz / (4 * n * n)) / den
    return [max(0.0, c - h), min(1.0, c + h)]


def gate_counts(pos_d, neg_d):
    """PLAN section 4: G1-G3 from the usable positive (lens) and negative (binary) D values."""
    kp, kn = sum(d <= T_D for d in pos_d), sum(d > T_D for d in neg_d)
    g1 = len(pos_d) >= GATE["min_controls"] and len(neg_d) >= GATE["min_controls"]
    g2 = bool(pos_d) and kp / len(pos_d) >= GATE["pos_frac"]
    g3 = bool(neg_d) and kn / len(neg_d) >= GATE["neg_frac"]
    return {
        "T": T_D,
        "rules": GATE,
        "usable_positive": len(pos_d),
        "usable_negative": len(neg_d),
        "positive_consistent": kp,
        "negative_different": kn,
        "positive_frac": kp / len(pos_d) if pos_d else None,
        "positive_wilson95": wilson(kp, len(pos_d)),
        "negative_frac": kn / len(neg_d) if neg_d else None,
        "negative_wilson95": wilson(kn, len(neg_d)),
        "G1": bool(g1),
        "G2": bool(g2),
        "G3": bool(g3),
        "passed": bool(g1 and g2 and g3),
    }


def cmd_validate(args):
    out = Path(args.out)
    res = analyse(out, {"control", "binary"})
    pos = [r for r in res if r["role"] == "control" and r.get("D") is not None]
    neg = [r for r in res if r["role"] == "binary" and r.get("D") is not None]
    gate = {
        "label": "derived (control spectra only; no rejected-role spectrum retrieved; see "
        "overlap_rejected_with_controls)",
        **gate_counts([r["D"] for r in pos], [r["D"] for r in neg]),
    }
    # review additions (reported, not gating): PLAN section 4's close-pair list, and rejected
    # pairs whose DESI TARGETIDs are also a control pair (their spectra were read as controls)
    gate["close_pairs_le_2arcsec"] = [
        {"role": r["role"], "name": r["name"], "sep_desi": r["sep_desi"], "D": r["D"]}
        for r in pos + neg
        if r["sep_desi"] <= 2.0
    ]
    P = json.load(open(out / "pairs.json"))["systems"]
    ctl = {
        frozenset((r["targetid_a"], r["targetid_b"])): r["name"]
        for r in P
        if r["role"] in ("control", "binary") and r["status"] == "pair"
    }
    gate["overlap_rejected_with_controls"] = [
        {"rejected": r["name"], "control": ctl[frozenset((r["targetid_a"], r["targetid_b"]))]}
        for r in P
        if r["role"] == "rejected"
        and r["status"] == "pair"
        and frozenset((r["targetid_a"], r["targetid_b"])) in ctl
    ]
    json.dump(gate, open(out / "gate.json", "w"), indent=1)
    (out / "gate.json").open("a").write("\n")
    json.dump(
        {"label": "derived; observed DESI spectra", "systems": res},
        open(out / "controls.json", "w"),
        indent=1,
        default=float,
    )
    (out / "controls.json").open("a").write("\n")
    print(json.dumps({k: v for k, v in gate.items() if k != "rules"}, indent=1))


def cmd_apply(args):
    out = Path(args.out)
    if not json.load(open(out / "gate.json"))["passed"]:
        raise SystemExit("gate not passed: targets are not read (PLAN section 4)")
    res = analyse(out, {"rejected"})
    json.dump(
        {
            "label": "derived; observed DESI spectra; classes are spectral only, "
            "never a dark-lens claim",
            "T": T_D,
            "systems": res,
        },
        open(out / "targets.json", "w"),
        indent=1,
        default=float,
    )
    (out / "targets.json").open("a").write("\n")
    for r in res:
        print(r["name"], r["status"], r.get("D"), r.get("dv_kms"), r.get("sep_desi"))


def cmd_sheet(args):
    """Contact sheet of the control pairs: the counted line with the best weaker-side S/N."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out)
    res = [r for r in json.load(open(out / "controls.json"))["systems"] if r.get("D") is not None]
    fig, axs = plt.subplots(4, 4, figsize=(14, 11))
    for ax, r in zip(axs.flat, res, strict=False):
        used = {k: v for k, v in r["lines"].items() if v.get("used")}
        line = max(used, key=lambda k: min(used[k]["snr_a"], used[k]["snr_b"]))
        d = load_pair(out / "spectra" / f"{r['role']}_{r['name'].replace('/', '_')}.csv.gz")
        z = 0.5 * (r["z_a"] + r["z_b"])
        a = line_profile(
            *(d[f"{x}_a"].astype(float) for x in ("wave", "flux", "ivar")), z, LINES[line]
        )
        b = line_profile(
            *(d[f"{x}_b"].astype(float) for x in ("wave", "flux", "ivar")), z, LINES[line]
        )
        ax.plot(a[0], a[1], lw=0.5, label="A")
        ax.plot(b[0], used[line]["s"] * b[1], lw=0.5, label=f"B x {used[line]['s']:.2f}")
        ax.set_title(
            f"{r['role']} {r['name'][:22]}\n{line} chi2nu {used[line]['chi2_nu']:.2f}; "
            f'D {r["D"]:.2f}; sep {r["sep_desi"]:.1f}"',
            fontsize=7,
        )
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=6)
    for ax in axs.flat[len(res) :]:
        ax.axis("off")
    fig.suptitle("W2 DESI pairs: control line profiles (observed; v in km/s)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "controls_sheet.png", dpi=70)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=["pairs", "fetch", "validate", "apply", "sheet"])
    ap.add_argument("--out", default="results/w2_desi_pairs")
    ap.add_argument("--cache", default="data/cache/w2_desi_pairs")
    ap.add_argument("--roles", default="control,binary")
    a = ap.parse_args(argv)
    stages = {"pairs": cmd_pairs, "fetch": cmd_fetch, "validate": cmd_validate}
    stages.update(apply=cmd_apply, sheet=cmd_sheet)
    stages[a.stage](a)


if __name__ == "__main__":
    main()

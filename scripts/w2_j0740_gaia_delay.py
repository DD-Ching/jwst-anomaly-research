"""W2 J0740+2926: can Gaia DR3 G epoch photometry identify a delayed replay of one source?

Pre-registered in ``results/w2_j0740_gaia_delay/PREREG.md`` (frozen at commit 4ef8f8c, before any
light curve of the pair was fetched). Estimator: PyCS3 3.0.4 spline technique (Millon et al. 2020,
JOSS 5, 2654), regdiff as a cross-check on the target. Controls: simulated same-source and
independent pairs on the actual Gaia epochs and errors (damped random walk, common per-transit
errors, A->B crosstalk), and a real pair of unrelated Gaia DR3 AGN.

    uv pip install -e ".[delays]"                  # PyCS3 3.0.4
    python scripts/w2_j0740_gaia_delay.py fetch   # TAP + DataLink, cached under $JWST_ANOMALY_DATA
    python scripts/w2_j0740_gaia_delay.py run     # controls, then the target; writes result.json
"""

from __future__ import annotations

import argparse
import io
import json
import logging
import re
import sys
import urllib.parse
import urllib.request
import warnings
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from jwst_anomaly import gaia_mulens as gm
from jwst_anomaly import paths

A_ID, B_ID = "878723302007625344", "878723302009111040"
RA0, DEC0 = 115.05598, 29.44640
# Frozen settings (PREREG.md); all are ASSUMPTIONs of the test design.
KNOTSTEP = 100.0  # d, PyCS3 spline knot step
DT_MAX = 400.0  # d, search range |dt| <= DT_MAX
ZERO_LAG = 5.0  # d, |dt| below this is not interpreted
STARTS = np.arange(-400.0, 400.1, 50.0)  # 17 initial delays
TOL = 20.0  # d, recovery tolerance
GROUP_LEN = 1.0  # d, transit groups used for the error scale
DRW_TAU, DRW_SIG = 300.0, 0.10  # d, mag
COMMON_SIG, XTALK, SAME_TRANSIT = 0.01, 0.01, 0.01  # mag, flux fraction A->B, d
C1_DELAYS = (-200.0, -100.0, 50.0, 100.0, 200.0, 300.0)
N_C1, N_C2, N_POWER = 40, 100, 40
GRID = np.arange(-400.0, 400.1, 5.0)
NIT_START, NIT_FINE, NIT_FIXED = 3, 10, 5  # PyCS3 opt_fine iterations (implementation detail)
# Deviation (implementation, after the first run aborted): PyCS3's opt_fine default bokeps = 10
# raised 'Knot spacing min < epsilon' on these sparse curves; use the spline fit() default 2.0,
# skip a failed start and record a failed fixed-delay fit as not evaluated (NaN).
BOKEPS = 2.0
SEED = 740
CONTROL_QUERY = (
    "SELECT TOP 2 g.source_id, g.ra, g.dec, g.phot_g_mean_mag, "
    f"DISTANCE(POINT(g.ra, g.dec), POINT({RA0}, {DEC0})) AS dist "
    "FROM gaiadr3.gaia_source AS g JOIN gaiadr3.vari_agn AS v ON v.source_id = g.source_id "
    f"WHERE 1 = CONTAINS(POINT(g.ra, g.dec), CIRCLE({RA0}, {DEC0}, 1.0)) "
    "AND g.phot_g_mean_mag BETWEEN 18.0 AND 19.8 "
    f"AND g.source_id NOT IN ({A_ID}, {B_ID}) ORDER BY dist"
)


def out_dir() -> Path:
    d = paths.repo_root() / "results" / "w2_j0740_gaia_delay"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cache_dir() -> Path:
    d = paths.data_root() / "raw" / "w2_j0740_gaia_delay"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ----------------------------------------------------------------------------- data


def fetch() -> dict:
    ctrl_path = cache_dir() / "control_pair.csv"
    if not ctrl_path.exists():
        ctrl_path.write_text(gm._tap_csv(CONTROL_QUERY))
    ctrl = [line.split(",") for line in ctrl_path.read_text().strip().splitlines()[1:]]
    ids = [A_ID, B_ID] + [r[0] for r in ctrl]
    todo = [i for i in ids if not (cache_dir() / f"{i}.csv").exists()]
    if todo:
        body = urllib.parse.urlencode(
            {
                "RETRIEVAL_TYPE": "EPOCH_PHOTOMETRY",
                "DATA_STRUCTURE": "DATAMODEL_STANDARD",
                "FORMAT": "CSV",
                "VALID_DATA": "false",
                "ID": ",".join(todo),
            }
        ).encode()
        with urllib.request.urlopen(gm.DATALINK, data=body, timeout=600) as resp:  # noqa: S310
            raw = resp.read()
        if zipfile.is_zipfile(io.BytesIO(raw)):
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                for name in z.namelist():
                    m = re.search(r"(\d{6,})\.csv$", name)
                    if m:
                        (cache_dir() / f"{m.group(1)}.csv").write_bytes(z.read(name))
        elif len(todo) == 1:
            (cache_dir() / f"{todo[0]}.csv").write_bytes(raw)
    return {"ids": ids, "control_rows": ctrl, "query": CONTROL_QUERY}


def load(sid: str) -> dict:
    """Q1 (DR3 rejection flags, as ``gaia_mulens.parse_epoch_csv``) and Q3 (unbinned)."""
    text = (cache_dir() / f"{sid}.csv").read_text()
    lc = gm.parse_epoch_csv(text)
    t = np.asarray(lc["time"], float) - gm.GAIA_T0
    return {
        "t": t,
        "m": np.asarray(lc["mag"], float),
        "e": np.asarray(lc["mag_err"], float),
        "n_rejected": int(lc.meta["n_rejected"]),
    }


def groups(t: np.ndarray) -> list[np.ndarray]:
    """Indices of transit groups no longer than GROUP_LEN (t sorted)."""
    out, start = [], 0
    for i in range(1, len(t) + 1):
        if i == len(t) or t[i] - t[start] > GROUP_LEN:
            out.append(np.arange(start, i))
            start = i
    return out


def error_scale(curves: list[dict]) -> dict:
    """Q2: s = max(1, sqrt(pooled chi2_nu about group means)), groups <= 1 d, n >= 2."""
    chi2, dof, ngr = 0.0, 0, 0
    for c in curves:
        for g in groups(c["t"]):
            if len(g) < 2:
                continue
            w = 1.0 / c["e"][g] ** 2
            mu = np.sum(w * c["m"][g]) / np.sum(w)
            chi2 += float(np.sum(((c["m"][g] - mu) / c["e"][g]) ** 2))
            dof += len(g) - 1
            ngr += 1
    chi2nu = chi2 / dof if dof else float("nan")
    return {"chi2_nu": chi2nu, "dof": dof, "groups": ngr, "s": max(1.0, float(np.sqrt(chi2nu)))}


# ----------------------------------------------------------------------------- PyCS3


def _quiet():
    warnings.filterwarnings("ignore")
    logging.disable(logging.CRITICAL)


def _lcs(a: dict, b: dict, dt: float):
    from pycs3.gen.lc_func import factory
    from pycs3.gen.polyml import addtolc

    la = factory(a["t"], a["m"], a["e"], telescopename="Gaia", object="A")
    lb = factory(b["t"], b["m"] - (np.median(b["m"]) - np.median(a["m"])), b["e"], object="B")
    lb.telescopename = "Gaia"
    addtolc(lb, nparams=2, autoseasonsgap=1e4)
    lb.timeshift = -dt  # B lags A by dt  <=>  B is moved dt earlier onto A
    return [la, lb]


def _dt(lcs) -> float:
    return float(lcs[0].timeshift - lcs[1].timeshift)


def r2_fixed(a: dict, b: dict, dt: float, nit: int = NIT_FIXED) -> float:
    from pycs3.spl.topopt import opt_fine

    lcs = _lcs(a, b, dt)
    try:
        sp = opt_fine(
            lcs, knotstep=KNOTSTEP, nit=nit, shifttime=False, bokeps=BOKEPS, verbose=False
        )
    except RuntimeError:
        return float("nan")
    return float(sp.lastr2nostab)


def r2_single(c: dict) -> float:
    from pycs3.gen.lc_func import factory
    from pycs3.spl.topopt import opt_fine

    lc = factory(c["t"], c["m"], c["e"], telescopename="Gaia", object="X")
    try:
        sp = opt_fine(
            [lc], knotstep=KNOTSTEP, nit=NIT_FIXED, shifttime=False, bokeps=BOKEPS, verbose=False
        )
    except RuntimeError:
        return float("nan")
    return float(sp.lastr2nostab)


def estimate(a: dict, b: dict) -> dict:
    """Spline technique from the 17 starts; the best start is refined; result kept if in range."""
    from pycs3.spl.topopt import opt_fine

    best = None
    for s in STARTS:
        lcs = _lcs(a, b, float(s))
        try:
            sp = opt_fine(
                lcs, knotstep=KNOTSTEP, nit=NIT_START, shifttime=True, bokeps=BOKEPS, verbose=False
            )
        except RuntimeError:
            continue
        d, r2 = _dt(lcs), float(sp.lastr2nostab)
        if abs(d) <= DT_MAX and (best is None or r2 < best[1]):
            best = (d, r2)
    if best is None:
        return {"dt_hat": None, "r2_joint": None}
    lcs = _lcs(a, b, best[0])
    try:
        sp = opt_fine(
            lcs, knotstep=KNOTSTEP, nit=NIT_FINE, shifttime=True, bokeps=BOKEPS, verbose=False
        )
        d, r2 = _dt(lcs), float(sp.lastr2nostab)
    except RuntimeError:
        d, r2 = best
    if abs(d) > DT_MAX or r2 > best[1]:
        d, r2 = best
    return {"dt_hat": d, "r2_joint": r2}


def regdiff(a: dict, b: dict) -> dict:
    from pycs3.regdiff.multiopt import opt_ts

    best = None
    for s in STARTS:
        lcs = _lcs(a, b, float(s))
        try:
            _, wtv = opt_ts(lcs, verbose=False)
        except RuntimeError:
            continue
        d = _dt(lcs)
        if abs(d) <= DT_MAX and (best is None or wtv < best[1]):
            best = (d, float(wtv))
    return {"dt_hat": best[0] if best else None, "minwtv": best[1] if best else None}


def analyse(a: dict, b: dict, dt_true: float | None = None, scan=None) -> dict:
    _quiet()
    sep = r2_single(a) + r2_single(b)
    out = estimate(a, b)
    out["r2_sep"] = sep
    d_hat = None if out["r2_joint"] is None else out["r2_joint"] - sep
    out["D_hat"] = None if d_hat is None or not np.isfinite(d_hat) else d_hat
    if dt_true is not None:
        out["D_true"] = r2_fixed(a, b, dt_true) - sep
    if scan is not None:
        out["scan"] = [r2_fixed(a, b, float(x)) - sep for x in scan]
    return out


# ----------------------------------------------------------------------------- simulations


def drw(t: np.ndarray, rng) -> np.ndarray:
    o = np.argsort(t)
    ts = t[o]
    x = np.empty(len(t))
    x[0] = rng.normal(0.0, DRW_SIG)
    for i in range(1, len(t)):
        q = np.exp(-(ts[i] - ts[i - 1]) / DRW_TAU)
        x[i] = q * x[i - 1] + rng.normal(0.0, DRW_SIG * np.sqrt(1.0 - q * q))
    out = np.empty(len(t))
    out[o] = x
    return out


def simulate(a: dict, b: dict, dt: float | None, rng) -> tuple[dict, dict]:
    """Same source (dt given) or independent (dt None) on the actual epochs and errors."""
    na = len(a["t"])
    if dt is None:
        sa, sb = drw(a["t"], rng), drw(b["t"], rng)
        sa_at_b = None
    else:
        s = drw(np.concatenate([a["t"], b["t"] - dt, b["t"]]), rng)
        sa, sb, sa_at_b = s[:na], s[na : na + len(b["t"])], s[na + len(b["t"]) :]
    if sa_at_b is None:  # A's intrinsic curve at B's epochs, for the crosstalk term
        sa_at_b = np.interp(b["t"], a["t"], sa)
    ma = 18.24 + sa + rng.normal(0.0, a["e"])
    mb_int = 19.57 + sb
    fa = 10 ** (-0.4 * (18.24 + sa_at_b))
    mb = -2.5 * np.log10(10 ** (-0.4 * mb_int) + XTALK * fa) + rng.normal(0.0, b["e"])
    for i, tb in enumerate(b["t"]):  # common per-transit term
        j = np.flatnonzero(np.abs(a["t"] - tb) < SAME_TRANSIT)
        if len(j):
            c = rng.normal(0.0, COMMON_SIG)
            ma[j[0]] += c
            mb[i] += c
    return {"t": a["t"], "m": ma, "e": a["e"]}, {"t": b["t"], "m": mb, "e": b["e"]}


def _job(args):
    kind, dt, seed, a, b, scan = args
    rng = np.random.default_rng(seed)
    sa, sb = simulate(a, b, dt, rng)
    res = analyse(sa, sb, dt_true=dt, scan=scan)
    res.update(kind=kind, dt_true=dt, seed=seed)
    return res


# ----------------------------------------------------------------------------- run


def intervals(mask: np.ndarray, grid: np.ndarray) -> list[list[float]]:
    out, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        if start is not None and (not m or i == len(mask) - 1):
            end = i if m else i - 1
            out.append([float(grid[start]), float(grid[end])])
            start = None
    return out


def run(procs: int = 4) -> dict:
    meta = fetch()
    a, b = load(A_ID), load(B_ID)
    c_ids = [r[0] for r in meta["control_rows"]]
    ca, cb = load(c_ids[0]), load(c_ids[1])
    for c in (ca, cb):
        o = np.argsort(c["t"])
        c.update(t=c["t"][o], m=c["m"][o], e=c["e"][o])
    o = np.argsort(a["t"]), np.argsort(b["t"])
    a.update(t=a["t"][o[0]], m=a["m"][o[0]], e=a["e"][o[0]])
    b.update(t=b["t"][o[1]], m=b["m"][o[1]], e=b["e"][o[1]])
    q2 = error_scale([a, b])
    q2c = error_scale([ca, cb])
    for c in (a, b):
        c["e"] = c["e"] * q2["s"]
    for c in (ca, cb):
        c["e"] = c["e"] * q2c["s"]
    sampling = {
        "n_A": len(a["t"]),
        "n_B": len(b["t"]),
        "rejected_A": a["n_rejected"],
        "rejected_B": b["n_rejected"],
        "groups_A": len(groups(a["t"])),
        "groups_B": len(groups(b["t"])),
        "span_d": float(max(a["t"].max(), b["t"].max()) - min(a["t"].min(), b["t"].min())),
        "median_err_A": float(np.median(a["e"])),
        "median_err_B": float(np.median(b["e"])),
        "same_transit_pairs": int(sum(np.any(np.abs(a["t"] - x) < SAME_TRANSIT) for x in b["t"])),
    }
    # Controls first (C1, C2 with a power scan on the first N_POWER, C3), then the target.
    jobs = []
    k = 0
    for dt in C1_DELAYS:
        for _ in range(N_C1):
            jobs.append(("C1", dt, SEED + k, a, b, None))
            k += 1
    for i in range(N_C2):
        jobs.append(("C2", None, SEED + k, a, b, GRID if i < N_POWER else None))
        k += 1
    with ProcessPoolExecutor(procs) as ex:
        sims = list(ex.map(_job, jobs, chunksize=4))
    c1 = [s for s in sims if s["kind"] == "C1"]
    c2 = [s for s in sims if s["kind"] == "C2"]
    d_all = np.array([s["D_true"] for s in c1])
    d_true = d_all[np.isfinite(d_all)]
    t95, t99 = float(np.percentile(d_true, 95)), float(np.percentile(d_true, 99))
    rec = {}
    for dt in C1_DELAYS:
        hats = np.array(
            [np.nan if s["dt_hat"] is None else s["dt_hat"] for s in c1 if s["dt_true"] == dt]
        )
        ok = np.abs(hats - dt) <= TOL
        rec[str(dt)] = {
            "fraction": float(np.mean(ok)),
            "median_dt_hat": float(np.nanmedian(hats)),
            "pass": bool(np.mean(ok) >= 0.68),
        }

    def flagged(s):
        return (
            s["D_hat"] is not None and s["D_hat"] <= t95 and ZERO_LAG <= abs(s["dt_hat"]) <= DT_MAX
        )

    fpr = float(np.mean([flagged(s) for s in c2]))
    power = [float(np.mean(np.array(s["scan"]) > t99)) for s in c2 if "scan" in s]
    c3 = analyse(ca, cb)
    target = analyse(a, b, scan=GRID)
    rd = regdiff(a, b)
    excl = np.array(target["scan"]) > t99
    lens = (GRID >= 20.0) & (GRID <= DT_MAX)
    c1_pass = any(v["pass"] for v in rec.values())
    c2_pass = fpr <= 0.05
    c3_pass = c3["D_hat"] is None or c3["D_hat"] > t95
    supports = (
        c1_pass
        and c2_pass
        and c3_pass
        and target["D_hat"] is not None
        and target["D_hat"] <= t95
        and ZERO_LAG <= abs(target["dt_hat"]) <= DT_MAX
        and rd["dt_hat"] is not None
        and abs(rd["dt_hat"] - target["dt_hat"]) <= TOL
    )
    if supports:
        verdict = "supports a same-source delay (lens candidate only)"
    elif excl.any():
        verdict = "excludes the listed delay intervals"
    else:
        verdict = "data insufficient"
    res = {
        "label": "derived from observed Gaia DR3 epoch photometry; controls simulated",
        "prereg": "results/w2_j0740_gaia_delay/PREREG.md (commit 4ef8f8c)",
        "pycs3": __import__("pycs3").__version__
        if hasattr(__import__("pycs3"), "__version__")
        else "3.0.4",
        "settings": {
            "knotstep": KNOTSTEP,
            "dt_max": DT_MAX,
            "zero_lag": ZERO_LAG,
            "tol": TOL,
            "drw_tau": DRW_TAU,
            "drw_sig": DRW_SIG,
            "common_sig": COMMON_SIG,
            "xtalk": XTALK,
            "n_c1": N_C1,
            "n_c2": N_C2,
            "n_power": N_POWER,
            "nit": [NIT_START, NIT_FINE, NIT_FIXED],
            "seed": SEED,
        },
        "sampling": sampling,
        "error_scale": q2,
        "control_pair": {"rows": meta["control_rows"], "query": CONTROL_QUERY, "error_scale": q2c},
        "thresholds": {"T95": t95, "T99": t99},
        "C1": rec,
        "C2": {
            "fpr": fpr,
            "pass": c2_pass,
            "dt_hat": [s["dt_hat"] for s in c2],
            "D_hat": [s["D_hat"] for s in c2],
            "power_fraction_excluded": {
                "median": float(np.median(power)),
                "min": float(np.min(power)),
                "max": float(np.max(power)),
            },
        },
        "C3": {**c3, "pass": c3_pass},
        "C1_D_true": d_true.tolist(),
        "C1_D_true_not_evaluated": int(np.sum(~np.isfinite(d_all))),
        "C2_D_hat_not_evaluated": int(sum(s["D_hat"] is None for s in c2)),
        "target": {
            "dt_hat": target["dt_hat"],
            "D_hat": target["D_hat"],
            "r2_joint": target["r2_joint"],
            "r2_sep": target["r2_sep"],
            "regdiff": rd,
            "scan_grid": GRID.tolist(),
            "scan_D": target["scan"],
            "excluded_intervals": intervals(excl, GRID),
            "R_lens_excluded": bool(excl[lens].all()),
            "R_all_excluded": bool(excl.all()),
            "fraction_R_lens_excluded": float(excl[lens].mean()),
        },
        "controls_pass": {"C1_any": c1_pass, "C2": c2_pass, "C3": c3_pass},
        "verdict": verdict,
    }
    curves = {
        k: {"t": v["t"].tolist(), "m": v["m"].tolist(), "e": v["e"].tolist()}
        for k, v in (("A", a), ("B", b), ("CA", ca), ("CB", cb))
    }
    (out_dir() / "result.json").write_text(json.dumps(res, indent=1) + "\n")
    (cache_dir() / "curves.json").write_text(json.dumps(curves))
    return res


def post_hoc() -> dict:
    """Post hoc (after the run, labelled): applies PREREG's zero-lag rule to the exclusion, counts
    the C2 delay aliases and rescales the target's D for the error-scale uncertainty (Q2 used
    14 dof). Uniform error rescaling by s'/s multiplies r2, hence D, by (s/s')^2; T99 is kept."""
    from scipy.stats import chi2

    r = json.loads((out_dir() / "result.json").read_text())
    g = np.array(r["target"]["scan_grid"])
    d = np.array(r["target"]["scan_D"])
    t99 = r["thresholds"]["T99"]
    lens = (g >= 20.0) & (g <= DT_MAX)
    interp = np.abs(g) >= ZERO_LAG
    q = r["error_scale"]
    lo = q["chi2_nu"] * q["dof"] / chi2.ppf(0.975, q["dof"])
    hi = q["chi2_nu"] * q["dof"] / chi2.ppf(0.025, q["dof"])
    rows = {}
    for name, s_new in (("adopted", q["s"]), ("s_95_low", np.sqrt(lo)), ("s_95_high", np.sqrt(hi))):
        ex = (d * (q["s"] / s_new) ** 2 > t99) & interp
        rows[name] = {
            "s": float(s_new),
            "excluded_intervals": intervals(ex, g),
            "fraction_grid_excluded": float(ex[interp].mean()),
            "fraction_R_lens_excluded": float(ex[lens].mean()),
        }
    h = np.array([x for x in r["C2"]["dt_hat"] if x is not None])
    res = {
        "label": "post hoc (computed after the frozen run; derived from result.json)",
        "zero_lag_rule": "PREREG section 2: |dt| < 5 d is not interpreted; the run's exclusion "
        "omitted this mask, so it is applied here",
        "error_scale_sensitivity": rows,
        "C2_dt_hat_counts": {
            "abs_ge_200": int(np.sum(np.abs(h) >= 200)),
            "abs_lt_100": int(np.sum(np.abs(h) < 100)),
            "n": int(len(h)),
        },
        "C2_D_hat_max": float(np.max([x for x in r["C2"]["D_hat"] if x is not None])),
        "C1_D_true_median": float(np.median(r["C1_D_true"])),
    }
    cur = json.loads((cache_dir() / "curves.json").read_text())
    res["last_season_minus_earlier_median_mag"] = {
        k: float(
            np.median(np.array(cur[k]["m"])[np.array(cur[k]["t"]) > 2600.0])
            - np.median(np.array(cur[k]["m"])[np.array(cur[k]["t"]) <= 2600.0])
        )
        for k in ("A", "B")
    }
    res["last_season_definition"] = "Gaia time > 2600 d (the last group of visits, 2638-2690)"
    (out_dir() / "post_hoc.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def sis() -> dict:
    """[der] SIS delay for the observed configuration (model prediction, not a measurement).

    Image radii from the separation and the Gaia G flux ratio (ASSUMPTION: the G ratio is the
    macro-magnification ratio); z_s = mean DESI redshift; Planck18; lens redshift unknown, so a
    grid.
    The brighter image is the minimum and leads, so dt > 0 (B lags) in this script's convention.
    """
    import astropy.units as u
    from astropy.constants import c
    from astropy.cosmology import Planck18 as cosmo

    sep, dg, zs = 2.64, 19.574408 - 18.240225, (0.9800 + 0.9779) / 2
    f = 10 ** (0.4 * dg)
    tp, tm = sep * f / (1 + f), sep / (1 + f)
    dth2 = ((tp * u.arcsec) ** 2 - (tm * u.arcsec) ** 2).to(u.rad**2).value
    rows = {}
    for zl in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7):
        dd = (
            cosmo.angular_diameter_distance(zl)
            * cosmo.angular_diameter_distance(zs)
            / cosmo.angular_diameter_distance_z1z2(zl, zs)
        )
        rows[str(zl)] = float(((1 + zl) / (2 * c) * dd * dth2).to(u.day).value)
    res = {
        "label": "model prediction (SIS), derived; not a measurement",
        "inputs": {"sep_arcsec": sep, "dG_mag": dg, "flux_ratio": f, "z_s": zs},
        "theta_plus_minus_arcsec": [tp, tm],
        "dt_days_by_z_lens": rows,
    }
    (out_dir() / "sis.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def plot() -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res = json.loads((out_dir() / "result.json").read_text())
    cur = json.loads((cache_dir() / "curves.json").read_text())
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    for k, col, dy in (("A", "C0", 0.0), ("B", "C3", None), ("CA", "C2", 0.6), ("CB", "C4", None)):
        c = cur[k]
        m = np.array(c["m"])
        if dy is None:  # align the partner's median on its reference curve
            ref = np.array(cur["A" if k == "B" else "CA"]["m"])
            dy = (0.0 if k == "B" else 0.6) + np.median(ref) - np.median(m)
        else:
            dy = dy + 0.0 * np.median(m)
        ax[0, 0].errorbar(c["t"], m + dy, c["e"], fmt="o", ms=3, color=col, label=k)
    ax[0, 0].invert_yaxis()
    ax[0, 0].set(xlabel="Gaia time (BJD(TCB) - 2455197.5, d)", ylabel="G (B and CB moved to A, CA)")
    ax[0, 0].legend(fontsize=8, ncol=4)
    ax[0, 0].set_title("J0740 A/B and control pair CA/CB (+0.6 mag)", fontsize=9)
    g = np.array(res["target"]["scan_grid"])
    ax[0, 1].plot(g, res["target"]["scan_D"], "k-", label="J0740 D(dt)")
    for name, v, ls in (
        ("T95", res["thresholds"]["T95"], "--"),
        ("T99", res["thresholds"]["T99"], ":"),
    ):
        ax[0, 1].axhline(v, color="r", ls=ls, label=name)
    if res["target"]["dt_hat"] is not None:
        ax[0, 1].axvline(res["target"]["dt_hat"], color="C0", lw=1, label="spline dt_hat")
    ax[0, 1].set(xlabel="dt (d; >0: B lags A)", ylabel="D = r2_joint - r2_A - r2_B")
    ax[0, 1].legend(fontsize=8)
    d1 = np.array(res["C1_D_true"])
    ax[1, 0].hist(d1, bins=30, alpha=0.6, label="C1 same source, D(dt_true)")
    d2 = np.array([x for x in res["C2"]["D_hat"] if x is not None])
    ax[1, 0].hist(d2, bins=30, alpha=0.6, label="C2 independent, D(dt_hat)")
    if res["target"]["D_hat"] is not None:
        ax[1, 0].axvline(res["target"]["D_hat"], color="k", label="J0740 D(dt_hat)")
    if res["C3"]["D_hat"] is not None:
        ax[1, 0].axvline(res["C3"]["D_hat"], color="C2", ls="--", label="C3 real pair")
    ax[1, 0].axvline(res["thresholds"]["T95"], color="r", ls="--")
    ax[1, 0].set(xlabel="D", ylabel="count")
    ax[1, 0].legend(fontsize=8)
    h2 = [x for x in res["C2"]["dt_hat"] if x is not None]
    ax[1, 1].hist(h2, bins=np.arange(-400, 401, 20), alpha=0.6, label="C2 dt_hat (independent)")
    for k, v in res["C1"].items():
        ax[1, 1].plot(float(k), v["fraction"] * 20, "r^")
    ax[1, 1].set(xlabel="dt (d)", ylabel="count;  red: C1 recovery fraction x 20")
    ax[1, 1].legend(fontsize=8)
    fig.suptitle(f"W2 J0740+2926 Gaia DR3 G delay test: {res['verdict']}", fontsize=10)
    fig.tight_layout()
    path = out_dir() / "contact_sheet.png"
    fig.savefig(path, dpi=80)
    return path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("cmd", choices=["fetch", "run", "sis", "plot", "post_hoc"])
    p.add_argument("--procs", type=int, default=4)
    args = p.parse_args(argv)
    if args.cmd == "fetch":
        print(json.dumps(fetch(), indent=1))
    elif args.cmd == "sis":
        print(json.dumps(sis(), indent=1))
    elif args.cmd == "plot":
        print(plot())
    elif args.cmd == "post_hoc":
        print(json.dumps(post_hoc(), indent=1))
    else:
        res = run(args.procs)
        print(
            json.dumps(
                {k: res[k] for k in ("thresholds", "C1", "controls_pass", "verdict")}, indent=1
            )
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""W3 TESS 120 s product-level injection trial (owner authorization 2026-10-26, after #561).

Plan: ``results/w3_kepler/tess120_inject/PLAN.md``. Injections of the n = 1 negative-mass point-lens
model (``exotic_sim``) into published SPOC 120 s PDCSAP light curves, each paired with its own
noise-free ideal (the #561 shape criterion on a gap-free cadence grid). No real-event search.

Stages (``python scripts/w3_tess120_inject.py <stage>``):

- ``select``: the star-sector list (rate-weighted draw from the #561 TIC sample) and its download
  order (``sample.json``); metadata only;
- ``run``: fetch each light curve (ledger), record its sha256, run the background and the injections
  for every scenario, delete the file (``per_ss/*.json``, ``files.json``);
- ``summarize``: paired counts per scenario and t_E bin, losses per step, bootstrap intervals, the
  frozen verdict, and the effect on #561's reach (``summary.json``).

Labels: injected events and their recoveries are simulations; the light curves are observed data
used only as a background; the efficiency is product level (after SPOC processing), not sky to
result.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import net_ledger as nl  # noqa: E402
import w3_shortcad_reach as ws  # noqa: E402

from jwst_anomaly import exotic_sim as es  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w3_kepler" / "tess120_inject"
DATA = Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data")) / "tess120_inject" / "raw"
CAP = 1_000_000_000  # owner limit for this round, bytes
DL_URL = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/{obs_id}_lc.fits"

SEED = 564
N_STAR_SECTORS = 200
N_INJ = 30  # injections per draw and scenario
SCENARIOS = {  # (|M| in M_sun, v in km/s): #561 TESS 120 s best mass and one lower mass per speed
    "v220_M5.6e-6": (5.623413251903491e-06, 220.0),
    "v220_M1e-7": (1e-07, 220.0),
    "v30_M3.2e-6": (3.162277660168379e-06, 30.0),
    "v30_M1e-8": (1e-08, 30.0),
}
U0_MAX = 2.65  # the #561 shape grid's u0 range
T_HALF = ws.T_HALF  # model window, +- t_E
BASE_WINDOW_D = (
    3.0  # running-median baseline window; the longest frozen-scenario dark phase is <= 1.2 d
)
NSIG = 5.0  # D and S also need >= NSIG robust sigma (the Kepler LOCAL_NSIG)
CONTIG = 1.5  # neighbouring cadences are contiguous if their spacing is <= CONTIG cadences
SHOULDER_MIN_CAD = 4  # shoulder window per side: max(dark-run duration, 4 cadences)
TE_BINS = [0.0, 2.5, 5.0, 15.0, 50.0, np.inf]  # t_E in cadences
N_BOOT = 2000
SPAN_CAP_D = 40.0  # a single-sector record longer than this is not used (none expected)


# -------------------------------------------------------------------------------------- selection


def weighted_draws(weights: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """n draws with replacement, probability proportional to weight: each draw then carries equal
    weight, as in the reach integral (a draw without replacement would cap the heaviest records)."""
    w = np.asarray(weights, float)
    return rng.choice(len(w), size=n, replace=True, p=w / w.sum())


def stage_select(_a) -> None:
    lg = nl.Ledger(OUT / "net_ledger.jsonl", CAP)
    base = (
        "FROM dbo.ObsPointing WHERE obs_collection='TESS' AND dataproduct_type='timeseries'"
        " AND provenance_name='SPOC' AND t_exptime=120 AND obs_id NOT LIKE '%-s____-s____-%'"
        f" AND target_name LIKE '%{ws.SAMPLE_SUFFIX['tess_120']}'"
    )
    recs = ws.tap(
        ws.CAOM, "SELECT obs_id, target_name, sequence_number, t_min, t_max " + base, lg, 8_000_000
    )
    stars = {int(r["id"]): r for r in ws._read_csv("tess_120.csv.gz")}
    rows = []
    for r in recs:
        tic, sector = int(r["target_name"]), int(r["sequence_number"])
        s = stars.get(tic)
        span = float(r["t_max"]) - float(r["t_min"])
        if s is None or not s["radius"] or not s["dist"] or sector == 1751 or span > SPAN_CAP_D:
            continue
        rows.append(
            {
                "obs_id": r["obs_id"],
                "tic": tic,
                "sector": sector,
                "radius": float(s["radius"]),
                "dist": float(s["dist"]),
                "span_d": round(span, 5),
            }
        )
    rows.sort(key=lambda r: r["obs_id"])
    w = np.array([r["span_d"] * r["dist"] ** 1.5 for r in rows])
    draws = weighted_draws(w, N_STAR_SECTORS, np.random.default_rng(SEED))
    pick: dict[int, dict] = {}  # download order: first draw of each record
    for k, i in enumerate(draws):
        pick.setdefault(int(i), rows[int(i)] | {"weight": float(f"{w[i]:.6g}"), "draws": []})
        pick[int(i)]["draws"].append(k)
    pick = list(pick.values())
    srt = np.sort(w)[::-1]
    cum = np.cumsum(srt) / srt.sum()
    out = {
        "label": "metadata (MAST CAOM) and #561 TIC values; weighted draw frozen before any flux",
        "plan": "results/w3_kepler/tess120_inject/PLAN.md",
        "pool_records": len(rows),
        "pool_weight_sum": float(f"{w.sum():.6g}"),
        "weight": "span_d * dist_pc**1.5 (the Level-0 crossing rate per star-sector, any |M|)",
        "seed": SEED,
        "draws": N_STAR_SECTORS,
        "unique_star_sectors": len(pick),
        "share_of_weight_in_top_frac": {
            f: float(f"{cum[max(0, int(round(f * len(srt))) - 1)]:.4g}")
            for f in (0.01, 0.05, 0.1, 0.25, 0.5)
        },
        "star_sectors": pick,
    }
    (OUT / "sample.json").write_text(json.dumps(out, indent=1) + "\n")
    print(len(rows), "records in the pool;", len(pick), "drawn; bytes", lg.used())


# -------------------------------------------------------------------------------------- the model


def draw_events(rng, n, m_msun, v_kms, r_sun, d_pc, t_lo, t_hi) -> list[dict]:
    """n umbra-crossing events (u0 < U0_MAX) of one population at one star: x from the crossing
    rate (proportional to R_E(x), i.e. sqrt(x(1 - x))), u0 uniform, t0 uniform over the sector
    span (gaps included). t_E and rho follow from |M|, x, v, R_* and D_S. Identifiability is
    not used in the draw, so it is counted once, in the paired ideal."""
    xg = ws.X_GRID
    pdf = np.sqrt(xg * (1 - xg))
    cdf = np.r_[0, np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(xg))]
    cdf /= cdf[-1]
    out = []
    for _ in range(n):
        x = float(np.interp(rng.random(), cdf, xg))
        r_e_pc = math.sqrt(ws.SCHW4_PC * m_msun * x * (1 - x) * d_pc)
        t_e_d = r_e_pc * ws.PC_KM / v_kms / ws.DAY_S
        rho = x * r_sun * ws.R_SUN_KM / (r_e_pc * ws.PC_KM)
        out.append(
            {
                "x": x,
                "t_E_d": t_e_d,
                "rho": rho,
                "u0": float(rng.uniform(0, U0_MAX)),
                "t0": float(rng.uniform(t_lo, t_hi)),
            }
        )
    return out


def amp_table(rho: float) -> tuple[np.ndarray, np.ndarray]:
    """A(beta) of the uniform-disk source (n = 1, sign -1), dense near the caustic edges."""
    beta = np.unique(
        np.r_[
            np.linspace(0, T_HALF + 0.5, 2201),
            2 + rho * np.linspace(-1.2, 1.2, 241),
            rho * np.linspace(0, 2, 41),
        ]
    )
    return beta, es.finite_source_magnification(beta, rho, 1.0, -1)


def model_flux(t_c: np.ndarray, ev: dict, cad_d: float, table=None) -> np.ndarray:
    """Exposure-averaged magnification at cadence centres t_c (contiguous exposures of cad_d),
    1 outside +- T_HALF t_E of t0."""
    f = np.ones(len(t_c))
    te, rho, u0, t0 = ev["t_E_d"], ev["rho"], ev["u0"], ev["t0"]
    lo, hi = t0 - T_HALF * te, t0 + T_HALF * te
    sel = (t_c + cad_d / 2 > lo) & (t_c - cad_d / 2 < hi)
    if not sel.any():
        return f
    beta, amp = table if table is not None else amp_table(rho)
    dt = max(min(cad_d / 20, rho * te / 10, te / 200), 1e-7)
    a, b = t_c[sel].min() - cad_d, t_c[sel].max() + cad_d
    tt = np.arange(a, b + dt, dt)
    at = np.interp(np.hypot(u0, (tt - t0) / te), beta, amp, right=1.0)
    cum = np.r_[0, np.cumsum(0.5 * (at[1:] + at[:-1]) * dt)]
    tc = t_c[sel]
    f[sel] = (np.interp(tc + cad_d / 2, tt, cum) - np.interp(tc - cad_d / 2, tt, cum)) / cad_d
    return f


def ideal_identifiable(ev: dict, cad_d: float, t_ref: float, table=None) -> bool:
    """The #561 criterion (D and S, noise-free, no window) on a gap-free cadence grid with the
    light curve's cadence phase (t_ref) over +- T_HALF t_E."""
    te = ev["t_E_d"]
    k0 = math.floor((ev["t0"] - T_HALF * te - t_ref) / cad_d) - 1
    k1 = math.ceil((ev["t0"] + T_HALF * te - t_ref) / cad_d) + 1
    t_c = t_ref + np.arange(k0, k1 + 1) * cad_d
    return ws.identifiable(model_flux(t_c, ev, cad_d, table))[1]


# ----------------------------------------------------------------------------------- the detector


def detect(t: np.ndarray, f: np.ndarray, sigma: float, cad_d: float) -> list[dict]:
    """One fixed detector on good cadences (t ascending, f relative flux). D: >= MIN_RUN contiguous
    cadences with f <= DIM_FLUX and 1 - f >= NSIG sigma. S: in a window of max(run duration,
    SHOULDER_MIN_CAD cadences) on each side, >= SHOULDER_N cadences with f >= SHOULDER_FLUX and
    f - 1 >= NSIG sigma. Returns every D run with its S result."""
    dim = (f <= ws.DIM_FLUX) & (1 - f >= NSIG * sigma)
    bright = (f >= ws.SHOULDER_FLUX) & (f - 1 >= NSIG * sigma)
    brk = np.r_[True, np.diff(t) > CONTIG * cad_d]  # a new contiguous block starts here
    out = []
    i, n = 0, len(f)
    while i < n:
        if not dim[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and dim[j + 1] and not brk[j + 1]:
            j += 1
        if j - i + 1 >= ws.MIN_RUN:
            w = max(t[j] - t[i] + cad_d, SHOULDER_MIN_CAD * cad_d)
            pre = bright & (t < t[i]) & (t >= t[i] - w)
            post = bright & (t > t[j]) & (t <= t[j] + w)
            out.append(
                {
                    "t_start": float(t[i]),
                    "t_end": float(t[j]),
                    "n": int(j - i + 1),
                    "S": bool(pre.sum() >= ws.SHOULDER_N and post.sum() >= ws.SHOULDER_N),
                }
            )
        i = j + 1
    return out


def running_median(f: np.ndarray, n_win: int) -> np.ndarray:
    import pandas as pd

    return (
        pd.Series(f)
        .rolling(n_win, center=True, min_periods=max(10, n_win // 20))
        .median()
        .to_numpy()
    )


def matched(dets: list[dict], ev: dict, need_s: bool = True) -> bool:
    """Post hoc pairing with the truth: a detection overlapping the umbra interval of the event."""
    half = math.sqrt(max(4 - ev["u0"] ** 2, 0)) * ev["t_E_d"]
    lo, hi = ev["t0"] - half, ev["t0"] + half
    return any(d["t_end"] >= lo and d["t_start"] <= hi and (d["S"] or not need_s) for d in dets)


# ---------------------------------------------------------------------------------------- the run


def read_lc(raw: bytes) -> dict:
    from astropy.io import fits

    with fits.open(io.BytesIO(raw)) as h:
        d = h[1].data
        cad_d = float(h[1].header.get("TIMEDEL", 120 / 86400))
        return {
            "time": np.array(d["TIME"], float),
            "flux": np.array(d["PDCSAP_FLUX"], float),
            "quality": np.array(d["QUALITY"], int),
            "cad_d": cad_d,
        }


def process(ss: dict, raw: bytes) -> dict:
    """Background and injections for one star-sector; returns per-event booleans and the
    background flags. The detector never sees the truth; it is used only in matched()."""
    lc = read_lc(raw)
    cad = lc["cad_d"]
    tfin = np.isfinite(lc["time"])
    t_all = lc["time"][tfin]
    good = tfin & np.isfinite(lc["flux"]) & (lc["quality"] == 0)  # our quality screen
    n_win = int(round(BASE_WINDOW_D / cad))

    def l2(flux_full: np.ndarray) -> list[dict]:
        t, fl = lc["time"][good], flux_full[good]
        base = running_median(fl, n_win)
        rel = fl / base
        ok = np.isfinite(rel)
        r = rel[ok] - 1
        sigma = 1.4826 * float(np.median(np.abs(r - np.median(r)))) if ok.any() else np.inf
        return detect(t[ok], rel[ok], sigma, cad)

    bg = l2(lc["flux"])
    res = {
        "obs_id": ss["obs_id"],
        "draws": ss["draws"],
        "n_cadences": int(tfin.sum()),
        "n_good": int(good.sum()),
        "background_D": len(bg),
        "background_DS": sum(d["S"] for d in bg),
        "scenarios": {},
    }
    t_lo, t_hi = float(np.nanmin(lc["time"])), float(np.nanmax(lc["time"]))
    t_ref = float(t_all[0])
    for si, (name, (m, v)) in enumerate(SCENARIOS.items()):
        evs = []
        for k in ss["draws"]:  # each draw of this record gets its own events
            rng = np.random.default_rng([SEED, k, si])
            for ev in draw_events(rng, N_INJ, m, v, ss["radius"], ss["dist"], t_lo, t_hi):
                evs.append(ev | {"draw": k})
        rows = []
        for ev in evs:
            table = amp_table(ev["rho"])
            ideal = ideal_identifiable(ev, cad, t_ref, table)
            f_all = model_flux(t_all, ev, cad, table)  # every cadence in the file (gaps are real)
            l1a = matched(detect(t_all, f_all, 0.0, cad), ev)
            f_full = np.ones(len(lc["time"]))
            f_full[tfin] = f_all
            l1b = matched(detect(lc["time"][good], f_full[good], 0.0, cad), ev)
            inj = lc["flux"] * f_full  # injected before our quality screen and baseline
            dets = l2(inj)
            rows.append(
                {
                    "draw": ev["draw"],
                    "t_E_cad": round(ev["t_E_d"] / cad, 4),
                    "rho": round(ev["rho"], 5),
                    "u0": round(ev["u0"], 4),
                    "ideal": bool(ideal),
                    "L1a": bool(l1a),
                    "L1b": bool(l1b),
                    "L2": bool(matched(dets, ev)),
                    "L2_D_only": bool(matched(dets, ev, need_s=False)),
                    "bg_overlap": bool(matched(bg, ev, need_s=False)),
                }
            )
        res["scenarios"][name] = rows
    return res


def _work(item):
    ss, path = item
    raw = Path(path).read_bytes()
    try:
        return process(ss, raw)
    finally:
        Path(path).unlink(missing_ok=True)


def stage_run(a) -> None:
    sample = json.loads((OUT / "sample.json").read_text())["star_sectors"]
    (OUT / "per_ss").mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    lg = nl.Ledger(OUT / "net_ledger.jsonl", CAP)
    files_path = OUT / "files.json"
    files = json.loads(files_path.read_text()) if files_path.exists() else {}
    todo = [s for s in sample if not (OUT / "per_ss" / f"{s['obs_id']}.json").exists()]

    def fetch(ss):
        b = nl.fetch(DL_URL.format(obs_id=ss["obs_id"]), lg, 2_500_000)
        p = DATA / f"{ss['obs_id']}_lc.fits"
        p.write_bytes(b)
        return ss, str(p), hashlib.sha256(b).hexdigest(), len(b)

    with (
        ThreadPoolExecutor(max_workers=6) as io_pool,
        ProcessPoolExecutor(max_workers=a.workers) as cpu,
    ):
        futs = []
        for ss, path, sha, nb in io_pool.map(fetch, todo):
            files[ss["obs_id"]] = {"sha256": sha, "bytes": nb}
            futs.append(cpu.submit(_work, (ss, path)))
        for fu in futs:
            r = fu.result()
            (OUT / "per_ss" / f"{r['obs_id']}.json").write_text(json.dumps(r) + "\n")
    files_path.write_text(json.dumps(dict(sorted(files.items())), indent=1) + "\n")
    print(len(todo), "processed; ledger bytes", lg.used())


# ---------------------------------------------------------------------------------------- summary


def verdict(n_ideal: int, ret: float, lo: float, hi: float, n_min: int = 30) -> str:
    """Frozen rule (PLAN section 4): retention = recovered / ideal over the same events."""
    if n_ideal < n_min:
        return "insufficient"
    if lo >= 0.5:
        return "retains at least half"
    if hi < 0.5:
        return "retains less than half"
    return "interval crosses one half (not a robust pass)"


def paired_counts(rows: list[dict]) -> dict:
    ideal = [r for r in rows if r["ideal"]]
    return {
        "events": len(rows),
        "ideal": len(ideal),
        "L1a": sum(r["L1a"] for r in ideal),
        "L1b": sum(r["L1b"] for r in ideal),
        "L2": sum(r["L2"] for r in ideal),
        "L2_not_ideal": sum(r["L2"] and not r["ideal"] for r in rows),
        "bg_overlap_in_ideal": sum(r["bg_overlap"] for r in ideal),
    }


def boot_retention(groups: list[list[dict]], rng, n_boot: int = N_BOOT) -> tuple[float, float]:
    """5-95 % of recovered/ideal over star-sectors resampled with replacement."""
    num = np.array([sum(r["ideal"] and r["L2"] for r in g) for g in groups], float)
    den = np.array([sum(r["ideal"] for r in g) for g in groups], float)
    vals = []
    for _ in range(n_boot):
        ix = rng.integers(0, len(groups), len(groups))
        d = den[ix].sum()
        if d > 0:
            vals.append(num[ix].sum() / d)
    return tuple(float(f"{q:.4g}") for q in np.percentile(vals, [5, 95])) if vals else (0.0, 1.0)


def stage_summarize(_a) -> None:
    per = [json.loads(p.read_text()) for p in sorted((OUT / "per_ss").glob("*.json"))]
    rng = np.random.default_rng(SEED)
    reach = json.loads((ws.OUT / "reach.json").read_text())
    masses = reach["masses_msun"]
    out = {
        "label": "simulation (injections) into observed TESS 120 s PDCSAP curves; product level",
        "plan": "results/w3_kepler/tess120_inject/PLAN.md",
        "star_sectors": len(per),
        "draws": sum(len(p["draws"]) for p in per),
        "background": {
            "D_runs_total": sum(p["background_D"] for p in per),
            "DS_flags_total": sum(p["background_DS"] for p in per),
            "star_sectors_with_DS_flag": sum(p["background_DS"] > 0 for p in per),
        },
        "scenarios": {},
    }
    for name, (m, v) in SCENARIOS.items():
        rows = [r for p in per for r in p["scenarios"][name]]
        by_draw: dict[int, list] = {}
        for r in rows:  # bootstrap unit: one draw of a star-sector
            by_draw.setdefault(r["draw"], []).append(r)
        groups = list(by_draw.values())
        c = paired_counts(rows)
        ret = c["L2"] / c["ideal"] if c["ideal"] else float("nan")
        lo, hi = boot_retention(groups, rng)
        rec = c | {"retention": float(f"{ret:.4g}"), "retention_p5_p95": [lo, hi]}
        rec["verdict"] = verdict(c["ideal"], ret, lo, hi)
        rec["ideal_share_of_events"] = (
            float(f"{c['ideal'] / c['events']:.4g}") if c["events"] else None
        )
        bins = {}
        for b0, b1 in zip(TE_BINS[:-1], TE_BINS[1:], strict=True):
            sub = [r for r in rows if b0 <= r["t_E_cad"] < b1]
            cc = paired_counts(sub)
            if cc["events"]:
                cc["retention"] = float(f"{cc['L2'] / cc['ideal']:.4g}") if cc["ideal"] else None
                bins[f"tE_cad_{b0:g}-{b1:g}"] = cc
        rec["by_tE_cadences"] = bins
        # effect on #561: this scenario's mass, v; same TIC sample and scaling as #561
        k = int(np.argmin([abs(math.log10(x / m)) for x in masses]))
        sname = "S1_halo_like_220" if v == 220.0 else "S2_disk_like_30"
        n1 = reach["products"]["tess_120"][sname]["N1_per_density"][k]
        nb = reach["products"]["kepler_lc_baseline"][sname]["N1_per_density"][k]
        rec["reach_561"] = {
            "mass_msun": masses[k],
            "N1_per_density_ideal": n1,
            "N1_per_density_x_retention": float(f"{n1 * ret:.4g}") if n1 else 0.0,
            "N1_x_retention_p5_p95": [float(f"{n1 * lo:.4g}"), float(f"{n1 * hi:.4g}")],
            "density_for_1_ideal": float(f"{1 / n1:.4g}") if n1 else None,
            "density_for_1_x_retention": float(f"{1 / (n1 * ret):.4g}") if n1 and ret > 0 else None,
            "baseline_N1_per_density": nb,
        }
        out["scenarios"][name] = rec
    (OUT / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["select", "run", "summarize"])
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    {"select": stage_select, "run": stage_run, "summarize": stage_summarize}[a.stage](a)


if __name__ == "__main__":
    main()

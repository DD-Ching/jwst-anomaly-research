"""W3 short-cadence physical reach: Kepler SC and TESS 20 s / 120 s against the #559 baseline.

Plan: ``results/w3_kepler/shortcad/PLAN.md``. Only observation metadata, catalogues and pure model
calculations: no target flux, pixels, injection or search.

Stages (``python scripts/w3_shortcad_reach.py <stage>``):

- ``inventory``: MAST CAOM records of Kepler SC and TESS SPOC 20 s / 120 s light curves, their
  coverage, and catalogue R_* and D_S (DR25 for Kepler, TIC v8.2 for TESS) (``inventory.json``,
  ``inputs/*.csv.gz``; network through the ledger);
- ``shape``: w(rho, t_E / cadence), the phase-averaged u0 half-width over which the exposure-
  integrated model light curve keeps an identifiable W3 shape (``shape.json``);
- ``reach``: expected events per unit density for one-|M|, one-speed populations, per product,
  and the investment gate (``reach.json``).

Labels: catalogue radii and distances are model values; coverage spans are upper bounds; the
model light curves are noise-free model predictions; everything else is derived.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import os
import sys
import urllib.parse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import net_ledger as nl  # noqa: E402
import w3_kepler as wk  # noqa: E402
import w3_kepler_physical as wp  # noqa: E402

from jwst_anomaly import exotic_sim as es  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w3_kepler" / "shortcad"
CAP = 150_000_000  # owner limit for this round, bytes
CAOM = "https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync"
TIC = "https://mast.stsci.edu/vo-tap/api/v0.1/tic/sync"
MAST_API = "https://mast.stsci.edu/api/v0/invoke"
SAMPLE_SUFFIX = {"tess_120": "37", "tess_20": "7"}  # TIC ID endings sampled (about 1 % / 10 %)

DAY_S = 86400.0
CADENCE_S = {"kepler_lc": wk.CAD * DAY_S, "kepler_sc": 58.84876, "tess_120": 120.0, "tess_20": 20.0}
SC_MONTH_D = 35.0  # upper bound on one Kepler SC month (Q0: SC_Q0_D)
SC_Q0_D = 10.0
DIM_FLUX = 10 ** (-0.4 * wk.wz.DIM_MAG)  # (D) >= 1 mag: the Kepler pre-screen depth (ASSUMPTION)
MIN_RUN = wk.MIN_RUN  # (D) consecutive dim exposures
SHOULDER_FLUX = 1 + wk.LOCAL_MIN  # (S) bright exposure: F >= 1.10, noise-free (ASSUMPTION)
SHOULDER_N = wk.LOCAL_NCAD  # (S) per side
SCHW4_PC = wp.SCHW4_AU / wp.PC_AU  # 4 G M_sun / c^2 in pc
PC_KM = wp.PC_AU * wp.AU_KM
R_SUN_KM = wp.R_SUN_AU * wp.AU_KM
MASSES = 10.0 ** np.arange(-9.0, 1.0001, 0.25)  # M_sun
GATE_FACTOR = 10.0
GATE_MIN_POINTS = 3  # contiguous grid points (>= 0.5 dex)


# ---------------------------------------------------------------------------------------- network


def _ledger() -> nl.Ledger:
    return nl.Ledger(OUT / "net_ledger.jsonl", CAP)


def tap(url: str, adql: str, lg: nl.Ledger, expect: int = 500_000) -> list[dict]:
    q = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql}
    )
    b = nl.fetch(url + "?" + q, lg, expect)
    return list(csv.DictReader(io.StringIO(b.decode())))


def _write_csv(name: str, rows: list[dict]) -> None:
    (OUT / "inputs").mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / "inputs" / name, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def _read_csv(name: str) -> list[dict]:
    with gzip.open(OUT / "inputs" / name, "rt") as fh:
        return list(csv.DictReader(fh))


def sc_months(obs_id: str) -> tuple[int, float]:
    """(SC months, upper-bound days) from a Kepler CAOM obs_id such as
    ``kplr008264588_sc_Q010000000000000000``: one digit per quarter Q0-Q17 giving the number of SC
    months in it (checked against product lists in ``inventory``)."""
    digits = obs_id.rsplit("_Q", 1)[1]
    n = [int(c) for c in digits]
    days = n[0] * SC_Q0_D + sum(n[1:]) * SC_MONTH_D
    return sum(n), days


def stage_inventory(_a) -> None:
    lg = _ledger()
    inv: dict = {
        "label": "observation metadata (MAST CAOM) and catalogue values; spans are upper bounds"
    }
    # Kepler SC
    sc = tap(
        CAOM,
        "SELECT obsid, obs_id, target_name FROM dbo.ObsPointing WHERE obs_collection='Kepler'"
        " AND dataproduct_type='timeseries' AND t_exptime < 100",
        lg,
        2_000_000,
    )
    for r in sc:
        r["months"], r["days_ub"] = sc_months(r["obs_id"])
    kids = sorted({int(r["target_name"].replace("kplr", "")) for r in sc})
    # check the quarter-string reading against the product lists of three stars
    check = []
    for r in sorted(sc, key=lambda r: -r["months"])[:1] + sorted(sc, key=lambda r: r["months"])[:2]:
        req = {"service": "Mast.Caom.Products", "params": {"obsid": r["obsid"]}, "format": "json"}
        b = nl.fetch(MAST_API + "?request=" + urllib.parse.quote(json.dumps(req)), lg, 3_000_000)
        files = {d["productFilename"] for d in json.loads(b)["data"]}
        check.append(
            {
                "obs_id": r["obs_id"],
                "months_decoded": r["months"],
                "slc_files": sum(f.endswith("_slc.fits") for f in files),
            }
        )
    rows = []
    for i in range(0, len(kids), 500):
        part = ",".join(str(k) for k in kids[i : i + 500])
        q = (
            f"select {wp.COLS} from keplerstellar where st_delivname='q1_q17_dr25_stellar'"
            f" and kepid in ({part})"
        )
        b = nl.fetch(wp.TAP + "?query=" + urllib.parse.quote(q) + "&format=csv", lg, 300_000)
        rows += list(csv.DictReader(io.StringIO(b.decode())))
    by_id = {int(r["kepid"]): r for r in rows if (r.get("kepid") or "").isdigit()}
    out = []
    for r in sc:
        k = int(r["target_name"].replace("kplr", ""))
        s = by_id.get(k, {})
        out.append(
            {
                "id": k,
                "months": r["months"],
                "T_days_ub": r["days_ub"],
                "radius": s.get("radius", ""),
                "dist": s.get("dist", ""),
            }
        )
    _write_csv("kepler_sc.csv.gz", out)
    inv["kepler_sc"] = {
        "records": len(sc),
        "stars": len(kids),
        "sc_months_total": int(sum(r["months"] for r in sc)),
        "quarter_string_check": check,
        "dr25_rows": len(by_id),
    }
    # TESS: totals server-side, per-record rows for 20 s (all) and 120 s (sample)
    base = (
        "FROM dbo.ObsPointing WHERE obs_collection='TESS' AND dataproduct_type='timeseries'"
        " AND provenance_name='SPOC' AND t_exptime={c}"
    )
    for name, c in (("tess_20", 20), ("tess_120", 120)):
        where = f" AND target_name LIKE '%{SAMPLE_SUFFIX[name]}'"
        tot = tap(
            CAOM,
            "SELECT COUNT(*) AS records, COUNT(DISTINCT target_name) AS stars, MIN(sequence_number)"
            " AS s_min, MAX(sequence_number) AS s_max, SUM(t_max - t_min) AS span_d "
            + base.format(c=c),
            lg,
        )[0]
        recs = tap(
            CAOM,
            "SELECT target_name, sequence_number, t_min, t_max " + base.format(c=c) + where,
            lg,
            8_000_000,
        )
        per: dict[str, dict[int, float]] = {}
        for r in recs:  # one span per star and sector (dedupe repeated records)
            d = per.setdefault(r["target_name"], {})
            s = int(r["sequence_number"])
            d[s] = max(d.get(s, 0.0), float(r["t_max"]) - float(r["t_min"]))
        ids = sorted(int(t) for t in per)
        cat = {}
        for i in range(0, len(ids), 400):
            part = ",".join(str(k) for k in ids[i : i + 400])
            for r in tap(
                TIC, f"SELECT ID, rad, d FROM dbo.CatalogRecord WHERE ID IN ({part})", lg, 400_000
            ):
                cat[int(r["id"])] = r
        rows_t = [
            {
                "id": k,
                "sectors": len(per[str(k)]),
                "T_days_ub": round(sum(per[str(k)].values()), 4),
                "radius": cat.get(k, {}).get("rad", ""),
                "dist": cat.get(k, {}).get("d", ""),
            }
            for k in ids
        ]
        _write_csv(f"{name}.csv.gz", rows_t)
        inv[name] = {
            "records_total": int(tot["records"]),
            "stars_total": int(tot["stars"]),
            "sectors": [int(tot["s_min"]), int(tot["s_max"])],
            "span_days_total_ub": float(tot["span_d"]),
            "sample": f"TIC ID ending in {SAMPLE_SUFFIX[name]}",
            "sample_records": len(recs),
            "sample_stars": len(ids),
            "sample_duplicate_records": len(recs) - sum(len(v) for v in per.values()),
            "tic_rows": len(cat),
        }
    inv["bytes_used_ledger"] = lg.used()
    (OUT / "inventory.json").write_text(json.dumps(inv, indent=1) + "\n")
    print(json.dumps(inv, indent=1))


# ------------------------------------------------------------------------------------------ shape

RHO_GRID = np.round(np.logspace(-3, math.log10(3.0), 25), 6)
TAU_GRID = np.round(np.logspace(-0.5, 3.0, 36), 6)  # t_E / cadence
U0_GRID = np.round(np.arange(0.0, 2.6501, 0.05), 4)
N_PHASE = 16
T_HALF = 5.0  # t_E; A(beta = 5) - 1 < 0.004


def identifiable(f: np.ndarray) -> tuple[bool, bool]:
    """(D, D and S) for one exposure series: D = >= MIN_RUN consecutive F <= DIM_FLUX; S = >=
    SHOULDER_N exposures with F >= SHOULDER_FLUX before the first such run and after it."""
    dim = f <= DIM_FLUX
    if dim.sum() < MIN_RUN:
        return False, False
    idx = np.flatnonzero(np.diff(np.r_[0, dim.astype(np.int8), 0]))
    starts, ends = idx[::2], idx[1::2]  # runs [start, end)
    long = np.flatnonzero(ends - starts >= MIN_RUN)
    if not len(long):
        return False, False
    i, j = starts[long[0]], ends[long[0]]
    bright = f >= SHOULDER_FLUX
    return True, bool(bright[:i].sum() >= SHOULDER_N and bright[j:].sum() >= SHOULDER_N)


def shape_for_rho(rho: float) -> tuple[np.ndarray, np.ndarray]:
    """w_D and w_DS on (TAU_GRID) for one rho: trapezoid over U0_GRID of the share of N_PHASE
    exposure phases with an identifiable shape."""
    beta = np.unique(
        np.r_[
            np.linspace(0, 7.5, 7501),
            2 + rho * np.linspace(-1.2, 1.2, 481),
            rho * np.linspace(0, 2, 81),
        ]
    )
    beta = beta[beta >= 0]
    amp = es.finite_source_magnification(beta, rho, 1.0, -1)
    dt = min(1e-3, rho / 20)
    t = np.arange(-T_HALF, T_HALF + dt / 2, dt)
    pd = np.zeros((len(U0_GRID), len(TAU_GRID)))
    pds = np.zeros_like(pd)
    for a, u0 in enumerate(U0_GRID):
        a_t = np.interp(np.hypot(u0, t), beta, amp)
        cum = np.r_[0, np.cumsum(0.5 * (a_t[1:] + a_t[:-1]) * dt)]
        for b, tau in enumerate(TAU_GRID):
            for ph in range(N_PHASE):
                edges = (
                    np.arange(math.floor(-T_HALF * tau) - 1, math.ceil(T_HALF * tau) + 1)
                    + ph / N_PHASE
                ) / tau
                edges = edges[(edges >= t[0]) & (edges <= t[-1])]
                f = np.diff(np.interp(edges, t, cum)) * tau
                d, ds = identifiable(f)
                pd[a, b] += d / N_PHASE
                pds[a, b] += ds / N_PHASE
    return np.trapezoid(pd, U0_GRID, axis=0), np.trapezoid(pds, U0_GRID, axis=0)


def stage_shape(a) -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        res = list(ex.map(shape_for_rho, RHO_GRID))
    out = {
        "label": "model prediction (noise-free, exposure-integrated, finite source, n = 1 sign -1)",
        "plan": "results/w3_kepler/shortcad/PLAN.md",
        "criterion": {
            "D": f">= {MIN_RUN} consecutive exposures with F <= {DIM_FLUX:.4f}",
            "S": f">= {SHOULDER_N} exposures with F >= {SHOULDER_FLUX} on each side of that run",
            "phases": N_PHASE,
            "u0_grid": [float(U0_GRID[0]), float(U0_GRID[-1]), 0.05],
        },
        "rho": RHO_GRID.tolist(),
        "tau_tE_over_cadence": TAU_GRID.tolist(),
        "w_D": [[round(float(v), 4) for v in r[0]] for r in res],
        "w_DS": [[round(float(v), 4) for v in r[1]] for r in res],
    }
    (OUT / "shape.json").write_text(json.dumps(out) + "\n")
    for rho, r in zip(RHO_GRID, res, strict=True):
        print(rho, np.round(r[1][::5], 3))


# ------------------------------------------------------------------------------------------ reach


def load_shape():
    s = json.loads((OUT / "shape.json").read_text())
    return np.log10(s["rho"]), np.log10(s["tau_tE_over_cadence"]), np.array(s["w_DS"])


def w_lookup(rho, tau, shape) -> np.ndarray:
    """Bilinear in (log rho, log tau); 0 outside the rho grid or below the tau grid, the last tau
    column above it (the shape no longer changes with finer sampling)."""
    lr, lt, w = shape
    x = np.log10(np.maximum(rho, 1e-30))
    y = np.clip(np.log10(np.maximum(tau, 1e-30)), None, lt[-1])
    ok = (x >= lr[0]) & (x <= lr[-1]) & (y >= lt[0])
    i = np.clip(np.searchsorted(lr, x) - 1, 0, len(lr) - 2)
    j = np.clip(np.searchsorted(lt, y) - 1, 0, len(lt) - 2)
    fx = np.clip((x - lr[i]) / (lr[i + 1] - lr[i]), 0, 1)
    fy = np.clip((y - lt[j]) / (lt[j + 1] - lt[j]), 0, 1)
    v = (
        w[i, j] * (1 - fx) * (1 - fy)
        + w[i + 1, j] * fx * (1 - fy)
        + w[i, j + 1] * (1 - fx) * fy
        + w[i + 1, j + 1] * fx * fy
    )
    return np.where(ok, v, 0.0)


X_GRID = np.unique(np.r_[np.logspace(-7, -1, 241), np.linspace(0.1, 1 - 1e-7, 361)])


def events_per_density(
    d_pc, r_sun, t_yr, m_msun, v_kms, cadence_s, shape=None, dl_max_pc=None
) -> float:
    """Sum_i T_i (v D_S / |M|) int_0^1 2 R_E(x) w(x) dx for a uniform lens mass density (M_sun /
    pc^3); w = 2 (every umbra crossing) when shape is None, else w(rho(x), t_E(x)/cadence)."""
    d = np.asarray(d_pc, float)[:, None]
    x = X_GRID[None, :]
    r_e_pc = np.sqrt(SCHW4_PC * m_msun * x * (1 - x) * d)
    if shape is None:
        w = 2.0
    else:
        t_e_s = r_e_pc * PC_KM / v_kms
        rho = x * np.asarray(r_sun, float)[:, None] * R_SUN_KM / (r_e_pc * PC_KM)
        w = w_lookup(rho, t_e_s / cadence_s, shape)
    if dl_max_pc is not None:  # post hoc sensitivity: lenses only within dl_max_pc of the observer
        w = w * (x * d <= dl_max_pc)
    integral = np.trapezoid(2 * r_e_pc * w, X_GRID, axis=1)
    v_pc_yr = v_kms * 365.25 * DAY_S / PC_KM
    return float((np.asarray(t_yr) * v_pc_yr * d[:, 0] / m_msun * integral).sum())


def event_properties(d_pc, r_sun, t_yr, m_msun, v_kms, cadence_s, shape) -> dict | None:
    """Rate-weighted 5/50/95 % of t_E (s), rho, D_L (pc) and x over the Level-1 (identifiable)
    events of one population."""
    d = np.asarray(d_pc, float)[:, None]
    x = X_GRID[None, :]
    r_e_pc = np.sqrt(SCHW4_PC * m_msun * x * (1 - x) * d)
    t_e_s = r_e_pc * PC_KM / v_kms
    rho = x * np.asarray(r_sun, float)[:, None] * R_SUN_KM / (r_e_pc * PC_KM)
    w = w_lookup(rho, t_e_s / cadence_s, shape)
    dx = np.gradient(X_GRID)[None, :]
    wt = (np.asarray(t_yr)[:, None] * d * 2 * r_e_pc * w * dx).ravel()
    if not wt.sum() > 0:
        return None

    def q(a):
        a = np.broadcast_to(a, w.shape).ravel()
        o = np.argsort(a)
        c = np.cumsum(wt[o]) / wt.sum()
        return [float(f"{a[o][np.searchsorted(c, p)]:.3g}") for p in (0.05, 0.5, 0.95)]

    return {"t_E_s": q(t_e_s), "rho": q(rho), "D_L_pc": q(x * d), "x": q(x)}


def gate(ratios, factor: float = GATE_FACTOR, min_points: int = GATE_MIN_POINTS) -> dict:
    """PLAN section 4: passes if >= min_points contiguous grid points have ratio >= factor
    (None = baseline 0 with a non-zero product reach, counted as >= factor)."""
    ok = [(r is None) or (r >= factor) for r in ratios]
    runs, cur = [], []
    for i, v in enumerate(ok):
        if v:
            cur.append(i)
        else:
            if cur:
                runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    runs = [r for r in runs if len(r) >= min_points]
    return {"passes": bool(runs), "intervals_idx": [[r[0], r[-1]] for r in runs]}


def products():
    base = wp.load_stellar()
    yield (
        "kepler_lc_baseline",
        {
            "d": [s["dist"] for s in base],
            "r": [s["radius"] for s in base],
            "t": [1459.49 / 365.25] * len(base),
            "cad": CADENCE_S["kepler_lc"],
            "scale": 1.0,
            "n_used": len(base),
        },
    )
    inv = json.loads((OUT / "inventory.json").read_text())
    for name, cad, tot in (
        ("kepler_sc", "kepler_sc", None),
        ("tess_20", "tess_20", "tess_20"),
        ("tess_120", "tess_120", "tess_120"),
    ):
        rows = _read_csv(f"{name}.csv.gz")
        use = [r for r in rows if r["radius"] and r["dist"]]
        scale = 1.0
        if tot and inv[tot]["sample_stars"] < inv[tot]["stars_total"]:
            scale = inv[tot]["stars_total"] / inv[tot]["sample_stars"]
        yield (
            name,
            {
                "d": [float(r["dist"]) for r in use],
                "r": [float(r["radius"]) for r in use],
                "t": [float(r["T_days_ub"]) / 365.25 for r in use],
                "cad": CADENCE_S[cad],
                "scale": scale,
                "n_used": len(use),
                "n_rows": len(rows),
            },
        )


def stage_reach(_a) -> None:
    shape = load_shape()
    prods = dict(products())
    res = {}
    rng = np.random.default_rng(559)
    for name, p in prods.items():
        rec = {
            "stars_used": p["n_used"],
            "stars_in_table": p.get("n_rows", p["n_used"]),
            "scale": p["scale"],
        }
        rec["cadence_s"] = p["cad"]
        rec["sum_T_yr_used"] = round(float(np.sum(p["t"])), 3)
        rec["median_D_pc"] = round(float(np.median(p["d"])), 1)
        rec["median_R_rsun"] = round(float(np.median(p["r"])), 3)
        for sname, v in wp.SPEEDS.items():
            l0, l1 = [], []
            for m in MASSES:
                l0.append(p["scale"] * events_per_density(p["d"], p["r"], p["t"], m, v, p["cad"]))
                l1.append(
                    p["scale"] * events_per_density(p["d"], p["r"], p["t"], m, v, p["cad"], shape)
                )
            rec[sname] = {
                "N0_per_density": [float(f"{x:.4g}") for x in l0],
                "N1_per_density": [float(f"{x:.4g}") for x in l1],
            }
            k_best = int(np.argmax(l1))
            rec[sname]["best_mass_msun"] = float(MASSES[k_best])
            rec[sname]["events_at_best_mass"] = event_properties(
                p["d"], p["r"], p["t"], MASSES[k_best], v, p["cad"], shape
            )
            if p["scale"] > 1:  # bootstrap over the sampled stars at the best Level-1 mass
                k = int(np.argmax(l1))
                n = len(p["d"])
                bs = []
                for _ in range(100):
                    ix = rng.integers(0, n, n)
                    bs.append(
                        p["scale"]
                        * events_per_density(
                            np.take(p["d"], ix),
                            np.take(p["r"], ix),
                            np.take(p["t"], ix),
                            MASSES[k],
                            v,
                            p["cad"],
                            shape,
                        )
                    )
                rec[sname]["N1_best_bootstrap_p5_p95"] = [
                    float(f"{q:.4g}") for q in np.percentile(bs, [5, 95])
                ]
        res[name] = rec
    base = res["kepler_lc_baseline"]
    for name in ("kepler_sc", "tess_20", "tess_120"):
        for sname in wp.SPEEDS:
            n1, b1, b0 = (
                res[name][sname]["N1_per_density"],
                base[sname]["N1_per_density"],
                base[sname]["N0_per_density"],
            )
            r1 = [
                None if (b == 0 and n > 0) else (n / b if b > 0 else 0.0)
                for n, b in zip(n1, b1, strict=True)
            ]
            res[name][sname]["ratio_vs_baseline_N1"] = [
                None if r is None else float(f"{r:.4g}") for r in r1
            ]
            res[name][sname]["ratio_vs_baseline_N0"] = [
                float(f"{n / b:.4g}") for n, b in zip(n1, b0, strict=True)
            ]
            g = gate(r1)
            g["intervals_msun"] = [
                [float(MASSES[i]), float(MASSES[j])] for i, j in g["intervals_idx"]
            ]
            res[name][sname]["gate"] = g
    out = {
        "label": "derived (optimistic reach per unit density; not an abundance, limit, exclusion)",
        "plan": "results/w3_kepler/shortcad/PLAN.md",
        "masses_msun": [float(f"{m:.4g}") for m in MASSES],
        "speeds_kms": wp.SPEEDS,
        "references_msun_pc3": wp.REF,
        "products": res,
    }
    (OUT / "reach.json").write_text(json.dumps(out, indent=1) + "\n")
    for name, rec in res.items():
        for sname in wp.SPEEDS:
            n1 = rec[sname]["N1_per_density"]
            k = int(np.argmax(n1))
            print(
                name,
                sname,
                "best M",
                MASSES[k],
                "N1",
                n1[k],
                "rho_for_1",
                1 / n1[k] if n1[k] else None,
                rec[sname].get("gate"),
            )


DL_MAX_PC = 1000.0  # post hoc sensitivity (RESULT section 3): uniform lens density only to 1 kpc


def _sens_one(item):
    name, p, shape = item
    out = {}
    for sname, v in wp.SPEEDS.items():
        vals = []
        for m in MASSES[::2]:
            n = events_per_density(p["d"], p["r"], p["t"], m, v, p["cad"], shape, DL_MAX_PC)
            vals.append(float(f"{p['scale'] * n:.4g}"))
        out[sname] = vals
    return name, out


def stage_sensitivity(a) -> None:
    """Post hoc: Level-1 reach with lenses only at D_L <= DL_MAX_PC, on the 0.5-dex mass grid."""
    shape = load_shape()
    items = [(n, p, shape) for n, p in products()]
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        res = dict(ex.map(_sens_one, items))
    base = res["kepler_lc_baseline"]
    out = {
        "label": "derived, post hoc sensitivity (not in PLAN.md): lenses only within D_L <= 1 kpc",
        "dl_max_pc": DL_MAX_PC,
        "masses_msun": [float(f"{m:.4g}") for m in MASSES[::2]],
        "N1_per_density": res,
        "ratio_vs_baseline_N1": {
            n: {
                s: [
                    None if (b == 0 and x > 0) else (float(f"{x / b:.4g}") if b > 0 else 0.0)
                    for x, b in zip(res[n][s], base[s], strict=True)
                ]
                for s in wp.SPEEDS
            }
            for n in res
            if n != "kepler_lc_baseline"
        },
    }
    (OUT / "sensitivity_dl1kpc.json").write_text(json.dumps(out, indent=1) + "\n")
    for n in res:
        for s in wp.SPEEDS:
            k = int(np.argmax(res[n][s]))
            print(n, s, "best M", MASSES[::2][k], "N1", res[n][s][k])


def stage_summary(_a) -> None:
    """Per-product distributions for question 1 (stars, coverage per star, R_*, D_S)."""
    inv = json.loads((OUT / "inventory.json").read_text())
    out = {"label": "observation metadata and catalogue values; T per star is an upper bound"}
    for name, p in products():
        rows = _read_csv(f"{name}.csv.gz") if name != "kepler_lc_baseline" else None
        rec = {
            "cadence_s": round(p["cad"], 3),
            "stars_with_R_and_D": p["n_used"],
            "stars_in_table": p.get("n_rows", p["n_used"]),
            "scale_to_all_stars": round(p["scale"], 4),
            "T_days_per_star_p5_50_95": wp._pct([t * 365.25 for t in p["t"]]),
            "sum_T_star_years_used_x_scale": float(f"{np.sum(p['t']) * p['scale']:.4g}"),
            "D_pc_p5_50_95": wp._pct(p["d"]),
            "R_rsun_p5_50_95": wp._pct(p["r"]),
            "sum_T_D15_pc15_yr_x_scale": float(
                f"{np.sum(np.asarray(p['t']) * np.asarray(p['d']) ** 1.5) * p['scale']:.4g}"
            ),
        }
        if rows:
            key = "months" if name == "kepler_sc" else "sectors"
            rec[f"{key}_per_star_p5_50_95"] = wp._pct([float(r[key]) for r in rows])
            rec["share_R_lt_0p6"] = round(float(np.mean(np.asarray(p["r"]) < 0.6)), 4)
            rec["share_R_gt_2"] = round(float(np.mean(np.asarray(p["r"]) > 2.0)), 4)
        out[name] = rec
    out["tess_note"] = {
        "sector_numbers": "1-107 and one 2026 segment labelled 1751 (about 7 d)",
        "totals": {k: inv[k] for k in ("tess_20", "tess_120")},
    }
    (OUT / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["inventory", "shape", "reach", "summary", "sensitivity"])
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    {
        "inventory": stage_inventory,
        "shape": stage_shape,
        "reach": stage_reach,
        "summary": stage_summary,
        "sensitivity": stage_sensitivity,
    }[a.stage](a)


if __name__ == "__main__":
    main()

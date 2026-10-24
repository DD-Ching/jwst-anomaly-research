"""W3 (negative-mass point lens, n = 1, finite source) pilot on Kepler DR25 long-cadence SAP flux.

Plan: ``results/w3_kepler/PLAN.md`` (frozen before any target flux was read). A thin data interface
around the existing W3 model (``exotic_sim``) and the ZTF W3 chain (``w3_ztf``: ``runs``, ``vet``
with the v23 "bounds" edge rule, ``edge_bounds_alt``, ``control_model``); no photometry, PSF or
optimiser is built here.

Stages (``python scripts/w3_kepler.py <stage>``):

- ``sample``: metadata-only selection and the dev / val / search split (``sample.json``);
- ``fetch --group dev|val|search``: DR25 ``_llc.fits`` files through the download ledger; only the
  needed columns are cached (data root, gitignored); sizes and sha256 go to ``files_<group>.json``;
- ``inject --phase dev|val``: W3 cells and ordinary controls on that phase's carriers, plus the
  uninjected screen of the same stars (``injections_<phase>.json``, ``real_<phase>.json``);
- ``search``: only after ``gate.json`` says the validation passed.

Labels: injections are simulated; light curves and their flags observed; outcomes derived.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import pickle
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import net_ledger as nl  # noqa: E402
import w3_ztf as wz  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w3_kepler"
CACHE = Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data")) / "cache" / "w3_kepler"
LEDGER_CAP = 2_000_000_000  # owner limit for this round, bytes
SALT = "w3-kepler-2026-10-24"
POOL_SHA256 = "013b8a0d0f464a9574a3e4e8a1dce0c1b306116b468c18f638efe83fa13c60b0"
GROUPS = {"dev": (0, 20), "val": (20, 40), "search": (40, 100)}
ARCHIVE = "https://archive.stsci.edu/missions/kepler/lightcurves"

# Quality: a cadence is excluded if any of these SAP_QUALITY bits is set: attitude tweak 1, safe
# mode 2, coarse point 4, Earth point 8, Argabrightening 16, desaturation 32, aperture cosmic ray
# 64, manual exclude 128, impulsive outlier 512, collateral cosmic ray 1024, bad calibration 16384,
# no data 65536.
QMASK = 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 | 512 | 1024 | 16384 | 65536
CAD = 0.02043359821692  # d, TIMEDEL of long cadence
NSUB = 21  # sub-samples per cadence for the exposure average
MIN_RUN = 3  # consecutive dim cadences
MIN_BRACKET = 20  # non-dim cadences on each side (anywhere in the light curve), as ZTF
SHOULDER_MIN = 0.10  # W3 shoulder: rel >= 1 + SHOULDER_MIN ...
SHOULDER_SNR = 5.0  # ... and (rel - 1) / sigma >= SHOULDER_SNR, on both sides
CELLS = [
    (rho, blend, te)
    for rho, blend in ((0.01, 1.0), (0.1, 1.0), (0.01, 0.7), (0.1, 0.7))
    for te in (0.1, 0.3, 1.0, 3.0, 10.0)
]
CONTROLS = ("eclipse", "single_eclipse", "fade", "pspl", "common", "quarter_step")
N_INJ = 100
SEEDS = {"dev": (3100, 3200), "val": (3300, 3400)}
GATE_REC = 0.50
GATE_CTRL = 0.10
U_MAX = 12.0  # beyond this impact parameter the W3 magnification is set to 1 (A(12) - 1 = 1e-4)


# ------------------------------------------------------------------------------------------ sample


def stage_sample(a) -> None:
    raw = Path(a.pool).read_bytes()
    if hashlib.sha256(raw).hexdigest() != POOL_SHA256:
        raise SystemExit("pool file differs from the pinned sha256")
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    pool = [r for r in rows if r["st_quarters"] == "1" * 17]
    rank = sorted(
        pool, key=lambda r: hashlib.sha256(f"{SALT}:{int(r['kepid'])}".encode()).hexdigest()
    )
    stars = []
    for i, r in enumerate(rank[:100]):
        g = next(k for k, (lo, hi) in GROUPS.items() if lo <= i < hi)
        stars.append(
            {
                "rank": i,
                "kepid": int(r["kepid"]),
                "group": g,
                **{k: float(r[k]) for k in ("kepmag", "teff", "logg", "ra", "dec")},
            }
        )
    out = {
        "label": "observed (metadata)",
        "pool_query": a.query,
        "pool_sha256": POOL_SHA256,
        "pool_rows": len(rows),
        "pool_all_q1_q17": len(pool),
        "rule": f"sort by sha256('{SALT}:<kepid>'); ranks 0-19 dev, 20-39 val, 40-99 search",
        "interface_test_star": {"rank": 100, "kepid": int(rank[100]["kepid"]), "target": False},
        "stars": stars,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "sample.json").write_text(json.dumps(out, indent=1) + "\n")
    print({g: sum(s["group"] == g for s in stars) for g in GROUPS})


def sample_stars(groups) -> list[dict]:
    stars = json.loads((OUT / "sample.json").read_text())["stars"]
    return [s for s in stars if s["group"] in groups]


def log_usage(stage: str, groups) -> None:
    rec = {"t": round(time.time(), 1), "stage": stage, "groups": sorted(groups)}
    with open(OUT / "source_usage.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


# ------------------------------------------------------------------------------------------- fetch

COLS = ("TIME", "CADENCENO", "SAP_FLUX", "SAP_FLUX_ERR", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR")
COLS += ("SAP_QUALITY", "MOM_CENTR1", "MOM_CENTR2", "SAP_BKG")


def _fetch_star(kepid: int) -> dict:
    from astropy.io import fits

    lg = nl.Ledger(OUT / "net_ledger.jsonl", LEDGER_CAP)
    kid = f"{kepid:09d}"
    base = f"{ARCHIVE}/{kid[:4]}/{kid}/"
    listing = nl.fetch(base, lg, 20000).decode()
    names = sorted(set(re.findall(r"kplr\d{9}-\d{13}_llc\.fits", listing)))
    parts, files = [], []
    for name in names:
        b = nl.fetch(base + name, lg, est_bytes=520_000)
        with fits.open(io.BytesIO(b)) as h:
            h0, h1, d = h[0].header, h[1].header, h[1].data
            q = int(h0["QUARTER"])
            files.append(
                {
                    "name": name,
                    "bytes": len(b),
                    "sha256": hashlib.sha256(b).hexdigest(),
                    "quarter": q,
                    "channel": int(h0["CHANNEL"]),
                    "data_rel": int(h0["DATA_REL"]),
                    "crowdsap": h1.get("CROWDSAP"),
                }
            )
            if not 1 <= q <= 17:
                continue  # Q0 is not used (plan)
            arr = {c: np.asarray(d[c], float) for c in COLS}
            arr["QUARTER"] = np.full(len(d), q, float)
            arr["CHANNEL"] = np.full(len(d), int(h0["CHANNEL"]), float)
            parts.append(arr)
    cat = {c: np.concatenate([p[c] for p in parts]) for c in parts[0]}
    o = np.argsort(cat["TIME"], kind="stable")
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / f"{kid}.npz", **{c: v[o] for c, v in cat.items()})
    return {"kepid": kepid, "files": files}


def stage_fetch(a) -> None:
    stars = sample_stars([a.group])
    log_usage(f"fetch {a.group}", [a.group])
    with ThreadPoolExecutor(8) as ex:
        info = list(ex.map(_fetch_star, [s["kepid"] for s in stars]))
    out = {"label": "observed", "stars": info}
    (OUT / f"files_{a.group}.json").write_text(json.dumps(out) + "\n")
    used = nl.Ledger(OUT / "net_ledger.jsonl", LEDGER_CAP).used()
    print(a.group, len(info), "stars; ledger", used)


def load_star(kepid: int) -> dict:
    with np.load(CACHE / f"{kepid:09d}.npz") as z:
        return {k: z[k] for k in z.files}


# ---------------------------------------------------------------------------------- thin interface


def good_mask(st: dict) -> np.ndarray:
    """Plan section 3: finite time and SAP flux, positive SAP error, no QMASK bit set."""
    qual = st["SAP_QUALITY"].astype(np.int64)
    return (
        np.isfinite(st["TIME"])
        & np.isfinite(st["SAP_FLUX"])
        & (st["SAP_FLUX_ERR"] > 0)
        & ((qual & QMASK) == 0)
    )


_ATAB: dict[float, tuple[np.ndarray, np.ndarray]] = {}


def a_table(rho: float) -> tuple[np.ndarray, np.ndarray]:
    """exotic_sim's finite-source W3 magnification A(u; rho) (n = 1, sign -1) on a dense u grid,
    interpolated (checked against direct evaluation in the tests)."""
    if rho not in _ATAB:
        from jwst_anomaly import exotic_sim as es

        u = np.unique(
            np.concatenate(
                [
                    np.linspace(0, 1.5, 301),
                    np.linspace(1.5, 2.5, 20001),
                    np.linspace(2.5, U_MAX, 951),
                ]
            )
        )
        _ATAB[rho] = (u, np.asarray(es.light_curve(u, 0.0, 1.0, 0.0, 1.0, -1, rho)))
    return _ATAB[rho]


def w3_flux_avg(t, t0, te, u0, rho, blend) -> np.ndarray:
    """Exposure-averaged W3 flux factor (blend * A + 1 - blend) over [t - CAD/2, t + CAD/2]."""
    t = np.asarray(t, float)
    f = np.ones(len(t))
    near = np.abs(t - t0) <= te * U_MAX + CAD
    if near.any():
        u_tab, a_tab = a_table(rho)
        sub = t[near, None] + np.linspace(-CAD / 2, CAD / 2, NSUB)[None, :]
        u = np.hypot(u0, (sub - t0) / te)
        a = np.where(u <= U_MAX, np.interp(u, u_tab, a_tab), 1.0)
        f[near] = (blend * a + 1.0 - blend).mean(1)
    return f


def dmag_flux_avg(model, t) -> np.ndarray:
    """Exposure-averaged flux factor of a w3_ztf.control_model-like dmag(t) closure."""
    t = np.asarray(t, float)
    sub = t[:, None] + np.linspace(-CAD / 2, CAD / 2, NSUB)[None, :]
    dm = np.asarray(model(sub.ravel()), float).reshape(sub.shape)
    return (10 ** (-0.4 * dm)).mean(1)


def inject_sap(sap, err, f, rng) -> tuple[np.ndarray, np.ndarray]:
    """Plan section 4: SAP_inj = SAP * F plus Gaussian noise, so that the variance is
    err^2 * max(F, 0.09) (the sky-limited floor of w3_microlensing.inject_w3); err_inj is the square
    root of that variance. Cadences with F == 1 are unchanged. No baseline is estimated first."""
    f = np.asarray(f, float)
    sap_i, err_i = np.array(sap, float), np.array(err, float)
    m = f != 1.0
    if m.any():
        var_t = err[m] ** 2 * np.maximum(f[m], wz.INJ_FLOOR)
        extra = np.sqrt(np.maximum(var_t - (err[m] * f[m]) ** 2, 0.0))
        sap_i[m] = sap[m] * f[m] + rng.normal(0.0, 1.0, m.sum()) * extra
        err_i[m] = np.sqrt(var_t)
    return sap_i, err_i


def normalise(st, sap, err, good) -> tuple[np.ndarray, np.ndarray]:
    """Per-quarter median of the good cadences (plan section 3); rel and sigma per good cadence. A
    quarter with fewer than 2 * MIN_BRACKET good cadences is left out."""
    rel, sig = np.full(len(sap), np.nan), np.full(len(sap), np.nan)
    for q in np.unique(st["QUARTER"]):
        m = (st["QUARTER"] == q) & good
        if m.sum() < 2 * MIN_BRACKET:
            continue
        med = np.median(sap[m])
        if med > 0:
            rel[m], sig[m] = sap[m] / med, err[m] / med
    return rel, sig


def series(st, rel, sig, frac) -> tuple[list, np.ndarray]:
    """w3_ztf series tuples (t, dm, nd, lim, q, k, base, frac) for the normalised good cadences.
    rel <= 0 is a non-detection (dim). lim is set so that w3_ztf's maglimit error formula returns
    the Kepler sigma of dm."""
    idx = np.flatnonzero(np.isfinite(rel) & np.isfinite(sig))
    r, sg = rel[idx], sig[idx]
    pos = r > 0
    rr = np.where(pos, r, 1.0)
    dm = np.where(pos, -2.5 * np.log10(rr), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        lim = np.where(pos, dm - 2.5 * np.log10(1.0857 * sg / rr / 0.2171), -2.5 * np.log10(5 * sg))
    ks = st["CADENCENO"][idx].astype(int)
    out = [
        (float(a), float(b), not bool(c), float(d), 0, int(k), 0.0, frac.get(int(k), 0.0))
        for a, b, c, d, k in zip(st["TIME"][idx], dm, pos, lim, ks, strict=True)
    ]
    return out, idx


def flagged_runs_kepler(s) -> list[tuple[int, int]]:
    """w3_ztf.flagged_runs without the ground-based night condition: >= MIN_RUN consecutive dim
    cadences and >= MIN_BRACKET non-dim cadences on each side."""
    d = wz.dim_state(s)
    out = []
    for i, j in wz.runs(d):
        if j - i + 1 >= MIN_RUN and (~d[:i]).sum() >= MIN_BRACKET:
            if (~d[j + 1 :]).sum() >= MIN_BRACKET:
                out.append((i, j))
    return out


def edge_bounds_kepler(s, quarter, run, span) -> dict:
    """edge_bounds (v23) per side on the cadences of the quarter holding that run edge only: an
    edge whose nearest normal cadence is in another quarter (separate aperture and normalisation)
    is not bounded."""
    i, j = run
    ar = np.arange(len(s))
    pre_idx = np.flatnonzero((quarter == quarter[i]) & (ar <= j)).tolist()
    post_idx = np.flatnonzero((quarter == quarter[j]) & (ar >= i)).tolist()
    pre = wz.edge_bounds_alt(
        [s[k] for k in pre_idx], (pre_idx.index(i), len(pre_idx) - 1), span, ()
    )
    post = wz.edge_bounds_alt([s[k] for k in post_idx], (0, post_idx.index(j)), span, ())
    cls = {"pre": pre["edge_pre_class"], "post": post["edge_post_class"]}
    if all(c == "fast" for c in cls.values()):
        run_cls = "fast"
    else:
        run_cls = "slow" if "slow" in cls.values() else "insufficient"
    return {
        "edge_class": run_cls,
        "edge_pre_class": cls["pre"],
        "edge_post_class": cls["post"],
        "edge_pre_d": pre["edge_pre_d"],
        "edge_post_d": post["edge_post_d"],
        "edge_pre_max_d": pre["edge_pre_max_d"],
        "edge_post_max_d": post["edge_post_max_d"],
    }


def shoulder(s, quarter, rel, sig, run, span) -> dict:
    """W3 shoulder test: on each side, within EDGE_FRAC * span of the run edge and in that edge's
    quarter, some cadence with rel >= 1 + SHOULDER_MIN and (rel - 1) / sigma >= SHOULDER_SNR.
    Both sides are required."""
    i, j = run
    t = np.array([x[0] for x in s])
    w = wz.EDGE_FRAC * span
    hit = (rel >= 1 + SHOULDER_MIN) & ((rel - 1) / sig >= SHOULDER_SNR)
    pre = hit & (quarter == quarter[i]) & (t < t[i]) & (t >= t[i] - w)
    post = hit & (quarter == quarter[j]) & (t > t[j]) & (t <= t[j] + w)
    ok = bool(pre.any() and post.any())
    return {"shoulder_pre": bool(pre.any()), "shoulder_post": bool(post.any()), "shoulder_ok": ok}


def screen_star(st, sap, err, frac, frac_override=None) -> list[dict]:
    """Quality mask, per-quarter normalisation, run finding and vetting of one (possibly injected)
    light curve. pass = common mode, bounded fast edges (Kepler quarter rule), one-off and W3
    shoulders. The achromatic test is not applicable (single band) and never counts as passed.
    pass_no_shoulder: the same without the shoulder test (secondary)."""
    good = good_mask(st)
    rel, sig = normalise(st, sap, err, good)
    fr = dict(frac)
    if frac_override:
        fr.update(frac_override)
    s, idx = series(st, rel, sig, fr)
    if len(s) < 2 * MIN_BRACKET + MIN_RUN:
        return []
    quarter = st["QUARTER"][idx]
    r_s, sg_s = rel[idx], sig[idx]
    out = []
    for run in flagged_runs_kepler(s):
        i, j = run
        v = wz.vet(s, run, None, None)  # common mode and one-off only; edges: edge_bounds_kepler
        span = max(s[j][0] - s[i][0], 1.0)
        e = edge_bounds_kepler(s, quarter, run, span)
        sh = shoulder(s, quarter, r_s, sg_s, run, span)
        rec = {
            "t_start": float(s[i][0]),
            "t_end": float(s[j][0]),
            "n": j - i + 1,
            "quarter_start": int(quarter[i]),
            "common_mode_frac": v["common_mode_frac"],
            "common_mode_ok": bool(v["common_mode_ok"]),
            "n_other_runs": v["n_other_runs"],
            "oneoff_ok": bool(v["oneoff_ok"]),
            "achromatic": "not tested (single band)",
            **e,
            **sh,
        }
        rec["sharp_ok"] = e["edge_class"] == "fast"
        base_ok = rec["common_mode_ok"] and rec["sharp_ok"] and rec["oneoff_ok"]
        rec["pass_no_shoulder"] = bool(base_ok)
        rec["pass"] = bool(base_ok and sh["shoulder_ok"])
        rec["insufficient"] = bool(
            e["edge_class"] == "insufficient" and rec["common_mode_ok"] and rec["oneoff_ok"]
        )
        out.append(rec)
    return out


def common_frac(stars: dict) -> dict:
    """Per cadence number: counts of the given uninjected stars with a good cadence there (den) and
    of those that are dim (num; rel <= 10^(-0.4 DIM_MAG), rel <= 0 included); per star, its own
    good and dim cadences, so that frac_for can leave it out."""
    num: dict[int, int] = {}
    den: dict[int, int] = {}
    per: dict[int, tuple[set, set]] = {}
    for kid, st in stars.items():
        good = good_mask(st)
        rel, _ = normalise(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], good)
        ok = np.isfinite(rel)
        ks = st["CADENCENO"][ok].astype(int)
        dim = rel[ok] <= 10 ** (-0.4 * wz.DIM_MAG)
        per[kid] = (set(ks.tolist()), set(ks[dim].tolist()))
        for k in ks.tolist():
            den[k] = den.get(k, 0) + 1
        for k in ks[dim].tolist():
            num[k] = num.get(k, 0) + 1
    return {"num": num, "den": den, "per": per}


def frac_for(cf: dict, kid: int) -> dict:
    """Leave-one-out dim fraction of the other stars per cadence (non-zero entries only)."""
    own_good, own_dim = cf["per"].get(kid, (set(), set()))
    out = {}
    for k, n in cf["num"].items():
        n -= k in own_dim
        d = cf["den"][k] - (k in own_good)
        if d > 0 and n > 0:
            out[k] = n / d
    return out


def pdc_in_run(st, rec) -> float | None:
    """Diagnostic, not a test: median PDCSAP flux inside a flagged run over the median of its
    quarter's good cadences (same observation as SAP, never an independent confirmation)."""
    good = good_mask(st) & np.isfinite(st["PDCSAP_FLUX"])
    t = st["TIME"]
    m = good & (t >= rec["t_start"]) & (t <= rec["t_end"])
    if not m.any():
        return None
    q = st["QUARTER"][m][0]
    med = np.median(st["PDCSAP_FLUX"][good & (st["QUARTER"] == q)])
    return float(np.median(st["PDCSAP_FLUX"][m]) / med) if med > 0 else None


# -------------------------------------------------------------------------------------- injections


def _box(t1: float, t2: float, depth: float = 1.5):
    def f(x):
        x = np.asarray(x)
        return np.where((x >= t1) & (x < t2), depth, 0.0)

    return f


def control_flux(kind, st, good, rng):
    """Exposure-averaged flux factor F(t) of an ordinary control and the common-mode override
    (common only). eclipse, fade and pspl reuse w3_ztf.control_model unchanged."""
    t = st["TIME"]
    tg = t[good]
    fin = np.isfinite(t)
    f = np.ones(len(t))
    if kind in ("eclipse", "fade", "pspl"):
        f[fin] = dmag_flux_avg(wz.control_model(kind, tg, rng), t[fin])
        return f, None
    if kind == "single_eclipse":
        t1, dur = rng.uniform(tg.min(), tg.max()), rng.uniform(0.1, 2.0)
        f[fin] = dmag_flux_avg(_box(t1, t1 + dur), t[fin])
        return f, None
    if kind == "quarter_step":
        t1 = rng.uniform(tg.min(), tg.max())
        q = st["QUARTER"][good][np.searchsorted(tg, t1).clip(0, len(tg) - 1)]
        t_end = t[(st["QUARTER"] == q) & fin].max() + CAD
        f[fin] = dmag_flux_avg(_box(t1, t_end), t[fin])
        return f, None
    if kind == "common":
        gi = np.flatnonzero(good)
        k0 = rng.integers(MIN_BRACKET, len(gi) - MIN_BRACKET - 5)
        f[gi[k0 : k0 + 5]] = 10 ** (-0.4 * 3.0)
        return f, {int(st["CADENCENO"][g]): 0.5 for g in gi[k0 : k0 + 5]}
    raise ValueError(kind)


_CTX: dict = {}


def _init(ctx_path: str) -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    with open(ctx_path, "rb") as fh:
        _CTX.update(pickle.load(fh))


def _one(job):
    """One injection; job = (kind, cell index or control name, j, seed)."""
    kind, c, j, seed = job
    rng = np.random.default_rng([seed, j])
    kids = _CTX["kids"]
    kid = kids[rng.integers(len(kids))]
    st, frac, base = _CTX["stars"][kid], _CTX["frac"][kid], _CTX["base"][kid]
    good = good_mask(st)
    tg = st["TIME"][good]
    over = None
    if kind == "w3":
        rho, blend, te = CELLS[c]
        t0, u0 = rng.uniform(tg.min(), tg.max()), rng.uniform(0, 2)
        f = w3_flux_avg(st["TIME"], t0, te, u0, rho, blend)
        half = te * np.sqrt(4 - u0**2)
    else:
        f, over = control_flux(c, st, good, rng)
    sap, err = inject_sap(st["SAP_FLUX"], st["SAP_FLUX_ERR"], f, rng)
    res = screen_star(st, sap, err, frac, over)
    new = [x for x in res if (x["t_start"], x["t_end"]) not in base]
    if kind == "w3":
        new = [x for x in new if x["t_start"] <= t0 + half and x["t_end"] >= t0 - half]
    row = {
        "kind": kind,
        "cell": c,
        "j": j,
        "kepid": int(kid),
        "flagged": bool(new),
        "pass": any(x["pass"] for x in new),
        "pass_no_shoulder": any(x["pass_no_shoulder"] for x in new),
        "insufficient": (not any(x["pass_no_shoulder"] for x in new))
        and any(x["insufficient"] for x in new),
    }
    if kind == "w3":
        row.update({"t0": float(t0), "u0": float(u0)})
        best = [x for x in new if x["pass"]] or [x for x in new if x["pass_no_shoulder"]] or new
        if best:
            keys = ("edge_class", "edge_pre_max_d", "edge_post_max_d", "edge_pre_d")
            keys += ("edge_post_d", "shoulder_ok", "oneoff_ok", "common_mode_ok")
            row.update({k: best[0][k] for k in keys})
    return row


def wilson(k: int, n: int) -> list[float]:
    if n == 0:
        return [0.0, 1.0]
    z, p = 1.96, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(float(c - h), 3), round(float(c + h), 3)]


def build_ctx(carrier_groups, frac_groups, tag: str) -> str:
    """Load the carriers and the common-mode comparison stars, screen the uninjected carriers."""
    need = sample_stars(set(carrier_groups) | set(frac_groups))
    stars = {s["kepid"]: load_star(s["kepid"]) for s in need}
    fk = {s["kepid"] for s in sample_stars(frac_groups)}
    cf = common_frac({k: v for k, v in stars.items() if k in fk})
    kids = [s["kepid"] for s in sample_stars(carrier_groups)]
    frac = {k: frac_for(cf, k) for k in kids}
    real = {}
    for k in kids:
        real[k] = screen_star(stars[k], stars[k]["SAP_FLUX"], stars[k]["SAP_FLUX_ERR"], frac[k])
        for x in real[k]:
            x["pdc_rel_in_run"] = pdc_in_run(stars[k], x)
    base = {k: {(x["t_start"], x["t_end"]) for x in real[k]} for k in kids}
    ctx = {"kids": kids, "stars": {k: stars[k] for k in kids}, "frac": frac, "base": base}
    ctx["real"] = real
    p = CACHE / f"ctx_{tag}.pkl"
    with open(p, "wb") as fh:
        pickle.dump(ctx, fh)
    return str(p)


def summarise(rows) -> dict:
    summ = {"w3": [], "controls": {}}
    for c, (rho, blend, te) in enumerate(CELLS):
        r = [x for x in rows if x["kind"] == "w3" and x["cell"] == c]
        k = sum(x["pass"] for x in r)
        summ["w3"].append(
            {
                "rho": rho,
                "blend": blend,
                "tE": te,
                "n": len(r),
                "flagged": sum(x["flagged"] for x in r),
                "recovered": k,
                "recovered_wilson95": wilson(k, len(r)),
                "recovered_no_shoulder": sum(x["pass_no_shoulder"] for x in r),
                "insufficient": sum(x["insufficient"] for x in r),
                "distinct_stars": len({x["kepid"] for x in r}),
            }
        )
    for kind in CONTROLS:
        r = [x for x in rows if x["kind"] == "ctl" and x["cell"] == kind]
        k = sum(x["pass"] for x in r)
        summ["controls"][kind] = {
            "n": len(r),
            "flagged": sum(x["flagged"] for x in r),
            "passed": k,
            "passed_wilson95": wilson(k, len(r)),
            "passed_no_shoulder": sum(x["pass_no_shoulder"] for x in r),
            "insufficient": sum(x["insufficient"] for x in r),
        }
    return summ


def gate_of(summ) -> dict:
    """Plan section 6: every cell recovered / n >= 0.50, every control family passed / n <= 0.10."""
    g = {
        "all_cells_ge_0.50": all(x["recovered"] / x["n"] >= GATE_REC for x in summ["w3"]),
        "all_controls_le_0.10": all(
            v["passed"] / v["n"] <= GATE_CTRL for v in summ["controls"].values()
        ),
    }
    g["pass"] = g["all_cells_ge_0.50"] and g["all_controls_le_0.10"]
    return g


def stage_inject(a) -> None:
    phase = a.phase
    frac_groups = ["dev"] if phase == "dev" else ["dev", "val"]
    log_usage(f"inject {phase}", set(frac_groups))
    ctx_path = build_ctx([phase], frac_groups, phase)
    with open(ctx_path, "rb") as fh:
        ctx = pickle.load(fh)
    real = {"label": "derived from observed light curves (uninjected)"}
    real["stars"] = {str(k): v for k, v in ctx["real"].items()}
    (OUT / f"real_{phase}.json").write_text(json.dumps(real, indent=1) + "\n")
    s_w3, s_ctl = SEEDS[phase]
    jobs = [("w3", c, j, s_w3 + c) for c in range(len(CELLS)) for j in range(N_INJ)]
    jobs += [("ctl", k, j, s_ctl + i) for i, k in enumerate(CONTROLS) for j in range(N_INJ)]
    t0 = time.time()
    with ProcessPoolExecutor(a.cores, initializer=_init, initargs=(ctx_path,)) as ex:
        rows = list(ex.map(_one, jobs, chunksize=8))
    el = time.time() - t0
    summ = summarise(rows)
    gate = gate_of(summ)
    out = {"label": "simulated", "phase": phase, "carriers": len(ctx["kids"])}
    out.update({"elapsed_s": round(el, 1), "summary": summ, "gate": gate, "rows": rows})
    (OUT / f"injections_{phase}.json").write_text(json.dumps(out) + "\n")
    if phase == "val":
        (OUT / "gate.json").write_text(json.dumps({"label": "derived", **gate}, indent=1) + "\n")
    print(json.dumps({"gate": gate, "elapsed_s": round(el, 1)}))


def stage_search(_a) -> None:
    gate = json.loads((OUT / "gate.json").read_text())
    if not gate["pass"]:
        raise SystemExit("validation gate not passed: the search group is not read")
    log_usage("search", {"dev", "val", "search"})
    with open(build_ctx(["search"], ["dev", "val", "search"], "search"), "rb") as fh:
        real = pickle.load(fh)["real"]
    order = [s["kepid"] for s in sample_stars(["search"])]
    surv = [
        {"kepid": k, "rank_in_search": i, **x}
        for i, k in enumerate(order)
        for x in real[k]
        if x["pass"]
    ]
    out = {"label": "derived from observed light curves", "survivors_in_order": surv}
    out["stars"] = {str(k): v for k, v in real.items()}
    (OUT / "search.json").write_text(json.dumps(out, indent=1) + "\n")
    n_runs = sum(len(v) for v in real.values())
    print(len(order), "stars;", n_runs, "flagged runs;", len(surv), "survivors")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("sample")
    p.add_argument("--pool", required=True)
    p.add_argument("--query", required=True)
    p.set_defaults(fn=stage_sample)
    p = sub.add_parser("fetch")
    p.add_argument("--group", choices=list(GROUPS), required=True)
    p.set_defaults(fn=stage_fetch)
    p = sub.add_parser("inject")
    p.add_argument("--phase", choices=["dev", "val"], required=True)
    p.add_argument("--cores", type=int, default=4)
    p.set_defaults(fn=stage_inject)
    p = sub.add_parser("search")
    p.set_defaults(fn=stage_search)
    a = ap.parse_args(argv)
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

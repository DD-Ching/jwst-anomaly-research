"""X3 G1: full W3 screen on unselected Gaia epoch photometry, calibrated on dev before held-out.

Unit file: docs/state/units/X3.md (G1). Builds on scripts/x3_gaps.py (G0: split, fetch, v0 fitter
flag, v1 astrometric veto). The G1 pipeline, applied identically to every light curve:

1. >= MIN_POINTS accepted G epochs (x3_gaps.MIN_POINTS).
2. v1 astrometric-quality veto (x3_gaps.v1_veto).
3. Pre-screen: >= PRE_N G epochs below the median by > PRE_SIG sigma (a dip must be sampled).
4. v0 fitter flag: dBIC(N1neg - best of PSPL/PAR) < w3_microlensing.P.flag_dbic.
5. Repeat-dip veto (no period search, so no alias or scan-law false periods): outside the fitted
   event window, >= REPEAT_N epochs as faint as half the fitted dip means a recurring dimmer.
6. Multi-band common-magnification test: BP and RP are fitted with the G-fitted A(t) and their own
   source and blend fluxes (F_b = fs_b A + fb_b), so blend colour is free. A band is "informative"
   when it has >= INFO_N epochs inside the event with |A - 1| >= 0.2 and fs_b is measured to better
   than 50 %. Verdict per light curve: "consistent" (every informative band fits, fs_b > 0),
   "inconsistent" (some informative band does not), or "insufficient" (no informative band).
   Only "inconsistent" vetoes.

Calibration sets (dev only; labels: light curves `observed`, injections `simulated`):
- dev real light curves already fetched in G0 (seeds 2026, 4242);
- W3 injections with independent blend fractions per band, PSPL injections, and dust-dimming
  injections (extinction ratios ASSUMED: A_BP/A_G = 1.27, A_RP/A_G = 0.75);
- selection scenarios on the same W3 injections: "floor" (epochs fainter than the per-band
  detection limit removed), "outlier" (isolated > 5 sigma epochs removed, as an aggressive stand-in
  for Gaia's outlier rejection);
- DR3 classified variables (ECL, YSO, LPV; `vari_classifier_result`), fetched by DataLink.

Subcommands: ``calib`` (BP/RP error scales on dev), ``refs`` (fetch reference classes),
``run --set {dev,inject,refs}``, ``summary``.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import io
import json
import math
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import exotic_sim as es
from jwst_anomaly import gaia_mulens as gm

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


x3 = _load("x3_gaps")
w3 = x3.w3
wg = x3.wg

# ---- frozen G1 screen parameters (docs/state/units/X3.md G1) ----
PRE_N, PRE_SIG = 2, 3.0
REPEAT_N, REPEAT_SIG = 2, 3.0
INFO_N = 3
CHI2_P_MIN = 1e-3
DET_LIMIT = {"G": 20.7, "BP": 20.9, "RP": 20.2}  # ASSUMPTION: per-transit detection limits (mag)
EXT_RATIO = {"G": 1.0, "BP": 1.27, "RP": 0.75}  # ASSUMPTION: A_band / A_G for dust
REF_CLASSES = ("ECL", "YSO", "LPV")
SEEDS_DEV = (2026, 4242)


def out_dir() -> Path:
    return x3.out_dir()


def results_dir() -> Path:
    return x3.results_dir()


# ----------------------------------------------------------------------------- data
def parse_bands(text: str) -> dict:
    """G (as in G0: magnitudes, rescaled errors) plus BP and RP fluxes (e-/s) with their errors."""
    rows = list(csv.DictReader(io.StringIO(text)))
    out = {}
    for band, tcol, fcol, ecol, rej in (
        ("BP", "bp_obs_time", "bp_flux", "bp_flux_error", "variability_flag_bp_reject"),
        ("RP", "rp_obs_time", "rp_flux", "rp_flux_error", "variability_flag_rp_reject"),
    ):
        keep = [
            r
            for r in rows
            if r[tcol] and r[fcol] and r[ecol] and float(r[ecol]) > 0 and r[rej] != "true"
        ]
        t = np.array([float(r[tcol]) for r in keep]) + gm.GAIA_T0
        f = np.array([float(r[fcol]) for r in keep])
        e = np.array([float(r[ecol]) for r in keep])
        out[band] = (t, f, e)
    return out


def load_source(sid: str, radec, survey) -> tuple[object, dict] | None:
    path = survey._lc_path(sid)
    if not path.exists():
        return None
    text = path.read_text()
    tab = gm.parse_epoch_csv(text)
    if len(tab) < x3.MIN_POINTS:
        return None
    lc = wg.gaia_light_curve(tab, sid, radec)
    return lc, parse_bands(text)


# ----------------------------------------------------------------------------- screen steps
def robust_sigma(x) -> float:
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def prescreen(lc) -> bool:
    med = np.median(lc.f)
    return int(np.sum(lc.f < med - PRE_SIG * lc.sf)) >= PRE_N


def event_window(p: dict) -> tuple[float, float]:
    u0 = min(abs(p["u0"]), 1.99)
    half = p["tE"] * max(math.sqrt(4.0 - u0 * u0), 1.0) * 1.5
    return p["t0"] - half, p["t0"] + half


def repeat_dip(lc, p: dict, fs: float, fb: float) -> bool:
    """True when the baseline outside the fitted event window holds a second dip."""
    a = w3.magnification("N1neg", lc, p)
    model = fs * a + fb
    base = fs + fb
    depth = base - np.min(model)
    lo, hi = event_window(p)
    out = (lc.t < lo) | (lc.t > hi)
    faint = lc.f < base - 0.5 * depth - REPEAT_SIG * lc.sf
    return int(np.sum(out & faint)) >= REPEAT_N


def band_test(t, f, e, p: dict, scale: float) -> dict:
    """Common-magnification fit of one band: F = fs A(t) + fb, linear in (fs, fb)."""
    if t.size < 5:
        return {"verdict": "insufficient", "n": int(t.size)}
    lc_b = w3.LightCurve(t, f, e * scale)
    a = w3.magnification("N1neg", lc_b, p)
    w = 1.0 / (e * scale) ** 2
    X = np.vstack([a, np.ones_like(a)]).T
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    fs, fb = cov @ (X.T @ (w * f))
    chi2 = float(np.sum(w * (f - fs * a - fb) ** 2))
    dof = t.size - 2
    n_in = int(np.sum(np.abs(a - 1.0) >= 0.2))
    sfs = math.sqrt(cov[0, 0])
    informative = n_in >= INFO_N and fs > 0 and sfs < 0.5 * abs(fs)
    from scipy.stats import chi2 as chi2_dist

    p_val = float(chi2_dist.sf(chi2, dof)) if dof > 0 else 0.0
    if not informative:
        verdict = "insufficient"
    elif p_val >= CHI2_P_MIN and fs > 0:
        verdict = "consistent"
    else:
        verdict = "inconsistent"
    return {
        "verdict": verdict,
        "n": int(t.size),
        "n_in": n_in,
        "fs": float(fs),
        "fb": float(fb),
        "sfs": sfs,
        "chi2_dof": chi2 / max(dof, 1),
        "p": p_val,
    }


def achromatic(bands: dict, p: dict, scales: dict) -> tuple[str, dict]:
    res = {b: band_test(*bands[b], p, scales.get(b, 1.0)) for b in ("BP", "RP")}
    v = [r["verdict"] for r in res.values()]
    if "inconsistent" in v:
        return "inconsistent", res
    if "consistent" in v:
        return "consistent", res
    return "insufficient", res


def screen_full(lc, bands: dict, vetoed: bool, scales: dict) -> dict:
    row = {"v1_vetoed": bool(vetoed), "pre": False, "fit_flag": False, "repeat": False}
    row.update({"achrom": "", "final": False, "dbic_N1neg": np.nan})
    if vetoed:
        return row
    row["pre"] = bool(prescreen(lc))
    if not row["pre"]:
        return row
    k = np.convolve(lc.f, np.ones(3) / 3, mode="same")
    t0g = float(lc.t[int(np.argmax(k))])
    models = ("PSPL", "PAR", "N1neg")
    res = w3.fit_event(lc, t0g, 30.0, 0.3, models=models, te_grid=(10.0, 100.0))
    best = min(res[m]["bic"] for m in w3.ORDINARY if m in res)
    row["dbic_N1neg"] = float(res["N1neg"]["bic"] - best)
    row["fit_flag"] = bool(row["dbic_N1neg"] < w3.P.flag_dbic)
    if not row["fit_flag"]:
        return row
    p = {k_: res["N1neg"][k_] for k_ in ("t0", "tE", "u0", "rho")}
    row.update({f"fit_{k_}": float(v) for k_, v in p.items()})
    row["repeat"] = bool(repeat_dip(lc, p, res["N1neg"]["fs"], res["N1neg"]["fb"]))
    if row["repeat"]:
        return row
    verdict, br = achromatic(bands, p, scales)
    row["achrom"] = verdict
    for b, r in br.items():
        row[f"{b}_verdict"] = r["verdict"]
        row[f"{b}_n_in"] = r.get("n_in", 0)
    row["final"] = verdict != "inconsistent"
    return row


# ----------------------------------------------------------------------------- injections
def inject_band(t, f, e, frac, prm, n=1.0, sign=-1, rho=None):
    base = float(np.median(f))
    rho_ = prm["rho"] if rho is None else rho
    inj = es.inject_light_curve(t, base, prm["t0"], prm["tE"], prm["u0"], n, sign, rho_, frac)
    fnew = np.asarray(inj["flux"], float)
    scale = np.sqrt(np.clip(fnew / base, 0.09, None))
    return fnew + (f - base) * scale, e * scale


def inject_dust(t, f, e, prm, band):
    base = float(np.median(f))
    dmag = prm["depth_g"] * EXT_RATIO[band] * np.exp(-0.5 * ((t - prm["t0"]) / prm["width"]) ** 2)
    fac = 10.0 ** (-0.4 * dmag)
    scale = np.sqrt(np.clip(fac, 0.09, None))
    return base * fac + (f - base) * scale, e * scale


def mag_of_flux(f, ref_flux, ref_mag):
    return ref_mag - 2.5 * np.log10(np.clip(f / ref_flux, 1e-9, None))


def apply_selection(kind: str, lc, bands, g_mean: float, bands_base):
    """Remove epochs as a selection scenario would (``floor`` or ``outlier``)."""
    if kind == "base":
        return lc, bands
    if kind == "floor":
        gmag = mag_of_flux(lc.f, np.median(lc.f), g_mean)
        lc2 = lc.subset(gmag <= DET_LIMIT["G"])
        b2 = {}
        for b, (t, f, e) in bands.items():
            f0 = np.median(bands_base[b][1]) if bands_base[b][1].size else 1.0
            m = mag_of_flux(f, f0, g_mean)  # same zero point as G: ASSUMPTION (colour ignored)
            keep = m <= DET_LIMIT[b]
            b2[b] = (t[keep], f[keep], e[keep])
        return lc2, b2
    if kind == "outlier":
        f = lc.f
        keep = np.ones(f.size, bool)
        for i in range(1, f.size - 1):
            nb = 0.5 * (f[i - 1] + f[i + 1])
            if abs(f[i] - nb) > 5.0 * lc.sf[i] and np.sign(f[i - 1] - f[i]) == np.sign(
                f[i + 1] - f[i]
            ):
                keep[i] = False
        return lc.subset(keep), bands
    raise ValueError(kind)


# ----------------------------------------------------------------------------- jobs
def _job(args):
    kind, sid, lc, bands, vetoed, scales, extra = args
    row = {"kind": kind, "event_id": sid, **{k: v for k, v in extra.items() if np.isscalar(v)}}
    try:
        row.update(screen_full(lc, bands, vetoed, scales))
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001
        row.update({"final": False, "error": repr(exc)[:200]})
    return row


def dev_sources() -> list[tuple[str, tuple, float]]:
    out = []
    for seed in SEEDS_DEV:
        tab = Table.read(out_dir() / f"g0_seed{seed}.ecsv")
        real = tab[tab["kind"] == "real"]
        ids = {str(e): float(g) for e, g in zip(real["event_id"], real["G_mean"], strict=True)}
        dev = x3.ids_table(19.0)
        pos = {
            str(s): (float(r), float(d))
            for s, r, d in zip(dev["source_id"], dev["ra"], dev["dec"], strict=True)
        }
        out.extend((s, pos[s], g) for s, g in ids.items())
    return out


def dev_veto() -> dict[str, bool]:
    veto = {}
    for seed in SEEDS_DEV:
        veto.update(x3.v1_veto(Table.read(out_dir() / f"attrs_seed{seed}.csv", format="ascii.csv")))
    return veto


def calib_scales() -> dict:
    path = results_dir() / "g1_bp_rp_error_scales.json"
    return json.loads(path.read_text())["scales"] if path.exists() else {"BP": 1.0, "RP": 1.0}


def run_calib() -> dict:
    """BP/RP error scale: median sqrt(chi2/dof) of a constant fit on dev light curves that are not
    vetoed and not pre-screened (quiet), per band; G uses D-061's published rescaling."""
    survey = x3.GapsDev()
    veto = dev_veto()
    vals = {"BP": [], "RP": []}
    for sid, radec, _g in dev_sources():
        if veto.get(sid, True):
            continue
        src = load_source(sid, radec, survey)
        if src is None or prescreen(src[0]):
            continue
        for b in ("BP", "RP"):
            t, f, e = src[1][b]
            if f.size >= 10:
                w = 1.0 / e**2
                mu = np.sum(w * f) / np.sum(w)
                vals[b].append(math.sqrt(np.sum(w * (f - mu) ** 2) / (f.size - 1)))
    rec = {
        "scales": {b: float(np.median(v)) for b, v in vals.items()},
        "n_sources": {b: len(v) for b, v in vals.items()},
        "provenance": "derived from observed dev light curves (quiet, not vetoed)",
    }
    (results_dir() / "g1_bp_rp_error_scales.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def build_jobs(which: str, seed: int) -> list:
    survey = x3.GapsDev()
    scales = calib_scales()
    rng = np.random.default_rng(seed)
    jobs = []
    if which in ("dev", "inject"):
        veto = dev_veto()
        for sid, radec, gmean in dev_sources():
            src = load_source(sid, radec, survey)
            if src is None:
                continue
            lc, bands = src
            vet = veto.get(sid, True)
            if which == "dev":
                jobs.append(("real", sid, lc, bands, vet, scales, {"G_mean": gmean}))
                continue
            # W3 with independent per-band source fractions (blend colour free)
            prm = {
                "tE": float(rng.choice(x3.TE_CELLS)),
                "rho": float(rng.choice(x3.RHO_CELLS)),
                "u0": float(rng.uniform(0.0, 2.0)),
                "t0": float(rng.uniform(lc.t.min(), lc.t.max())),
                "G_mean": gmean,
            }
            fr = {b: float(rng.uniform(0.5, 1.0)) for b in ("G", "BP", "RP")}
            prm.update({f"fs_{b}": v for b, v in fr.items()})
            base_g = float(np.median(lc.f))
            lcw = w3.inject_w3(
                lc, base_g, fr["G"] * base_g, prm["t0"], prm["tE"], prm["u0"], prm["rho"]
            )
            bw = {}
            for b in ("BP", "RP"):
                t, f, e = bands[b]
                bw[b] = (t, *inject_band(t, f, e, fr[b], prm)) if f.size else (t, f, e)
            for sel in ("base", "floor", "outlier"):
                lcs, bs = apply_selection(sel, lcw, bw, gmean, bands)
                jobs.append((f"W3_{sel}", sid, lcs, bs, vet, scales, prm))
            # PSPL control
            pp = {**prm, "u0": float(rng.uniform(0.0, 1.0)), "rho": 0.0}
            lcp = w3.inject_w3(
                lc, base_g, fr["G"] * base_g, pp["t0"], pp["tE"], pp["u0"], 0.0, 1.0, 1
            )
            bp_ = {}
            for b in ("BP", "RP"):
                t, f, e = bands[b]
                bp_[b] = (t, *inject_band(t, f, e, fr[b], pp, 1.0, 1, 0.0)) if f.size else (t, f, e)
            jobs.append(("PSPL", sid, lcp, bp_, vet, scales, pp))
            # dust-dimming negative control (Gaussian in time, chromatic)
            dp = {
                "t0": float(rng.uniform(lc.t.min(), lc.t.max())),
                "width": float(rng.choice((10.0, 30.0, 100.0))),
                "depth_g": float(rng.uniform(0.3, 1.5)),
                "G_mean": gmean,
            }
            fg, eg = inject_dust(lc.t, lc.f, lc.sf, dp, "G")
            lcd = w3.LightCurve(lc.t, fg, eg, lc.ra, lc.dec, lc.event_id)
            bd = {}
            for b in ("BP", "RP"):
                t, f, e = bands[b]
                bd[b] = (t, *inject_dust(t, f, e, dp, b)) if f.size else (t, f, e)
            jobs.append(("dust", sid, lcd, bd, vet, scales, dp))
    elif which == "refs":
        rsv = RefSurvey()
        for cls in REF_CLASSES:
            ids = Table.read(out_dir() / f"ref_{cls}.csv", format="ascii.csv")
            for r in ids:
                sid = str(r["source_id"])
                src = load_source(sid, (float(r["ra"]), float(r["dec"])), rsv)
                if src is None:
                    continue
                vet = ref_veto(r)
                jobs.append(
                    (
                        f"ref_{cls}",
                        sid,
                        src[0],
                        src[1],
                        vet,
                        scales,
                        {"G_mean": float(r["phot_g_mean_mag"])},
                    )
                )
    return jobs


def write_rows(rows: list[dict], path: Path, meta: dict) -> Table:
    """Rows of screen results -> ECSV with concrete column dtypes."""
    keys = sorted({k for r in rows for k in r})
    tab = Table()
    for k in keys:
        vals = [r.get(k) for r in rows]
        if all(isinstance(v, (bool, np.bool_)) or v is None for v in vals):
            tab[k] = np.array([bool(v) if v is not None else False for v in vals])
        elif all(isinstance(v, (int, float, np.floating, np.integer)) or v is None for v in vals):
            tab[k] = np.array([float(v) if v is not None else np.nan for v in vals])
        else:
            tab[k] = np.array(["-" if v in (None, "") else str(v) for v in vals])
    tab.meta.update(meta)
    tab.write(path, overwrite=True)
    return tab


def run_set(which: str, seed: int, procs: int) -> Path:
    jobs = build_jobs(which, seed)
    t1 = time.time()
    with Pool(procs) as pool:
        rows = list(pool.imap_unordered(_job, jobs, 4))
    meta = {
        "source": f"scripts/x3_g1.py run --set {which}",
        "provenance": "derived (observed light curves; simulated where kind is not real/ref)",
        "seed": seed,
        "wall_time_s": round(time.time() - t1, 1),
    }
    path = out_dir() / f"g1_{which}.ecsv"
    tab = write_rows(rows, path, meta)
    w3.write_ecsv_gz(tab, results_dir() / f"g1_{which}.ecsv.gz")
    print(f"wrote {path}: {len(tab)} rows in {time.time() - t1:.0f} s")
    return path


# ----------------------------------------------------------------------------- references
# References are drawn from the same sky (GAPS level-5 pixels, same scanning law) and from the dev
# digit only, so no held-out light curve is read. One query per pixel's source_id range (indexed);
# a global query on vari_classifier_result timed out (HTTP 408).
REF_QUERY = (
    "SELECT v.source_id, v.best_class_name, g.ra, g.dec, g.phot_g_mean_mag, g.ruwe, "
    "g.ipd_frac_multi_peak FROM gaiadr3.vari_classifier_result AS v JOIN gaiadr3.gaia_source AS g "
    "ON v.source_id = g.source_id WHERE v.source_id BETWEEN {lo} AND {hi} "
    "AND MOD(v.source_id, 10) = 0 AND v.best_class_name IN ('ECL', 'YSO', 'LPV') "
    "AND v.best_class_score > 0.5 AND g.phot_g_mean_mag BETWEEN 13 AND 19 "
    "AND g.has_epoch_photometry = 'True'"
)


class RefSurvey(gm.GaiaDR3Microlensing):
    name = "gaia-dr3-x3-refs"

    def __init__(self):
        d = x3.paths.data_root() / "raw" / "gaia_dr3_x3_refs"
        d.mkdir(parents=True, exist_ok=True)
        super().__init__(cache=d)


def ref_veto(r) -> bool:
    ruwe = (
        float(r["ruwe"])
        if r["ruwe"] not in (None, "", "--") and not np.ma.is_masked(r["ruwe"])
        else np.nan
    )
    ipd = (
        float(r["ipd_frac_multi_peak"]) if not np.ma.is_masked(r["ipd_frac_multi_peak"]) else np.nan
    )
    return bool(np.isnan(ruwe) or ruwe > 1.4 or ipd > 10)


# The M31 field holds few YSOs and LPVs in the dev digit (1 and 36 at G 13-19), so a second
# field around Orion (star-forming, YSO-rich) is added. Its scanning cadence differs from M31's:
# mimic rates from it are indicative, not a calibration of the GAPS cadence.
REF_FIELDS = {"M31": (x3.M31, x3.GAPS_RADIUS_DEG), "Orion": ((83.8, -5.4), 5.0)}


def _cone_pixels(center, radius) -> list[int]:
    import astropy.units as u
    from astropy.coordinates import Latitude, Longitude
    from astropy_healpix import HEALPix

    hp = HEALPix(nside=2**x3.HPX_LEVEL, order="nested")
    pix = hp.cone_search_lonlat(
        Longitude(center[0] * u.deg), Latitude(center[1] * u.deg), radius * u.deg
    )
    return sorted(int(p) for p in pix)


def ref_table() -> Table:
    path = out_dir() / "ref_candidates_fields.csv"
    if not path.exists():
        shift = 35 + 2 * (12 - x3.HPX_LEVEL)
        parts, header = [], None
        for field, (center, radius) in REF_FIELDS.items():
            for p in _cone_pixels(center, radius):
                q = REF_QUERY.format(lo=p << shift, hi=((p + 1) << shift) - 1)
                lines = gm._tap_csv(q).strip().splitlines()
                header = lines[0] + ",field"
                parts.extend(f"{ln},{field}" for ln in lines[1:])
        gm._write_atomic(path, ("\n".join([header, *parts]) + "\n").encode())
    return Table.read(path, format="ascii.csv")


def run_refs(n_per: int, seed: int = 11) -> None:
    allr = ref_table()
    rng = np.random.default_rng(seed)
    for cls in REF_CLASSES:
        sub = allr[np.array([str(c) == cls for c in allr["best_class_name"]])]
        if len(sub) > n_per:
            sub = sub[np.sort(rng.choice(len(sub), n_per, replace=False))]
        sub.write(out_dir() / f"ref_{cls}.csv", format="ascii.csv", overwrite=True)
        n = RefSurvey().prefetch([str(s) for s in sub["source_id"]])
        n_all = int(np.sum(allr["best_class_name"] == cls))
        print(f"{cls}: {len(sub)} of {n_all} listed, {n} fetched")


# ----------------------------------------------------------------------------- selection bound
# Injections keep the host's catalogue RUWE, so they cannot show that a real event would not itself
# be vetoed. Bound for one mechanism: a blended source dims, so the photocentre of source + blend
# (separation d) moves toward the blend; the extra AL scatter inflates RUWE.
#   x_i = d * b_i, b_i = f_b / (f_s A_i + f_b);  AL projection: variance factor 1/2 (random scans).
#   RUWE_new^2 = RUWE_0^2 + var(x_AL) / sigma_obs^2,  sigma_obs = pos_error * sqrt(N_obs_AL / 2).
# ASSUMPTIONS: flux-weighted photocentre (an upper bound on the shift once the pair is resolved),
# d log-uniform in [10, 1000] mas, f_s as injected, sigma_obs from each host's catalogue errors and
# scaled by 0.5 and 2 for the sensitivity range. Not modelled (stated, not bounded): astrometric
# microlensing of the source image itself (centroid shifts ~ theta_E), Gaia's on-board window
# allocation in crowded fields, and any reprocessing of DR3 photometry.
SEL_QUERY = (
    "SELECT source_id, ra_error, dec_error, astrometric_n_good_obs_al FROM gaiadr3.gaia_source "
    "WHERE source_id IN ({ids})"
)
D_RANGE_MAS = (10.0, 1000.0)
SIGMA_SCALES = (0.5, 1.0, 2.0)


def _sel_attrs(ids: list[str]) -> dict:
    path = out_dir() / "g1_sel_attrs.csv"
    if not path.exists():
        lines, header = [], None
        for k in range(0, len(ids), 500):
            rows = gm._tap_csv(SEL_QUERY.format(ids=",".join(ids[k : k + 500]))).strip()
            rows = rows.splitlines()
            header = rows[0]
            lines.extend(rows[1:])
        gm._write_atomic(path, ("\n".join([header, *lines]) + "\n").encode())
    t = Table.read(path, format="ascii.csv")
    out = {}
    for r in t:
        if np.ma.is_masked(r["ra_error"]) or np.ma.is_masked(r["astrometric_n_good_obs_al"]):
            continue
        pos = max(float(r["ra_error"]), float(r["dec_error"]))
        out[str(r["source_id"])] = pos * math.sqrt(float(r["astrometric_n_good_obs_al"]) / 2.0)
    return out


def selection_bound(seed: int = 5, n_d: int = 50) -> dict:
    tab = Table.read(out_dir() / "g1_inject.ecsv")
    kind = _strcol(tab, "kind")
    m = (kind == "W3_base") & np.asarray(tab["final"], bool)
    sub = tab[m]
    ids = sorted({str(e) for e in sub["event_id"]})
    sig = _sel_attrs(ids)
    ruwe0 = {}
    for sd in SEEDS_DEV:
        at = Table.read(out_dir() / f"attrs_seed{sd}.csv", format="ascii.csv")
        for s_, r_ in zip(at["source_id"], at["ruwe"], strict=True):
            if not np.ma.is_masked(r_):
                ruwe0[str(s_)] = float(r_)
    survey = x3.GapsDev()
    pos = dict((s_, (r_, d_)) for s_, r_, d_ in dev_sources_pos())
    rng = np.random.default_rng(seed)
    logd = rng.uniform(*np.log10(D_RANGE_MAS), size=n_d)
    d = 10.0**logd
    kept = {sc: [] for sc in SIGMA_SCALES}
    kept_band = {sc: {"10-100": [], "100-1000": []} for sc in SIGMA_SCALES}
    dmean_g, over19 = [], 0
    n_used = 0
    for r in sub:
        sid = str(r["event_id"])
        if sid not in sig or sid not in ruwe0:
            continue
        src = load_source(sid, pos[sid], survey)
        if src is None:
            continue
        lc = src[0]
        p = {k: float(r[k]) for k in ("t0", "tE", "u0", "rho")}
        a = w3.magnification("N1neg", lc, p)
        fs = float(r["fs_G"])
        fb = 1.0 - fs
        b = fb / (fs * a + fb)
        var_al = 0.5 * np.var(b) * d**2  # mas^2, one value per separation draw
        for sc in SIGMA_SCALES:
            ruwe = np.sqrt(ruwe0[sid] ** 2 + var_al / (sc * sig[sid]) ** 2)
            ok = ruwe <= 1.4
            kept[sc].append(float(ok.mean()))
            kept_band[sc]["10-100"].append(float(ok[d < 100].mean()) if (d < 100).any() else 1.0)
            kept_band[sc]["100-1000"].append(
                float(ok[d >= 100].mean()) if (d >= 100).any() else 1.0
            )
        dg = -2.5 * math.log10(float(np.mean(fs * a + fb)))
        dmean_g.append(dg)
        over19 += int(float(r["G_mean"]) <= 19.0 < float(r["G_mean"]) + dg)
        n_used += 1
    rec = {
        "label": "derived bound from simulated injections into observed dev light curves",
        "assumptions": {
            "photocentre": "flux-weighted (upper bound on shift)",
            "blend_separation_mas": list(D_RANGE_MAS),
            "separation_prior": "log-uniform",
            "sigma_obs": "max(ra_error, dec_error) * sqrt(astrometric_n_good_obs_al / 2)",
            "sigma_scales": list(SIGMA_SCALES),
        },
        "n_injections_passing_final": int(m.sum()),
        "n_used": n_used,
        "retained_fraction_mean": {f"{sc:g}": float(np.mean(v)) for sc, v in kept.items()},
        "retained_by_separation": {
            f"{sc:g}": {k: float(np.mean(v)) for k, v in bd.items()} for sc, bd in kept_band.items()
        },
        "mean_G_shift_mag": {
            "median": float(np.median(dmean_g)),
            "p95": float(np.percentile(dmean_g, 95)),
        },
        "n_crossing_G19_cut": over19,
        "not_bounded": [
            "astrometric microlensing of the source image (centroid shift ~ theta_E)",
            "on-board window allocation and dead time in crowded fields",
            "DR3 photometric reprocessing beyond the floor/outlier scenarios",
        ],
    }
    (results_dir() / "g1_selection_bound.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def dev_sources_pos() -> list[tuple[str, float, float]]:
    return [(s, r, d) for s, (r, d), _g in dev_sources()]


# ----------------------------------------------------------------------------- summary
def wilson(k: int, n: int, z: float = 1.96) -> list:
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [max(0.0, c - h), min(1.0, c + h)]


def _strcol(tab, col) -> np.ndarray:
    """String column with masked (empty on read-back) entries as ""."""
    if col not in tab.colnames:
        return np.array([""] * len(tab))
    c = tab[col]
    c = c.filled("") if hasattr(c, "filled") else c
    return np.array(["" if str(x) in ("-", "--") else str(x) for x in c])


def summary() -> dict:
    out = {
        "frozen_parameters": {
            "PRE_N": PRE_N,
            "PRE_SIG": PRE_SIG,
            "REPEAT_N": REPEAT_N,
            "REPEAT_SIG": REPEAT_SIG,
            "INFO_N": INFO_N,
            "CHI2_P_MIN": CHI2_P_MIN,
            "flag_dbic": w3.P.flag_dbic,
        }
    }
    for which in ("dev", "inject", "refs"):
        path = out_dir() / f"g1_{which}.ecsv"
        if not path.exists():
            continue
        tab = Table.read(path)
        err = _strcol(tab, "error")
        for kind in sorted(set(str(k) for k in tab["kind"])):
            m = (np.array([str(k) for k in tab["kind"]]) == kind) & (err == "")
            n = int(m.sum())
            stages = {}
            for col in ("v1_vetoed", "pre", "fit_flag", "repeat", "final"):
                if col in tab.colnames:
                    stages[col] = int(np.sum(np.asarray(tab[col], bool)[m]))
            k_final = stages.get("final", 0)
            ach = list(_strcol(tab, "achrom")[m & np.asarray(tab["final"], bool)])
            rec = {
                "n": n,
                "n_errors": int(
                    np.sum((np.array([str(k) for k in tab["kind"]]) == kind) & (err != ""))
                ),
                "stage_counts": stages,
                "final_pass": k_final,
                "final_rate": k_final / n if n else None,
                "final_rate_wilson95": wilson(k_final, n),
                "final_by_colour": {v: ach.count(v) for v in ("consistent", "insufficient")},
            }
            if kind.startswith("W3"):
                te = np.asarray(tab["tE"], float)[m]
                fin = np.asarray(tab["final"], bool)[m]
                rec["eff_by_tE"] = {
                    f"{t:g}": [int(fin[te == t].sum()), int((te == t).sum())] for t in x3.TE_CELLS
                }
                sel = (te >= 30) & (te <= 300)
                rec["eff_te30_300"] = [int(fin[sel].sum()), int(sel.sum())]
                rec["eff_te30_300_wilson95"] = wilson(int(fin[sel].sum()), int(sel.sum()))
            out[f"{which}:{kind}"] = rec
    real = out.get("dev:real")
    if real:
        cnt = json.loads((results_dir() / "sample_counts_g19.json").read_text())
        n_ho = cnt["n_heldout"]
        lo, hi = real["final_rate_wilson95"]
        out["expected_heldout_flags"] = {
            "n_heldout": n_ho,
            "point": real["final_rate"] * n_ho,
            "wilson95": [lo * n_ho, hi * n_ho],
        }
    path = results_dir() / "g1_summary.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=("calib", "refs", "run", "selection", "summary"))
    ap.add_argument("--set", dest="which", choices=("dev", "inject", "refs"), default="dev")
    ap.add_argument("--seed", type=int, default=77)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--n-ref", type=int, default=200)
    a = ap.parse_args(argv)
    if a.cmd == "calib":
        print(json.dumps(run_calib(), indent=1))
    elif a.cmd == "refs":
        run_refs(a.n_ref)
    elif a.cmd == "run":
        run_set(a.which, a.seed, a.procs)
    elif a.cmd == "selection":
        print(json.dumps(selection_bound(), indent=1))
    else:
        print(json.dumps(summary(), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

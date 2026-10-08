"""W3 in the Gaia DR3 microlensing candidates: can the published selection see a W3 event? (D-061)

Hypothesis: the Gaia DR3 candidates (``gaia_mulens``) were selected before any PSPL-shaped cut that
excludes W3 events, so fitting them could limit W3 where the Mróz OGLE samples cannot (D-057).
Test first, fit second: inject W3 events (n = 1, ε < 0, as in D-057) into real Gaia cadences and
pass them through an emulation of the published Sample A selection (Wyrzykowski et al. 2023,
Appendix C). Only if a useful fraction passes is a W3 fit of the candidates informative.

Subcommands (outputs under ``$JWST_ANOMALY_DATA/derived/w3_gaia/``; the fit table is also written
gzipped to ``results/w3_gaia/``):

- ``fetch``: event table, positions and epoch photometry of all 363 candidates (DataLink, cached).
- ``fit``: ordinary (PSPL, FSPL, PAR) and exotic (N1neg, E2pos, E2neg) fits of every G light curve
  with the D-057 fitter (``w3_microlensing``), plus the published-style Level 0 PSPL (no blend), and
  the emulated Sample A selection of each real candidate (the audit).
- ``inject``: W3 events and PSPL controls on PSPL-subtracted Gaia light curves; fitter flags and the
  emulated selection.

G errors are rescaled as in the paper (Eq. 9–10) before every fit. Exotic physics is a hypothesis;
a better exotic fit is an anomaly to vet, never a discovery.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import minimize
from scipy.stats import skew

from jwst_anomaly import gaia_mulens as gm
from jwst_anomaly import paths, schema

_DIR = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("w3_microlensing", _DIR / "w3_microlensing.py")
w3 = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("w3_microlensing", w3)
_spec.loader.exec_module(w3)


@dataclass(frozen=True)
class SelParams:
    """Sample A cuts of Wyrzykowski et al. 2023 (Table C.1, membership score in Appendix C) are
    marked "published"; every other value is an ASSUMPTION of this emulation."""

    score_min: float = 0.85  # published
    u0_max: float = 1.0  # published, Level 0
    te_min: float = 20.0  # published
    te_max: float = 300.0  # published
    g0_max: float = 19.0  # published, Level 0 baseline
    amp_bright: float = 0.25  # published: max − min (mag) when median < 17
    amp_faint: float = 0.5  # published: when median > 17
    n_event_min: int = 3  # published: Extractor points in the event > 3
    duration_min: float = 135.0  # published: days between the Extractor's outlying points
    max_sigma_min: float = 50.0  # published: Extractor "strength"
    par_chi2_dof_max: float = 1.9  # published: Level 2 (parallax) χ²/dof
    par_pi_max: float = 3.5  # published: |π_E,N|, |π_E,E|
    # published: ``paczynski0_tmax`` window of the stricter skew–Abbe cut. The archive's tmax is
    # BJD − 2455197.5 (≈ 1400–2830 for these events), where this window would never apply; read as
    # JD − 2450000 it is 2014-06-19 … 2015-06-19, the first mission year (ASSUMPTION).
    abbe_window: tuple = (6824.5, 7189.5)
    abbe_window_zero: float = 2450000.0  # ASSUMPTION: zero point of ``abbe_window`` (see above)
    # ASSUMPTIONS (the Extractor's definitions are not published):
    event_nsigma: float = 3.0  # an "outlying point" is brighter than the median by > 3σ
    abbe_log: str = "log10"  # the paper writes "log"; ``ln`` is also reported


S = SelParams()
INJ_TE = (10.0, 30.0, 100.0, 300.0)  # days (ASSUMPTION: grid; < 20 d fails the published cut)
INJ_RHO = (0.01, 0.1)


def out_dir() -> Path:
    d = paths.data_root() / "derived" / "w3_gaia"
    d.mkdir(parents=True, exist_ok=True)
    return d


def results_dir() -> Path:
    d = paths.repo_root() / "results" / "w3_gaia"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ----------------------------------------------------------------------------- light curves


def rescale_g_errors(mag, err) -> np.ndarray:
    """Wyrzykowski et al. 2023 Eq. 9–10: σ² + (√30 · 10^(0.17 max(G, 13.5) − 5.1))²."""
    g = np.maximum(np.asarray(mag, float), 13.5)
    s_exp = math.sqrt(30.0) * 10.0 ** (0.17 * g - 5.1)
    return np.sqrt(np.asarray(err, float) ** 2 + s_exp**2)


def gaia_light_curve(tab: Table, event_id: str, radec=None) -> w3.LightCurve:
    ra, dec = radec if radec is not None else (None, None)
    err = rescale_g_errors(tab["mag"], tab["mag_err"])
    return w3.LightCurve.from_mag(tab["time"], tab["mag"], err, ra=ra, dec=dec, event_id=event_id)


# ----------------------------------------------------------------------------- selection


def level0_fit(lc: w3.LightCurve, t0: float, te: float, u0: float) -> dict:
    """PSPL with no blending (F = F_s·A), the published Level 0 model; F_s from least squares."""

    def chi2(x):
        if not (-1.0 <= x[1] <= 3.7):
            return 1e30
        a = w3.pspl(w3.straight_beta(lc.t, x[0], 10.0 ** x[1], abs(x[2])))
        fs = np.sum(lc.w * a * lc.f) / np.sum(lc.w * a * a)
        return float(np.sum(lc.w * (lc.f - fs * a) ** 2))

    best = None
    for te0 in (te, 30.0, 100.0):
        for u in (u0, 0.3, 0.8):
            r = minimize(
                chi2,
                [t0, math.log10(max(te0, 0.2)), u],
                method="Nelder-Mead",
                options={"maxfev": 1500, "xatol": 1e-5, "fatol": 1e-4},
            )
            if best is None or r.fun < best.fun:
                best = r
    t0b, lte, ub = best.x
    a = w3.pspl(w3.straight_beta(lc.t, t0b, 10.0**lte, abs(ub)))
    fs = float(np.sum(lc.w * a * lc.f) / np.sum(lc.w * a * a))
    return {
        "t0": float(t0b),
        "tE": float(10.0**lte),
        "u0": float(abs(ub)),
        "g0": float(w3.flux_to_mag(fs)) if fs > 0 else np.inf,
        "chi2_dof": float(best.fun / max(lc.t.size - 4, 1)),
    }


def abbe(x) -> float:
    """Abbe value Σ(x_{i+1} − x_i)² / (2 Σ(x_i − x̄)²) · n / (n − 1)."""
    x = np.asarray(x, float)
    n = x.size
    den = np.sum((x - x.mean()) ** 2)
    if n < 3 or den <= 0:
        return np.nan
    return float(np.sum(np.diff(x) ** 2) / (2.0 * den) * n / (n - 1))


def membership_score(amp: float, chi2_dof: float, l0: dict, t_lo: float, t_hi: float, par=None):
    """Appendix C score; error bars are ignored (ASSUMPTION: point estimates)."""
    s = 1.0
    if amp < 0.25:
        s *= 0.6
    elif amp < 0.5:
        s *= 0.9
    if chi2_dof > 10:
        s *= 0.5
    elif chi2_dof >= 3:
        s *= 0.9
    if not t_lo <= l0["t0"] <= t_hi:
        s *= 0.9
    if l0["u0"] > 1.5:
        s *= 0.7
    if l0["u0"] < 0.001:
        s *= 0.9
    if l0["g0"] > 20.0:  # "baseline below 20 mag" read as fainter than 20 (ASSUMPTION)
        s *= 0.8
    if not 2.0 <= l0["tE"] <= 500.0:
        s *= 0.2
    if par is not None:
        big = (abs(par["pi_E_N"]) > 1) + (abs(par["pi_E_E"]) > 1)
        s *= {0: 1.0, 1: 0.95, 2: 0.9}[int(big)]
    return s


# Extractor cuts: computed with ASSUMED definitions but left out of ``selected``, because they fail
# 126 of the 163 real Sample A events (the audit, D-061); ``selected_ext`` includes them.
EXTRACTOR = ("n_event", "duration", "max_sigma")


def sample_a_selection(lc: w3.LightCurve, l0: dict, par: dict | None) -> dict:
    """Emulated Sample A cuts on one G light curve. Not emulated (no G-only equivalent): the BP−RP
    colour cut and the RP-robustness cut (injections keep the base star's colour: pass), the
    u0/te error cuts (pass), the visual inspection, and (in ``selected``) the Extractor cuts.
    ``par`` is this fitter's parallax fit with blending (the published Level 2 has none:
    ASSUMPTION)."""
    mag = w3.flux_to_mag(np.maximum(lc.f, 1e-6))
    smag = lc.sf / lc.f * 2.5 / math.log(10.0)
    med = float(np.median(mag))
    amp = float(mag.max() - mag.min())  # published "max − min"; a W3 dimming counts too
    sk = float(skew(mag))
    nsig = (med - mag) / smag
    ev = nsig > S.event_nsigma
    duration = float(np.ptp(lc.t[ev])) if ev.sum() >= 2 else 0.0
    wtime = l0["t0"] - S.abbe_window_zero
    ab = abbe(mag)
    lin = math.log10 if S.abbe_log == "log10" else math.log
    slope, icpt = (1.2, -0.84) if S.abbe_window[0] <= wtime <= S.abbe_window[1] else (0.8, -0.6)

    def abbe_ok(log):
        return bool(sk < 0 and ab > 0 and log(ab) < slope * log(-sk) + icpt)

    par_chi = par["chi2"] / par["dof"] if par else l0["chi2_dof"]
    score = membership_score(amp, l0["chi2_dof"], l0, lc.t.min(), lc.t.max(), par)
    passed = {
        "score": score >= S.score_min,
        "u0": l0["u0"] < S.u0_max,
        "te": S.te_min < l0["tE"] < S.te_max,
        "g0": l0["g0"] < S.g0_max,
        "skew": sk < 0,
        "amp": amp > (S.amp_bright if med < 17 else S.amp_faint),
        "t_first": l0["t0"] - l0["tE"] > lc.t.min(),
        "n_event": int(ev.sum()) > S.n_event_min,
        "duration": duration > S.duration_min,
        "max_sigma": float(nsig.max()) > S.max_sigma_min,
        "par_chi2": par_chi < S.par_chi2_dof_max,
        "par_pi": par is None
        or (abs(par["pi_E_N"]) < S.par_pi_max and abs(par["pi_E_E"]) < S.par_pi_max),
        "abbe": abbe_ok(lin),
    }
    return {
        "passed": passed,
        "selected": all(v for k, v in passed.items() if k not in EXTRACTOR),
        "selected_ln": all(
            v for k, v in {**passed, "abbe": abbe_ok(math.log)}.items() if k not in EXTRACTOR
        ),
        "selected_ext": all(passed.values()),
        "score": score,
        "skew": sk,
        "abbe": ab,
        "amp": amp,
        "n_event": int(ev.sum()),
        "duration": duration,
        "max_sigma": float(nsig.max()),
    }


# ----------------------------------------------------------------------------- fit


SEL_STATS = ("score", "skew", "abbe", "amp", "n_event", "duration", "max_sigma")


def _fails(tab: Table) -> dict[str, int]:
    """Number of rows failing each emulated cut (``sel_*`` columns)."""
    return {c[4:]: int(np.sum(~tab[c].astype(bool))) for c in tab.colnames if c.startswith("sel_")}


no_error = w3.no_error


def _fit_worker(job):
    ev, tab, radec = job
    t1 = time.time()
    row = {"event_id": ev["event_id"], "method": ev["method"]}
    try:
        if tab is None:
            raise KeyError(f"no epoch photometry for {ev['event_id']}")
        lc = gaia_light_curve(tab, ev["event_id"], radec)
        row["n_points"] = int(lc.t.size)
        res = w3.fit_event(lc, ev["t0_pub"], ev["tE_pub"], max(ev["u0_pub"], 0.01))
        row.update(w3.summarise(ev["event_id"], res, lc.t.size))
        l0 = level0_fit(lc, ev["t0_pub"], ev["tE_pub"], max(ev["u0_pub"], 0.01))
        sel = sample_a_selection(lc, l0, res.get("PAR"))
        row.update({f"L0_{k}": v for k, v in l0.items()})
        row.update({f"sel_{k}": v for k, v in sel["passed"].items()})
        row.update(
            selected=sel["selected"],
            selected_ln=sel["selected_ln"],
            selected_ext=sel["selected_ext"],
            **{k: sel[k] for k in SEL_STATS},
        )
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001 — one bad light curve must not stop the run
        row["error"] = repr(exc)[:200]
    row["seconds"] = time.time() - t1
    return row


def _table(rows: list[dict], fill_str=("event_id", "method", "error", "best_ordinary")) -> Table:
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "event_id", k))
    return Table({k: [r.get(k, "" if k in fill_str else np.nan) for r in rows] for k in keys})


def run_fit(procs: int, limit: int | None = None) -> Path:
    survey = gm.GaiaDR3Microlensing()
    ev = survey.events()
    pos = survey.radec()
    if limit:
        ev = ev[:limit]
    survey.prefetch(ev["event_id"])

    def _lc(eid):
        try:
            return survey.light_curve(eid)
        except KeyError:  # recorded as this event's error row
            return None

    cols = ("event_id", "method", "t0_pub", "tE_pub", "u0_pub")
    jobs = [
        (
            {c: (r[c].item() if hasattr(r[c], "item") else r[c]) for c in cols},
            _lc(r["event_id"]),
            pos.get(r["event_id"]),
        )
        for r in ev
    ]
    t1 = time.time()
    with Pool(procs) as pool:
        rows = list(pool.imap_unordered(_fit_worker, jobs, 2))
    tab = _table(rows)
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"scripts/w3_gaia.py fit on {survey.name} ({gm.REFERENCE})",
        fit_params=json.dumps(asdict(w3.P)),
        sel_params=json.dumps(asdict(S)),
        wall_time_s=round(time.time() - t1, 1),
    )
    # a --limit smoke test must not replace the full table that ``inject`` and ``summary`` read
    name = "fits_gaia_dr3.ecsv" if limit is None else f"fits_gaia_dr3_limit{limit}.ecsv"
    path = out_dir() / name
    tab.write(path, overwrite=True)
    if limit is None:
        w3.write_ecsv_gz(tab, results_dir() / "fits_gaia_dr3.ecsv.gz")
    print(f"wrote {path}: {len(tab)} events, {time.time() - t1:.0f} s wall")
    return path


# ----------------------------------------------------------------------------- inject


def _inject_worker(job):
    kind, base_lc, base_flux, fs, prm = job
    if kind == "W3":
        lc = w3.inject_w3(base_lc, base_flux, fs, prm["t0"], prm["tE"], prm["u0"], prm["rho"])
    else:
        lc = w3.inject_w3(base_lc, base_flux, fs, prm["t0"], prm["tE"], prm["u0"], 0.0, 1.0, 1)
    k = np.convolve(lc.f, np.ones(3) / 3, mode="same")
    t0g = float(lc.t[int(np.argmax(k))])
    row = {"kind": kind, "event_id": base_lc.event_id, **prm}
    try:
        models = ("PSPL", "PAR", "N1neg") if kind == "W3" else ("PSPL", "PAR")
        res = w3.fit_event(lc, t0g, 30.0, 0.3, models=models, te_grid=(10.0, 100.0))
        l0 = level0_fit(lc, t0g, 30.0, 0.3)
        sel = sample_a_selection(lc, l0, res.get("PAR"))
        row.update({k: sel[k] for k in ("selected", "selected_ln", "selected_ext")})
        row.update({f"sel_{k}": v for k, v in sel["passed"].items()})
        if kind == "W3":
            best = min(res[m]["bic"] for m in w3.ORDINARY if m in res)
            row["dbic_N1neg"] = res["N1neg"]["bic"] - best
            row["flagged"] = bool(row["dbic_N1neg"] < w3.P.flag_dbic)
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001
        row["error"] = repr(exc)[:200]
    return row


def run_inject(procs: int, per_cell: int, seed: int = 61) -> Path:
    """Bases: candidates whose PSPL fit is good (χ²/dof ≤ 1.5, ≥ 30 points); the fitted event is
    subtracted (``w3.baseline_from_fit``) and a simulated one injected at a uniform t0 in the
    DR3 window."""
    survey = gm.GaiaDR3Microlensing()
    fits = Table.read(out_dir() / "fits_gaia_dr3.ecsv")
    pos = survey.radec()
    ok = (
        no_error(fits)
        & (fits["PSPL_chi2"] / fits["PSPL_dof"] <= 1.5)
        & (fits["PSPL_fs"] > 0)
        & (fits["n_points"] >= 30)
    )
    bases = {}
    for r in fits[ok]:
        eid = str(r["event_id"])
        lc = gaia_light_curve(survey.light_curve(eid), eid, pos.get(eid))
        fit = {k: float(r[f"PSPL_{k}"]) for k in ("t0", "tE", "u0", "fs", "fb")}
        bases[eid] = w3.baseline_from_fit(lc, fit)
    ids = sorted(bases)
    rng = np.random.default_rng(seed)
    jobs = []
    for te in INJ_TE:
        for kind, rhos in (("W3", INJ_RHO), ("PSPL", (0.0,))):
            for rho in rhos:
                for _ in range(per_cell):
                    eid = ids[int(rng.integers(len(ids)))]
                    blc, bflux, fs = bases[eid]
                    prm = {
                        "tE": te,
                        "rho": rho,
                        "u0": float(rng.uniform(0.0, 2.0 if kind == "W3" else 1.0)),
                        "t0": float(rng.uniform(gm.T_FIRST, gm.T_LAST)),
                        "G_base": float(w3.flux_to_mag(bflux)),
                        "fs_frac": float(min(fs / bflux, 1.0)),
                    }
                    jobs.append((kind, blc, bflux, fs, prm))
    t1 = time.time()
    with Pool(procs) as pool:
        rows = list(pool.imap_unordered(_inject_worker, jobs, 2))
    tab = _table(rows, fill_str=("kind", "event_id", "error"))
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/w3_gaia.py inject: exotic_sim.inject_light_curve (simulated) on "
        f"PSPL-subtracted Gaia DR3 light curves of {len(ids)} candidates",
        sel_params=json.dumps(asdict(S)),
        seed=seed,
        wall_time_s=round(time.time() - t1, 1),
    )
    path = out_dir() / "injections_gaia_dr3.ecsv"
    tab.write(path, overwrite=True)
    w3.write_ecsv_gz(tab, results_dir() / "injections_gaia_dr3.ecsv.gz")
    print(f"wrote {path}: {len(tab)} rows, {time.time() - t1:.0f} s")
    return path


def summary() -> dict:
    out = {}
    fits = Table.read(out_dir() / "fits_gaia_dr3.ecsv")
    ok = no_error(fits)
    groups = (("A", np.isin(fits["method"], ["A", "A+B"])), ("B_only", fits["method"] == "B"))
    for name, m in groups:
        sub = fits[ok & m]
        out[f"audit_{name}"] = {
            "n": len(sub),
            "selected": int(np.sum(sub["selected"])),
            "selected_ln": int(np.sum(sub["selected_ln"])),
            "selected_ext": int(np.sum(sub["selected_ext"])),
            "fail": _fails(sub),
        }
    out["fit_errors"] = int(np.sum(~ok))
    out["flags"] = sorted(
        (str(r["event_id"]), round(float(r["dbic_min"]), 1), str(r["method"]))
        for r in fits[ok & (fits["dbic_min"] < w3.P.flag_dbic)]
    )
    path = out_dir() / "injections_gaia_dr3.ecsv"
    if path.exists():
        inj = Table.read(path)
        inj = inj[no_error(inj)]
        cells = {}
        for kind in ("W3", "PSPL"):
            for te in INJ_TE:
                for rho in INJ_RHO if kind == "W3" else (0.0,):
                    c = inj[(inj["kind"] == kind) & (inj["tE"] == te) & (inj["rho"] == rho)]
                    d = {
                        "n": len(c),
                        "selected": int(np.sum(c["selected"])),
                        "selected_ln": int(np.sum(c["selected_ln"])),
                        "selected_ext": int(np.sum(c["selected_ext"])),
                    }
                    if kind == "W3":
                        d["flagged"] = int(np.sum(c["flagged"]))
                    cells[f"{kind} tE={te:g} rho={rho:g}"] = d
        out["injections"] = cells
        out["W3_fail"] = _fails(inj[inj["kind"] == "W3"])
        out["PSPL_fail"] = _fails(inj[inj["kind"] == "PSPL"])
    print(json.dumps(out, indent=1))
    return out


def write_manifest() -> Path:
    """``data/manifests/gaia_dr3_mulens.ecsv``: URI, sha256 and size of every cached input file."""
    import urllib.parse

    import jwst_anomaly
    from jwst_anomaly import photometry

    def tap(query):
        return f"{gm.TAP_SYNC}?" + urllib.parse.urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
        )

    survey = gm.GaiaDR3Microlensing()
    ids = list(survey.events()["event_id"])
    files = {
        survey._cache / "vari_microlensing.csv": tap(gm.EVENTS_QUERY),
        survey._cache / "positions.csv": tap(
            "SELECT source_id, ra, dec FROM gaiadr3.gaia_source WHERE source_id IN "
            "(<the vari_microlensing source_ids>)"
        ),
    }
    for e in ids:
        files[survey._lc_path(e)] = f"{gm.DATALINK}?" + urllib.parse.urlencode(
            {
                "RETRIEVAL_TYPE": "EPOCH_PHOTOMETRY",
                "DATA_STRUCTURE": "INDIVIDUAL",
                "FORMAT": "CSV",
                "VALID_DATA": "false",
                "ID": e,
            }
        )
    files = {f: u for f, u in files.items() if f.exists()}
    tab = Table(
        {
            "uri": list(files.values()),
            "file": [f.name for f in files],
            "sha256": [photometry.sha256_file(f) for f in files],
            "size": [f.stat().st_size for f in files],
        }
    )
    tab.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=f"{gm.TAP_SYNC} ({gm.EVENTS_QUERY}; gaia_source positions) and {gm.DATALINK} "
        f"EPOCH_PHOTOMETRY INDIVIDUAL CSV; {gm.REFERENCE}; method from {gm.PAPER_SRC} "
        f"(sha256 {gm.PAPER_SHA})",
        retrieved=time.strftime("%Y-%m-%d", time.gmtime()),
        pipeline_version=jwst_anomaly.__version__,
    )
    path = paths.repo_root() / "data" / "manifests" / "gaia_dr3_mulens.ecsv"
    tab.write(path, overwrite=True)
    print(f"wrote {path}: {len(tab)} files")
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("fetch", help="event table, positions and epoch photometry (cached)")
    f = sub.add_parser("fit", help="fit every candidate and emulate the Sample A selection")
    f.add_argument("--procs", type=int, default=4)
    f.add_argument("--limit", type=int)
    i = sub.add_parser("inject", help="W3 and PSPL injections through the emulated selection")
    i.add_argument("--procs", type=int, default=4)
    i.add_argument("--per-cell", type=int, default=40)
    sub.add_parser("summary", help="audit, flags and injection pass rates")
    sub.add_parser("manifest", help="write data/manifests/gaia_dr3_mulens.ecsv")
    a = ap.parse_args(argv)
    if a.cmd == "fetch":
        s = gm.GaiaDR3Microlensing()
        print(len(s.events()), len(s.radec()), "fetched", s.prefetch())
    elif a.cmd == "fit":
        run_fit(a.procs, a.limit)
    elif a.cmd == "inject":
        run_inject(a.procs, a.per_cell)
    elif a.cmd == "manifest":
        write_manifest()
    else:
        summary()
    return 0


if __name__ == "__main__":
    sys.exit(main())

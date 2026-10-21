"""S3 (D-069, D-072 revisit): hybrid copies in SN Refsdal light curves.

Hypothesis (hypothesis only): image j carries a faint copy of the SN at lag delta relative to
its own signal, F_j(t) = m F_1(t - tau) + a F_1(t - tau - delta) + c, r = a / m. Ordinary:
microlensing, shared difference-imaging covariance, template interpolation error, and the slow
SN 1987A-like light curve making short-lag copies collinear with the main term.

Chain: CDS J/ApJ/948/93 table3 (Kelly et al. 2023; observed HST difference-image fluxes) -> S1 as
the template per filter, linear interpolation, zero outside its sampled span (ASSUMPTION) -> tau
grid fit of (m, c) (derived) -> delta grid fit of a with tau fixed -> r(delta), sigma_r (derived).
Frozen plan: results/s3_hybrid/refsdal/PLAN.md.

    python scripts/s3_refsdal.py run
"""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import numpy as np
import requests

TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
TABLE = "J/ApJ/948/93/table3"
TARGETS = ("S2", "S3", "S4", "SX")
FILTERS = ("F125W", "F160W")
TAU_GRID = {"SX": np.arange(250.0, 501.0, 1.0)}  # other images: -60..60 d
TAU_DEFAULT = np.arange(-60.0, 61.0, 1.0)
DELTA_GRID = np.concatenate([np.arange(-400.0, -29.0, 2.0), np.arange(30.0, 401.0, 2.0)])
FLAG_SNR = 5.0  # ASSUMPTION (PLAN.md)
OUT = Path(__file__).resolve().parents[1] / "results" / "s3_hybrid" / "refsdal"


def fetch() -> dict:
    """Stream table3 from TAPVizieR; returns {(obj, filt): (mjd, flux, err)} sorted by time."""
    q = f'SELECT Obj, Filt, MJD, Flux, e_Flux FROM "{TABLE}"'
    r = requests.get(
        TAP, params={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}, timeout=120
    )
    r.raise_for_status()
    rows = np.genfromtxt(
        io.StringIO(r.text), delimiter=",", names=True, dtype=None, encoding="utf-8"
    )
    out = {}
    for obj in np.unique(rows["Obj"]):
        for filt in np.unique(rows["Filt"]):
            s = (rows["Obj"] == obj) & (rows["Filt"] == filt)
            o = np.argsort(rows["MJD"][s])
            key = (str(obj).strip('"'), str(filt).strip('"'))
            out[key] = tuple(np.asarray(rows[c][s][o], float) for c in ("MJD", "Flux", "e_Flux"))
    return out


def template(t1: np.ndarray, f1: np.ndarray):
    """Linear interpolation of the S1 curve; 0 outside its sampled span (ASSUMPTION)."""
    return lambda t: np.interp(t, t1, f1, left=0.0, right=0.0)


def wls(design: np.ndarray, y: np.ndarray, w: np.ndarray):
    """Weighted least squares: coefficients, their covariance and chi^2."""
    a = design * w[:, None]
    cov = np.linalg.pinv(a.T @ a)
    beta = cov @ (a.T @ (y * w))
    chi2 = float(np.sum(((y - design @ beta) * w) ** 2))
    return beta, cov, chi2


def fit_tau(tpl, t, f, e, grid):
    """Best tau on the grid for f = m tpl(t - tau) + c."""
    w = 1.0 / e
    chi = []
    for tau in grid:
        d = np.column_stack([tpl(t - tau), np.ones_like(t)])
        chi.append(wls(d, f, w)[2])
    i = int(np.argmin(chi))
    return float(grid[i]), np.asarray(chi)


def scan_delta(tpl, t, f, e, tau):
    """r(delta) and sigma_r with tau fixed; m, a, c fitted jointly at each delta."""
    w = 1.0 / e
    r, sr = [], []
    for dl in DELTA_GRID:
        d = np.column_stack([tpl(t - tau), tpl(t - tau - dl), np.ones_like(t)])
        beta, cov, _ = wls(d, f, w)
        m, a = beta[0], beta[1]
        if not np.any(d[:, 1]):  # copy term outside the template span: not measurable
            r.append(np.nan)
            sr.append(np.nan)
            continue
        # delta-method error on a / m
        g = np.array([-a / m**2, 1.0 / m, 0.0])
        r.append(a / m)
        sr.append(float(np.sqrt(g @ cov @ g)))
    return np.asarray(r), np.asarray(sr)


def _frac(cond, values):
    """Fraction of measurable (finite) lags where cond holds."""
    ok = np.isfinite(values)
    return float(np.mean(cond[ok])) if ok.any() else float("nan")


def flagged(r1, s1, r2, s2, snr_min=FLAG_SNR):
    """Frozen rule (PLAN.md), per lag: r/sigma >= snr_min in both filters AND |r1 - r2| < 2 sigma.

    sigma is the combined error hypot(s1, s2); NaN lags (copy term unmeasurable) never flag.
    """
    with np.errstate(invalid="ignore"):
        both = (r1 / s1 >= snr_min) & (r2 / s2 >= snr_min)
        consistent = np.abs(r1 - r2) < 2.0 * np.hypot(s1, s2)
    return np.nan_to_num(both & consistent, nan=False).astype(bool)


def run(data: dict) -> dict:
    res = {"targets": {}, "gates": {}}
    scans = {}
    for filt in FILTERS:
        t1, f1, _ = data[("S1", filt)]
        tpl = template(t1, f1)
        for obj in TARGETS:
            t, f, e = data[(obj, filt)]
            tau, _ = fit_tau(tpl, t, f, e, TAU_GRID.get(obj, TAU_DEFAULT))
            main = np.column_stack([tpl(t - tau), np.ones_like(t)])
            chi2_main = wls(main, f, 1.0 / e)[2]
            r, sr = scan_delta(tpl, t, f, e, tau)
            snr = r / sr
            k = int(np.nanargmax(snr))
            scans[(obj, filt)] = (r, sr)
            res["targets"][f"{obj}_{filt}"] = {
                "tau_d": tau,
                "max_snr": float(snr[k]),
                "delta_at_max_d": float(DELTA_GRID[k]),
                "r_at_max": float(r[k]),
                "n_lags_measurable": int(np.isfinite(r).sum()),
                "r95_max": float(np.nanmax(r + 1.645 * sr)),
                "r95_median": float(np.nanmedian(r + 1.645 * sr)),
                "n_epochs": int(t.size),
                "chi2_dof_main": float(chi2_main / (t.size - 3)),  # post hoc diagnostic
            }
    # G1: known delays
    g1 = all(abs(res["targets"][f"SX_{fl}"]["tau_d"] - 376.0) <= 30 for fl in FILTERS) and all(
        abs(res["targets"][f"{o}_{fl}"]["tau_d"]) <= 30
        for o in ("S2", "S3", "S4")
        for fl in FILTERS
    )
    res["gates"]["G1"] = bool(g1)
    # G2: inject r = 0.2 at delta = +150 d into SX (F125W) using the fitted main term
    t1, f1, _ = data[("S1", "F125W")]
    tpl = template(t1, f1)
    t, f, e = data[("SX", "F125W")]
    tau = res["targets"]["SX_F125W"]["tau_d"]
    d = np.column_stack([tpl(t - tau), np.ones_like(t)])
    m = wls(d, f, 1.0 / e)[0][0]
    finj = f + 0.2 * m * tpl(t - tau - 150.0)
    r, sr = scan_delta(tpl, t, finj, e, tau)
    k = int(np.argmin(abs(DELTA_GRID - 150.0)))
    res["gates"]["G2"] = {
        "r_rec": float(r[k]),
        "sigma": float(sr[k]),
        "pass": bool(abs(r[k] - 0.2) < 0.1),
    }
    # frozen decision rule: flag only if both filters >= 5 sigma with consistent r
    flags = []
    for obj in TARGETS:
        hit = flagged(*scans[(obj, FILTERS[0])], *scans[(obj, FILTERS[1])])
        flags += [{"target": obj, "delta_d": float(dl)} for dl in DELTA_GRID[hit]]
    res["flags"] = flags
    # post hoc (added after the first run): the other filter at the F125W peak lag, and the F125W
    # peak significance with sigma rescaled by sqrt(chi2/dof) of the main-term fit
    post = {}
    for obj in TARGETS:
        (ra, sa), (rb, sb) = scans[(obj, FILTERS[0])], scans[(obj, FILTERS[1])]
        k = int(np.nanargmax(ra / sa))
        c = res["targets"][f"{obj}_{FILTERS[0]}"]["chi2_dof_main"]
        post[obj] = {
            "delta_d": float(DELTA_GRID[k]),
            "snr_F125W": float(ra[k] / sa[k]),
            "snr_F125W_rescaled": float(ra[k] / sa[k] / np.sqrt(c)),
            "snr_F160W_same_delta": float(rb[k] / sb[k]),
            "r_F125W": float(ra[k]),
            "r_F160W": float(rb[k]),
            "frac_lags_r_positive_F125W_delta_gt_0": _frac(
                ra[DELTA_GRID > 0] > 0, ra[DELTA_GRID > 0]
            ),
            "frac_lags_r_negative_F125W_delta_lt_0": _frac(
                ra[DELTA_GRID < 0] < 0, ra[DELTA_GRID < 0]
            ),
        }
    res["post_hoc_F125W_peaks"] = post
    return res


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["run"])
    p.parse_args()
    res = run(fetch())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()

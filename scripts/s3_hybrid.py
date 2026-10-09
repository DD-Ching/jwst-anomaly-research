"""S3 (D-069 round-1 survivor): hybrid images in COSMOGRAIL lensed-quasar light curves.

Hypothesis (A1 P2, docs/hypotheses/round-1/): a link near the lens plane joins the source-side leg
of image i to the observer-side leg of image j, so image j carries a faint copy of the quasar at
lag delta_ij = s_i - s_j relative to its own signal. Ordinary expectation: no copy at the model
lags; quasar red noise correlates the curves at every lag, which the off-model-lag null measures.

Lags (model_prediction, ASSUMPTION: singular isothermal sphere, Euclidean leg split with angular
diameter distances). For an SIS the geometric delays of the two images are equal, so the measured
delay Delta (trailing T minus leading L) is all potential delay. With the observer leg
o = D_l theta^2 / 2 and the source leg s = D_l^2 (theta - beta D_s / D_l)^2 / (2 D_ls), both times
(1 + z_l) / c, and a fraction f of the potential delay put on the source leg:
    s_L - s_T = -Delta * g,   g = D_ls / D_s + f,   f in [0, 1].
The split f is unobservable (D-069 B-on-A1), so the screen scans the whole window f in [0, 1] and
pays the trials factor through a null built from windows of the same width at off-model lags.

Chain: CDS J/A+A/640/A105 R-band curves (observed) -> fluxes normalised per image -> image j
modelled as B-spline microlensing + m * F_i(t - tau_ji) + a * F_i(t - tau_ji - delta) with F_i the
other image's curve interpolated inside a season (derived) -> r = a / m at each lag -> window
maximum against the off-model-window null -> injection-recovery of r through the same fit.

    python scripts/s3_hybrid.py run      # fetch (cached), validate, screen, inject, sensitivity
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import requests
from astropy.cosmology import Planck18
from scipy.interpolate import BSpline

from jwst_anomaly import paths

CDS = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/640/A105/"
OUT = Path(__file__).resolve().parents[1] / "results" / "s3_hybrid"


@dataclass(frozen=True)
class System:
    """One COSMOGRAIL double: published delay and redshifts (Millon et al. 2020, Tables 1 and 4)."""

    name: str
    file: str
    dt_ab: float  # observed (published) Delta t_AB = t_A - t_B in days; negative: A leads
    dt_err: float
    z_l: float
    z_s: float


# Doubles of Millon+2020 (arXiv:2002.05736) whose delay is not flagged "uncertain" and is > 3 sigma.
SYSTEMS = [
    System("HE0047-1756", "lcab/HE0047_Euler.dat", -10.4, 3.5, 0.407, 1.678),
    System("J0158-4325", "lcab/J0158_Euler.dat", -22.7, 3.6, 0.317, 1.29),
    System("J1226-0006", "lcab/J1226_Euler.dat", 33.7, 2.7, 0.517, 1.123),
    System("J1335+0118", "lcab/J1335_Euler.dat", -56.0, 6.1, 0.44, 1.570),
    System("J1455+1447", "lcab/J1455_Euler.dat", -47.2, 7.8, 0.42, 1.424),
    System("J1515+1511", "lcab/J1515_Euler.dat", -210.2, 5.7, 0.742, 2.054),
    System("J1620+1203", "lcab/J1620_Euler.dat", -171.5, 8.7, 0.398, 1.158),
]


@dataclass(frozen=True)
class Params:
    max_gap: float = (
        40.0  # days; template interpolation only across gaps shorter than this (ASSUMPTION)
    )
    knot_spacing: float = (
        3000.0  # days; microlensing B-spline knots; shorter spacings fail the delay validation
    )
    lag_step: float = 2.0  # days; lag grid
    null_span: float = 3.0  # null windows cover |lag| <= null_span * window span (ASSUMPTION)
    min_points: int = 40  # overlapping epochs needed for a fit
    wing: float = (
        10.0  # days; |lag| below this is the main-term wing, excluded everywhere (ASSUMPTION)
    )
    inject_r: float = 0.05  # injected copy / main flux ratio
    quantile: float = 0.95  # one-sided level for the window statistic and the limit (ASSUMPTION)


def fetch(system: System, cache: Path) -> Path:
    path = cache / system.file
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(CDS + system.file, timeout=60)
        r.raise_for_status()
        tmp = path.with_suffix(".part")
        tmp.write_text(r.text)
        tmp.replace(path)
    return path


def read_ab(path: Path) -> dict[str, np.ndarray]:
    """Byte-by-byte lcab format; magnitudes to fluxes normalised to their median."""
    t, ma, ea, mb, eb = [], [], [], [], []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        t.append(float(line[0:11]))
        ma.append(float(line[12:21]))
        ea.append(float(line[22:29]))
        mb.append(float(line[30:39]))
        eb.append(float(line[40:47]))
    t = np.array(t)
    out = {"t": t}
    for key, m, e in (("A", ma, ea), ("B", mb, eb)):
        f = 10 ** (-0.4 * np.array(m))
        med = np.median(f)
        out[key] = f / med
        out["e" + key] = np.array(e) * 0.4 * np.log(10) * f / med
    order = np.argsort(t)
    return {k: v[order] for k, v in out.items()}


def lag_window(system: System) -> dict:
    """Model lags (model_prediction) of both hybrid directions for f in [0, 1]."""
    d_s = Planck18.angular_diameter_distance(system.z_s).value
    d_ls = Planck18.angular_diameter_distance(system.z_l, system.z_s).value
    q = d_ls / d_s
    delta = abs(system.dt_ab)
    lead = "A" if system.dt_ab < 0 else "B"
    trail = "B" if lead == "A" else "A"
    # delta_ij relative to image j's own signal; positive = later.
    return {
        "dls_over_ds": q,
        "delay": delta,
        "lead": lead,
        "trail": trail,
        # copy at T of the light that travelled L's source leg: earlier than T by Delta * g
        "into_trail": (-delta * (q + 1.0), -delta * q),
        # copy at L of T's source leg: later than L by Delta * g
        "into_lead": (delta * q, delta * (q + 1.0)),
    }


def interp_template(
    t_src: np.ndarray, f_src: np.ndarray, t_eval: np.ndarray, max_gap: float
) -> np.ndarray:
    """Linear interpolation of a light curve, NaN outside its span or across gaps > max_gap."""
    out = np.interp(t_eval, t_src, f_src, left=np.nan, right=np.nan)
    k = np.searchsorted(t_src, t_eval)
    k = np.clip(k, 1, len(t_src) - 1)
    gap = t_src[k] - t_src[k - 1]
    out[gap > max_gap] = np.nan
    exact = np.isin(t_eval, t_src)
    out[exact] = np.interp(t_eval[exact], t_src, f_src)
    return out


def spline_basis(t: np.ndarray, spacing: float) -> np.ndarray:
    lo, hi = t.min(), t.max()
    n = max(int(np.ceil((hi - lo) / spacing)), 1)
    inner = np.linspace(lo, hi, n + 1)
    knots = np.r_[[lo] * 3, inner, [hi] * 3]
    return BSpline.design_matrix(np.clip(t, lo, hi), knots, 3).toarray()


def fit_copy(
    t: np.ndarray,
    f_j: np.ndarray,
    e_j: np.ndarray,
    t_i: np.ndarray,
    f_i: np.ndarray,
    tau_ji: float,
    lag: float | None,
    p: Params,
) -> tuple[float, float, int]:
    """Weighted fit of image j; returns (r = a / m, m, n). lag=None fits the main term only."""
    main = interp_template(t_i, f_i, t - tau_ji, p.max_gap)
    cols = [main]
    if lag is not None:
        cols.append(interp_template(t_i, f_i, t - tau_ji - lag, p.max_gap))
    good = np.all(np.isfinite(cols), axis=0)
    if good.sum() < p.min_points:
        return np.nan, np.nan, int(good.sum())
    basis = spline_basis(t[good], p.knot_spacing)
    x = np.column_stack([basis] + [c[good] - np.mean(c[good]) for c in cols])
    w = 1.0 / e_j[good]
    coef, *_ = np.linalg.lstsq(x * w[:, None], f_j[good] * w, rcond=None)
    m = coef[basis.shape[1]]
    a = coef[basis.shape[1] + 1] if lag is not None else 0.0
    return a / m, m, int(good.sum())


def scan(lc: dict, src: str, dst: str, tau_ji: float, lags: np.ndarray, p: Params, f_j=None):
    f_j = lc[dst] if f_j is None else f_j
    return np.array(
        [fit_copy(lc["t"], f_j, lc["e" + dst], lc["t"], lc[src], tau_ji, lag, p)[0] for lag in lags]
    )


def window_stat(lags: np.ndarray, r: np.ndarray, lo: float, hi: float, wing: float = 0.0) -> float:
    """Window maximum of r, ignoring lags on the main-term wing (|lag| < wing)."""
    sel = (lags >= lo) & (lags <= hi) & (np.abs(lags) >= wing) & np.isfinite(r)
    return float(np.max(r[sel])) if sel.any() else np.nan


def validate_delay(lc: dict, system: System, p: Params) -> dict:
    """Known case: the main-term fit quality peaks near the published delay."""
    taus = np.arange(-300.0, 300.0 + p.lag_step, p.lag_step)
    chi = []
    for tau in taus:
        main = interp_template(lc["t"], lc["A"], lc["t"] - tau, p.max_gap)
        good = np.isfinite(main)
        if good.sum() < p.min_points:
            chi.append(np.nan)
            continue
        basis = spline_basis(lc["t"][good], p.knot_spacing)
        x = np.column_stack([basis, main[good]])
        w = 1.0 / lc["eB"][good]
        coef, *_ = np.linalg.lstsq(x * w[:, None], lc["B"][good] * w, rcond=None)
        chi.append(float(np.mean(((x @ coef - lc["B"][good]) * w) ** 2)))
    chi = np.array(chi)
    best = float(taus[np.nanargmin(chi)])
    # tau = t_B - t_A = -dt_ab in the CDS convention
    return {"best_tau_BA": best, "published_tau_BA": -system.dt_ab, "dt_err": system.dt_err}


def run_system(system: System, cache: Path, p: Params, rng: np.random.Generator) -> dict:
    lc = read_ab(fetch(system, cache))
    win = lag_window(system)
    lead, trail, delay = win["lead"], win["trail"], win["delay"]
    res = {"system": asdict(system), "n_epochs": int(len(lc["t"])), "window": win}
    v = validate_delay(lc, system, p)
    miss = abs(v["best_tau_BA"] - v["published_tau_BA"])
    # The sign must match and the delay must differ from 0, or a delay-free fit would pass.
    v["passed"] = bool(
        miss <= 3 * system.dt_err + 2 * p.lag_step
        and np.sign(v["best_tau_BA"]) == np.sign(v["published_tau_BA"])
        and abs(v["best_tau_BA"]) > 2 * p.lag_step
    )
    res["validation"] = v
    out = {}
    if not v["passed"]:  # the template fit cannot find the known delay: no screen on this system
        res["directions"] = out
        return res
    for name, (src, dst, tau) in {
        "into_trail": (lead, trail, delay),
        "into_lead": (trail, lead, -delay),
    }.items():
        lo, hi = win[name]
        span = hi - lo
        reach = p.null_span * max(span, delay)
        lags = np.arange(-reach - span, reach + span + p.lag_step, p.lag_step)
        r = scan(lc, src, dst, tau, lags, p)
        # Null: same-width windows off the model window; the wing is excluded from both statistics.
        # The windows overlap (start every 2 lag steps), so the p-value is approximate (D-072).
        starts = np.arange(-reach - span, reach, p.lag_step * 2)
        null = [
            window_stat(lags, r, s, s + span, p.wing) for s in starts if s + span < lo or s > hi
        ]
        null = np.array([v for v in null if np.isfinite(v)])
        obs = window_stat(lags, r, lo, hi, p.wing)
        ok_null = len(null) >= 10 and np.isfinite(obs)
        thresh = float(np.quantile(null, p.quantile)) if ok_null else np.nan
        # Injection through the screen statistic: a copy of amplitude inject_r at a random lag in
        # the window, then the window maximum on the same lag grid, minus the maximum without it.
        m = fit_copy(lc["t"], lc[dst], lc["e" + dst], lc["t"], lc[src], tau, None, p)[1]
        sel = (lags >= lo) & (lags <= hi) & (np.abs(lags) >= p.wing)
        rec = []
        for _ in range(8 if ok_null else 0):
            lag = rng.uniform(lo, hi)
            copy = interp_template(lc["t"], lc[src], lc["t"] - tau - lag, p.max_gap)
            ok = np.isfinite(copy)
            sub = {k: v[ok] for k, v in lc.items()}
            f_inj = sub[dst] + p.inject_r * m * (copy[ok] - 1.0)
            r_inj = np.array(
                [
                    fit_copy(sub["t"], f_inj, sub["e" + dst], lc["t"], lc[src], tau, x, p)[0]
                    for x in lags[sel]
                ]
            )
            r_0 = np.array(
                [
                    fit_copy(sub["t"], sub[dst], sub["e" + dst], lc["t"], lc[src], tau, x, p)[0]
                    for x in lags[sel]
                ]
            )
            rec.append((np.nanmax(r_inj) - np.nanmax(r_0)) / p.inject_r)
        eff = float(np.nanmedian(rec)) if rec else np.nan
        # Sensitivity, not a limit (D-072): a copy r adds ~ eff * r to the window maximum; the
        # one-sided 95 % bound (ASSUMPTION) uses the null spread; eff > 1 is degeneracy, capped.
        if ok_null and eff > 0.5:
            limit = float(max(obs, np.median(null)) + (thresh - np.median(null))) / min(eff, 1.0)
        else:
            limit = np.nan
        out[name] = {
            "window_days": [lo, hi],
            "observed_window_max_r": obs,
            "null_n": int(len(null)),
            "null_median": float(np.median(null)) if len(null) else np.nan,
            "null_q95": thresh,
            "p_value": float((np.sum(null >= obs) + 1) / (len(null) + 1)) if ok_null else np.nan,
            "injection_efficiency": eff,
            "r95_sensitivity": limit,
            "lags": lags.tolist(),
            "r": [None if not np.isfinite(v) else round(float(v), 5) for v in r],
        }
    res["directions"] = out
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--seed", type=int, default=20261009)
    args = ap.parse_args(argv)
    p = Params()
    cache = paths.data_root() / "cosmograil_xix"
    rng = np.random.default_rng(args.seed)
    OUT.mkdir(parents=True, exist_ok=True)
    results = [run_system(s, cache, p, rng) for s in SYSTEMS]
    summary = {"params": asdict(p), "source": CDS, "provenance": "derived", "systems": []}
    for res in results:
        row = {
            "name": res["system"]["name"],
            "validation": res["validation"],
            "dls_over_ds": res["window"]["dls_over_ds"],
        }
        for k, d in res["directions"].items():
            row[k] = {kk: d[kk] for kk in d if kk not in ("lags", "r")}
        summary["systems"].append(row)
        print(json.dumps(row))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    (OUT / "scans.json").write_text(json.dumps(results) + "\n")


if __name__ == "__main__":
    main()

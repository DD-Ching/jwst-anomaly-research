"""E1 "causal event network": pair-count excess versus time lag between localized events.

Hypothesis E1 (owner idea 4; docs/hypotheses/round-1/summary.md "Correction to the causal-network
idea", B-on-A1.md section 7): time-stamped, localized events from *different* directions
(separation well beyond the localization errors) might be statistically dependent at time lags that
no ordinary path explains. Dependence, if found, is an anomaly, never evidence of anything exotic:
observer-side common causes (instrument uptime and exposure, Sun constraints, alert chains, shared
triggers, catalogue duplicates, background events) come first.

The statistic is a cross-correlation function in lag: for each catalogue pair (A, B) the number of
pairs in lag bins, split by angular separation into ``same`` (sep <= n_sigma * sigma_comb) and
``wide`` (sep > n_sigma * sigma_comb and > min_wide_deg). Events without a localization (GW from
the GWTC CSV) enter an ``all`` class only.

Two nulls. Both permute times among the events of the same catalogue and calendar year and keep
each event's declination and hour angle (sidereal scrambling), so Earth-fixed exposure such as the
CHIME transit or the IceCube zenith acceptance is kept:

- ``perm``: the permutation only. The multiset of times is unchanged, so the lag histogram summed
  over separation is identical to the data: this null tests only the coupling of separation and
  lag (is the wide fraction at lag tau unusual?).
- ``jit``: the permutation, then each time is shifted independently. For GBM the shift is
  k * Fermi orbit + U(-orbit_slop, +orbit_slop) with |shift| <= jitter, which keeps the orbital
  phase (hence SAA and orbit-locked uptime); for the other catalogues it is U(-jitter, +jitter).
  This emulates uptime and the long-term rate from each catalogue's own empirical times and
  destroys clustering below ~jitter, so it is the null for lag excess. (A jitstrap with
  replacement was rejected: duplicated draws make self-pairs at lags < 2 * jitter.)

All thresholds are ASSUMPTIONs and live in :class:`Params`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy.table import Table

from .schema import Provenance

DAY = 86400.0
SIDEREAL_DAY = 86164.0905
#: Fermi orbital period (s), ~95.6 min (Meegan+2009, ApJ 702, 791). ASSUMPTION: constant over the
#: mission.
FERMI_ORBIT = 5736.0
#: Catalogue codes used in tables and channel names.
CATS = ("GBM", "ICECAT", "GW", "CHIME")


@dataclass(frozen=True)
class Params:
    # --- lag bins (s); pre-registered by the task brief (ASSUMPTION) ---
    lag_edges: tuple[float, ...] = (0.0, 10.0, 100.0, 3600.0, DAY, 7 * DAY)
    # --- separation classes (ASSUMPTION) ---
    n_sigma: float = 3.0  # same-direction radius in units of the combined 1-sigma error
    min_wide_deg: float = (
        0.1  # a wide pair is also > 0.1 deg (the geometric limit of B-on-A1 section 7)
    )
    # --- localization (ASSUMPTION) ---
    gbm_sys_deg: float = (
        3.7  # GBM systematic, 68 % core (Connaughton+2015, arXiv:1411.2685), as in D-071
    )
    gbm_external_deg: float = (
        0.5  # error_radius below this = position from another instrument (0 in catalogue)
    )
    gbm_external_sigma_deg: float = 0.05  # 1-sigma used for externally localized GBM bursts
    icecat_90_to_sigma: float = 2.146  # 2-D Gaussian: r90 = 2.146 sigma
    # --- jit null (ASSUMPTION) ---
    jitter_s: float = 3 * DAY
    # shift quantum per catalogue; 0 = continuous. Hour angle is kept by construction, so the
    # ground instruments need no sidereal quantum.
    quantum_s: dict = field(default_factory=lambda: {"GBM": FERMI_ORBIT})
    orbit_slop_s: float = 300.0  # GBM: orbit phase kept to +-5 min (~5 % of the orbit)
    # --- vetting windows (ASSUMPTION): half-widths around k * period ---
    orbit_k: tuple[int, ...] = tuple(range(1, 16))
    orbit_half_s: float = 120.0
    day_k: tuple[int, ...] = tuple(range(1, 8))
    day_half_s: float = 600.0
    # --- injection (ASSUMPTION): anchors within this many days of the moved event ---
    inject_local_days: float = 30.0


# ----------------------------------------------------------------------------------------- events


def _events(name, cat, mjd, ra, dec, sigma, source: str) -> Table:
    t = Table(
        {
            "name": np.asarray(name, dtype=str),
            "cat": np.full(len(mjd), cat),
            "mjd": np.asarray(mjd, dtype=float),
            "ra": np.asarray(ra, dtype=float),
            "dec": np.asarray(dec, dtype=float),
            "sigma": np.asarray(sigma, dtype=float),
        }
    )
    t["year"] = mjd_year(t["mjd"])
    t.sort("mjd")
    t.meta = {"provenance": str(Provenance.DERIVED), "source": source}
    return t


def mjd_year(mjd) -> np.ndarray:
    """Calendar year (UTC) of each MJD."""
    from astropy.time import Time

    return (
        np.asarray(
            Time(np.asarray(mjd, dtype=float), format="mjd").datetime64.astype("datetime64[Y]")
        ).astype(int)
        + 1970
    )


def gbm_events(cat: Table, p: Params | None = None) -> Table:
    """HEASARC fermigbrst rows -> events. 1-sigma = stat (+) systematic; below gbm_external_deg
    (position from another instrument, usually error_radius = 0) it is max(error_radius, 0.05°)."""
    p = p or Params()
    er = np.asarray(cat["error_radius"], dtype=float)
    sig = np.where(
        er < p.gbm_external_deg,
        np.maximum(er, p.gbm_external_sigma_deg),
        np.hypot(er, p.gbm_sys_deg),
    )
    return _events(
        cat["trigger_name"], "GBM", cat["trigger_time"], cat["ra"], cat["dec"], sig, "fermigbrst"
    )


def icecat_events(df, p: Params | None = None, drop_cr_veto: bool = True) -> Table:
    """ICECAT-1 rows (pandas) -> events. 1-sigma = mean of the four 90 % errors / 2.146."""
    p = p or Params()
    if drop_cr_veto:
        df = df[~df["CR_VETO"].astype(str).str.upper().eq("TRUE")]
    e90 = df[["RA_ERR_PLUS", "RA_ERR_MINUS", "DEC_ERR_PLUS", "DEC_ERR_MINUS"]].to_numpy(float)
    sig = np.nanmean(np.abs(e90), axis=1) / p.icecat_90_to_sigma
    return _events(df["NAME"], "ICECAT", df["EVENTMJD"], df["RA"], df["DEC"], sig, "ICECAT-1 v4")


def gw_events(df) -> Table:
    """GWTC CSV rows (pandas) -> events with no localization (sky maps are not in the CSV)."""
    from astropy.time import Time

    mjd = Time(df["gps"].to_numpy(float), format="gps").utc.mjd
    n = len(df)
    return _events(
        df["name"],
        "GW",
        mjd,
        np.full(n, np.nan),
        np.full(n, np.nan),
        np.full(n, np.nan),
        "GWTC CSV",
    )


def chime_events(df, one_per_source: bool = True) -> Table:
    """CHIME/FRB Catalog 2 rows (pandas) -> events. Sub-bursts (sub_num > 0) are dropped; with
    ``one_per_source`` each repeater contributes only its first burst (repeaters are one node)."""
    df = df[df["sub_num"] == 0].sort_values("mjd_400")
    if one_per_source:
        key = df["repeater_name"].fillna(df["tns_name"]).astype(str)
        df = df[~key.duplicated()]
    sig = np.maximum(df["ra_err"].to_numpy(float), df["dec_err"].to_numpy(float))
    return _events(
        df["tns_name"], "CHIME", df["mjd_400"], df["ra"], df["dec"], sig, "CHIME/FRB Catalog 2"
    )


# ----------------------------------------------------------------------------------------- geometry


def gmst_deg(mjd) -> np.ndarray:
    """Greenwich mean sidereal angle (deg); only differences matter here."""
    return np.mod(280.46061837 + 360.98564736629 * (np.asarray(mjd, dtype=float) - 51544.5), 360.0)


def sep_deg(ra1, dec1, ra2, dec2) -> np.ndarray:
    """Great-circle separation (deg), haversine form."""
    r1, d1, r2, d2 = (np.radians(np.asarray(x, dtype=float)) for x in (ra1, dec1, ra2, dec2))
    a = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))


# ----------------------------------------------------------------------------------------- pairs


def pairs_within(ta: np.ndarray, tb: np.ndarray, max_lag_s: float, same: bool):
    """Index pairs (i, j) with |ta[i] - tb[j]| <= max_lag (times in MJD, any order). For ``same``
    (ta is tb) only i < j in time order is returned, each unordered pair once."""
    ob = np.argsort(tb, kind="stable")
    sb = tb[ob]
    w = max_lag_s / DAY
    lo = np.searchsorted(sb, ta - w, side="left")
    hi = np.searchsorted(sb, ta + w, side="right")
    if same:
        oa_rank = np.empty(len(ta), dtype=np.int64)
        oa_rank[ob] = np.arange(len(ta))
        lo = np.maximum(lo, oa_rank + 1)
    n = np.maximum(hi - lo, 0)
    i = np.repeat(np.arange(len(ta)), n)
    start = np.repeat(lo, n)
    off = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
    j = ob[start + off]
    return i, j


def sep_class(sep, sig_a, sig_b, p: Params) -> np.ndarray:
    """0 = same, 1 = wide, -1 = between (fails only the 0.1 deg floor), 2 = unlocalized."""
    sc = np.hypot(sig_a, sig_b)
    cls = np.where(
        sep <= p.n_sigma * sc, 0, np.where(sep > np.maximum(p.n_sigma * sc, p.min_wide_deg), 1, -1)
    )
    return np.where(np.isfinite(sc) & np.isfinite(sep), cls, 2)


@dataclass
class Sample:
    """Arrays of one catalogue (possibly scrambled)."""

    cat: str
    mjd: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    sigma: np.ndarray
    year: np.ndarray

    @classmethod
    def from_table(cls, t: Table) -> Sample:
        return cls(
            str(t["cat"][0]),
            np.asarray(t["mjd"], float),
            np.asarray(t["ra"], float),
            np.asarray(t["dec"], float),
            np.asarray(t["sigma"], float),
            np.asarray(t["year"], int),
        )

    def with_times(self, new_mjd: np.ndarray) -> Sample:
        """Same events at new times; Dec and hour angle kept (RA follows sidereal time)."""
        ra = np.mod(self.ra + gmst_deg(new_mjd) - gmst_deg(self.mjd), 360.0)
        return Sample(self.cat, new_mjd, ra, self.dec, self.sigma, self.year)


def windows_for(p: Params) -> list[tuple[str, list[tuple[float, float]]]]:
    """Lag windows: the pre-registered bins, then (vetting) orbit / day multiples."""
    e = p.lag_edges
    w = [(f"lag{i}", [(e[i], e[i + 1])]) for i in range(len(e) - 1)]
    w.append(
        (
            "orbit_k",
            [
                (k * FERMI_ORBIT - p.orbit_half_s, k * FERMI_ORBIT + p.orbit_half_s)
                for k in p.orbit_k
            ],
        )
    )
    w.append(
        (
            "sidereal_k",
            [(k * SIDEREAL_DAY - p.day_half_s, k * SIDEREAL_DAY + p.day_half_s) for k in p.day_k],
        )
    )
    w.append(("solar_k", [(k * DAY - p.day_half_s, k * DAY + p.day_half_s) for k in p.day_k]))
    return w


def lag_labels(p: Params) -> list[str]:
    def f(s):
        if s < 3600:
            return f"{s:g}s"
        if s < DAY:
            return f"{s / 3600:g}h"
        return f"{s / DAY:g}d"

    e = p.lag_edges
    return [f"{f(e[i])}-{f(e[i + 1])}" for i in range(len(e) - 1)]


def count_channel(a: Sample, b: Sample, p: Params, same: bool, windows=None) -> np.ndarray:
    """Pair counts, shape (n_windows, 3): classes (same, wide, all = every pair)."""
    windows = windows_for(p) if windows is None else windows
    max_lag = max(hi for _, iv in windows for _, hi in iv)
    i, j = pairs_within(a.mjd, b.mjd, max_lag, same)
    lag = np.abs(b.mjd[j] - a.mjd[i]) * DAY
    sep = sep_deg(a.ra[i], a.dec[i], b.ra[j], b.dec[j])
    cls = sep_class(sep, a.sigma[i], b.sigma[j], p)
    out = np.zeros((len(windows), 3), dtype=np.int64)
    c0, c1 = cls == 0, cls == 1
    for k, (_, ivs) in enumerate(windows):
        if len(ivs) == 1:
            lo, hi = ivs[0]
            m = (lag >= lo) & (lag < hi) if lo > 0 else (lag <= hi)
        else:  # union of disjoint intervals: one searchsorted over the sorted edges
            edges = np.ravel(ivs)
            m = (np.searchsorted(edges, lag, side="right") % 2) == 1
        out[k, 0] = np.count_nonzero(m & c0)
        out[k, 1] = np.count_nonzero(m & c1)
        out[k, 2] = np.count_nonzero(m)
    return out


# ----------------------------------------------------------------------------------------- nulls


def scramble_perm(s: Sample, rng: np.random.Generator) -> Sample:
    """Permute times among events of the same calendar year."""
    new = s.mjd.copy()
    for y in np.unique(s.year):
        idx = np.flatnonzero(s.year == y)
        new[idx] = s.mjd[rng.permutation(idx)]
    return s.with_times(new)


def scramble_jit(s: Sample, rng: np.random.Generator, p: Params) -> Sample:
    """Per-year permutation, then an independent shift per event (orbit phase kept for GBM)."""
    perm = scramble_perm(s, rng)
    q = p.quantum_s.get(s.cat, 0.0)
    if q > 0:
        kmax = int(p.jitter_s // q)
        shift = rng.integers(-kmax, kmax + 1, size=len(s.mjd)) * q
        shift = shift + rng.uniform(-p.orbit_slop_s, p.orbit_slop_s, size=len(s.mjd))
    else:
        shift = rng.uniform(-p.jitter_s, p.jitter_s, size=len(s.mjd))
    return s.with_times(perm.mjd + shift / DAY)


def empirical_p(obs: np.ndarray, null: np.ndarray) -> np.ndarray:
    """Upper-tail p per cell, (1 + #(null >= obs)) / (N + 1); scramble axis first."""
    return (1.0 + (null >= obs[None]).sum(axis=0)) / (null.shape[0] + 1.0)


def pooled_cell_p(obs: np.ndarray, null: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Per-cell upper-tail p for every row of the pooled ensemble (the N scrambles plus the
    observation as the last row), with one formula for all rows: (#rows >= value) / (N + 1).
    The observation is exchangeable with the scrambles under H0, so the observed row and the
    scramble rows are treated identically (review finding 1 on PR #112)."""
    allv = np.vstack([null[:, mask], obs[mask][None]])
    n = allv.shape[0]
    srt = np.sort(allv, axis=0)
    out = np.empty(allv.shape, dtype=float)
    for c in range(allv.shape[1]):
        out[:, c] = (n - np.searchsorted(srt[:, c], allv[:, c], side="left")) / n
    return out


def global_p(obs: np.ndarray, null: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    """Trials-corrected p over the ``mask`` cells: the observed minimum per-cell p, ranked among the
    same minimum for every row of the pooled ensemble (scrambles + observation). Uniform under H0;
    its floor is 1 / (N + 1)."""
    pv = pooled_cell_p(obs, null, mask)
    mins = pv.min(axis=1)
    return float(mins[-1]), float(np.count_nonzero(mins <= mins[-1]) / len(mins))


def empirical_floor(n_scrambles: int) -> float:
    """Smallest per-cell p an ensemble of n scrambles can give, 1 / (n + 1)."""
    return 1.0 / (n_scrambles + 1.0)


def reachable(n_scrambles: int, alpha: float, n_cells: int = 1) -> bool:
    """Can an empirical ensemble of this size cross a (Bonferroni) trials-corrected threshold
    alpha over n_cells? With 1e4 scrambles 3 sigma (1.35e-3) is reachable for one cell, but no
    family-wise 5 sigma claim (2.9e-7) is: that needs the analytic tail (:func:`analytic_p`)."""
    return empirical_floor(n_scrambles) * n_cells <= alpha


#: ASSUMPTION: below this null mean the cell is Poisson; above it Gaussian with the null sd.
POISSON_MEAN_MAX = 30.0


def analytic_p(obs: np.ndarray, null: np.ndarray) -> np.ndarray:
    """Upper-tail p per cell from a fitted null: Poisson(mean) when the null mean is below
    POISSON_MEAN_MAX, else Gaussian(mean, sd of the ensemble), which keeps over-dispersion. Unlike
    the empirical p it is not floored at 1 / (n + 1), so tails beyond the ensemble are reachable;
    it is an extrapolation (model_prediction) and is quoted beside the empirical p."""
    from scipy import stats

    mu, sd = null.mean(axis=0), null.std(axis=0)
    obs = np.asarray(obs, dtype=float)
    pois = stats.poisson.sf(obs - 1, np.maximum(mu, 1e-12))
    with np.errstate(divide="ignore", invalid="ignore"):
        gaus = stats.norm.sf(np.where(sd > 0, (obs - mu) / sd, 0.0))
    return np.where(mu < POISSON_MEAN_MAX, pois, np.where(sd > 0, gaus, 1.0))


def z_score(obs, null) -> np.ndarray:
    mu, sd = null.mean(axis=0), null.std(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(sd > 0, (obs - mu) / sd, 0.0)


def upper_limit_95(obs: np.ndarray, null: np.ndarray) -> np.ndarray:
    """95 % upper limit on the mean number of additional (dependent) pairs per cell.

    - Null mean < POISSON_MEAN_MAX: the classical Poisson upper limit on the observed count with no
      background subtraction, 0.5 * chi2.ppf(0.95, 2 obs + 2) (3.00 for obs = 0). Conservative.
    - Otherwise: s95 = obs - q05(null), i.e. P(null + s <= obs) < 5 %, floored at the 1.645-sigma
      sensitivity of the ensemble when the observation fluctuates low."""
    from scipy import stats

    obs = np.asarray(obs, dtype=float)
    q05 = np.quantile(null, 0.05, axis=0)
    gauss = np.maximum(obs - q05, 1.645 * null.std(axis=0))
    pois = 0.5 * stats.chi2.ppf(0.95, 2 * obs + 2)
    return np.where(null.mean(axis=0) < POISSON_MEAN_MAX, pois, gauss)


# -----------------------------------------------------------------------------------------
# injection


def inject_pairs(
    a: Sample,
    b: Sample,
    n: int,
    lag_lo: float,
    lag_hi: float,
    p: Params,
    same: bool,
    rng: np.random.Generator,
) -> Sample:
    """Return a copy of ``b`` in which up to ``n`` events are moved to t_A + lag (lag log-uniform
    in [lag_lo, lag_hi], random sign) of distinct ``a`` anchors, keeping their own RA/Dec. The
    anchor is drawn among the ``a`` events within ``p.inject_local_days`` of the moved event, so
    events stay in their own observing era, and the pair must be wide (or unlocalized). These are
    synthetic dependent wide-separation pairs (provenance: simulated)."""
    mjd = b.mjd.copy()
    order = np.argsort(a.mjd)
    ta = a.mjd[order]
    anchors: set[int] = set()  # a events used as anchors (never moved when a is b)
    moved: set[int] = set()  # b events moved (never used as anchors when a is b)
    for j in rng.permutation(len(b.mjd)):
        if len(moved) >= n:
            break
        if same and j in anchors:
            continue
        lo = np.searchsorted(ta, b.mjd[j] - p.inject_local_days)
        hi = np.searchsorted(ta, b.mjd[j] + p.inject_local_days)
        for i in order[lo + rng.permutation(hi - lo)][:20]:
            if i in anchors or (same and (i == j or i in moved)):
                continue
            sep = sep_deg(a.ra[i], a.dec[i], b.ra[j], b.dec[j])
            cls = sep_class(np.atleast_1d(sep), a.sigma[i : i + 1], b.sigma[j : j + 1], p)[0]
            if cls in (1, 2):
                lag = np.exp(rng.uniform(np.log(max(lag_lo, 1e-3)), np.log(lag_hi)))
                mjd[j] = a.mjd[i] + rng.choice((-1.0, 1.0)) * lag / DAY
                anchors.add(int(i))
                moved.add(int(j))
                break
    return Sample(b.cat, mjd, b.ra, b.dec, b.sigma, b.year)

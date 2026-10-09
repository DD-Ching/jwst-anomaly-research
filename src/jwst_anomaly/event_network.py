"""E1 "causal event network": pair-count excess versus time lag between localized events.

Hypothesis E1 (owner idea 4; docs/hypotheses/round-1/summary.md "Correction to the causal-network idea",
B-on-A1.md section 7): time-stamped, localized events from *different* directions (separation well beyond the
localization errors) might be statistically dependent at time lags that no ordinary path explains. Dependence,
if found, is an anomaly, never evidence of anything exotic: observer-side common causes (instrument uptime and
exposure, Sun constraints, alert chains, shared triggers, catalogue duplicates, background events) come first.

The statistic is a cross-correlation function in lag: for each catalogue pair (A, B) the number of pairs in lag
bins, split by angular separation into ``same`` (sep <= n_sigma * sigma_comb) and ``wide`` (sep > n_sigma *
sigma_comb and > min_wide_deg). Events without a localization (GW from the GWTC CSV) enter an ``all`` class only.

Two nulls (both keep each event's declination and hour angle, i.e. sidereal scrambling, so Earth-fixed
exposure such as CHIME transit or IceCube zenith acceptance is kept):

- ``perm``: times are permuted among the events of the same catalogue and calendar year. The multiset of times
  is unchanged, so the lag histogram summed over separation is identical to the data: this null tests only the
  coupling of separation and lag (is the wide fraction at lag tau unusual?).
- ``boot``: each time is redrawn from the same catalogue's empirical times of that year (with replacement) and
  shifted by k * quantum, |k * quantum| <= jitter (catalogue-specific quantum: the Fermi orbit for GBM so that
  orbital phase, hence SAA and Earth-occultation uptime, is kept; the sidereal day for ground instruments).
  It keeps the long-term rate and uptime and tests temporal clustering itself, so it is the null for lag excess.

All thresholds are ASSUMPTIONs and live in :class:`Params`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy.table import Table

from .schema import Provenance

DAY = 86400.0
SIDEREAL_DAY = 86164.0905
#: Fermi orbital period (s), ~95.6 min (Meegan+2009, ApJ 702, 791). ASSUMPTION: constant over the mission.
FERMI_ORBIT = 5736.0
#: Catalogue codes used in tables and channel names.
CATS = ("GBM", "ICECAT", "GW", "CHIME")


@dataclass(frozen=True)
class Params:
    # --- lag bins (s); pre-registered by the task brief (ASSUMPTION) ---
    lag_edges: tuple[float, ...] = (0.0, 10.0, 100.0, 3600.0, DAY, 7 * DAY)
    # --- separation classes (ASSUMPTION) ---
    n_sigma: float = 3.0  # same-direction radius in units of the combined 1-sigma error
    min_wide_deg: float = 0.1  # a wide pair is also > 0.1 deg (the geometric limit of B-on-A1 section 7)
    # --- localization (ASSUMPTION) ---
    gbm_sys_deg: float = 3.7  # GBM systematic, 68 % core (Connaughton+2015, arXiv:1411.2685), as in D-071
    gbm_external_deg: float = 0.5  # error_radius below this = position from another instrument (0 in catalogue)
    gbm_external_sigma_deg: float = 0.05  # 1-sigma used for externally localized GBM bursts
    icecat_90_to_sigma: float = 2.146  # 2-D Gaussian: r90 = 2.146 sigma
    # --- boot null (ASSUMPTION) ---
    jitter_s: float = 3 * DAY
    quantum_s: dict = field(
        default_factory=lambda: {"GBM": FERMI_ORBIT, "ICECAT": SIDEREAL_DAY, "CHIME": SIDEREAL_DAY, "GW": 0.0}
    )
    # --- vetting windows (ASSUMPTION): half-widths around k * period ---
    orbit_k: tuple[int, ...] = tuple(range(1, 16))
    orbit_half_s: float = 120.0
    day_k: tuple[int, ...] = tuple(range(1, 8))
    day_half_s: float = 600.0


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

    return np.asarray(Time(np.asarray(mjd, dtype=float), format="mjd").datetime64.astype("datetime64[Y]")).astype(
        int
    ) + 1970


def gbm_events(cat: Table, p: Params = Params()) -> Table:
    """HEASARC fermigbrst rows -> events. 1-sigma = stat (+) systematic, or a tiny radius when error_radius
    is 0 (position from another instrument)."""
    er = np.asarray(cat["error_radius"], dtype=float)
    sig = np.where(er < p.gbm_external_deg, p.gbm_external_sigma_deg, np.hypot(er, p.gbm_sys_deg))
    return _events(cat["trigger_name"], "GBM", cat["trigger_time"], cat["ra"], cat["dec"], sig, "fermigbrst")


def icecat_events(df, p: Params = Params(), drop_cr_veto: bool = True) -> Table:
    """ICECAT-1 rows (pandas) -> events. 1-sigma = mean of the four 90 % errors / 2.146."""
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
    return _events(df["name"], "GW", mjd, np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan), "GWTC CSV")


def chime_events(df, one_per_source: bool = True) -> Table:
    """CHIME/FRB Catalog 2 rows (pandas) -> events. Sub-bursts (sub_num > 0) are dropped; with
    ``one_per_source`` each repeater contributes only its first burst (repeaters are one node)."""
    df = df[df["sub_num"] == 0].sort_values("mjd_400")
    if one_per_source:
        key = df["repeater_name"].fillna(df["tns_name"]).astype(str)
        df = df[~key.duplicated()]
    sig = np.maximum(df["ra_err"].to_numpy(float), df["dec_err"].to_numpy(float))
    return _events(df["tns_name"], "CHIME", df["mjd_400"], df["ra"], df["dec"], sig, "CHIME/FRB Catalog 2")


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
    """Index pairs (i, j) with |ta[i] - tb[j]| <= max_lag (times in MJD, any order). For ``same`` (ta is
    tb) only i < j in time order is returned, each unordered pair once."""
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
    """0 = same, 1 = wide, -1 = between (not wide only by the 0.1 deg floor), 2 = unlocalized (all)."""
    sc = np.hypot(sig_a, sig_b)
    cls = np.where(sep <= p.n_sigma * sc, 0, np.where(sep > np.maximum(p.n_sigma * sc, p.min_wide_deg), 1, -1))
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
        """Same events at new times; declination and hour angle kept (RA follows the sidereal angle)."""
        ra = np.mod(self.ra + gmst_deg(new_mjd) - gmst_deg(self.mjd), 360.0)
        return Sample(self.cat, new_mjd, ra, self.dec, self.sigma, self.year)


def windows_for(p: Params) -> list[tuple[str, list[tuple[float, float]]]]:
    """Lag windows: the pre-registered bins, then (vetting only) unions around orbital / day multiples."""
    e = p.lag_edges
    w = [(f"lag{i}", [(e[i], e[i + 1])]) for i in range(len(e) - 1)]
    w.append(("orbit_k", [(k * FERMI_ORBIT - p.orbit_half_s, k * FERMI_ORBIT + p.orbit_half_s) for k in p.orbit_k]))
    w.append(("sidereal_k", [(k * SIDEREAL_DAY - p.day_half_s, k * SIDEREAL_DAY + p.day_half_s) for k in p.day_k]))
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
    """Pair counts, shape (n_windows, 3) for classes (same, wide, all). ``all`` counts every pair."""
    windows = windows_for(p) if windows is None else windows
    max_lag = max(hi for _, iv in windows for _, hi in iv)
    i, j = pairs_within(a.mjd, b.mjd, max_lag, same)
    lag = np.abs(b.mjd[j] - a.mjd[i]) * DAY
    sep = sep_deg(a.ra[i], a.dec[i], b.ra[j], b.dec[j])
    cls = sep_class(sep, a.sigma[i], b.sigma[j], p)
    out = np.zeros((len(windows), 3), dtype=np.int64)
    for k, (_, ivs) in enumerate(windows):
        m = np.zeros(len(lag), bool)
        for lo, hi in ivs:
            m |= (lag >= lo) & (lag < hi) if lo > 0 else (lag <= hi)
        out[k, 0] = np.count_nonzero(m & (cls == 0))
        out[k, 1] = np.count_nonzero(m & (cls == 1))
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


def scramble_boot(s: Sample, rng: np.random.Generator, p: Params) -> Sample:
    """Redraw each time from the same-year empirical times (with replacement), shifted by k * quantum."""
    q = p.quantum_s.get(s.cat, 0.0)
    new = np.empty_like(s.mjd)
    for y in np.unique(s.year):
        idx = np.flatnonzero(s.year == y)
        new[idx] = s.mjd[rng.choice(idx, size=len(idx), replace=True)]
    if q > 0:
        kmax = int(p.jitter_s // q)
        shift = rng.integers(-kmax, kmax + 1, size=len(new)) * q
    else:
        shift = rng.uniform(-p.jitter_s, p.jitter_s, size=len(new))
    return s.with_times(new + shift / DAY)


def empirical_p(obs: np.ndarray, null: np.ndarray) -> np.ndarray:
    """Upper-tail p per cell: (1 + #(null >= obs)) / (N + 1). ``null`` has the scramble axis first."""
    return (1.0 + (null >= obs[None]).sum(axis=0)) / (null.shape[0] + 1.0)


def global_p(obs: np.ndarray, null: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    """Trials-corrected p: the observed minimum per-cell p over ``mask`` cells, compared with the same
    minimum computed for every scramble against the ensemble (look-elsewhere over all tested cells)."""
    nn = null.shape[0]
    flat = null[:, mask]
    o = obs[mask]
    p_obs = (1.0 + (flat >= o[None]).sum(axis=0)) / (nn + 1.0)
    # per-scramble p: fraction of the ensemble >= its value (ties count as >=)
    srt = np.sort(flat, axis=0)
    p_null = np.empty_like(flat, dtype=float)
    for c in range(flat.shape[1]):
        col = srt[:, c]
        p_null[:, c] = (nn - np.searchsorted(col, flat[:, c], side="left")) / nn
    min_obs = p_obs.min()
    min_null = p_null.min(axis=1)
    return float(min_obs), float((1.0 + np.count_nonzero(min_null <= min_obs)) / (nn + 1.0))


def z_score(obs, null) -> np.ndarray:
    mu, sd = null.mean(axis=0), null.std(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(sd > 0, (obs - mu) / sd, 0.0)


def upper_limit_95(obs: np.ndarray, null: np.ndarray) -> np.ndarray:
    """Neyman-style 95 % upper limit on additional pairs s: P(null + s <= obs) < 5 % -> s95 = obs - q05(null),
    floored at the 1.645-sigma sensitivity when the observation fluctuates low."""
    q05 = np.quantile(null, 0.05, axis=0)
    floor = 1.645 * null.std(axis=0)
    return np.maximum(obs - q05, floor)


# ----------------------------------------------------------------------------------------- injection


def inject_pairs(
    a: Sample, b: Sample, n: int, lag_lo: float, lag_hi: float, p: Params, same: bool, rng: np.random.Generator
) -> Sample:
    """Return a copy of ``b`` in which ``n`` events are moved to t_A + lag (lag log-uniform in [lag_lo,
    lag_hi], random sign) of distinct random ``a`` events, keeping their own RA/Dec, accepted only when the
    pair is wide (or unlocalized). These are synthetic dependent wide-separation pairs (provenance: simulated)."""
    mjd = b.mjd.copy()
    ib = rng.permutation(len(b.mjd))
    ia = rng.permutation(len(a.mjd))
    done, k = 0, 0
    used_a: set[int] = set()
    for j in ib:
        if done >= n or k >= len(ia):
            break
        while k < len(ia):
            i = ia[k]
            k += 1
            if same and (i == j or j in used_a):
                continue
            sep = sep_deg(a.ra[i], a.dec[i], b.ra[j], b.dec[j])
            if sep_class(np.atleast_1d(sep), a.sigma[i : i + 1], b.sigma[j : j + 1], p)[0] in (1, 2):
                lag = np.exp(rng.uniform(np.log(max(lag_lo, 1e-3)), np.log(lag_hi)))
                mjd[j] = a.mjd[i] + rng.choice((-1.0, 1.0)) * lag / DAY
                used_a.add(i)
                done += 1
                break
    return Sample(b.cat, mjd, b.ra, b.dec, b.sigma, b.year)

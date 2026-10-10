"""Neutrino Frontier helpers (owner brief 2026-10-10, D-077): ghost-pair statistics for NF-H04.

NF-H04 (hypothesis): a fraction R_g of astrophysical neutrino tracks has a delayed copy at an
unrelated sky position. These helpers count wide-separation pairs by lag (plain and
signalness-weighted), scramble times with a cyclic jitter for long lags, and inject synthetic ghosts
(provenance: simulated). Pair geometry and the short-lag nulls are reused from
:mod:`jwst_anomaly.event_network` (D-074).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import event_network as en

DAY = en.DAY


@dataclass(frozen=True)
class GhostParams:
    # lag bins (s): D-074's five bins plus 7-30 d and 30-180 d (ASSUMPTION, N-B round-1 review)
    lag_edges: tuple[float, ...] = (0.0, 10.0, 100.0, 3600.0, DAY, 7 * DAY, 30 * DAY, 180 * DAY)
    # bins whose null is the cyclic long jitter instead of D-074 jit (lag >= 7 d)
    long_from_bin: int = 5
    # cyclic jitter half-width for the long bins (ASSUMPTION: alert rate stationary on ~1 yr)
    long_jitter_days: float = 365.0


def n_bins(gp: GhostParams) -> int:
    return len(gp.lag_edges) - 1


def wide_pair_stats(
    s: en.Sample,
    w: np.ndarray,
    gp: GhostParams,
    p: en.Params | None = None,
    group: np.ndarray | None = None,
) -> np.ndarray:
    """(n_bins, 2): wide-pair count and sum of w_i w_j per lag bin, unordered pairs of one catalogue
    (bin 0 includes lag 0). Wide = D-074 separation class 1. Pairs with equal ``group`` (one IceCube
    readout) are not counted."""
    p = p or en.Params()
    e = np.asarray(gp.lag_edges)
    i, j = en.pairs_within(s.mjd, s.mjd, e[-1], True)
    lag = np.abs(s.mjd[j] - s.mjd[i]) * DAY
    sep = en.sep_deg(s.ra[i], s.dec[i], s.ra[j], s.dec[j])
    wide = (en.sep_class(sep, s.sigma[i], s.sigma[j], p) == 1) & (lag <= e[-1])
    if group is not None:
        wide &= group[i] != group[j]
    k = np.clip(np.searchsorted(e, lag[wide], side="right") - 1, 0, n_bins(gp) - 1)
    out = np.zeros((n_bins(gp), 2))
    out[:, 0] = np.bincount(k, minlength=n_bins(gp))
    out[:, 1] = np.bincount(k, weights=(w[i] * w[j])[wide], minlength=n_bins(gp))
    return out


def scramble_cyclic(s: en.Sample, rng: np.random.Generator, half_width_days: float) -> en.Sample:
    """Independent uniform time shift per event within +-half_width, wrapped cyclically inside the
    catalogue span (no edge loss); Dec and hour angle kept (``Sample.with_times``)."""
    t0, t1 = s.mjd.min(), s.mjd.max()
    span = t1 - t0
    new = t0 + np.mod(s.mjd - t0 + rng.uniform(-half_width_days, half_width_days, len(s.mjd)), span)
    return s.with_times(new)


def inject_ghosts(
    s: en.Sample,
    w: np.ndarray,
    r_g: float,
    lag_lo_s: float,
    lag_hi_s: float,
    rng: np.random.Generator,
) -> tuple[en.Sample, np.ndarray]:
    """Copy of ``s`` plus ghosts (provenance: simulated): each event is a parent with probability
    r_g * w (its signalness); its ghost comes Delta t later (log-uniform in the bin; a ghost beyond
    the catalogue end is dropped), at RA uniform and Dec / error drawn from the catalogue, with the
    parent's weight."""
    parent = rng.random(len(s.mjd)) < np.clip(r_g * w, 0, 1)
    k = int(parent.sum())
    lag = np.exp(rng.uniform(np.log(max(lag_lo_s, 1e-3)), np.log(lag_hi_s), k)) / DAY
    t = s.mjd[parent] + lag
    keep = t <= s.mjd.max()
    t = t[keep]
    n = len(t)
    draw = rng.integers(0, len(s.mjd), n)
    ra = rng.uniform(0, 360, n)
    sample = en.Sample(
        s.cat,
        np.concatenate([s.mjd, t]),
        np.concatenate([s.ra, ra]),
        np.concatenate([s.dec, s.dec[draw]]),
        np.concatenate([s.sigma, s.sigma[draw]]),  # Dec and error from the same event
        np.concatenate([s.year, en.mjd_year(t)]),
    )
    return sample, np.concatenate([w, w[parent][keep]])


# --- IceTracks-DR2 (E-NF1b, D-079) -------------------------------------------------------------

#: columns of the IceTracks-DR2 ``events/<season>_exp.tab`` files (doi:10.7910/DVN/MMIIZA)
DR2_COLUMNS = (
    "run",
    "event",
    "subevent",
    "mjd",
    "log10e",
    "angerr",
    "ra",
    "dec",
    "azimuth",
    "zenith",
)


def read_icetracks_tab(lines) -> np.ndarray:
    """Float array from the lines of one IceTracks-DR2 events or uptime ``.tab`` file as Dataverse
    serves it (a ``#`` header line, then whitespace columns, each row in double quotes)."""
    rows = [
        ln.replace('"', "").split()
        for ln in lines
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    return np.array(rows, dtype=float)


def dedupe_events(ev: np.ndarray) -> np.ndarray:
    """Time-sorted events, one per (run, event, subevent): the seasons overlap by weeks."""
    key = np.round(ev[:, :3]).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    ev = ev[np.sort(idx)]
    return ev[np.argsort(ev[:, 3], kind="stable")]


def in_uptime(mjd: np.ndarray, start: np.ndarray, stop: np.ndarray) -> np.ndarray:
    """True where ``mjd`` lies in a good-run interval ``[start, stop]`` (intervals may overlap)."""
    o = np.argsort(start)
    s, e = start[o], np.maximum.accumulate(stop[o])
    k = np.searchsorted(s, mjd, side="right") - 1
    return (k >= 0) & (mjd <= e[np.clip(k, 0, None)])


def jitter_uptime(
    s: en.Sample,
    start: np.ndarray,
    stop: np.ndarray,
    half_width_days: float,
    rng: np.random.Generator,
    max_tries: int = 100,
) -> en.Sample:
    """Uptime-aware jitter null: an independent uniform shift of +-half_width per event, redrawn
    until the new time is inside a good run (an event that never lands keeps its time); Dec and hour
    angle kept (``Sample.with_times``)."""
    new = s.mjd.copy()
    todo = np.arange(len(new))
    for _ in range(max_tries):
        if not len(todo):
            break
        cand = s.mjd[todo] + rng.uniform(-half_width_days, half_width_days, len(todo))
        ok = in_uptime(cand, start, stop)
        new[todo[ok]] = cand[ok]
        todo = todo[~ok]
    return s.with_times(new)


def wide_pair_counts(
    s: en.Sample, group: np.ndarray, lag_edges, p: en.Params | None = None
) -> np.ndarray:
    """Wide-pair counts per lag bin (``wide_pair_stats`` column 0); pairs within one ``group``
    (one IceCube (run, event): split or coincident muons of one readout) are not counted."""
    gp = GhostParams(lag_edges=tuple(lag_edges))
    return wide_pair_stats(s, np.ones(len(s.mjd)), gp, p, group)[:, 0]


def union_days(start: np.ndarray, stop: np.ndarray) -> float:
    """Total length of the union of the intervals ``[start, stop]`` (overlapping good runs once)."""
    o = np.argsort(start)
    s, e = start[o], stop[o]
    total, cur_s, cur_e = 0.0, s[0], e[0]
    for a, b in zip(s[1:], e[1:], strict=True):
        if a > cur_e:
            total += cur_e - cur_s
            cur_s, cur_e = a, b
        else:
            cur_e = max(cur_e, b)
    return float(total + cur_e - cur_s)

"""S1 "burst twins": wide-separation, any-delay light-curve matches between Fermi GBM bursts.

Hypothesis S1 (docs/hypotheses/round-1/summary.md; B-on-A2.md section 5): two bursts far apart on the sky whose
light curves match after a time shift. This module holds the numerical pieces of the screen; the driver is
``scripts/s1_twins.py``. A high match score is an anomaly, never evidence of new physics.

Inputs are the GBM burst-catalogue "bcat" products (``glg_bcat_all_bn*.fit``): rmfit's time-resolved,
background-subtracted, deconvolved photon-flux spectra per detector (HDU 1 ``PHTCNTS``/``PHTERRS``, time-major,
8 CTIME channels) on the time bins of HDU 2 ``TIMEBIN``. They are reduced to three matched energy bands.

Thresholds are ASSUMPTIONs and live in :class:`Params`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal, stats

#: Matched energy bands (keV). A detector channel belongs to the band containing its centre energy, so 8-channel
#: CTIME and 128-channel TTE bcat files map the same way (CTIME: channels 1-2 / 3-4 / 5-6). ASSUMPTION.
BANDS_KEV = ((10.0, 50.0), (50.0, 300.0), (300.0, 1000.0))
BAND_NAMES = ("b1_10_50", "b2_50_300", "b3_300_1000")


@dataclass(frozen=True)
class Params:
    # --- light-curve reduction (ASSUMPTIONs) ---
    base_dt: float = 0.064  # s; bcat TTE-binned resolution
    max_bins: int = 256  # per stored light curve; dt = base_dt * 2**k with the smallest k that fits
    pad_frac: float = 0.25  # window = [t90_start - pad, t90_start + t90 + pad], pad = pad_frac*T90 + pad_s
    pad_s: float = 2.0
    bad_abs: float = 1e6  # |PHTCNTS| above this (or non-finite) is a failed bin
    # --- pulses (ASSUMPTION): a pulse is a peak of the summed-band S/N curve with this prominence ---
    pulse_prominence_sigma: float = 5.0
    pulse_min_sep_bins: int = 2
    min_pulses: int = 2
    # --- positions (ASSUMPTION): GBM systematic 3.7 deg 68 % core (Connaughton+2015, arXiv:1411.2685) ---
    sys_deg: float = 3.7
    sys_tail_deg: float = 14.0  # ~10 % of bursts; reported as a secondary, stricter cut
    pos_nsigma: float = 3.0
    # --- matching ---
    max_dk: int = 1  # compare only pairs whose stored resolutions differ by <= 2**max_dk
    stretch_grid: tuple[float, ...] = (0.5, 0.71, 1.41, 2.0)  # secondary s grid (s = 1 is primary)
    # --- vetting (ASSUMPTIONs) ---
    chi2_p_min: float = 1e-3  # a twin must be noise-consistent with a scaled, shifted copy
    retrigger_days: float = 1.0  # pairs closer in time than this are duplicate / re-trigger candidates


# ----------------------------------------------------------------------------------------------- reduction


def window_for(
    t90_start: float, t90: float, p: Params = Params(), min_dt: float | None = None
) -> tuple[float, float, float]:
    """Return (t_lo, t_hi, dt) of the stored window relative to the trigger.

    dt = base_dt * 2**k with the smallest k such that the window fits in max_bins and dt >= min_dt (the native
    bin width of the source data inside the window, if given).
    """
    t90 = max(float(t90), p.base_dt)
    pad = p.pad_frac * t90 + p.pad_s
    lo, hi = t90_start - pad, t90_start + t90 + pad
    k = max(0, int(np.ceil(np.log2((hi - lo) / (p.max_bins * p.base_dt)))))
    if min_dt is not None and min_dt > p.base_dt:
        k = max(k, int(np.ceil(np.log2(min_dt / p.base_dt) - 0.05)))
    dt = p.base_dt * 2**k
    n = int(np.ceil((hi - lo) / dt))
    return lo, lo + n * dt, dt


def rebin(tb: np.ndarray, flux: np.ndarray, var: np.ndarray, lo: float, hi: float, dt: float):
    """Rebin rates on bins ``tb`` (n, 2) onto the regular grid [lo, hi) of width dt by exact overlap.

    flux_k = sum_i o_ik f_i / c_k and var_k = sum_i o_ik^2 var_i / c_k^2, with o_ik the overlap (s) of source
    bin i with target bin k and c_k = sum_i o_ik. Target bins with coverage < 50 % are NaN.
    """
    n = int(round((hi - lo) / dt))
    a, b = tb[:, 0], tb[:, 1]
    keep = (b > lo) & (a < hi) & np.isfinite(flux) & np.isfinite(var)
    a, b, flux, var = a[keep], b[keep], flux[keep], var[keep]
    cov = np.zeros(n)
    f = np.zeros(n)
    v = np.zeros(n)
    if a.size:
        k0 = np.clip(np.floor((a - lo) / dt).astype(int), 0, n - 1)
        k1 = np.clip(np.floor((b - lo) / dt).astype(int), 0, n - 1)
        for j in range(int((k1 - k0).max()) + 1):
            k = k0 + j
            ok = k <= k1
            o = np.minimum(b, lo + (k + 1) * dt) - np.maximum(a, lo + k * dt)
            ok &= o > 0
            np.add.at(cov, k[ok], o[ok])
            np.add.at(f, k[ok], (o * flux)[ok])
            np.add.at(v, k[ok], (o * o * var)[ok])
    with np.errstate(invalid="ignore", divide="ignore"):
        good = cov >= 0.5 * dt
        fo = np.where(good, f / cov, np.nan)
        eo = np.where(good, np.sqrt(v) / cov, np.nan)
    return fo, eo


def reduce_bcat(hdul, t90_start: float, t90: float, p: Params = Params()) -> dict:
    """Reduce an open bcat HDUList to three-band photon-flux light curves (ph cm^-2 s^-1) in the S1 window.

    Detectors flagged INCLUDED are combined by inverse-variance weighting per bin and band.
    """
    tb = np.asarray(hdul[2].data["TIMEBIN"], float)
    n = len(tb)
    lo0, hi0, _ = window_for(t90_start, t90, p)
    inwin = (tb[:, 1] > lo0) & (tb[:, 0] < hi0)
    native = float(np.median(tb[inwin, 1] - tb[inwin, 0])) if inwin.any() else p.base_dt
    lo, hi, dt = window_for(t90_start, t90, p, min_dt=native)
    num = np.zeros((3, n))
    den = np.zeros((3, n))
    dets = []
    for row in hdul[1].data:
        if str(row["DETSTAT"]).strip().upper() != "INCLUDED" or not str(row["DETNAM"]).startswith("NAI"):
            continue
        e = np.asarray(row["E_EDGES"], float)
        nch = e.size - 1
        if nch < 3 or np.asarray(row["PHTCNTS"]).size != n * nch:
            continue
        dE = np.diff(e)
        centre = np.sqrt(e[:-1] * np.maximum(e[1:], 1e-3))
        cts = np.asarray(row["PHTCNTS"], float).reshape(n, nch)
        err = np.asarray(row["PHTERRS"], float).reshape(n, nch)
        bad = ~np.isfinite(cts) | ~np.isfinite(err) | (np.abs(cts) > p.bad_abs) | (err <= 0) | (err > p.bad_abs)
        cts = np.where(bad, 0.0, cts)
        err = np.where(bad, np.inf, err)
        dets.append(str(row["DETNAM"]).strip())
        for b, (elo, ehi) in enumerate(BANDS_KEV):
            ch = np.flatnonzero((centre >= elo) & (centre < ehi))
            f = (cts[:, ch] * dE[ch]).sum(1)
            v = ((err[:, ch] * dE[ch]) ** 2).sum(1)
            with np.errstate(invalid="ignore", divide="ignore"):
                wgt = np.where(np.isfinite(v) & (v > 0), 1.0 / v, 0.0)
            num[b] += wgt * f
            den[b] += wgt
    if not dets:
        raise ValueError("no usable INCLUDED NaI detector")
    flux = np.full((3, int(round((hi - lo) / dt))), np.nan)
    errs = np.full_like(flux, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        fb = np.where(den > 0, num / den, np.nan)
        vb = np.where(den > 0, 1.0 / den, np.nan)
    for b in range(3):
        flux[b], errs[b] = rebin(tb, fb[b], vb[b], lo, hi, dt)
    return {"t_lo": lo, "dt": dt, "flux": flux, "err": errs, "dets": ",".join(dets)}


# ------------------------------------------------------------------------------------- compact encoding


def encode(flux: np.ndarray, err: np.ndarray) -> tuple[list[float], list[str], list[str]]:
    """Encode (3, N) flux/err as integer tenths of a per-band unit (the 20th-percentile error).

    Missing bins are written as ``nan``. Returns (units, flux_strings, err_strings).
    """
    units, fs, es = [], [], []
    for b in range(flux.shape[0]):
        e = err[b][np.isfinite(err[b])]
        u = float(np.percentile(e, 20)) if e.size else 1.0
        u = u if u > 0 else 1.0
        units.append(float(f"{u:.5g}"))

        def q(a, u=u):
            return " ".join("nan" if not np.isfinite(x) else str(int(round(10 * x / u))) for x in a)

        fs.append(q(flux[b]))
        es.append(q(err[b]))
    return units, fs, es


def decode(units, fs, es) -> tuple[np.ndarray, np.ndarray]:
    flux = np.array([np.array(s.split(), float) * u / 10 for u, s in zip(units, fs, strict=True)])
    err = np.array([np.array(s.split(), float) * u / 10 for u, s in zip(units, es, strict=True)])
    return flux, err


def clean(flux: np.ndarray, err: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Missing bins -> flux 0 with the band's median error (keeps arrays aligned)."""
    f = np.where(np.isfinite(flux), flux, 0.0)
    med = np.nanmedian(np.where(np.isfinite(err), err, np.nan), axis=1, keepdims=True)
    med = np.where(np.isfinite(med), med, 1.0)
    e = np.where(np.isfinite(err) & np.isfinite(flux), err, med)
    return f, e


# ----------------------------------------------------------------------------------------------- pulses


def find_pulses(flux: np.ndarray, err: np.ndarray, p: Params = Params()) -> np.ndarray:
    """Indices of resolved pulses: peaks of the summed-band S/N curve with prominence >= N sigma."""
    tot = flux.sum(0)
    sig = np.sqrt((err**2).sum(0))
    unit = np.median(sig) if np.all(np.isfinite(sig)) else np.nanmedian(sig)
    snr = tot / unit
    peaks, _ = signal.find_peaks(snr, prominence=p.pulse_prominence_sigma, distance=p.pulse_min_sep_bins)
    return peaks


def pulse_segments(flux: np.ndarray, err: np.ndarray, p: Params = Params()) -> list[tuple[int, int]]:
    """Split the window at the minima between consecutive pulses; returns [start, stop) per pulse segment.

    The first segment starts at 0 and the last ends at N, so the segments tile the window.
    """
    peaks = find_pulses(flux, err, p)
    n = flux.shape[1]
    if len(peaks) < 2:
        return [(0, n)]
    tot = flux.sum(0)
    cuts = [0]
    for a, b in zip(peaks[:-1], peaks[1:], strict=True):
        cuts.append(int(a + np.argmin(tot[a : b + 1])))
    cuts.append(n)
    return [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)]


def pulse_shuffle(flux: np.ndarray, err: np.ndarray, rng: np.random.Generator, p: Params = Params()):
    """Surrogate light curve: the burst's pulse segments in a random order (all bands moved together).

    Keeps the pulse count, each pulse's shape and spectrum, and the total power; destroys the pulse order.
    For a burst whose order is unchanged by the draw, a random different permutation is used when possible.
    """
    seg = pulse_segments(flux, err, p)
    if len(seg) < 2:
        return flux.copy(), err.copy()
    order = rng.permutation(len(seg))
    if np.all(order == np.arange(len(seg))):
        order = np.roll(order, 1)
    fo = np.concatenate([flux[:, a:b] for a, b in (seg[i] for i in order)], axis=1)
    eo = np.concatenate([err[:, a:b] for a, b in (seg[i] for i in order)], axis=1)
    return fo, eo


# --------------------------------------------------------------------------------------------- positions


def angsep_deg(ra1, dec1, ra2, dec2):
    """Great-circle separation in degrees (vectorised; haversine form)."""
    r1, d1, r2, d2 = (np.radians(np.asarray(x, float)) for x in (ra1, dec1, ra2, dec2))
    s = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(s, 0, 1))))


def position_inconsistent(sep, err1, err2, sys_deg: float = Params.sys_deg, nsig: float = Params.pos_nsigma):
    """True where sep > nsig * sqrt(err1^2 + err2^2 + 2 sys^2): the two bursts cannot share a position.

    ``err`` is the catalogue 1-sigma statistical radius (deg). The GBM systematic applies to each burst, so it
    enters the pair variance twice. Bursts localised by another instrument (err = 0) still get the systematic,
    which only makes the cut stricter.
    """
    sigma = np.sqrt(np.asarray(err1, float) ** 2 + np.asarray(err2, float) ** 2 + 2 * sys_deg**2)
    return np.asarray(sep, float) > nsig * sigma


# --------------------------------------------------------------------------------------------- statistic


def rebin_factor(flux: np.ndarray, err: np.ndarray, m: int):
    """Sum ``m`` adjacent bins (rates averaged); trailing partial bin dropped."""
    if m == 1:
        return flux, err
    n = flux.shape[1] // m
    f = flux[:, : n * m].reshape(flux.shape[0], n, m).mean(2)
    e = np.sqrt((err[:, : n * m] ** 2).reshape(err.shape[0], n, m).sum(2)) / m
    return f, e


def xcorr_max(x: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    """Maximum over lag of the multi-band normalised cross-correlation of (B, N) and (B, M) arrays.

    rho(tau) = sum_b sum_t x_b(t) y_b(t + tau) / (||x|| ||y||), with zero padding; rho = 1 for y a scaled,
    shifted copy of x in every band (matched bands; band ratios must agree). Returns (rho_max, lag) with lag
    in bins of y relative to x (positive: y later).
    """
    nx, ny = x.shape[1], y.shape[1]
    L = 1 << int(np.ceil(np.log2(nx + ny)))
    fx = np.fft.rfft(x, L, axis=1)
    fy = np.fft.rfft(y, L, axis=1)
    cc = np.fft.irfft((np.conj(fx) * fy).sum(0), L)
    norm = np.sqrt((x * x).sum() * (y * y).sum())
    if norm <= 0:
        return 0.0, 0
    k = int(np.argmax(cc))
    lag = k if k < L - nx else k - L
    return float(cc[k] / norm), lag


def stretch(flux: np.ndarray, err: np.ndarray, s: float):
    """Time-stretch a (B, N) light curve by factor s (linear interpolation on bin centres; same dt)."""
    n = flux.shape[1]
    m = max(2, int(round(n * s)))
    t_new = (np.arange(m) + 0.5) / s - 0.5
    f = np.array([np.interp(t_new, np.arange(n), fb, left=0, right=0) for fb in flux])
    e = np.array([np.interp(t_new, np.arange(n), eb) for eb in err]) / np.sqrt(max(s, 1e-9))
    return f, e


def twin_chi2(x, ex, y, ey, lag: int) -> tuple[float, int, float]:
    """Noise-consistency of y with a scaled copy of x at ``lag``: chi2, dof and p-value.

    The scale a minimises sum (y - a x)^2 / (ey^2 + a^2 ex^2) (iterated twice). Only bins where either curve
    exceeds 2 sigma in the summed band enter, so long background stretches do not dilute the test.
    """
    nx, ny = x.shape[1], y.shape[1]
    lo, hi = max(0, lag), min(ny, nx + lag)
    if hi - lo < 3:
        return np.inf, 0, 0.0
    ys, eys = y[:, lo:hi], ey[:, lo:hi]
    xs, exs = x[:, lo - lag : hi - lag], ex[:, lo - lag : hi - lag]
    sel = (xs.sum(0) > 2 * np.sqrt((exs**2).sum(0))) | (ys.sum(0) > 2 * np.sqrt((eys**2).sum(0)))
    if sel.sum() < 3:
        return np.inf, 0, 0.0
    xs, exs, ys, eys = xs[:, sel], exs[:, sel], ys[:, sel], eys[:, sel]
    a = (xs * ys).sum() / max((xs * xs).sum(), 1e-30)
    for _ in range(2):
        w = 1.0 / (eys**2 + a * a * exs**2)
        a = (w * xs * ys).sum() / max((w * xs * xs).sum(), 1e-30)
    w = 1.0 / (eys**2 + a * a * exs**2)
    chi2 = float((w * (ys - a * xs) ** 2).sum())
    dof = int(xs.size - 1)
    return chi2, dof, float(stats.chi2.sf(chi2, dof))


# --------------------------------------------------------------------------------------------- injection


def inject_twin(
    flux_a: np.ndarray,
    err_a: np.ndarray,
    err_b: np.ndarray,
    ratio: float,
    rng: np.random.Generator,
):
    """A synthetic twin of burst A placed in burst B's slot at flux ratio ``ratio`` (<= 1).

    The copy is ratio * A's light curve (A's own noise scaled down) plus Gaussian noise that restores B's
    background level and the Poisson source variance of a fainter copy:
      var_add = max(u_B^2 - ratio^2 u_A^2, 0) + ratio (1 - ratio) max(err_A^2 - u_A^2, 0),
    with u the per-band 20th-percentile (background) error. Returns (flux, err) on A's grid (simulated).
    """
    u_a = np.percentile(err_a, 20, axis=1, keepdims=True)
    u_b = np.percentile(err_b, 20, axis=1, keepdims=True)
    var_add = np.maximum(u_b**2 - ratio**2 * u_a**2, 0) + ratio * (1 - ratio) * np.maximum(err_a**2 - u_a**2, 0)
    f = ratio * flux_a + rng.normal(size=flux_a.shape) * np.sqrt(var_add)
    e = np.sqrt(ratio**2 * err_a**2 + var_add)
    return f, e


def snr_total(flux: np.ndarray, err: np.ndarray) -> float:
    """Peak S/N of the summed-band curve (in units of its median error)."""
    tot = flux.sum(0)
    sig = np.sqrt((err**2).sum(0))
    return float(tot.max() / np.median(sig))


def poisson_upper_limit(k: int, cl: float = 0.95) -> float:
    """Classical one-sided Poisson upper limit on the mean for k observed events (k=0 -> 3.00 at 95 %)."""
    return float(stats.chi2.ppf(cl, 2 * (k + 1)) / 2)

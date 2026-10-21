"""S1 burst twins on CHIME/FRB Catalog 2: adapter and frozen rules (validation of <= 30 events).

This is the thin adapter of the S1 line to CHIME dynamic spectra. It turns one Catalog 2 Stokes-I
dynamic spectrum (``<TNS>_stokesi_dynamic_spectrum.h5``, doi:10.11570/25.0066; format facts from
issue #153 and arXiv:2601.09399 section 3.3) into four matched 100-MHz sub-band light curves in Jy
with per-band errors, and then reuses ``jwst_anomaly.burst_twins`` for the pulse count
(``find_pulses``), the statistic (``xcorr_max``), the surrogate (``pulse_shuffle``) and the D-071
chi2 (``twin_chi2``, secondary). The only new rule is the issue-#153 fix of the chi2 test (union of
both emission regions, one non-negative gain per sub-band that must be detected where either member
has signal). Frozen plan: ``results/s1_chime/PLAN.md``. Nothing here is a search, and a match is an
anomaly, never evidence of new physics.

Provenance labels: dynamic spectra and catalogue values are observed; sub-band profiles, rho, chi2
and gains are derived; injected copies are simulated; every threshold is an ASSUMPTION in
:class:`Params`.
"""

from __future__ import annotations

import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from jwst_anomaly import burst_twins as bt  # noqa: E402

#: CHIME/FRB time resolution (s), from the file attribute ``res_time`` (0.98304 ms).
RES_TIME_S = 0.00098304
#: Native channel width (MHz) of the Catalog 2 dynamic spectra (``res_freq``).
RES_FREQ_MHZ = 0.0244140625


@dataclass(frozen=True)
class Params:
    # --- reduction (ASSUMPTIONs) ---
    # four matched sub-bands (MHz, [lo, hi)); the top edge includes the last channel (800.2 MHz)
    bands_mhz: tuple[tuple[float, float], ...] = (
        (400.0, 500.0),
        (500.0, 600.0),
        (600.0, 700.0),
        (700.0, 800.3),
    )
    min_valid_chan_frac: float = 0.10  # a sub-band with fewer valid channels is "missing"
    pad_bins: int = 10  # crop = union of the file's pulse_emission_region +- pad (issue #153)
    min_off_bins: int = 20  # off-pulse bins needed for the noise estimate
    # --- pulses: D-071 eligibility via burst_twins.find_pulses (5 sigma prominence, >= 2) ---
    min_pulses: int = 2
    # --- prescreen: rho* = this quantile of the per-catalogue maximum over surrogates (D-071) ---
    n_surrogates: int = 100
    rho_star_q: float = 0.95
    # --- chi2 twin test (issue-#153 fixed rule; ASSUMPTIONs) ---
    chi2_p_min: float = 1e-3  # as D-071
    chi2_sys_frac: float = 0.10  # as D-071
    band_signal_snr: float = 5.0  # a band "has signal" if either member peaks >= this in the window
    gain_det_sigma: float = 3.0  # each band with signal needs a_b / sigma(a_b) >= this
    # sub-bin lag refinement: A is Fourier-shifted by delta in this grid (bins) and the best chi2 is
    # kept (1 extra fitted parameter). FRB pulses are ~1 bin wide, so an integer lag alone fails a
    # bright exact copy at a 0.3-bin offset (found on the synthetic test before any data were read).
    subbin_grid: tuple[float, ...] = tuple(np.round(np.linspace(-0.5, 0.5, 11), 2))
    # --- injections (simulated; ASSUMPTIONs) ---
    ratios: tuple[float, ...] = (1.0, 0.5, 0.3)
    shift_max_bins: float = 10.0  # copy offset from the host's region centre, U(-s, s) (bins)
    beam_u_max: float = 0.5  # G_beam: radial offset / FWHM(600 MHz), U(0, u_max) per member
    path_tau600_ms: float = 1.0  # P_path: extra scattering time at 600 MHz, nu^-4
    path_ddm: float = 300.0  # P_path: DM difference (pc cm^-3) -> intra-channel smearing, nu^-3
    # --- transferability verdict (ASSUMPTIONs, frozen in PLAN.md) ---
    t1_eff_min: float = 0.5  # recovery at ratio 1, G_cal set
    t2_chance_pass_max: float = 0.01  # fraction of different-source pairs passing the chain
    t3_quantile: float = (
        0.95  # median injected rho (ratio 1, G_cal) must exceed this chance quantile
    )
    t4_repeater_pass_max: float = 0.05  # same-source pass fraction above this = "confused"
    min_eligible: int = 10  # fewer eligible bursts: "not testable"


DEFAULT = Params()
GAIN_SETS = ("G0", "G_cal", "G_beam", "P_path")


def band_centres(p: Params = DEFAULT) -> np.ndarray:
    return np.array([(a + min(b, 800.0)) / 2 for a, b in p.bands_mhz])


# -------------------------------------------------------------------------------- selection


def select_sample(cat, files: set[str], reps: int = 3, per_rep: int = 5, nonrep: int = 15):
    """Frozen sample (non-outcome criteria only), from the Catalog 2 table (pandas DataFrame).

    Per event: n_comp = number of table rows (fitburst sub-bursts), snr = max snr_fitb over rows.
    Clean = excluded_flag == 0 and sidelobe_flag == 0 and n_comp >= 2 and a dynamic-spectrum file
    in the listing. Repeaters: the ``reps`` repeaters with the most clean events (ties: name), each
    with its ``per_rep`` highest-snr clean events. Non-repeaters: the ``nonrep`` highest-snr clean
    events without a repeater_name. Ties on snr break by TNS name.
    """
    import pandas as pd

    g = cat.groupby("tns_name").agg(
        n_comp=("sub_num", "size"),
        repeater_name=("repeater_name", "first"),
        snr=("snr_fitb", "max"),
        excluded_flag=("excluded_flag", "max"),
        sidelobe_flag=("sidelobe_flag", "max"),
    )
    g["has_file"] = [f"{n}_stokesi_dynamic_spectrum.h5" in files for n in g.index]
    clean = g[(g.n_comp >= 2) & (g.excluded_flag == 0) & (g.sidelobe_flag == 0) & g.has_file]
    clean = clean.reset_index().sort_values(["snr", "tns_name"], ascending=[False, True])
    rep = clean[clean.repeater_name.notna() & (clean.repeater_name != "")]
    counts = rep.groupby("repeater_name").size().reset_index(name="n")
    counts = counts.sort_values(["n", "repeater_name"], ascending=[False, True])
    chosen = []
    for name in counts.repeater_name.iloc[:reps]:
        chosen.append(rep[rep.repeater_name == name].iloc[:per_rep].assign(group="repeater"))
    non = clean[clean.repeater_name.isna() | (clean.repeater_name == "")]
    chosen.append(non.iloc[:nonrep].assign(group="non_repeater"))
    out = pd.concat(chosen, ignore_index=True)
    out["source"] = [
        r if isinstance(r, str) and r else n
        for r, n in zip(out.repeater_name, out.tns_name, strict=True)
    ]
    return out


# -------------------------------------------------------------------------------- reduction


def reduce_arrays(
    data,
    model,
    flag,
    good,
    cal_good,
    fcf,
    cal_spec,
    freqs,
    region,
    p: Params = DEFAULT,
    chan_sigma=None,
):
    """Reduce one dynamic spectrum to sub-band light curves (Jy) with errors.

    ``data``/``model``/``flag`` are (n_freq, n_time); ``good``, ``cal_good``, ``fcf`` (Jy per data
    unit), ``cal_spec`` and ``freqs`` (MHz) are (n_freq,); ``region`` is the file's
    ``pulse_emission_region`` ([[t0, t1], ...], bins). Per band b and bin t:
    S_b(t) = K_b * mean over valid channels of data_f(t) / sigma_f (equal weight per channel in S/N
    units), with sigma_f the channel noise in data units (``chan_sigma``; if None, the 1.4826 MAD of
    the channel's off-pulse data in this window) and K_b = median over the band's valid channels of
    fcf * sigma_f (Jy per unit S/N); valid = good & cal_good & finite positive fcf & finite positive
    sigma_f & flag. Fixes before any outcome: v1 took the unweighted mean of fcf * data (a few
    channels with a near-zero calibration response, huge fcf, dominated the noise); v2 weighted by
    1 / (fcf sigma_f)^2 (a few channels with a tiny fcf then carried ~45 % of the weight and the
    band S/N of a 230-sigma burst fell to ~2).
    The off-pulse median (outside the crop) is subtracted from the data (not from the model); the
    error sigma_b is 1.4826 * MAD of the data over the off-pulse bins (radiometer noise; ASSUMPTION:
    no self-noise). A band with < min_valid_chan_frac valid channels is missing: flux, model and err
    are 0 there and ``band_ok`` is False. ``cal_shape`` is ln(mean cal_spec over the band's valid
    channels) minus its band mean (the measured CHIME response shape for that file's calibrator).
    """
    data = np.asarray(data, np.float64)
    model = np.asarray(model, np.float64)
    flag = np.asarray(flag, bool)
    nf, nt = data.shape
    ch_ok = (
        np.asarray(good, bool)
        & np.asarray(cal_good, bool)
        & np.isfinite(fcf)
        & (np.asarray(fcf) > 0)
    )
    reg = np.asarray(region, int).reshape(-1, 2)
    r0, r1 = int(reg[:, 0].min()), int(reg[:, 1].max())
    c0, c1 = max(0, r0 - p.pad_bins), min(nt, r1 + p.pad_bins + 1)
    off = np.ones(nt, bool)
    off[c0:c1] = False
    if chan_sigma is None:
        dm = np.where(flag[:, off], data[:, off], np.nan)
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(dm, axis=1)
            chan_sigma = 1.4826 * np.nanmedian(np.abs(dm - med[:, None]), axis=1)
    chan_sigma = np.asarray(chan_sigma, np.float64)
    ch_ok &= np.isfinite(chan_sigma) & (chan_sigma > 0)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        inv = np.where(ch_ok, 1.0 / chan_sigma, 0.0)
        jy = np.asarray(fcf, np.float64) * chan_sigma
    nb = len(p.bands_mhz)
    flux = np.zeros((nb, nt))
    mod = np.zeros((nb, nt))
    err = np.zeros((nb, nt))
    band_ok = np.zeros(nb, bool)
    cal_shape = np.full(nb, np.nan)
    nvalid = np.zeros(nb, int)
    for b, (lo, hi) in enumerate(p.bands_mhz):
        inb = (freqs >= lo) & (freqs < hi)
        ch = inb & ch_ok
        nvalid[b] = int(ch.sum())
        if inb.sum() == 0 or ch.sum() < p.min_valid_chan_frac * inb.sum():
            continue
        w = flag[ch].astype(np.float64)  # (n_ch, nt)
        cnt = w.sum(0)
        f = inv[ch][:, None]
        kb = float(np.median(jy[ch]))
        with np.errstate(invalid="ignore", divide="ignore"):
            s = np.where(cnt > 0, kb * (w * f * np.nan_to_num(data[ch])).sum(0) / cnt, np.nan)
            m = np.where(cnt > 0, kb * (w * f * np.nan_to_num(model[ch])).sum(0) / cnt, np.nan)
        okt = np.isfinite(s) & (cnt >= p.min_valid_chan_frac * inb.sum())
        if (off & okt).sum() < p.min_off_bins:
            continue
        base = float(np.median(s[off & okt]))
        sig = float(1.4826 * np.median(np.abs(s[off & okt] - base)))
        if not np.isfinite(sig) or sig <= 0:
            continue
        flux[b] = np.where(okt, s - base, 0.0)
        mod[b] = np.where(okt, m, 0.0)
        err[b] = sig
        band_ok[b] = True
        cal_shape[b] = float(np.log(np.mean(cal_spec[ch])))
    if band_ok.any():
        cal_shape[band_ok] -= np.nanmean(cal_shape[band_ok])
    return {
        "flux": flux,
        "model": mod,
        "err": err,
        "band_ok": band_ok,
        "crop": (c0, c1),
        "region": (r0, r1),
        "cal_shape": cal_shape,
        "n_valid_chan": nvalid,
        "n_time": nt,
    }


def reduce_file(path, p: Params = DEFAULT) -> dict:
    """Open a Catalog 2 dynamic-spectrum HDF5 file and reduce it (see :func:`reduce_arrays`).

    Some files have no ``calibration`` group (e.g. FRB20200314A). There is then no finite flux
    conversion factor, so under rule 3 of the plan no channel is valid and every band is missing
    (``has_calibration`` False); the file is ineligible and cannot host an injection.
    """
    import h5py

    with h5py.File(path, "r") as f:
        nf = f["data"].shape[0]
        has_cal = "calibration" in f and "flux_conversion_factor" in f["calibration"]
        nan = np.full(nf, np.nan)
        chan_sigma = None  # per-channel off-pulse MAD of this window (reduce_arrays)
        red = reduce_arrays(
            f["data"][...],
            f["model"][...],
            f["flag"][...],
            f["good_freq"][...],
            f["calibration/good_freq"][...] if has_cal else np.zeros(nf, bool),
            f["calibration/flux_conversion_factor"][...] if has_cal else nan,
            f["calibration/spectrum"][...] if has_cal else nan,
            f["index_map/freqs"][...],
            np.asarray(f.attrs["pulse_emission_region"]),
            p,
            chan_sigma,
        )
        red["has_calibration"] = bool(has_cal)
        red["tns_name"] = str(f.attrs["tns_name"])
        red["dm_incoherent"] = float(f.attrs["dm_incoherent"])
        red["res_time"] = float(f.attrs["res_time"])
        red["beam_number"] = int(f.attrs["beam_number"])
        red["cal_source"] = str(f["calibration"].attrs.get("source_name", "")) if has_cal else ""
    return red


def crop(red: dict, key: str = "flux"):
    c0, c1 = red["crop"]
    return red[key][:, c0:c1], red["err"][:, c0:c1]


def n_pulses(red: dict, p: Params = DEFAULT) -> int:
    f, e = crop(red)
    return int(len(bt.find_pulses(f, e)))


def is_eligible(red: dict, p: Params = DEFAULT) -> bool:
    return bool(red["band_ok"].any()) and n_pulses(red, p) >= p.min_pulses


# -------------------------------------------------------------------------------- statistic


def pair_rho(ra: dict, rb: dict, fa=None, fb=None) -> tuple[float, int]:
    """burst_twins.xcorr_max on the two crops; returns (rho, lag) with the lag in full-window bins
    (B index = A index + lag). ``fa``/``fb`` override the crop fluxes (surrogates)."""
    xa = crop(ra)[0] if fa is None else fa
    xb = crop(rb)[0] if fb is None else fb
    rho, lag_c = bt.xcorr_max(xa, xb)
    return rho, int(lag_c + rb["crop"][0] - ra["crop"][0])


def chi2_union(ra: dict, rb: dict, lag: int, p: Params = DEFAULT) -> dict:
    """:func:`chi2_union_at` at ``lag`` + delta for delta in ``subbin_grid`` (A Fourier-shifted by
    delta); the delta with the largest p is kept and dof is reduced by 1 for it (if the grid has
    more than one value). Returns that result with ``delta``."""
    grid = p.subbin_grid if p.subbin_grid else (0.0,)
    best = None
    for d in grid:
        if d == 0:
            rs = ra
        else:
            rs = dict(ra, flux=fshift(ra["flux"], d, ra["n_time"]))
        c = chi2_union_at(rs, rb, lag, p)
        if c["testable"] and len(grid) > 1:
            c["dof"] -= 1
            c["p"] = float(stats.chi2.sf(c["chi2"], c["dof"])) if c["dof"] > 0 else 0.0
            c["passed"] = bool(c["gains_ok"] and c["p"] >= p.chi2_p_min)
        c["delta"] = float(d)
        if best is None or (c["p"], -abs(d)) > (best["p"], -abs(best["delta"])):
            best = c
    return best


def chi2_union_at(ra: dict, rb: dict, lag: int, p: Params = DEFAULT) -> dict:
    """Issue-#153 fixed chi2 twin test: B = a_b * A shifted by ``lag``, one gain per sub-band.

    The window is the union of A's crop (moved by ``lag``) and B's crop, in B's frame; it must lie
    inside both full windows, else the pair is "untestable" (fails). Per usable band (band_ok in
    both) a_b >= 0 minimises sum w (y - a x)^2 with w = 1 / (ey^2 + a^2 ex^2 + (sys a x)^2)
    (3 iterations, clipped at 0). chi2 sums over bands and bins, dof = n_bins * n_bands - n_bands.
    A band "has
    signal" if either member's peak S/N in the window is >= band_signal_snr; every such band needs
    a_b / sigma(a_b) >= gain_det_sigma, sigma(a_b) = (sum w x^2)^-1/2. Pass = testable and at least
    one band has signal and every signal band's gain is detected and p >= chi2_p_min.
    """
    x, ex, y, ey = ra["flux"], ra["err"], rb["flux"], rb["err"]
    na, nbt = x.shape[1], y.shape[1]
    lo = min(ra["crop"][0] + lag, rb["crop"][0])
    hi = max(ra["crop"][1] + lag, rb["crop"][1])
    out = {
        "chi2": np.inf,
        "dof": 0,
        "p": 0.0,
        "gains": [np.nan] * x.shape[0],
        "gain_snr": [np.nan] * x.shape[0],
    }
    out.update({"testable": False, "n_signal_bands": 0, "gains_ok": False, "passed": False})
    if lo < 0 or hi > nbt or lo - lag < 0 or hi - lag > na or hi - lo < 3:
        return out
    use = ra["band_ok"] & rb["band_ok"]
    if not use.any():
        return out
    chi2, nbin, nfit = 0.0, 0, 0
    sig_ok = True
    nsig = 0
    for b in np.flatnonzero(use):
        xs, ys = x[b, lo - lag : hi - lag], y[b, lo:hi]
        exs, eys = ex[b, lo - lag : hi - lag], ey[b, lo:hi]
        a = max(0.0, float((xs * ys).sum() / max((xs * xs).sum(), 1e-300)))
        for _ in range(3):
            w = 1.0 / (eys**2 + a * a * exs**2 + (p.chi2_sys_frac * a * xs) ** 2)
            a = max(0.0, float((w * xs * ys).sum() / max((w * xs * xs).sum(), 1e-300)))
        w = 1.0 / (eys**2 + a * a * exs**2 + (p.chi2_sys_frac * a * xs) ** 2)
        chi2 += float((w * (ys - a * xs) ** 2).sum())
        nbin += xs.size
        nfit += 1
        sa = 1.0 / np.sqrt(max(float((w * xs * xs).sum()), 1e-300))
        out["gains"][b] = a
        out["gain_snr"][b] = a / sa
        has = max(float(np.max(xs / exs)), float(np.max(ys / eys))) >= p.band_signal_snr
        if has:
            nsig += 1
            sig_ok &= a / sa >= p.gain_det_sigma
    dof = nbin - nfit
    pv = float(stats.chi2.sf(chi2, dof)) if dof > 0 else 0.0
    out.update({"chi2": chi2, "dof": int(dof), "p": pv, "testable": True, "n_signal_bands": nsig})
    out["gains_ok"] = bool(nsig > 0 and sig_ok)
    out["passed"] = bool(out["gains_ok"] and pv >= p.chi2_p_min)
    return out


def chi2_d071(ra: dict, rb: dict, lag: int) -> tuple[float, int, float]:
    """Secondary (reported, not deciding): burst_twins.twin_chi2 unchanged on the full windows."""
    use = ra["band_ok"] & rb["band_ok"]
    if not use.any():
        return np.inf, 0, 0.0
    return bt.twin_chi2(ra["flux"][use], ra["err"][use], rb["flux"][use], rb["err"][use], lag)


# -------------------------------------------------------------------------------- null


def surrogate_rho_max(reds: list[dict], pairs: list[tuple[int, int]], rng) -> float:
    """Maximum rho over ``pairs`` in one surrogate catalogue: every burst's crop pulse-shuffled
    (burst_twins.pulse_shuffle, all bands moved together)."""
    sur = [bt.pulse_shuffle(*crop(r), rng)[0] for r in reds]
    return max(pair_rho(reds[i], reds[j], sur[i], sur[j])[0] for i, j in pairs)


# -------------------------------------------------------------------------------- injections


def fshift(x: np.ndarray, s: float, n_out: int) -> np.ndarray:
    """Shift (B, N) curves by s bins (Fourier, zero-padded so nothing wraps) onto n_out bins."""
    nb, n = x.shape
    L = 1 << int(np.ceil(np.log2(n + n_out + 2 * int(abs(s)) + 4)))
    buf = np.zeros((nb, L))
    buf[:, :n] = x
    k = np.fft.rfftfreq(L)
    y = np.fft.irfft(np.fft.rfft(buf, axis=1) * np.exp(-2j * np.pi * k * s), L, axis=1)
    return y[:, :n_out]


def path_broaden(x: np.ndarray, tau600_ms: float, ddm: float, p: Params = DEFAULT) -> np.ndarray:
    """Extra line-of-sight broadening of a copy (simulated), per sub-band at its centre frequency:
    a unit-area one-sided exponential with tau = tau600 (nu/600 MHz)^-4 (scattering), then a
    unit-area boxcar of width 8.3 us * ddm * 0.0244 MHz * (nu/GHz)^-3 (intra-channel DM smearing of
    a DM difference ddm). Widths below one bin leave the band unchanged by that kernel."""
    out = np.array(x, float)
    for b, nu in enumerate(band_centres(p)):
        tau = tau600_ms * 1e-3 * (nu / 600.0) ** -4 / RES_TIME_S
        if tau > 0.05:
            t = np.arange(int(np.ceil(10 * tau)) + 1)
            k = np.exp(-t / tau)
            out[b] = np.convolve(out[b], k / k.sum())[: out.shape[1]]
        w = 8.3e-6 * ddm * RES_FREQ_MHZ * (nu / 1000.0) ** -3 / RES_TIME_S
        m = int(round(w))
        if m >= 2:
            k = np.ones(m) / m
            out[b] = np.convolve(out[b], k, mode="same")
    return out


def draw_gains(kind: str, rng, cal_shapes: np.ndarray, p: Params = DEFAULT) -> np.ndarray:
    """Per-sub-band gain of the copy relative to the original (simulated), geometric mean 1.

    G0: 1. G_cal and P_path: exp(l_X - l_Y), l = the measured cal_shape of two distinct sample files
    X, Y (rows of ``cal_shapes`` without NaN). G_beam: Gaussian formed beam with FWHM ~ 1/nu,
    ln g_b = -4 ln2 (u_B^2 - u_A^2) (nu_b / 600 MHz)^2, u_A, u_B ~ U(0, beam_u_max).
    """
    nb = len(p.bands_mhz)
    if kind == "G0":
        return np.ones(nb)
    if kind in ("G_cal", "P_path"):
        ok = np.flatnonzero(np.isfinite(cal_shapes).all(1))
        i, j = rng.choice(ok, 2, replace=False)
        lg = cal_shapes[i] - cal_shapes[j]
    elif kind == "G_beam":
        ua, ub = rng.uniform(0, p.beam_u_max, 2)
        lg = -4 * np.log(2) * (ub**2 - ua**2) * (band_centres(p) / 600.0) ** 2
    else:
        raise ValueError(kind)
    return np.exp(lg - lg.mean())


def make_copy(
    ra: dict,
    rb: dict,
    ratio: float,
    gains: np.ndarray,
    shift: float,
    path: bool,
    p: Params = DEFAULT,
) -> dict:
    """Injected twin of A in B's noise (simulated): B' = ratio * g_b * shift(model_A) + (data_B -
    model_B). The copy's emission region is A's moved so its centre lands at B's region centre plus
    ``shift`` bins; B's own burst is removed with B's fitburst model, so the noise (and B's model
    misfit) is B's real one. Errors are B's off-pulse sigma. Bands missing in A or B are missing."""
    ca = 0.5 * (ra["region"][0] + ra["region"][1])
    cb = 0.5 * (rb["region"][0] + rb["region"][1])
    s = cb - ca + shift
    tmpl = ra["model"]
    if path:
        tmpl = path_broaden(tmpl, p.path_tau600_ms, p.path_ddm, p)
    nbt = rb["n_time"]
    cp = fshift(tmpl, s, nbt)
    resid = rb["flux"] - rb["model"]
    ok = ra["band_ok"] & rb["band_ok"]
    flux = np.where(ok[:, None], ratio * gains[:, None] * cp + resid, 0.0)
    err = np.where(ok[:, None], rb["err"], 0.0)
    r0, r1 = int(round(ra["region"][0] + s)), int(round(ra["region"][1] + s))
    c0, c1 = max(0, r0 - p.pad_bins), min(nbt, r1 + p.pad_bins + 1)
    return {
        "flux": flux,
        "err": err,
        "model": np.zeros_like(flux),
        "band_ok": ok,
        "crop": (c0, c1),
        "region": (r0, r1),
        "n_time": nbt,
        "true_lag": s,
    }


def run_chain(ra: dict, rb: dict, rho_star: float, p: Params = DEFAULT) -> dict:
    """The whole frozen chain on one pair: eligibility of both, rho > rho*, chi2_union pass."""
    ea, eb = is_eligible(ra, p), is_eligible(rb, p)
    rho, lag = pair_rho(ra, rb)
    c = chi2_union(ra, rb, lag, p)
    flag = bool(ea and eb and rho > rho_star)
    return {
        "elig_a": ea,
        "elig_b": eb,
        "rho": rho,
        "lag": lag,
        "flag": flag,
        "chi2": c,
        "passed": bool(flag and c["passed"]),
    }


# -------------------------------------------------------------------------------- verdict


def verdict(s: dict, p: Params = DEFAULT) -> dict:
    """Frozen transferability rule (PLAN.md). ``s`` needs: n_eligible, eff_gcal_r1 (recovery at
    ratio 1, G_cal), n_chance_pairs, n_chance_pass, rho_inj_gcal_r1_median, rho_chance_q (the t3
    quantile of chance rho), n_rep_pairs, n_rep_pass. Each test is True/False/None (None = not
    computable).

    carries_over = T1 and T2 and T3 (all True). T4 is reported separately: True means same-source,
    different-event pairs pass the chain more often than t4_repeater_pass_max ("confused")."""
    if s["n_eligible"] < p.min_eligible or s["n_chance_pairs"] == 0:
        return {
            "testable": False,
            "T1": None,
            "T2": None,
            "T3": None,
            "T4": None,
            "carries_over": False,
        }
    t1 = bool(s["eff_gcal_r1"] >= p.t1_eff_min)
    t2 = bool(s["n_chance_pass"] / s["n_chance_pairs"] <= p.t2_chance_pass_max)
    t3 = bool(s["rho_inj_gcal_r1_median"] > s["rho_chance_q"])
    t4 = (
        None
        if s["n_rep_pairs"] == 0
        else bool(s["n_rep_pass"] / s["n_rep_pairs"] > p.t4_repeater_pass_max)
    )
    return {
        "testable": True,
        "T1": t1,
        "T2": t2,
        "T3": t3,
        "T4": t4,
        "carries_over": bool(t1 and t2 and t3),
    }

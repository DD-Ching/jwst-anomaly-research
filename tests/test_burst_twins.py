"""Offline tests for the S1 burst-twin screen (jwst_anomaly.burst_twins)."""

import numpy as np
import pytest
from astropy.io import fits

from jwst_anomaly import burst_twins as bt


def pulse(t, t0, rise, decay, amp):
    """Norris-like FRED pulse."""
    x = np.zeros_like(t)
    m = t > t0
    x[m] = amp * np.exp(-rise / (t[m] - t0) - (t[m] - t0) / decay)
    return x


def burst(n=200, pulses=((40, 1, 8, 10), (110, 1, 12, 6)), noise=0.3, seed=0, ratio=(1.0, 0.6)):
    rng = np.random.default_rng(seed)
    t = np.arange(n, dtype=float)
    s = sum(pulse(t, *p) for p in pulses)
    flux = np.array([r * s for r in ratio]) + rng.normal(0, noise, (len(ratio), n))
    err = np.full_like(flux, noise)
    return flux, err


# ------------------------------------------------------------------------------------------
# statistic


def test_xcorr_identical_copy_with_shift_is_one():
    f, _ = burst(noise=0.0)
    g = np.zeros((2, 260))
    g[:, 37 : 37 + 200] = 3.0 * f  # scaled and delayed by 37 bins
    rho, lag = bt.xcorr_max(f, g)
    assert rho == pytest.approx(1.0, abs=1e-9)
    assert lag == 37


def test_xcorr_band_ratio_mismatch_lowers_rho():
    f, _ = burst(noise=0.0, ratio=(1.0, 0.6))
    g, _ = burst(noise=0.0, ratio=(0.6, 1.0))  # same time profile, different spectrum
    assert bt.xcorr_max(f, g)[0] < 0.95


def test_xcorr_different_pulse_pattern_lower_than_twin():
    f, _ = burst(seed=1)
    twin, _ = burst(seed=2)
    other, _ = burst(pulses=((30, 1, 5, 10), (60, 1, 20, 3), (150, 2, 4, 8)), seed=3)
    assert bt.xcorr_max(f, twin)[0] > 0.9
    assert bt.xcorr_max(f, other)[0] < bt.xcorr_max(f, twin)[0] - 0.1


def test_twin_chi2_accepts_noisy_copy_rejects_different():
    f, e = burst(seed=4)
    twin, et = burst(seed=5)
    other, eo = burst(pulses=((40, 1, 8, 10), (110, 1, 12, 20)), seed=6)
    _, lag = bt.xcorr_max(f, twin)
    assert bt.twin_chi2(f, e, twin, et, lag)[2] > 1e-3
    _, lag = bt.xcorr_max(f, other)
    assert bt.twin_chi2(f, e, other, eo, lag)[2] < 1e-3


def test_stretch_recovers_time_dilated_copy():
    t = np.arange(400, dtype=float)
    a = np.array([pulse(t, 40, 1, 8, 10) + pulse(t, 90, 1, 6, 6)])
    b = np.array([pulse(t, 80, 2, 16, 10) + pulse(t, 180, 2, 12, 6)])  # same burst at s = 2
    ys, _ = bt.stretch(b, np.ones_like(b), 0.5)
    assert bt.xcorr_max(a, ys)[0] > bt.xcorr_max(a, b)[0]
    assert bt.xcorr_max(a, ys)[0] > 0.97


def test_rebin_factor_conserves_mean_rate():
    f = np.arange(12, dtype=float).reshape(1, 12)
    g, e = bt.rebin_factor(f, np.ones_like(f), 4)
    assert g.tolist() == [[1.5, 5.5, 9.5]]
    assert e[0, 0] == pytest.approx(0.5)


def test_rebin_overlap_wide_and_narrow_sources():
    tb = np.array([[0.0, 0.256], [0.256, 0.512], [0.512, 0.576], [0.576, 0.640]])
    flux = np.array([1.0, 2.0, 3.0, 5.0])
    fo, eo = bt.rebin(tb, flux, np.ones(4), 0.0, 0.64, 0.128)
    assert fo[:4].tolist() == [1.0, 1.0, 2.0, 2.0]
    assert fo[4] == pytest.approx(4.0)
    assert eo[0] == pytest.approx(1.0)  # a wide source bin keeps its own error


# ---------------------------------------------------------------------------------------------
# pulses


def test_pulse_count_two_vs_one_and_noise_robust():
    f2, e2 = burst()
    assert len(bt.find_pulses(f2, e2)) == 2
    f1, e1 = burst(pulses=((40, 1, 20, 10),))
    assert len(bt.find_pulses(f1, e1)) == 1
    # a bright single pulse with Poisson-like errors must not split into many pulses
    t = np.arange(300, dtype=float)
    s = pulse(t, 50, 2, 60, 500)
    rng = np.random.default_rng(7)
    err = np.sqrt(1 + s)[None]
    f = s[None] + rng.normal(size=(1, 300)) * err
    assert len(bt.find_pulses(f, err)) == 1


def test_pulse_shuffle_keeps_pulse_count_and_power():
    f, e = burst(pulses=((30, 1, 6, 10), (90, 1, 10, 4), (150, 1, 5, 7)), noise=0.2)
    g, _ = bt.pulse_shuffle(f, e, np.random.default_rng(0))
    assert g.shape == f.shape
    assert (g * g).sum() == pytest.approx((f * f).sum())
    assert len(bt.find_pulses(g, e)) == len(bt.find_pulses(f, e))
    assert not np.allclose(g, f)


# ------------------------------------------------------------------------------------------
# positions


def test_position_cut():
    sys = bt.Params().sys_deg
    thr = 3 * np.sqrt(2 * sys**2 + 2**2 + 3**2)
    assert bt.position_inconsistent(thr + 0.1, 2.0, 3.0)
    assert not bt.position_inconsistent(thr - 0.1, 2.0, 3.0)
    # catalogue err = 0 (localised by another instrument) still gets the systematic
    assert not bt.position_inconsistent(10.0, 0.0, 0.0)
    assert bt.angsep_deg(0, 0, 90, 0) == pytest.approx(90)
    assert bt.angsep_deg(10, 89, 190, 89) == pytest.approx(2)


# ------------------------------------------------------------------------------------------
# injection


def test_inject_twin_noise_level_and_recovery():
    f, e = burst(noise=0.0, seed=8)
    e = np.full_like(f, 0.2)
    eb = np.full_like(e, 0.5)  # host slot is noisier
    rng = np.random.default_rng(9)
    fa, ea, fb, ebo = bt.inject_twin(f, e, eb, 1.0, rng)
    assert np.median(ebo) == pytest.approx(0.5, rel=0.05)
    assert np.std(fb - f) == pytest.approx(0.5, rel=0.15)
    assert np.std(fa - f) == pytest.approx(0.2, rel=0.15)
    # independent noise realisations in the two members
    assert abs(np.corrcoef((fa - f).ravel(), (fb - f).ravel())[0, 1]) < 0.1
    assert bt.xcorr_max(fa, fb)[0] > 0.8
    _, _, weak, _ = bt.inject_twin(f, e, eb, 0.1, rng)
    assert bt.xcorr_max(fa, weak)[0] < bt.xcorr_max(fa, fb)[0]


def test_template_suppresses_noise():
    rng = np.random.default_rng(3)
    n = rng.normal(size=(1, 2000))
    assert np.std(bt.template(n)) < 0.6


def test_poisson_upper_limit():
    assert bt.poisson_upper_limit(0) == pytest.approx(2.9957, abs=1e-3)
    assert bt.poisson_upper_limit(1) == pytest.approx(4.7439, abs=1e-3)


# ------------------------------------------------------------------------------------- reduce /
# encode


def fake_bcat(n=400, dt=0.064):
    tb = np.c_[np.arange(n) * dt - 5.0, np.arange(1, n + 1) * dt - 5.0]
    t = tb.mean(1)
    sig = pulse(t, 0.0, 0.1, 1.0, 20) + pulse(t, 4.0, 0.1, 1.0, 10)
    ftot = np.c_[sig, np.full(n, 0.5)]
    fb = np.c_[0.6 * sig, np.full(n, 0.3)]
    ftot[7] = [-9.9e36, 0.0]  # fill value
    cols = [
        fits.Column("TIMEBIN", "2D", array=tb),
        fits.Column("PHTFLUX", "2E", array=ftot.astype(np.float32)),
        fits.Column("PHTFLUXB", "2E", array=fb.astype(np.float32)),
    ]
    return fits.HDUList(
        [fits.PrimaryHDU(), fits.BinTableHDU.from_columns([]), fits.BinTableHDU.from_columns(cols)]
    )


def test_reduce_bcat_and_encode_roundtrip():
    h = fake_bcat()
    red = bt.reduce_bcat(h, t90_start=0.0, t90=6.0)
    assert red["flux"].shape[0] == 2
    assert red["flux"].shape[1] <= bt.Params().max_bins
    assert np.nanmax(red["flux"][0]) == pytest.approx(0.6 * np.nanmax(red["flux"].sum(0)), rel=0.05)
    units, fs, es = bt.encode(red["flux"], red["err"])
    f2, e2 = bt.decode(units, fs, es)
    ok = np.isfinite(red["flux"])
    assert np.allclose(f2[ok], red["flux"][ok], atol=max(units) / 10)

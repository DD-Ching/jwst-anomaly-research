"""Offline tests for the S1 CHIME/FRB Catalog 2 adapter and rules (scripts/s1_chime_lib.py)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import s1_chime_lib as lib  # noqa: E402

P = lib.Params()


def gauss(t, t0, w, a):
    return a * np.exp(-0.5 * ((t - t0) / w) ** 2)


def red_from(flux, err, crop, region=None, model=None, band_ok=None, cal_shape=None):
    """A reduced burst dict as reduce_arrays returns it (for rule tests)."""
    nb, nt = flux.shape
    return {
        "flux": np.asarray(flux, float),
        "model": np.zeros_like(flux) if model is None else np.asarray(model, float),
        "err": np.broadcast_to(np.asarray(err, float).reshape(-1, 1), (nb, nt)).copy(),
        "band_ok": np.ones(nb, bool) if band_ok is None else np.asarray(band_ok, bool),
        "crop": crop,
        "region": region
        if region is not None
        else (crop[0] + P.pad_bins, crop[1] - P.pad_bins - 1),
        "n_time": nt,
        "cal_shape": np.zeros(nb) if cal_shape is None else cal_shape,
    }


def double(nt=160, t1=70.0, t2=80.0, amps=(30.0, 18.0), spec=(1.0, 0.8, 0.6, 0.3)):
    t = np.arange(nt, dtype=float)
    s = gauss(t, t1, 1.2, amps[0]) + gauss(t, t2, 1.5, amps[1])
    return np.array([c * s for c in spec])


# ------------------------------------------------------------------------------- reduction


def test_reduce_arrays_recovers_band_flux_noise_and_missing_band():
    rng = np.random.default_rng(1)
    nf, nt = 1024, 160
    freqs = np.linspace(400.2, 800.2, nf, endpoint=False)
    t = np.arange(nt)
    prof = gauss(t, 80, 1.5, 1.0)
    spec = np.where(freqs < 600, 2.0, 1.0)
    data = spec[:, None] * prof[None, :] + rng.normal(0, 0.5, (nf, nt))
    model = spec[:, None] * prof[None, :]
    flag = np.ones((nf, nt), bool)
    good = freqs < 700  # top band entirely bad -> missing
    fcf = np.full(nf, 3.0)  # Jy per data unit
    cal = np.exp(0.1 * (freqs - 600) / 200)
    red = lib.reduce_arrays(
        data, model, flag, good, np.ones(nf, bool), fcf, cal, freqs, [[78, 82]], P
    )
    assert red["band_ok"].tolist() == [True, True, True, False]
    assert red["crop"] == (68, 93)
    np.testing.assert_allclose(red["flux"][0, 80], 6.0, atol=0.6)  # 2 * 3 Jy
    np.testing.assert_allclose(red["flux"][2, 80], 3.0, atol=0.6)
    np.testing.assert_allclose(red["model"][0, 80], 6.0, rtol=1e-6)
    nch = (freqs >= 400) & (freqs < 500)
    np.testing.assert_allclose(red["err"][0, 0], 3.0 * 0.5 / np.sqrt(nch.sum()), rtol=0.25)
    assert np.all(red["flux"][3] == 0) and np.all(red["err"][3] == 0)
    assert np.nansum(red["cal_shape"][:3]) == pytest.approx(0.0, abs=1e-12)
    assert np.isnan(red["cal_shape"][3])


def test_reduce_arrays_weights_keep_huge_fcf_channels_from_dominating():
    # v1 bug: an unweighted mean let a few channels with a near-zero calibration response (huge
    # flux conversion factor) dominate the band noise
    rng = np.random.default_rng(11)
    nf, nt = 1024, 160
    freqs = np.linspace(400.2, 800.2, nf, endpoint=False)
    t = np.arange(nt)
    prof = gauss(t, 80, 1.5, 1.0)
    data = prof[None, :] + rng.normal(0, 0.5, (nf, nt))
    fcf = np.full(nf, 3.0)
    fcf[[3, 10, 300, 600, 900]] = 1e5
    args = (data, prof[None, :] * np.ones((nf, 1)), np.ones((nf, nt), bool), np.ones(nf, bool))
    red = lib.reduce_arrays(*args, np.ones(nf, bool), fcf, np.ones(nf), freqs, [[78, 82]], P)
    red_s = lib.reduce_arrays(
        *args, np.ones(nf, bool), fcf, np.ones(nf), freqs, [[78, 82]], P, np.full(nf, 0.5)
    )
    nch = ((freqs >= 400) & (freqs < 500)).sum()
    for r in (red, red_s):
        assert r["err"][0, 0] < 3.0 * 0.5 / np.sqrt(nch) * 1.3
        np.testing.assert_allclose(r["flux"][:, 80], 3.0, atol=0.5)


# ------------------------------------------------------------------------------- eligibility


def test_eligibility_needs_two_resolved_pulses():
    rng = np.random.default_rng(2)
    two = double() + rng.normal(0, 1, (4, 160))
    one = double(amps=(30.0, 0.0)) + rng.normal(0, 1, (4, 160))
    assert lib.is_eligible(red_from(two, np.ones(4), (55, 95)))
    assert not lib.is_eligible(red_from(one, np.ones(4), (55, 95)))


# ------------------------------------------------------------------------------- statistic


def test_pair_rho_lag_is_in_full_window_frame():
    a = double()
    b = np.zeros((4, 200))
    b[:, 37:197] = 0.5 * a  # B index = A index + 37
    ra = red_from(a, np.ones(4), (55, 95))
    rb = red_from(b, np.ones(4), (92, 132))
    rho, lag = lib.pair_rho(ra, rb)
    assert rho == pytest.approx(1.0, abs=1e-9)
    assert lag == 37


def test_chi2_union_passes_scaled_copy_with_per_band_gains():
    rng = np.random.default_rng(3)
    a = double()
    g = np.array([1.3, 1.0, 0.8, 0.6])
    b = np.zeros((4, 180))
    b[:, 10:170] = 0.5 * g[:, None] * a
    ra = red_from(a + rng.normal(0, 1, a.shape), np.ones(4), (55, 95))
    rb = red_from(b + rng.normal(0, 1, b.shape), np.ones(4), (65, 105))
    rho, lag = lib.pair_rho(ra, rb)
    assert lag == 10
    c = lib.chi2_union(ra, rb, lag, P)
    assert c["testable"] and c["gains_ok"] and c["passed"]
    np.testing.assert_allclose(np.array(c["gains"])[:3], 0.5 * g[:3], rtol=0.15)


def test_chi2_union_subbin_refinement_needed_for_narrow_bright_pulses():
    rng = np.random.default_rng(8)
    a = double()
    b = lib.fshift(a, 10.3, 180)
    ra = red_from(a + rng.normal(0, 1, a.shape), np.ones(4), (55, 95))
    rb = red_from(b + rng.normal(0, 1, b.shape), np.ones(4), (65, 105))
    rho, lag = lib.pair_rho(ra, rb)
    assert lag == 10
    c = lib.chi2_union(ra, rb, lag, P)
    assert c["passed"] and c["delta"] == pytest.approx(0.3)
    c0 = lib.chi2_union(ra, rb, lag, lib.Params(subbin_grid=(0.0,)))
    assert not c0["passed"]


def test_chi2_union_rejects_single_pulse_matching_one_pulse_of_a_double():
    # issue #153 failure: a single pulse slid onto one pulse of a double must not pass
    rng = np.random.default_rng(4)
    a = double(amps=(30.0, 30.0), t1=60.0, t2=90.0)
    single = double(amps=(30.0, 0.0), t1=60.0, t2=90.0)
    ra = red_from(a + rng.normal(0, 1, a.shape), np.ones(4), (45, 105))
    rb = red_from(single + rng.normal(0, 1, a.shape), np.ones(4), (45, 75), region=(55, 64))
    rho, lag = lib.pair_rho(ra, rb)
    c = lib.chi2_union(ra, rb, lag, P)
    assert c["testable"]
    assert not c["passed"]


def test_chi2_union_gain_must_be_detected_in_signal_bands():
    # B has signal only in band 0; A has strong signal in band 3 that B lacks -> gain 0 there, fails
    rng = np.random.default_rng(5)
    a = double(spec=(1.0, 0.0, 0.0, 1.0))
    b = double(spec=(1.0, 0.0, 0.0, 0.0))
    ra = red_from(a + rng.normal(0, 1, a.shape), np.ones(4), (55, 95))
    rb = red_from(b + rng.normal(0, 1, b.shape), np.ones(4), (55, 95))
    c = lib.chi2_union(ra, rb, 0, P)
    assert c["testable"] and not c["gains_ok"] and not c["passed"]


def test_chi2_union_untestable_when_union_leaves_a_window():
    a = double()
    ra = red_from(a, np.ones(4), (55, 95))
    rb = red_from(a, np.ones(4), (55, 95))
    c = lib.chi2_union(ra, rb, 120, P)  # A's crop moved to 175-215 > 160
    assert not c["testable"] and not c["passed"]


# ------------------------------------------------------------------------------- injections


def test_draw_gains_sets_and_geometric_mean():
    rng = np.random.default_rng(6)
    cal = np.array([[0.1, 0.0, -0.05, -0.05], [-0.2, 0.1, 0.1, 0.0], [np.nan, 0, 0, 0]])
    assert np.all(lib.draw_gains("G0", rng, cal) == 1)
    for kind in ("G_cal", "G_beam", "P_path"):
        g = lib.draw_gains(kind, rng, cal)
        assert np.mean(np.log(g)) == pytest.approx(0.0, abs=1e-12)
    g = lib.draw_gains("G_cal", rng, cal)
    d = cal[0] - cal[1]
    assert np.allclose(np.log(g), d - d.mean()) or np.allclose(np.log(g), -d + d.mean())
    with pytest.raises(ValueError):
        lib.draw_gains("nope", rng, cal)


def test_path_broaden_keeps_area_and_broadens_low_band_most():
    x = np.zeros((4, 200))
    x[:, 50] = 1.0
    y = lib.path_broaden(x, 1.0, 300.0)
    np.testing.assert_allclose(y.sum(1), 1.0, rtol=1e-6)
    widths = [np.sum(row > 0.01) for row in y]
    assert widths[0] > widths[3]


def test_fshift_moves_by_fractional_bins():
    t = np.arange(100, dtype=float)
    x = gauss(t, 30.0, 2.0, 1.0)[None, :]
    y = lib.fshift(x, 12.5, 120)
    assert np.sum(np.arange(120) * y[0]) / y[0].sum() == pytest.approx(42.5, abs=1e-3)


def test_make_copy_in_quiet_host_is_recovered_by_the_chain():
    rng = np.random.default_rng(7)
    m = double()
    ra = red_from(m + rng.normal(0, 1, m.shape), np.ones(4), (55, 95), region=(65, 84), model=m)
    host = np.zeros((4, 170))
    rb = red_from(host + rng.normal(0, 1, host.shape), np.ones(4), (60, 100), region=(70, 89))
    cp = lib.make_copy(ra, rb, 1.0, np.ones(4), 3.3, False, P)
    out = lib.run_chain(ra, cp, 0.5, P)
    assert out["elig_a"] and out["elig_b"]
    assert abs(out["lag"] - cp["true_lag"]) <= 1
    assert out["rho"] > 0.9 and out["passed"]


# ------------------------------------------------------------------------------- verdict


def _vin(**kw):
    base = {
        "n_eligible": 20,
        "eff_gcal_r1": 0.7,
        "n_chance_pairs": 200,
        "n_chance_pass": 1,
        "rho_inj_gcal_r1_median": 0.9,
        "rho_chance_q": 0.6,
        "n_rep_pairs": 20,
        "n_rep_pass": 0,
    }
    base.update(kw)
    return base


def test_verdict_rules():
    v = lib.verdict(_vin())
    assert v["testable"] and v["carries_over"] and v["T4"] is False
    assert not lib.verdict(_vin(eff_gcal_r1=0.49))["carries_over"]
    assert not lib.verdict(_vin(n_chance_pass=3))["carries_over"]  # 1.5 % > 1 %
    assert lib.verdict(_vin(n_chance_pass=2))["T2"]  # exactly 1 %
    assert not lib.verdict(_vin(rho_inj_gcal_r1_median=0.6))["T3"]  # must exceed, not equal
    assert lib.verdict(_vin(n_rep_pass=2))["T4"] is True  # 10 % > 5 %
    assert lib.verdict(_vin(n_rep_pairs=0))["T4"] is None
    nt = lib.verdict(_vin(n_eligible=9))
    assert not nt["testable"] and not nt["carries_over"]


# ------------------------------------------------------------------------------- selection


def test_select_sample_rules():
    rows = []
    for k in range(4):  # repeater R1: 4 clean events
        rows += [
            dict(
                tns_name=f"R1e{k}",
                sub_num=s,
                repeater_name="R1",
                snr_fitb=10 + k,
                excluded_flag=0,
                sidelobe_flag=0,
            )
            for s in range(2)
        ]
    for k in range(2):  # repeater R2: 2 clean events, one single-component
        rows += [
            dict(
                tns_name=f"R2e{k}",
                sub_num=s,
                repeater_name="R2",
                snr_fitb=50,
                excluded_flag=0,
                sidelobe_flag=0,
            )
            for s in range(2)
        ]
    rows += [
        dict(
            tns_name="R2single",
            sub_num=0,
            repeater_name="R2",
            snr_fitb=99,
            excluded_flag=0,
            sidelobe_flag=0,
        )
    ]
    rows += [
        dict(
            tns_name="N1",
            sub_num=s,
            repeater_name=np.nan,
            snr_fitb=30,
            excluded_flag=0,
            sidelobe_flag=0,
        )
        for s in range(3)
    ]
    rows += [
        dict(
            tns_name="N2",
            sub_num=s,
            repeater_name=np.nan,
            snr_fitb=40,
            excluded_flag=1,
            sidelobe_flag=0,
        )
        for s in range(2)
    ]
    rows += [
        dict(
            tns_name="N3",
            sub_num=s,
            repeater_name=np.nan,
            snr_fitb=20,
            excluded_flag=0,
            sidelobe_flag=0,
        )
        for s in range(2)
    ]
    rows += [
        dict(
            tns_name="N4",
            sub_num=s,
            repeater_name=np.nan,
            snr_fitb=60,
            excluded_flag=0,
            sidelobe_flag=0,
        )
        for s in range(2)
    ]
    cat = pd.DataFrame(rows)
    files = {
        f"{n}_stokesi_dynamic_spectrum.h5" for n in cat.tns_name if n != "N4"
    }  # N4 has no file
    s = lib.select_sample(cat, files, reps=1, per_rep=3, nonrep=2)
    assert list(s.tns_name) == ["R1e3", "R1e2", "R1e1", "N1", "N3"]
    assert list(s.source) == ["R1", "R1", "R1", "N1", "N3"]

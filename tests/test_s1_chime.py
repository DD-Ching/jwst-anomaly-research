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
    # K_b * mean(1 / sigma_f) with estimated sigma_f: a ~1 % scale bias, shared by data and model
    np.testing.assert_allclose(red["model"][0, 80], 6.0, rtol=0.05)
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


def test_select_delta_prefers_a_passing_delta_over_higher_p_with_failed_gains():
    # /code-review failure case: d = 0.2 (p 0.02, gains detected) must win over d = -0.3 (p 0.05,
    # gains not detected), so the pair passes
    good = {"testable": True, "gains_ok": True, "p": 0.02, "delta": 0.2, "passed": True}
    bad = {"testable": True, "gains_ok": False, "p": 0.05, "delta": -0.3, "passed": False}
    for order in ([good, bad], [bad, good]):
        c = lib.select_delta(order)
        assert c["delta"] == 0.2 and c["passed"]
    # no delta with detected gains: best p overall is reported, and it fails
    c = lib.select_delta([bad, dict(bad, p=0.01, delta=0.1)])
    assert c["delta"] == -0.3 and not c["passed"]
    # untestable deltas never win over testable ones with detected gains
    unt = {"testable": False, "gains_ok": False, "p": 0.0, "delta": 0.0, "passed": False}
    assert lib.select_delta([unt, good])["delta"] == 0.2
    # ties go to the smaller |delta|
    assert lib.select_delta([dict(good, delta=-0.4), dict(good, delta=0.1)])["delta"] == 0.1


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


def test_posthoc_members_and_path_kernels():
    rng = np.random.default_rng(9)
    m = double()
    ra = red_from(m, np.full(4, 0.5), (55, 95), region=(65, 84), model=m)
    a2 = lib.noisy_template(ra, rng)
    assert np.std(a2["flux"] - m) == pytest.approx(0.5, rel=0.1)
    host = np.zeros((4, 170))
    rb = red_from(host + 100.0, np.ones(4), (60, 100), region=(70, 89), model=host)
    cp = lib.make_copy(ra, rb, 1.0, np.ones(4), 0.0, False, P, rng)
    assert abs(np.median(cp["flux"])) < 0.2  # Gaussian host noise, not the host residual (100)
    k = lib.path_kernel_bins(1.0, 300.0)
    assert k["tau_bins"][0] > k["tau_bins"][3]
    assert max(k["dm_box_bins"]) < 2  # the DM-smearing boxcar of the P_path set is not applied


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


# ------------------------------------------------------------------------------- fetch window


def test_bounded_map_keeps_at_most_window_files_on_disk(tmp_path):
    import threading
    import time

    import s1_chime as drv

    lock = threading.Lock()
    state = {"running": 0, "max_running": 0}

    def fake_download(k):
        with lock:
            state["running"] += 1
            state["max_running"] = max(state["max_running"], state["running"])
        time.sleep(0.01)
        f = tmp_path / f"f{k}.h5"
        f.write_bytes(b"x")
        with lock:
            state["running"] -= 1
        return k, f

    seen, max_on_disk = [], 0
    for k, f in drv.bounded_map(fake_download, range(12), window=3):
        time.sleep(0.02)  # slow "reduction": downloads must not run ahead
        max_on_disk = max(max_on_disk, len(list(tmp_path.iterdir())))
        f.unlink()
        seen.append(k)
    assert seen == list(range(12))  # submission order
    assert max_on_disk <= 3
    assert state["max_running"] <= 3


def _fake_fetch(tmp_path, fail_at=None):
    import time

    def fetch(job, cancel=None, touched=None):
        k = job
        f = tmp_path / f"f{k}.h5"
        if touched is not None:
            touched.add(f)
        f.write_bytes(b"x")  # (partial) transfer on disk
        time.sleep(0.01)
        if k == fail_at:
            raise RuntimeError("download failed")
        return f"f{k}.h5", {"dest": str(f), "downloaded": True}, None

    return fetch


def test_fetch_and_reduce_deletes_everything_when_the_reducer_raises(tmp_path):
    import s1_chime as drv

    done = []

    def handle(res):
        if len(done) == 2:
            raise ValueError("reader bug")
        done.append(res[0])
        Path(res[1]["dest"]).unlink()

    with pytest.raises(ValueError):
        drv.fetch_and_reduce(range(10), _fake_fetch(tmp_path), handle, window=3)
    assert done == ["f0.h5", "f1.h5"]
    assert list(tmp_path.iterdir()) == []


def test_bounded_map_cleanup_alone_deletes_unhandled_results(tmp_path):
    import s1_chime as drv

    gen = drv.bounded_map(_fake_fetch(tmp_path), range(10), 3, drv._drop)
    with pytest.raises(ValueError):
        try:
            for k, res in enumerate(gen):
                if k == 2:
                    raise ValueError("reader bug")
                Path(res[1]["dest"]).unlink()
        except BaseException:
            gen.close()
            raise
    assert list(tmp_path.iterdir()) == []


def test_bounded_map_cleans_the_head_future_on_interrupt_while_waiting(tmp_path):
    # /code-review round 3, finding 1: an interrupt during the head future's result() wait must
    # still clean that future's file (v2 had popped it from pending before waiting)
    import os
    import signal
    import time

    import s1_chime as drv

    def fetch(k):
        f = tmp_path / f"f{k}.h5"
        f.write_bytes(b"x")
        if k == 0:
            os.kill(os.getpid(), signal.SIGINT)  # main thread is blocked in result()
            time.sleep(0.3)
        return f"f{k}.h5", {"dest": str(f), "downloaded": True}, None

    gen = drv.bounded_map(fetch, range(4), 2, drv._drop)
    with pytest.raises(KeyboardInterrupt):
        for res in gen:
            Path(res[1]["dest"]).unlink()
    assert list(tmp_path.iterdir()) == []


def test_fetch_and_reduce_deletes_everything_when_a_download_raises(tmp_path):
    import s1_chime as drv

    def handle(res):
        Path(res[1]["dest"]).unlink()

    with pytest.raises(RuntimeError):
        drv.fetch_and_reduce(range(10), _fake_fetch(tmp_path, fail_at=4), handle, window=3)
    assert list(tmp_path.iterdir()) == []


def test_fetch_and_reduce_keep_raw_keeps_files(tmp_path):
    import s1_chime as drv

    def handle(res):
        raise ValueError("reader bug")

    with pytest.raises(ValueError):
        drv.fetch_and_reduce(range(5), _fake_fetch(tmp_path), handle, keep_raw=True, window=2)
    assert len(list(tmp_path.iterdir())) >= 1


def _fake_curl(fail_names=()):
    def curl(dest, url, cancel):
        if any(url.endswith(n) for n in fail_names):
            Path(dest).write_bytes(b"partial")
            return 22, "HTTP 503"
        Path(dest).write_bytes(b"full file")
        return 0, ""

    return curl


def test_fetch_and_reduce_keeps_files_it_did_not_download(tmp_path, monkeypatch):
    # finding 2: on an exception only files this run downloaded or started writing are deleted
    import hashlib

    import s1_chime as drv

    monkeypatch.setattr(drv, "_curl", _fake_curl())
    raw = tmp_path / "raw"
    raw.mkdir()
    reused = raw / "A.h5"  # kept by an earlier --keep-raw run, sha-verified against the log
    reused.write_bytes(b"earlier run")
    other = raw / "other_run.h5"  # not ours at all
    other.write_bytes(b"x")
    log = tmp_path / "downloads.log"
    h = hashlib.sha256(b"earlier run").hexdigest()
    log.write_text(f"2026-10-21T00:00:00Z\t{drv.BASE}A.h5\t11\t{h}\tearlier\n")
    jobs = [(n, raw, log) for n in ("A.h5", "B.h5")]

    def handle(res):
        raise ValueError("reader bug")

    with pytest.raises(ValueError):
        drv.fetch_and_reduce(jobs, drv._fetch_one, handle, window=2)
    assert sorted(f.name for f in raw.iterdir()) == ["A.h5", "other_run.h5"]


def test_fetch_and_reduce_cleanup_is_prompt_while_a_download_backs_off(tmp_path, monkeypatch):
    # finding 3: the cancel event stops the 10 / 30 s back-off of a retrying download
    import time

    import s1_chime as drv

    monkeypatch.setattr(drv, "_curl", _fake_curl(fail_names=("B.h5",)))
    raw = tmp_path / "raw"
    raw.mkdir()
    log = tmp_path / "downloads.log"
    log.write_text("")
    jobs = [(n, raw, log) for n in ("A.h5", "B.h5")]

    def handle(res):
        time.sleep(0.2)  # B has failed once and is in its 10 s back-off
        raise ValueError("reader bug")

    t0 = time.time()
    with pytest.raises(ValueError):
        drv.fetch_and_reduce(jobs, drv._fetch_one, handle, window=2)
    assert time.time() - t0 < 3.0
    assert list(raw.iterdir()) == []


def test_run_cancellable_terminates_a_running_transfer():
    import threading
    import time

    import s1_chime as drv

    ev = threading.Event()
    threading.Timer(0.3, ev.set).start()
    t0 = time.time()
    rc, err = drv._run_cancellable(["sleep", "30"], ev, poll_s=0.1)
    assert rc is None and err == "cancelled" and time.time() - t0 < 3.0


def test_fetch_one_deletes_partial_file_after_failed_attempts(tmp_path, monkeypatch):
    import s1_chime as drv

    calls = []

    def fake_curl(dest, url, cancel):
        Path(dest).write_bytes(b"partial")
        calls.append(dest)
        return 22, "HTTP 503"

    monkeypatch.setattr(drv, "_curl", fake_curl)
    monkeypatch.setattr(drv, "RETRY_SLEEP_S", (0.0, 0.0))
    log = tmp_path / "downloads.log"
    log.write_text("")
    raw = tmp_path / "raw"
    raw.mkdir()
    name, info, err = drv._fetch_one(("FRBX_stokesi_dynamic_spectrum.h5", raw, log))
    assert len(calls) == 3 and info is None and "curl failed" in err
    assert list(raw.iterdir()) == [] and log.read_text() == ""

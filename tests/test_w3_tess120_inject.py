"""Offline tests of scripts/w3_tess120_inject.py (results/w3_kepler/tess120_inject/PLAN.md)."""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import w3_shortcad_reach as ws  # noqa: E402
import w3_tess120_inject as wi  # noqa: E402

from jwst_anomaly import exotic_sim as es  # noqa: E402

CAD = 120 / 86400


def test_weighted_draws_are_proportional():
    rng = np.random.default_rng(1)
    d = wi.weighted_draws(np.array([1.0, 3.0]), 40000, rng)
    assert abs((d == 1).mean() - 0.75) < 0.01


def test_draw_events_couples_tE_and_rho_to_the_lens():
    rng = np.random.default_rng(2)
    m, v, r, d = 1e-6, 220.0, 1.0, 500.0
    for ev in wi.draw_events(rng, 50, m, v, r, d, 0.0, 27.0):
        r_e_pc = math.sqrt(ws.SCHW4_PC * m * ev["x"] * (1 - ev["x"]) * d)
        assert abs(ev["t_E_d"] * 86400 * v / (r_e_pc * ws.PC_KM) - 1) < 1e-9
        assert abs(ev["rho"] * r_e_pc * ws.PC_KM / (ev["x"] * r * ws.R_SUN_KM) - 1) < 1e-9
        assert 0 <= ev["u0"] <= wi.U0_MAX and 0 <= ev["t0"] <= 27.0


def test_model_flux_matches_point_sampling_for_long_events():
    ev = {"t_E_d": 2.0, "rho": 0.05, "u0": 0.5, "t0": 10.0}
    t = 10.0 + np.array([-3.0, -1.0, 0.0, 1.5, 2.5])
    f = wi.model_flux(t, ev, CAD)
    a = es.finite_source_magnification(np.hypot(0.5, (t - 10.0) / 2.0), 0.05, 1.0, -1)
    assert np.allclose(f, a, rtol=2e-3, atol=2e-3)


def test_ideal_identifiable_for_a_resolved_small_rho_event():
    ev = {"t_E_d": 50 * CAD, "rho": 0.05, "u0": 0.3, "t0": 5.0}
    assert wi.ideal_identifiable(ev, CAD, 0.0)
    short = ev | {"t_E_d": 0.5 * CAD}
    assert not wi.ideal_identifiable(short, CAD, 0.0)


def test_detect_and_match():
    t = np.arange(100) * CAD
    f = np.ones(100)
    f[45:50] = 0.05  # dark run
    f[41:45] = 1.3  # shoulders
    f[50:54] = 1.3
    dets = wi.detect(t, f, 0.0, CAD)
    assert len(dets) == 1 and dets[0]["S"] and dets[0]["n"] == 5
    ev = {"t0": t[47], "t_E_d": 2 * CAD, "u0": 0.0}
    assert wi.matched(dets, ev)
    assert not wi.matched(dets, ev | {"t0": t[90]})
    f2 = f.copy()
    f2[50:54] = 1.0  # no post shoulder
    assert not wi.detect(t, f2, 0.0, CAD)[0]["S"]
    tg = np.delete(t, 47)  # a gap splits the dark run into 2 + 2
    assert wi.detect(tg, np.delete(f, 47), 0.0, CAD) == []


def test_verdict():
    assert wi.verdict(10, 0.9, 0.8, 1.0) == "insufficient"
    assert wi.verdict(100, 0.7, 0.55, 0.8) == "retains at least half"
    assert wi.verdict(100, 0.3, 0.2, 0.45) == "retains less than half"
    assert wi.verdict(100, 0.55, 0.4, 0.7).startswith("interval crosses")


def test_lc_url_matches_the_mast_redirect():
    oid = "tess2024030031500-s0075-0000000298053137-0270-s"
    assert wi.lc_url(oid) == (
        "https://stpubdata.s3.us-east-1.amazonaws.com/tess/public/tid/s0075/0000/0002/9805/3137/"
        + oid
        + "_lc.fits"
    )


def test_process_on_a_flat_noise_free_light_curve_pairs_the_levels():
    import io

    from astropy.io import fits

    n = int(10 / CAD)
    cols = [
        fits.Column(name="TIME", format="D", array=2000.0 + np.arange(n) * CAD),
        fits.Column(name="PDCSAP_FLUX", format="E", array=np.full(n, 1000.0)),
        fits.Column(name="QUALITY", format="J", array=np.zeros(n, int)),
    ]
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header["TIMEDEL"] = CAD
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
    ss = {"obs_id": "synthetic", "draws": [0], "radius": 0.3, "dist": 50.0}
    res = wi.process(ss, buf.getvalue())
    assert res["background_D"] == 0 and res["n_good"] == n
    rows = [r for v in res["scenarios"].values() for r in v]
    assert len(rows) == wi.N_INJ * len(wi.SCENARIOS)
    assert all(r["L1a"] == r["L1b"] for r in rows)  # no flags, no gaps
    short = [r for r in rows if r["t_E_cad"] < 30]  # far wings inside the 3 d median window
    assert short and all(r["L2"] == r["L1b"] for r in short)  # flat, noise-free: L2 = L1b
    assert any(r["ideal"] for r in rows)

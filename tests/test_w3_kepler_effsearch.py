"""Offline tests of scripts/w3_kepler_effsearch.py (results/w3_kepler/effsearch/PLAN.md)."""

import io
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w3_kepler as wk  # noqa: E402
import w3_kepler_effsearch as es  # noqa: E402


def _star(n_q=2, n=400):
    t = np.concatenate([q * 100.0 + np.arange(n) * wk.CAD for q in range(n_q)])
    return {
        "TIME": t,
        "CADENCENO": np.arange(len(t), dtype=float),
        "SAP_FLUX": np.full(len(t), 1000.0),
        "SAP_FLUX_ERR": np.full(len(t), 0.1),
        "PDCSAP_FLUX": np.full(len(t), 1000.0),
        "SAP_QUALITY": np.zeros(len(t)),
        "QUARTER": np.repeat(np.arange(1, n_q + 1, dtype=float), n),
    }


def _event(st, t0, te, u0, rho=0.01, blend=1.0, base=frozenset(), seed=1):
    rng = np.random.default_rng(seed)
    f = wk.w3_flux_avg(st["TIME"], t0, te, u0, rho, blend)
    sap, err = wk.inject_sap(st["SAP_FLUX"], st["SAP_FLUX_ERR"], f, rng)
    res = wk.screen_star(st, sap, err, {}, local=True)
    h = te * np.sqrt(4 - u0**2)
    t_fin = st["TIME"][np.isfinite(st["TIME"])]
    nid = es.n_dim_ideal(t_fin, t0, te, u0, rho, blend)
    return es.classify(st, sap, err, res, set(base), f, t0, h, nid), res


def test_columns_skips_q0():
    from astropy.io import fits

    def llc(q):
        cols = [fits.Column(name=c, format="D", array=np.arange(3.0)) for c in wk.COLS]
        hdu0 = fits.PrimaryHDU()
        hdu0.header["QUARTER"], hdu0.header["CHANNEL"] = q, 7
        buf = io.BytesIO()
        fits.HDUList([hdu0, fits.BinTableHDU.from_columns(cols)]).writeto(buf)
        return buf.getvalue()

    assert es.columns(llc(0)) is None
    arr = es.columns(llc(5))
    assert set(arr) == set(wk.COLS) | {"QUARTER", "CHANNEL"}
    assert arr["QUARTER"].tolist() == [5.0] * 3 and arr["CHANNEL"].tolist() == [7.0] * 3


def test_seed_streams_reproducible_and_distinct():
    def draw(rank, c, rep):
        return np.random.default_rng([es.SEED, rank, c, rep]).uniform(size=2).tolist()

    assert draw(3, 4, 0) == draw(3, 4, 0)
    assert len({tuple(draw(r, c, p)) for r in range(3) for c in range(3) for p in range(2)}) == 18


def test_n_dim_ideal_long_and_short():
    t = np.arange(0, 20, wk.CAD)
    n_long = es.n_dim_ideal(t, 10.0, 1.0, 0.0, 0.01, 1.0)
    assert abs(n_long - 4.0 / wk.CAD) < 6  # umbra |t - t0| < 2 t_E
    assert es.n_dim_ideal(t, 10.0, 0.1, 1.9999, 0.01, 1.0) < wk.MIN_RUN


def test_classify_recovered_event():
    ev, _ = _event(_star(), 4.0, 0.3, 0.0)
    assert ev["recovered"] and ev["stage"] == "recovered" and ev["n_matched"] >= 1
    assert not ev["spurious"]


def test_classify_gap_quality_short():
    st = _star()
    ev, _ = _event(st, 50.0, 0.3, 0.0)  # centre between the two quarters
    assert ev["stage"] == "gap" and ev["n_exist"] == 0
    st_q = _star()
    st_q["SAP_QUALITY"][np.abs(st_q["TIME"] - 4.0) < 1.0] = 1  # attitude tweak bit
    ev, _ = _event(st_q, 4.0, 0.3, 0.0)
    assert ev["stage"] == "quality" and ev["n_exist"] >= 3 and ev["n_good"] < 3
    ev, _ = _event(st, 4.0, 0.1, 1.9999)
    assert ev["stage"] == "short"


def test_existing_run_is_not_a_recovery():
    st = _star()
    _, res = _event(st, 4.0, 0.3, 0.0)
    base = {(x["t_start"], x["t_end"]) for x in res}
    ev, _ = _event(st, 4.0, 0.3, 0.0, base=base)
    assert not ev["recovered"] and ev["n_matched"] == 0 and ev["stage"] == "not_flagged"


def test_volume_uses_sampling_window_and_clusters():
    v = es.volume(np.array([4.0, 4.0, 4.0, 4.0]), np.array([2, 0, 1, 1]))
    assert v["effective_star_years"] == 8.0 and v["recovery_rate"] == 0.5
    assert v["n_eff"] == 8.0  # Kish: 2 (sum T)^2 / sum T^2
    lo, hi = es.wilson_p(0.5, 8.0)
    assert v["effective_mc95"] == [round(16 * lo, 1), round(16 * hi, 1)]
    assert v["per_star_mean_yr"] == 2.0
    assert abs(v["per_star_se_yr"] - np.std([4, 0, 2, 2], ddof=1) / 2) < 1e-4
    assert v["stars_k"] == [1, 2, 1]


def test_limit_and_multi_event_share():
    lim = es.limit(1000.0, 4.0)
    assert abs(lim["rate95_per_star_year"] - 0.002996) < 1e-6
    assert abs(lim["lambda_per_star"] - 0.012) < 1e-4
    assert 0.005 < lim["share_hosts_ge2"] < 0.007  # about lambda / 2
    assert es.limit(0.0, 4.0) is None


def test_lengthened_existing_run_is_not_matched():
    st = _star()
    t = st["TIME"]
    box = (t >= 3.0) & (t < 4.2)
    st["SAP_FLUX"] = np.where(box, 100.0, 1000.0)  # an existing deep dip of ~59 cadences
    base_runs = wk.screen_star(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], {}, local=True)
    base = {(x["t_start"], x["t_end"]) for x in base_runs}
    te, u0 = 0.1, 1.0
    h = te * np.sqrt(4 - u0**2)
    ev, res = _event(st, 4.1 + h, te, u0, base=base)  # the umbra starts 0.1 d before the dip ends
    new = [x for x in res if (x["t_start"], x["t_end"]) not in base]
    assert len(new) == 1 and new[0]["t_start"] < 3.1 and new[0]["t_end"] > 4.2  # one merged run
    assert ev["n_matched"] == 0 and not ev["recovered"]


def test_box_dip_stops_at_the_shoulder_test():
    st = _star()
    t = st["TIME"]
    f = np.where((t >= 4.0) & (t < 4.6), 0.1, 1.0)  # fast edges, no shoulders
    sap, err = wk.inject_sap(st["SAP_FLUX"], st["SAP_FLUX_ERR"], f, np.random.default_rng(2))
    res = wk.screen_star(st, sap, err, {}, local=True)
    ev = es.classify(st, sap, err, res, set(), f, 4.3, 0.3, 30)
    assert ev["n_matched"] == 1 and ev["stage"] == "shoulder" and ev["sub"] == "no"


def test_unmatched_passing_run_is_spurious():
    st = _star()
    f = np.ones(len(st["TIME"]))
    fake = {"t_start": 104.0, "t_end": 104.5, "pass_local": True}
    for k in es.ORDER:
        fake[k] = True
    ev = es.classify(st, st["SAP_FLUX"], st["SAP_FLUX_ERR"], [fake], set(), f, 4.0, 0.3, 30)
    assert ev["spurious"] and not ev["recovered"] and ev["n_matched"] == 0

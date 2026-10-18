from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w3_ztf as wz  # noqa: E402

Q = (1, 1, 2)
COV = {Q: {"dimfrac": {}}}


def mk(n=80, nd_at=(), dim_at=None):
    dim_at = dim_at or {}
    return [(float(i), dim_at.get(i, 0.0), i in nd_at, 3.0, Q, i, 16.0) for i in range(n)]


def test_runs_and_flag():
    assert wz.runs(np.array([0, 1, 1, 0, 1], bool)) == [(1, 2), (4, 4)]
    s = mk(nd_at={40, 41, 42})
    assert wz.flagged_runs(s) == [(40, 42)]
    assert wz.flagged_runs(mk(nd_at={40, 41})) == []  # too short
    assert wz.flagged_runs(mk(nd_at={5, 6, 7})) == []  # not bracketed by 20 baseline epochs


def test_vet_rules():
    s = mk(nd_at={40, 41, 42})
    v = wz.vet(s, (40, 42), [], COV)
    assert v["pass"]
    # gradual recovery: intermediate detections after the run
    s2 = mk(nd_at={40, 41, 42}, dim_at={43: 0.6, 44: 0.5})
    assert not wz.vet(s2, (40, 42), [], COV)["sharp_ok"]
    # a second dim run elsewhere
    s3 = mk(nd_at={40, 41, 42, 60, 61, 62})
    assert not wz.vet(s3, (40, 42), [], COV)["oneoff_ok"]
    # common mode
    assert not wz.vet(s, (40, 42), [], COV, {40: 0.5, 41: 0.5, 42: 0.5})["common_mode_ok"]
    # other filter bright inside the run
    o = [(41.0, 0.0, False, 3.0, Q, 41, 16.0)]
    assert not wz.vet(s, (40, 42), o, COV)["achromatic_ok"]


def test_apply_dmag_limit():
    s = mk(n=3)
    out = wz.apply_dmag(s, np.array([0.5, 5.0, np.inf]))
    assert out[0][1] == 0.5 and not out[0][2]
    assert out[1][2] and out[2][2]  # beyond the limit (3 mag) or zero flux -> non-detection


def test_exposure_match_straddles_old_bin_edges():
    # one exposure seen at slightly different hmjd by two stars across a 0.0005-d bin edge
    et = np.array([100.000255, 101.50005])
    idx = wz.match_exposures([100.00024, 100.00027, 101.5, 101.50009, 103.0], et)
    assert list(idx) == [0, 0, 1, 1, -1]


def test_night_boundary_keeps_one_night():
    # three epochs within one Palomar night (MJD fractions 0.15-0.45) count as one night
    s = mk(n=60, nd_at={30, 31, 32})
    s = [
        (58000.15 + 0.15 * (i - 30), *x[1:]) if 30 <= i <= 32 else (57900.0 + i, *x[1:])
        for i, x in enumerate(s)
    ]
    s.sort(key=lambda x: x[0])
    assert wz.flagged_runs(s) == []


def _toy(n_exp=60, n_comp=10):
    """One quadrant, filter r: n_exp exposures 1 d apart; n_comp comparison stars and one target."""
    from astropy.table import Table

    obsjd = 2459000.5 + np.arange(n_exp, dtype=float)
    expo = Table(
        {
            "field": [1] * n_exp,
            "ccdid": [1] * n_exp,
            "qid": [1] * n_exp,
            "filtercode": ["zr"] * n_exp,
            "obsjd": obsjd,
            "expid": np.arange(n_exp) + 100,
            "infobits": [0] * n_exp,
            "maglimit": [20.5] * n_exp,
        }
    )
    et = wz.exposure_times(expo, 10.0, 60.0)

    def star(times, m=16.0):
        k = len(times)
        return {
            "filt": 2,
            "quad": (1, 0),
            "ra": 10.0,
            "dec": 60.0,
            "t": np.asarray(times, float),
            "m": np.full(k, m),
            "e": np.full(k, 0.02),
            "f": np.zeros(k, np.int64),
        }

    rows = {i: star(et) for i in range(n_comp)}
    return expo, et, rows, star


def test_v3_duplicates_count_once_and_dimfrac_in_range():
    expo, et, rows, star = _toy()
    # comparison star 0 has every detection twice (two catalogue matches in one exposure)
    rows[0] = star(np.repeat(et, 2))
    cov = wz.coverage(rows, set(range(10)), set(), expo)
    c = cov[(1, 0, 2)]
    assert c["clean"].max() <= c["inplay"].max() == 10


def test_v3_target_not_its_own_comparison_and_span():
    expo, et, rows, star = _toy()
    # target 50 is also a comparison star, missing at exposures 30-32;
    # star 9 only exists after exposure 40
    tgt = np.delete(et, [30, 31, 32])
    rows[50] = star(tgt)
    rows[9] = star(et[41:])
    cov = wz.coverage(rows, set(range(10)) | {50}, {50}, expo)
    s = wz.series({"members": [50]}, rows, cov, 2)
    nd = [x for x in s if x[2]]
    assert len(nd) == 3
    # dimfrac at the target gap: target excluded, star 9 (not yet in play) not missing
    assert all(x[7] == 0.0 for x in nd)


def test_v3_infobits_and_state_priority():
    expo, et, rows, star = _toy()
    expo["infobits"][5] = 4
    # star 0: a dim and a clean detection in exposure 7; the clean one wins, with its magnitude
    t0 = np.append(et, et[7])
    r0 = star(t0)
    r0["m"][-1] = 18.0
    rows[0] = r0
    cov = wz.coverage(rows, set(range(10)), set(), expo)
    c = cov[(1, 0, 2)]
    assert not c["good"][5] and c["good"][6]
    st, med, mag = c["states"][0]
    assert st[7] == 1 and mag[7] == 16.0


def test_v3_no_nd_outside_target_span():
    expo, et, rows, star = _toy()
    rows[50] = star(et[10:])  # the target's light curve starts at exposure 10
    cov = wz.coverage(rows, set(range(10)), set(), expo)
    s = wz.series({"members": [50]}, rows, cov, 2)
    assert not any(x[2] for x in s) and min(x[0] for x in s) > et[10]


def test_apply_dmag_noise_grows_with_depth():
    s = [(float(i), 0.0, False, 4.0, Q, i, 16.0, 0.0) for i in range(400)]
    rng = np.random.default_rng(0)
    out = wz.apply_dmag(s, np.full(400, 3.0), rng=rng)
    d = np.array([x[1] for x in out])
    assert abs(d.mean() - 3.0) < 0.02 and 0.05 < d.std() < 0.2  # 1 mag above the limit: sigma ~0.13

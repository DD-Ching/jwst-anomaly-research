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


def _corners(expo, ra, dec, half_deg):
    """Square footprints of half-size ``half_deg`` around (ra, dec), corners in mixed order."""
    c = np.cos(np.radians(dec))
    pts = [(-1, -1), (1, 1), (1, -1), (-1, 1)]
    for i, (sx, sy) in enumerate(pts, 1):
        expo[f"ra{i}"] = np.full(len(expo), ra + sx * half_deg / c)
        expo[f"dec{i}"] = np.full(len(expo), dec + sy * half_deg)


def test_v4_footprint_test():
    from astropy.table import Table

    expo = Table({"expid": [1, 2]})
    _corners(expo, 359.99, 20.0, 0.05)  # wraps through RA 0
    assert wz.in_quad(expo, 359.99, 20.0).all()
    assert not wz.in_quad(expo, 359.99, 20.0 + 0.05 - 5 / 3600).any()  # within the 10" margin
    assert not wz.in_quad(expo, 0.2, 20.0).any()
    assert wz.in_quad(Table({"expid": [1]}), 0.0, 0.0).all()  # no corners: pre-v4 behaviour


def test_v4_no_nd_where_star_is_off_the_quadrant():
    # the target is outside the quadrant at exposures 30-32 (pointing moved it to a neighbour)
    expo, et, rows, star = _toy()
    _corners(expo, 10.0, 60.0, 0.2)
    for k in (30, 31, 32):
        for i in range(1, 5):
            expo[f"dec{i}"][k] += 0.5
    rows[50] = star(np.delete(et, [30, 31, 32]))
    cov = wz.coverage(rows, set(range(10)), set(), expo)
    s = wz.series({"members": [50]}, rows, cov, 2)
    assert not any(x[2] for x in s) and len(s) == 57 - 2  # no ND; first/last exposure out of span


def test_v4_objectid_quadrant_ambiguity_and_dedupe():
    assert wz.oid_quad(826215400000661) == (826, 59, 2)
    assert wz.oid_quad(1856105100003845) == (1856, 16, 1)
    assert not wz.ambiguous([826215400000661, 826115400000374, 1857208200005132])
    assert wz.ambiguous([826215400000661, 826215400000999])
    a = (1.0, 0.0, True, 3.0, Q, 7, 16.0, 0.0)
    b = (1.0001, 0.2, False, 3.0, Q, 7, 16.0, 0.0)
    assert wz.dedupe([a, b]) == [b] and wz.dedupe([b, a]) == [b]


def test_v4_footprint_great_circle_edges_at_high_dec():
    from astropy.table import Table

    c = dict(ra1=99.0, dec1=65.3, ra2=101.0, dec2=65.3, ra3=101.0, dec3=64.3, ra4=99.0, dec4=64.3)
    expo = Table({k: [v] for k, v in c.items()})
    # the northern edge is a great circle bowing ~12" poleward of dec 65.3 at RA 100
    assert wz.in_quad(expo, 100.0, 65.3 + 5 / 3600, 0.0).all()
    assert not wz.in_quad(expo, 100.0, 65.3 + 20 / 3600, 0.0).any()
    expo["dec2"][0] = np.nan
    assert not wz.in_quad(expo, 100.0, 65.0, 0.0).any()  # a non-finite corner counts as outside


def test_v8_exposure_mid_time_uses_exptime():
    from astropy.table import Table

    # a 299-s private-survey exposure: the detection sits at the mid-time, 149.5 s after the start
    obsjd = 2460065.7368519
    expo = Table({"obsjd": [obsjd], "exptime": [299.0]})
    et = wz.exposure_times(expo, 184.0, 33.0)
    t_det = wz.exposure_times(
        Table({"obsjd": [obsjd + 134.5 / 86400.0]}), 184.0, 33.0
    )  # 30-s rule + 134.5 s
    assert wz.match_exposures(t_det, et)[0] == 0
    assert (
        wz.match_exposures(t_det, wz.exposure_times(Table({"obsjd": [obsjd]}), 184.0, 33.0))[0]
        == -1
    )


def test_v10_psfcat_match_unique_missing_ambiguous():
    ra, dec = np.array([10.0, 10.0 + 2.0 / 3600]), np.array([60.0, 60.0])
    c = np.cos(np.radians(60.0))
    idx = wz.psfcat_match(ra, dec, [10.0, 10.0 + 1.0 / 3600 / c, 11.0], [60.0, 60.0, 60.0])
    assert idx[0] == -2  # second source 1" away: ambiguous
    assert idx[2] == -1  # nothing within 1.5": missing
    one = wz.psfcat_match(ra[:1], dec[:1], [10.0 + 1.0 / 3600 / c], [60.0])
    assert one[0] == 0


def test_v10_merge_psfcat_fills_without_overwriting_or_inventing_nd():
    expo, et, rows, star = _toy()
    gap = [30, 31, 32, 33, 34]
    for i in range(10):  # comparison stars: DR24 lacks the gap exposures
        rows[i] = star(np.delete(et, gap))
    tgt = star(np.delete(et, [30, 32, 33, 34]))  # the target keeps its DR24 detection at 31
    tgt["m"][30] = 16.3
    rows[50] = tgt
    src = (np.array([10.0]), np.array([60.0]), np.array([16.0]), np.array([0.02]))
    cats = {}
    for k in gap[:-1]:  # exposure 34: no source at all (missing for every star)
        cats[(1, 0, 2, int(expo["expid"][k]))] = (*src, np.array([4 if k == 33 else 0]))
    cats[(1, 0, 2, int(expo["expid"][34]))] = tuple(np.array([]) for _ in range(5))
    log = wz.merge_psfcat(rows, cats, expo)
    assert log[(1, 0, 2)]["dr24_kept"] == 1 and log[(1, 0, 2)]["missing"] == 11
    cov = wz.coverage(rows, set(range(10)), {50}, expo)
    s = wz.series({"members": [50]}, rows, cov, 2)
    ids = {x[5]: x for x in s}
    e = [int(expo["expid"][k]) for k in gap]
    assert e[0] in ids and e[2] in ids and abs(ids[e[1]][1] - 0.3) < 1e-6  # DR24 value kept
    assert e[3] not in ids and e[4] not in ids  # flagged: unknown; missing: unknown, not ND
    assert not any(x[2] for x in s)
    wz.FILL.clear()


def test_v13_edge_time_uses_time_not_count():
    # 20 dim epochs (days 40-59): dense intermediates at the edges pass, one far out fails
    wz.EDGE_RULE, wz.EDGE_FRAC = "time", 0.25
    try:
        dim = {i: 0.0 for i in range(40, 60)}
        s = mk(n=120, nd_at=set(dim))
        s = [(*x, 0.0) for x in s]  # dimfrac column
        edge = {37: 0.6, 38: 0.7, 39: 0.8, 60: 0.8, 61: 0.6, 62: 0.5}  # 6 intermediates within 3 d
        s_ok = [(x[0], edge.get(i, x[1]), x[2], 6.0, *x[4:]) for i, x in enumerate(s)]
        v = wz.vet(s_ok, (40, 59), [], COV)
        assert v["n_intermediate"] == 6 and v["sharp_ok"] and v["edge_post_d"] == 3.0
        far = {**edge, 70: 0.6}  # 11 d after a 19 d run: > 0.25 x 19 d
        s_bad = [(x[0], far.get(i, x[1]), x[2], 6.0, *x[4:]) for i, x in enumerate(s)]
        assert not wz.vet(s_bad, (40, 59), [], COV)["sharp_ok"]
        # an insignificant intermediate (within 2 sigma of 0.3 at a shallow limit) does not count
        noisy = [(x[0], 0.35, x[2], 0.5, *x[4:]) if i == 70 else x for i, x in enumerate(s_ok)]
        assert wz.vet(noisy, (40, 59), [], COV)["edge_post_d"] == 3.0
    finally:
        wz.EDGE_RULE = "count"


def test_v23_edge_bounds_gap_is_insufficient_not_fast():
    # run on days 40-59 (span 19 d, threshold 4.75 d); baseline detections elsewhere
    wz.EDGE_RULE, wz.EDGE_FRAC = "bounds", 0.25
    try:
        dim = set(range(40, 60))
        base = [(*x, 0.0) for x in mk(n=120, nd_at=dim)]
        base = [(x[0], x[1], x[2], 6.0, *x[4:]) for x in base]
        # dense sampling, no intermediate: both edges bounded by 1 d -> fast
        v = wz.vet(base, (40, 59), [], COV)
        assert v["edge_class"] == "fast" and v["sharp_ok"] and v["edge_post_max_d"] == 1.0
        assert not v["insufficient"] and v["pass"]
        # a 10-d gap after the run (days 60-69 absent), no intermediate: insufficient, not fast
        gap = [x for i, x in enumerate(base) if not 60 <= i < 70]
        v = wz.vet(gap, (40, 59), [], COV)
        assert v["edge_post_d"] == 0.0 and v["edge_post_max_d"] == 11.0
        assert v["edge_post_class"] == "insufficient" and v["edge_pre_class"] == "fast"
        assert v["edge_class"] == "insufficient" and not v["sharp_ok"]
        assert v["insufficient"] and not v["pass"]
        # the same gap under v13 records 0 d and passes
        wz.EDGE_RULE = "time"
        assert wz.vet(gap, (40, 59), [], COV)["sharp_ok"]
        wz.EDGE_RULE = "bounds"
        # a significant intermediate 8 d after the run: slow (L > threshold), as in v13
        slow = [(x[0], 0.6, *x[2:]) if i == 67 else x for i, x in enumerate(base)]
        v = wz.vet(slow, (40, 59), [], COV)
        assert v["edge_post_class"] == "slow" and v["edge_class"] == "slow"
        assert not v["insufficient"]
        # an ambiguous epoch (0.3 <= dm < 0.3 + 2 sigma) is neither normal nor intermediate
        amb = [(x[0], 0.35, x[2], 0.5, *x[4:]) if i == 60 else x for i, x in enumerate(base)]
        v = wz.vet(amb, (40, 59), [], COV)
        assert v["edge_post_d"] == 0.0 and v["edge_post_max_d"] == 2.0 and v["sharp_ok"]
    finally:
        wz.EDGE_RULE = "count"


def test_altdata_edge_bounds_alt():
    # run on days 40-59 (span 19 d, threshold 4.75 d); a 10-d r gap after the run (days 60-69)
    wz.EDGE_FRAC = 0.25
    dim = set(range(40, 60))
    base = [(*x, 0.0) for x in mk(n=120, nd_at=dim)]
    base = [(x[0], x[1], x[2], 6.0, *x[4:]) for x in base]
    gap = [x for i, x in enumerate(base) if not 60 <= i < 70]
    slow = [(x[0], 0.6, *x[2:]) if i == 67 else x for i, x in enumerate(base)]
    # no alt point: identical to edge_bounds
    for s in (base, gap, slow):
        assert wz.edge_bounds_alt(s, (40, 59), 19.0, []) == wz.edge_bounds(s, (40, 59), 19.0)
    # a normal alt point (S/N 10) 3 d after the run fills the gap: fast
    v = wz.edge_bounds_alt(gap, (40, 59), 19.0, [(62.0, 0.05, 0.1)])
    assert (
        v["edge_post_class"] == "fast" and v["edge_post_max_d"] == 3.0 and v["edge_class"] == "fast"
    )
    # the same point at S/N < 3, or 6 d after the run (> 4.75 d), does not decide the side
    v = wz.edge_bounds_alt(gap, (40, 59), 19.0, [(62.0, 0.05, 0.5)])
    assert v["edge_post_class"] == "insufficient"
    v = wz.edge_bounds_alt(gap, (40, 59), 19.0, [(65.0, 0.05, 0.1)])
    assert v["edge_post_class"] == "insufficient" and v["edge_post_max_d"] == 6.0
    # a significant alt intermediate 8 d after the run: slow; an ambiguous one does not decide
    v = wz.edge_bounds_alt(gap, (40, 59), 19.0, [(67.0, 0.6, 0.05)])
    assert v["edge_post_class"] == "slow" and v["edge_post_d"] == 8.0
    v = wz.edge_bounds_alt(gap, (40, 59), 19.0, [(67.0, 0.35, 0.05)])
    assert v["edge_post_class"] == "insufficient"
    # alt points inside the run are never used, and an alt normal point nearer than a significant
    # intermediate does not bound the edge (as for r epochs in edge_bounds)
    v = wz.edge_bounds_alt(gap, (40, 59), 19.0, [(50.0, 0.0, 0.01)])
    assert v == wz.edge_bounds(gap, (40, 59), 19.0)
    v = wz.edge_bounds_alt(gap, (40, 59), 19.0, [(61.0, 0.0, 0.05), (63.0, 0.7, 0.05)])
    assert v["edge_post_d"] == 4.0 and v["edge_post_max_d"] == 11.0
    assert v["edge_post_class"] == "insufficient"
    # an alt point without a finite positive sigma is not used
    v = wz.edge_bounds_alt(gap, (40, 59), 19.0, [(67.0, 0.6, 0.0), (62.0, 0.05, np.nan)])
    assert v == wz.edge_bounds(gap, (40, 59), 19.0)
    # pre side: a 10-d r gap before the run (days 30-39); an alt normal point 2 d before decides it,
    # a significant alt intermediate 7 d before makes it slow
    gap_pre = [x for i, x in enumerate(base) if not 30 <= i < 40]
    assert wz.edge_bounds_alt(gap_pre, (30, 49), 19.0, [])["edge_pre_class"] == "insufficient"
    v = wz.edge_bounds_alt(gap_pre, (30, 49), 19.0, [(38.0, 0.0, 0.05)])
    assert v["edge_pre_class"] == "fast" and v["edge_pre_max_d"] == 2.0
    v = wz.edge_bounds_alt(gap_pre, (30, 49), 19.0, [(33.0, 0.6, 0.05)])
    assert v["edge_pre_class"] == "slow" and v["edge_pre_d"] == 7.0


def test_joint_extent_and_edge_bounds_joint():
    # r: run on days 40-59 (NDs, limit 6), baseline elsewhere, no r epochs on days 60-69
    wz.EDGE_FRAC = 0.25
    base = [(*x, 0.0) for x in mk(n=120, nd_at=set(range(40, 60)))]
    base = [(x[0], x[1], x[2], 6.0, *x[4:]) for x in base]
    gap = [x for i, x in enumerate(base) if not 60 <= i < 70]
    run = (40, 59)

    def g(t, dm, nd=False, lim=6.0):
        return (float(t), dm, nd, lim, (1, 1, 1), 10_000 + int(t * 10), 16.0, 0.0)

    # g dim at 61 and 63, g normal at 64: the edge moves to 63, span 23 d, U 1 d -> fast;
    # with the g points only as edge points (#459) U is 5 d > 4.75 d -> insufficient
    other = [g(61, 1.5), g(63, 2.0), g(64, 0.0)]
    assert wz.joint_extent(gap, run, other) == (40.0, 63.0)
    v = wz.edge_bounds_joint(gap, run, other)
    assert v["joint_span"] == 23.0 and v["edge_post_max_d"] == 1.0 and v["edge_class"] == "fast"
    b = wz.edge_bounds_alt(gap, run, 19.0, wz.alt_detections(other))
    assert b["edge_post_class"] == "insufficient"
    # an intermediate or ambiguous g detection stops the extension, and a significant one counts
    # towards a slow edge
    other = [g(61, 0.6), g(62, 2.0), g(68, 0.0)]
    assert wz.joint_extent(gap, run, other) == (40.0, 59.0)
    v = wz.edge_bounds_joint(gap, run, other)
    assert v["edge_post_d"] == 2.0 and v["edge_post_class"] == "insufficient"
    assert wz.joint_extent(gap, run, [g(61, 0.35, lim=0.5), g(62, 2.0)])[1] == 59.0
    # a non-detection extends only when its limit implies dm >= 1; a shallow one is unknown
    assert wz.joint_extent(gap, run, [g(61, 0.0, True, 3.0)])[1] == 61.0
    assert wz.joint_extent(gap, run, [g(61, 0.0, True, 0.5), g(62, 1.5)])[1] == 62.0
    assert wz.joint_extent(gap, run, [g(61, 0.0, True, 0.5)])[1] == 59.0
    # at equal times a stop comes first; the extension never passes an r detection (day 70)
    assert wz.joint_extent(gap, run, [g(61, 0.0), g(61, 2.0)])[1] == 59.0
    assert wz.joint_extent(gap, run, [g(71, 2.0)])[1] == 59.0
    # pre side: a 10-d r gap before the run; g dim at 37 extends it, g normal at 36 bounds it
    gap_pre = [x for i, x in enumerate(base) if not 30 <= i < 40]
    v = wz.edge_bounds_joint(gap_pre, (30, 49), [g(36, 0.0), g(37, 1.2)])
    assert v["joint_lo"] == 37.0 and v["edge_pre_max_d"] == 1.0 and v["edge_pre_class"] == "fast"
    # pre side: a stop first at equal times, a shallow non-detection unknown, a deep one dim
    assert wz.joint_extent(gap_pre, (30, 49), [g(38, 0.0), g(38, 2.0)])[0] == 40.0
    assert wz.joint_extent(gap_pre, (30, 49), [g(38, 0.0, True, 0.5), g(37, 1.5)])[0] == 37.0
    assert wz.joint_extent(gap_pre, (30, 49), [g(38, 0.0, True, 0.5)])[0] == 40.0
    assert wz.joint_extent(gap_pre, (30, 49), [g(38, 0.0, True, 3.0), g(36, 0.0)])[0] == 38.0
    # no other band: identical to bounds
    v = wz.edge_bounds_joint(gap, run, [])
    w = wz.edge_bounds(gap, run, 19.0)
    assert {k: v[k] for k in w} == w


def test_v15_configure_twice_does_not_stack():
    saved = (wz.LC_PIX, wz.LC_PART, wz.OBJ_PART, wz.CACHE, wz.SAMPLE, wz.VDIR, wz.VERSION, wz.N_MAX)
    try:
        wz.configure(3837, "v15/3837")
        wz.configure(183, "v15/183")
        assert "Npix=183/" in wz.LC_PART and "Npix=45/" in wz.OBJ_PART
        assert wz.CACHE == wz._BASE[2] / "183"
        wz.configure(3839, "v5")  # back to the v1-v3 partition: defaults restored
        assert "Npix=3839/" in wz.LC_PART and "Npix=959/" in wz.OBJ_PART and wz.N_MAX == 1000
    finally:
        (wz.LC_PIX, wz.LC_PART, wz.OBJ_PART, wz.CACHE, wz.SAMPLE, wz.VDIR, wz.VERSION, wz.N_MAX) = (
            saved
        )

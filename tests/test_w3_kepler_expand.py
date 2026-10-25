"""Offline tests of scripts/w3_kepler_expand.py (results/w3_kepler/expand/PLAN.md)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w3_kepler as wk  # noqa: E402
import w3_kepler_expand as we  # noqa: E402


def _star(seed, dip=None, n_q=2, n=400):
    rng = np.random.default_rng(seed)
    t = np.concatenate([q * (n * wk.CAD + 0.5) + np.arange(n) * wk.CAD for q in range(n_q)])
    sap = 1000 + rng.normal(0, 0.1, len(t))
    if dip is not None:
        sap[dip[0] : dip[1]] = 100.0
    return {
        "TIME": t,
        "CADENCENO": np.arange(len(t), dtype=float),
        "SAP_FLUX": sap,
        "SAP_FLUX_ERR": np.full(len(t), 0.1),
        "PDCSAP_FLUX": sap.copy(),
        "SAP_QUALITY": np.zeros(len(t)),
        "QUARTER": np.repeat(np.arange(1, n_q + 1, dtype=float), n),
    }


def test_frac_loo_matches_pilot_frac_for():
    sts = {1: _star(1, (100, 110)), 2: _star(2, (100, 120)), 3: _star(3, (105, 108)), 4: _star(4)}
    cf = wk.common_frac(sts)
    den = np.zeros(80000, np.int32)
    num = np.zeros(80000, np.int32)
    for st in sts.values():
        g, d = we.own_cadences(st)
        np.add.at(den, g, 1)
        np.add.at(num, d, 1)
    for k, st in sts.items():
        assert we.frac_loo(den, num, st) == wk.frac_for(cf, k)
    # a star missing some cadences: equal on its own cadences (the only ones the screen looks up)
    short = _star(5, (100, 110))
    short = {c: v[:600] for c, v in short.items()}
    cf2 = wk.common_frac({**sts, 5: short})
    den2, num2 = den.copy(), num.copy()
    g, d = we.own_cadences(short)
    np.add.at(den2, g, 1)
    np.add.at(num2, d, 1)
    own = set(short["CADENCENO"].astype(int).tolist())
    ref = {c: v for c, v in wk.frac_for(cf2, 5).items() if c in own}
    assert we.frac_loo(den2, num2, short) == ref


def test_layers_order_and_reasons():
    def run(c, s, o, loc):
        return {
            "common_mode_ok": c,
            "sharp_ok": s,
            "oneoff_ok": o,
            "local_ok": loc,
            "pass_local": c and s and o and loc,
        }

    res = [
        {
            "kepid": 10,
            "runs": [run(True, True, True, True)],
            "screened": True,
            "applicable": True,
            "effective_days": 1.0,
            "dim_cadences": 5,
        },
        {
            "kepid": 11,
            "runs": [run(True, False, True, False)],
            "screened": True,
            "applicable": True,
            "effective_days": 1.0,
            "dim_cadences": 3,
        },
        {
            "kepid": 12,
            "runs": [],
            "screened": True,
            "applicable": False,
            "effective_days": 1.0,
            "dim_cadences": 0,
        },
    ]
    lay = we.layers(res, {11})
    assert lay["L2b_flagged_runs"] == 2 and lay["L2b_flagged_sources"] == 2
    assert lay["L3_sequential"]["after_sharp_ok"]["runs"] == 1
    assert lay["L4_survivor_sources_in_order"] == [10] and lay["L4_vet_list"] == [10]
    assert lay["flagged_source_reasons"]["11"]["first_failed"] == "sharp_ok"
    assert lay["flagged_source_reasons"]["11"]["in_kirk2016_eb_catalogue"]
    assert lay["L1_applicable"] == 2 and lay["L2a_stars_with_dim_cadence"] == 2

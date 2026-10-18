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
    return [(float(i), dim_at.get(i, 0.0), i in nd_at, 3.0, Q, i) for i in range(n)]


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
    o = [(41.0, 0.0, False, 3.0, Q, 41)]
    assert not wz.vet(s, (40, 42), o, COV)["achromatic_ok"]


def test_apply_dmag_limit():
    s = mk(n=3)
    out = wz.apply_dmag(s, np.array([0.5, 5.0, np.inf]))
    assert out[0][1] == 0.5 and not out[0][2]
    assert out[1][2] and out[2][2]  # beyond the limit (3 mag) or zero flux -> non-detection

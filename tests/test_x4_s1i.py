"""Offline tests for scripts/x4_s1i_dynesty.py (X4 S1i): the search-distribution mapping, the prior
transform, the run status rule and the frozen decision. No lalsuite or dynesty needed."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

_S = Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _S / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


s = _load("x4_s1i_dynesty")


def test_loglike_is_strictly_increasing_in_ff_and_keeps_the_maximum():
    ff = np.linspace(0.0, 1.0, 1001)
    ll = np.array([s.loglike_from_ff(v) for v in ff])
    assert np.all(np.diff(ll) > 0)
    assert s.loglike_from_ff(1.0) == 0.0
    assert s.loglike_from_ff(-0.5) == s.loglike_from_ff(0.0)  # invalid points share the floor


def test_prior_transform_maps_the_unit_cube_onto_the_box():
    assert np.allclose(s.prior_transform(np.zeros(s.NDIM)), s.LO)
    assert np.allclose(s.prior_transform(np.ones(s.NDIM)), s.HI)
    assert s.PERIODIC == [4, 7, 9]
    assert np.allclose(s.HI[s.PERIODIC], 2 * np.pi) and np.allclose(s.LO[s.PERIODIC], 0.0)


def test_status_rule_and_targets():
    assert s.counted({"stop": "dlogz"}) and s.counted({"stop": "maxcall"})
    assert not s.counted({"stop": "deadline"})
    assert abs(s.FF_TARGET - (0.9135105524586078 - 2e-4)) < 1e-15
    assert s.EPS_PASS == 5e-4


def test_decide_outcomes(tmp_path, monkeypatch):
    monkeypatch.setattr(s, "OUT", tmp_path)
    monkeypatch.setattr(s, "verify", lambda x: {"stub": True})

    def write(name, rows):
        (tmp_path / name).write_text(json.dumps({"rows": rows}))

    good = {"stop": "maxcall", "ncall": 1}
    write(
        "s1i_validate.json",
        [{**good, "name": "V7", "eps_final": 1e-6}, {**good, "name": "V10", "eps_final": 4e-4}],
    )
    assert s.decide()["validation_pass"]
    write(
        "s1i_r1.json",
        [
            {**good, "seed": 31, "FF_pad64": 0.9134, "x": []},
            {**good, "seed": 32, "FF_pad64": 0.95, "x": []},
        ],
    )
    d = s.decide()
    assert d["outcome"].startswith("both runs recovered") and 32 in d["verified"]
    write(
        "s1i_r1.json",
        [
            {**good, "seed": 31, "FF_pad64": 0.8923, "x": []},
            {**good, "seed": 32, "FF_pad64": 0.9134, "x": []},
        ],
    )
    assert s.decide()["outcome"] == "method did not pass this test (1 of 2 runs reached the target)"
    write(
        "s1i_r1.json",
        [
            {**good, "seed": 31, "FF_pad64": 0.95, "x": []},
            {"stop": "deadline", "ncall": 1, "seed": 32, "FF_pad64": 0.95, "x": []},
        ],
    )
    assert s.decide()["outcome"] == "unfinished"
    write(
        "s1i_validate.json",
        [{**good, "name": "V7", "eps_final": 1e-6}, {**good, "name": "V10", "eps_final": 6e-4}],
    )
    assert not s.decide()["validation_pass"]

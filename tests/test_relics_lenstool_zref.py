"""Offline checks of the RELICS Lenstool v1 z_ref test (PLAN.md H_U rule)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import relics_lenstool_zref as rz  # noqa: E402


def test_dls_ds_reproduces_a3192_zref():
    # the Abell 3192 factor fixes z_ref, so the prediction at its own z_lens returns it
    assert rz.dls_ds(0.425, rz.Z_REF) == pytest.approx(0.6918, abs=1e-3)
    assert rz.dls_ds(0.425, 1.0) < rz.dls_ds(0.425, 3.0) < 1.0


def test_h_u_rule(monkeypatch):
    pred = rz.dls_ds(0.412, rz.Z_REF)
    ok = rz.h_u({"macs0553": pred + 0.09})
    assert ok["rows"]["macs0553"]["consistent"] and not ok["rejected"]
    bad = rz.h_u({"macs0553": pred + 0.11})
    assert bad["rejected"]
    # one inconsistent cluster rejects H_U even if another agrees
    mixed = rz.h_u({"macs0553": pred, "cl0152": 1.045})
    assert mixed["rejected"]


def test_stored_factors_cover_five_clusters():
    f = rz.stored_factors()
    assert set(f) == {"rxc0232", "rxc2211", "macs0553", "cl0152", "a3192"}
    for c in f:
        assert f"{c}-relics-lenstool" in rz.lc.MODELS

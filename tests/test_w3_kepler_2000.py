"""Offline tests of scripts/w3_kepler_2000.py (results/w3_kepler/s2000/PLAN.md)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w3_kepler as wk  # noqa: E402
import w3_kepler_2000 as s2  # noqa: E402


def _star(n_q=2, n=400, depth_at=None, ch=7):
    t = np.concatenate([q * 100.0 + np.arange(n) * wk.CAD for q in range(n_q)])
    sap = np.full(len(t), 1000.0)
    if depth_at is not None:
        for t1, t2 in depth_at:
            sap[(t >= t1) & (t < t2)] = 100.0
    return {
        "TIME": t,
        "CADENCENO": np.arange(len(t), dtype=float),
        "SAP_FLUX": sap,
        "SAP_FLUX_ERR": np.full(len(t), 0.1),
        "PDCSAP_FLUX": sap.copy(),
        "SAP_QUALITY": np.zeros(len(t)),
        "QUARTER": np.repeat(np.arange(1, n_q + 1, dtype=float), n),
        "CHANNEL": np.full(len(t), float(ch)),
        "MOM_CENTR1": np.full(len(t), 500.0),
        "MOM_CENTR2": np.full(len(t), 600.0),
    }


def test_select_excludes_used_and_is_deterministic():
    pool = list(range(1000, 1100))
    used = {1000, 1001, 1050}
    a = s2.select(pool, used, n=20)
    assert a == s2.select(list(reversed(pool)), used, n=20)
    assert len(a) == 20 and not set(a) & used


def test_q0_names_are_skipped():
    assert s2.is_q0("kplr010064130-2009131105131_llc.fits")
    assert not s2.is_q0("kplr010064130-2009166043257_llc.fits")


def test_joint_limit_classical_and_wider():
    t = 2000.0
    g_big = s2.joint_limit(t, 0.8, 1e6)
    assert abs(g_big * t * 0.8 / -np.log(0.05) - 1) < 1e-3  # exact efficiency -> 2.996 / E
    g = s2.joint_limit(t, 0.8, 1000)
    assert g > g_big and g < -np.log(0.05) / (t * 0.75)  # wider, but < 2.996 / E at a low edge
    assert s2.joint_limit(t, 0.0, 1000) is None


def test_vet_common_mode_is_artifact():
    run_span = (4.0, 4.5)
    st = _star(depth_at=[run_span])
    others = {i: _star(depth_at=[run_span]) for i in range(3)}
    run = {"t_start": 4.0, "t_end": 4.49, "quarter_start": 1, "pdc_rel_in_run": 0.1}
    d = s2.vet_diagnostics(st, run, others, in_kirk=False)
    assert d["V3_common"] and s2.vet_class(d) == "artifact_supported"


def test_vet_repeating_dips_are_ordinary():
    dips = [(1.0, 1.2), (3.0, 3.2), (6.0, 6.2)]
    st = _star(depth_at=dips + [(4.0, 4.5)])
    run = {"t_start": 4.0, "t_end": 4.49, "quarter_start": 1, "pdc_rel_in_run": 0.1}
    d = s2.vet_diagnostics(st, run, {1: _star()}, in_kirk=False)
    assert d["V5_other_dips_ge_0.3mag"] >= 2 and d["V5_repeating"]
    assert s2.vet_class(d) == "ordinary_explanation_supported"


def test_vet_class_needs_pixels_for_confirmation():
    d = {k: False for k in ("V1_quality", "V3_common", "V5_repeating", "V6_flare_like")}
    d["V7_baseline_step"] = False
    assert s2.vet_class(d) == "undetermined"
    assert s2.vet_class(d, "on_target") == "real_dimming_confirmed_unexplained_by_V5_V7"
    assert s2.vet_class(d, "off_target") == "artifact_supported"

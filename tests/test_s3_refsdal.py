import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location(
    "s3_refsdal", Path(__file__).parents[1] / "scripts" / "s3_refsdal.py"
)
s3r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s3r)


def test_flagged_rule_needs_both_filters_and_consistency():
    r1 = np.array([0.5, 0.5, 0.5, np.nan])
    s1 = np.array([0.05, 0.05, 0.05, 0.05])
    r2 = np.array([0.5, 0.1, 0.9, 0.5])
    s2 = np.array([0.05, 0.05, 0.05, 0.05])
    # lag 0: both 10 sigma, equal -> flag; lag 1: second 2 sigma -> no;
    # lag 2: inconsistent -> no; lag 3: nan -> no
    assert s3r.flagged(r1, s1, r2, s2).tolist() == [True, False, False, False]


def test_scan_recovers_injected_copy():
    t = np.arange(0.0, 900.0, 5.0)
    f1 = np.exp(-0.5 * ((t - 300.0) / 60.0) ** 2)
    tpl = s3r.template(t, f1)
    tj = np.arange(0.0, 900.0, 7.0)
    fj = 2.0 * tpl(tj - 20.0) + 0.3 * 2.0 * tpl(tj - 20.0 - 150.0)
    e = np.full_like(tj, 0.01)
    tau, _ = s3r.fit_tau(tpl, tj, 2.0 * tpl(tj - 20.0), e, np.arange(-60.0, 61.0, 1.0))
    assert abs(tau - 20.0) <= 2
    r, _ = s3r.scan_delta(tpl, tj, fj, e, 20.0)
    k = int(np.argmin(abs(s3r.DELTA_GRID - 150.0)))
    assert abs(r[k] - 0.3) < 0.02

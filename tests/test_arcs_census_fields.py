import importlib.util
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_P = Path(__file__).resolve().parents[1] / "scripts" / "arcs_census_fields.py"
_spec = importlib.util.spec_from_file_location("arcs_census_fields", _P)
acf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(acf)


def test_i2d_uri_from_catalogue_name():
    u = acf.i2d_uri("jw01208-o002_t001_nircam_clear-f200w_cat.ecsv")
    assert u == (
        "s3://stpubdata/jwst/public/jw01208/L3/t/o002/jw01208-o002_t001_nircam_clear-f200w_i2d.fits"
    )
    with pytest.raises(ValueError):
        acf.i2d_uri("something_else.ecsv")


def test_strong_background_selection():
    t = Table(
        {
            "strong_shear": [True, True, False, True],
            "population": ["background", "unknown", "background", "background"],
        }
    )
    assert len(acf.strong_background(t)) == 2


def test_mixture_limits():
    p_cp, p_sn = [0.5] * 4, [0.0] * 4
    assert acf.mixture_p(p_cp, p_sn, 0.0, 1) == 0.0
    assert acf.mixture_p(p_cp, p_sn, 1.0, 4) == pytest.approx(0.5**4)
    # more clump pairs can only raise P(≥ k) when p_cp > p_sn
    ps = [acf.mixture_p([0.2] * 24, [0.02] * 24, f, 4) for f in np.linspace(0, 1, 6)]
    assert np.all(np.diff(ps) > 0)


def test_plan_outcome_threshold():
    assert acf.plan_outcome(0.05).startswith("S1063 excess explained")
    assert acf.plan_outcome(0.0499).startswith("S1063 excess remains")
    assert acf.plan_outcome(float("nan")).startswith("undetermined")


def test_summarize_counts_only_gate_passed():
    t = Table(
        {
            "gate_passed": [True, True, False, True],
            "deblend_class": ["two-knot", "single", "two-knot", "single"],
        }
    )
    s = acf.summarize(t)
    assert (s["n_gate_passed"], s["n_two_knot_gate_passed"]) == (3, 1)
    assert s["f_two_knot"] == pytest.approx(1 / 3)
    assert s["interpretable_alone"] is False

"""Offline tests for scripts/a2744_uncover.py (the WLS family source; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_uncover as au  # noqa: E402


def test_wls_source_equal_weights_is_mean():
    beta = np.array([[1.0, 2.0], [3.0, -2.0], [-1.0, 0.0]])
    jac = np.broadcast_to(np.eye(2), (3, 2, 2))
    assert np.allclose(au.wls_source(beta, jac), beta.mean(0))


def test_wls_source_follows_the_most_magnified_image():
    # image 0 magnified 10x isotropically (A = 0.1 I): its source dominates by 100:1 in weight
    beta = np.array([[0.0, 0.0], [1.0, 1.0]])
    jac = np.array([0.1 * np.eye(2), np.eye(2)])
    got = au.wls_source(beta, jac)
    assert np.allclose(got, [1.0 / 101.0, 1.0 / 101.0])


def test_wls_source_minimises_image_plane_offsets():
    rng = np.random.default_rng(1)
    beta = rng.normal(size=(4, 2))
    jac = rng.normal(size=(4, 2, 2)) + 2 * np.eye(2)
    best = au.wls_source(beta, jac)
    minv = np.linalg.inv(jac)

    def cost(b):
        return sum(np.sum((minv[i] @ (b - beta[i])) ** 2) for i in range(4))

    for d in rng.normal(scale=0.05, size=(20, 2)):
        assert cost(best) <= cost(best + d)


def test_uncover_system_ids():
    assert au.uncover_system("4.11") == "4.1"
    assert au.uncover_system("80.63") == "80.6"
    assert au.uncover_system("83.1") == "83"  # single-knot system, image 1
    assert au.uncover_system(repr(float("83.3"))) == "83"


def _xm_tables():
    from astropy.table import Table

    dec0 = -30.4
    step = 1.0 / 3600.0 / np.cos(np.deg2rad(dec0))  # 1" in RA
    berg = Table(
        {
            "system": ["1.1", "1.1", "2.1", "4.1", "4.1"],
            "ra": [3.58, 3.58 + 10 * step, 3.59, 3.57, 3.57 + 20 * step],
            "dec": [dec0] * 5,
            "magnification": [5.0, -3.0, 2.0, 40.0, 3.0],
            "image_class": ["observed", "missing", "missing", "observed", "observed"],
            "forced_class": ["", "absent", "absent", "", ""],
            "pred_snr": [np.nan, 20.0, 12.0, np.nan, np.nan],
            "flux_ratio": [np.nan, 0.05, 0.1, np.nan, np.nan],
        }
    )
    berg.meta["forced"] = {"halo_filter_arcsec": au.HALO_FILTER_ARCSEC}
    unc = Table(
        {
            "system": ["1.1", "1.1", "2.1", "3.1", "3.1", "4.1", "4.1"],
            "ra": [
                3.58,
                3.58 + 10.8 * step,
                3.59 + 2 * step,
                3.6,
                3.61,
                3.57 + 2 * step,
                3.57 + 20 * step,
            ],
            "dec": [dec0] * 7,
            "magnification": [4.0, -2.5, 2.2, 7.0, 3.0, -50.0, 3.0],
            "image_class": [
                "observed",
                "missing",
                "missing",
                "observed",
                "faint",
                "missing",
                "observed",
            ],
            "forced_class": ["", "recovered", "", "", "", "absent", ""],
            "pred_snr": [np.nan, 15.0, np.nan, np.nan, np.nan, 900.0, np.nan],
            "flux_ratio": [np.nan, 0.9, np.nan, np.nan, np.nan, 0.01, np.nan],
        }
    )
    return {"bergamini23": berg, "uncover2": unc}


def test_crossmodel_matches_extras_of_both_models():
    tables = _xm_tables()
    xm = au.crossmodel(tables)
    b = xm[xm["model"] == "bergamini23"]
    assert list(b["system"]) == ["1.1", "2.1"]
    assert np.allclose(b["other_sep_arcsec"], [0.8, 2.0], atol=1e-3)
    assert list(b["both"]) == [True, False]
    assert list(b["both_loose"]) == [True, True]
    u = xm[xm["model"] == "uncover2"]
    assert list(u["system"]) == ["1.1", "2.1", "3.1", "4.1"]
    assert list(u["other_forced_class"][:2]) == ["absent", "absent"]
    assert list(u["other_image_class"][2:]) == ["none", "none"]  # no Bergamini extra
    assert np.isnan(u["other_sep_arcsec"][2])
    assert np.isclose(u["other_observed_sep_arcsec"][3], 2.0, atol=1e-3)


def test_crossmodel_verdicts():
    import json

    tables = _xm_tables()
    s = au.crossmodel_summary(au.crossmodel(tables), tables)
    json.dumps(s, allow_nan=False)  # NaNs become null
    assert s["bergamini23"]["extra_forced"]["absent"] == {"rows": 2, "both": 1, "both_loose": 2}
    assert s["bergamini23"]["forced_all"] == {"absent": 2}
    assert s["uncover2"]["extra_forced"]["not_forced"]["rows"] == 2
    assert [(r["model"], r["system"], r["d041"]) for r in s["unmet"]] == [
        ("bergamini23", "1.1", "the other model's image is recovered (position offset; ordinary)"),
        ("bergamini23", "2.1", "both models (loose tier)"),
        ("uncover2", "4.1", "the other model's catalogued image, misplaced (ordinary)"),
    ]
    assert s["unmet_in_both_models"] == 1


def test_check_forced_table_refuses_other_filters():
    import pytest

    berg = _xm_tables()["bergamini23"]
    au.check_forced_table(berg)
    berg.meta["forced"] = {"halo_filter_arcsec": 0.6}
    with pytest.raises(SystemExit):
        au.check_forced_table(berg)

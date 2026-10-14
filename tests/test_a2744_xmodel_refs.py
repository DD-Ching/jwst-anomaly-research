"""Offline tests for scripts/a2744_xmodel_refs.py (no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_xmodel_refs as xr  # noqa: E402


def test_galaxy_of_knot():
    assert xr.galaxy("3.2") == "3"
    assert xr.galaxy("700.1") == "700"


def test_predictions_per_sibling_and_pair():
    # siblings: flux 10 at |mu| 50 and flux 30 at |mu| 50 -> pair source flux 40/100 = 0.4
    p = xr.predictions({"a": (10.0, 50.0), "b": (30.0, 50.0)}, mu_t=-5.0, e_t=0.5)
    assert np.isclose(p["from"]["a"], 10 / 50 * 5 / 0.5)
    assert np.isclose(p["from"]["b"], 30 / 50 * 5 / 0.5)
    assert np.isclose(p["pair_flux"], 2.0)
    assert np.isclose(p["pair"], 4.0)


def test_predictions_skip_unusable_siblings():
    p = xr.predictions({"a": (-1.0, 50.0), "b": (8.0, np.nan), "c": (4.0, 2.0)}, 1.0, 1.0)
    assert set(p["from"]) == {"c"}
    assert np.isclose(p["pair"], 2.0)
    empty = xr.predictions({"a": (-1.0, 50.0)}, 1.0, 1.0)
    assert empty["from"] == {} and np.isnan(empty["pair"])
    assert np.isnan(xr.predictions({"c": (4.0, 2.0)}, 1.0, np.nan)["pair"])


def _pred_table():
    return Table(
        rows=[
            ("33.1", 3.584699, -30.403125, -44.5, "observed", "33.1a", 0.1, "", 1.0),
            ("33.1", 3.600557, -30.395284, 2.94, "no_flux_ref", "", 56.8, "no_reference", 1.0),
            ("33.1", 3.584980, -30.403520, -32.1, "no_flux_ref", "", 1.7, "no_reference", 1.0),
        ],
        names=[
            "system",
            "ra",
            "dec",
            "magnification",
            "image_class",
            "matched_image",
            "sep_image_arcsec",
            "forced_class",
            "search_arcsec",
        ],
    )


def test_targets_are_the_both_model_rows_only():
    pred = _pred_table()
    xm = Table(
        rows=[
            ("bergamini23", "33.1", 3.600557, -30.395284, True),
            ("bergamini23", "33.1", 3.584980, -30.403520, False),
        ],
        names=["model", "system", "ra", "dec", "both_loose"],
    )
    t = xr.targets(pred, xm, "bergamini23", "33.1")
    assert len(t) == 1 and np.isclose(t["magnification"][0], 2.94)
    assert len(xr.targets(pred, xm, "uncover2", "33.1")) == 0


def test_sibling_mu_uses_observed_rows():
    assert xr.sibling_mu(_pred_table(), "33.1") == {"33.1a": 44.5}


def test_knot_match_skips_own_system():
    arcs = Table(
        rows=[("3.1c", "3.1", 3.5766219, -30.4017919), ("3.2a", "3.2", 3.5892122, -30.3938244)],
        names=["image_id", "system", "ra", "dec"],
    )
    kid, sep = xr.knot_match(arcs, "3.2", 3.57667, -30.401646)
    assert kid == "3.1c" and 0.4 < sep < 0.7
    assert xr.knot_match(arcs[1:], "3.2", 3.57667, -30.401646)[0] == ""

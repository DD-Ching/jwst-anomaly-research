"""Frozen decision rule of results/dimming/nexus_o014/PLAN.md (summarise.explanation)."""

import importlib.util
from pathlib import Path

import numpy as np

_P = Path(__file__).resolve().parents[1] / "results/dimming/nexus_o014/summarise.py"
_spec = importlib.util.spec_from_file_location("nexus_o014_summarise", _P)
sm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sm)

NE = 4


def _row(flux=(1.0, 1.0, 1.0, 1.0), vanish=False, survives=False, cut=False):
    r = {f"F200W_e{k}_fflux": f for k, f in enumerate(flux)}
    r.update(
        forced_vanish=vanish,
        survives_artefact_tests=survives,
        cut_edge=cut,
        cut_low_weight=False,
        cut_spike=False,
    )
    return r


def test_survivor_wins():
    assert sm.explanation(_row(vanish=True, survives=True), NE) == "survivor"


def test_present_in_forced_photometry_is_catalogue_depth():
    assert sm.explanation(_row(), NE) == "catalogue_depth_or_deblend"


def test_uncovered_other_epochs_is_footprint():
    nan = np.nan
    assert sm.explanation(_row(flux=(nan, nan, 1.0, nan)), NE) == "footprint_or_image_quality"
    assert sm.explanation(_row(vanish=True, cut=True), NE) == "footprint_or_image_quality"


def test_vanish_removed_by_a_later_test():
    assert sm.explanation(_row(vanish=True), NE) == "other_forced_test"

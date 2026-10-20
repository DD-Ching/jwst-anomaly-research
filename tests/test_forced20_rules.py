"""Frozen decision rule of results/dimming/forced20/PLAN.md (summarise.explanation, per-field
detection epoch)."""

import importlib.util
from pathlib import Path

import numpy as np

_P = Path(__file__).resolve().parents[1] / "results/dimming/forced20/summarise.py"
_spec = importlib.util.spec_from_file_location("forced20_summarise", _P)
sm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sm)

NE = 3  # Abell 2744 epochs; detection epoch 1


def _row(flux=(1.0, 1.0, 1.0), vanish=False, survives=False, cut=False):
    r = {f"F200W_e{k}_fflux": f for k, f in enumerate(flux)}
    r.update(
        forced_vanish=vanish,
        survives_artefact_tests=survives,
        cut_edge=cut,
        cut_low_weight=False,
        cut_spike=False,
    )
    return r


def test_detection_epochs():
    assert sm.DET_EPOCH == {"macs0416": 2, "abell2744": 1}


def test_survivor_wins():
    assert sm.explanation(_row(vanish=True, survives=True), NE, 1) == "survivor"


def test_present_in_forced_photometry_is_catalogue_depth():
    assert sm.explanation(_row(), NE, 1) == "catalogue_depth_or_deblend"


def test_only_the_detection_epoch_finite_is_footprint():
    nan = np.nan
    assert sm.explanation(_row(flux=(nan, 1.0, nan)), NE, 1) == "footprint_or_image_quality"
    # the detection epoch's own flux does not count as "other epoch covered"
    assert sm.explanation(_row(flux=(nan, nan, 1.0)), NE, 1) == "catalogue_depth_or_deblend"
    assert sm.explanation(_row(flux=(nan, nan, 1.0)), NE, 2) == "footprint_or_image_quality"
    assert sm.explanation(_row(vanish=True, cut=True), NE, 1) == "footprint_or_image_quality"


def test_vanish_removed_by_a_later_test():
    assert sm.explanation(_row(vanish=True), NE, 1) == "other_forced_test"

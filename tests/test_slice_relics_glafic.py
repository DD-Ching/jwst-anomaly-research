"""Offline checks of the SLICE x RELICS GLAFIC rule functions (PLAN.md rules 7-9)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import slice_relics_glafic as sg  # noqa: E402


def test_amplitude_class_bands():
    assert sg.amplitude_class(1.0) == "near_1"
    assert sg.amplitude_class(0.85) == "near_1"
    assert sg.amplitude_class(1.15) == "near_1"
    assert sg.amplitude_class(0.66) == "lenstool_like"
    assert sg.amplitude_class(0.75) == "lenstool_like"
    assert sg.amplitude_class(0.80) == "intermediate"
    assert sg.amplitude_class(1.3) == "intermediate"


def test_gate_p_needs_both_ratios():
    ok = {"kappa_map": {"median_ratio": 1.0}, "magnification_z2": {"median_ratio": 1.04}}
    bad = {"kappa_map": {"median_ratio": 1.0}, "magnification_z2": {"median_ratio": 1.06}}
    assert sg.gate_p(ok)
    assert not sg.gate_p(bad)


def test_two_model_statement():
    assert sg.two_model_11(False, 5.0, True) == "glafic_not_usable"
    assert sg.two_model_11(True, 3.5, True) == "missed_by_both"
    assert sg.two_model_11(True, 2.9, True) == "lenstool_only"
    assert sg.two_model_11(True, 3.5, False) == "glafic_only"
    assert sg.two_model_11(True, 2.9, False) == "neither"


def test_models_cover_srh_clusters():
    import slice_relics_heldout as srh

    assert set(sg.MODELS) <= set(srh.CLUSTERS)
    for model in sg.MODELS.values():
        assert model in sg.lc.MODELS

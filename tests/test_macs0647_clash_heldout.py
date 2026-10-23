"""Offline checks of the MACS0647 CLASH held-out rules (PLAN.md rules 2-4) and inputs."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import macs0647_clash_heldout as ho  # noqa: E402


def test_load_images_drops_uncertain_and_singletons():
    t = ho.load_images()
    ids = set(t["image_id"])
    assert not {"8.3", "16.3", "16.4"} & ids
    assert len(t) == 69  # 72 minus the three uncertain rows; every family keeps >= 2 images
    fam = np.asarray(t["system"], str)
    assert all(np.sum(fam == f) >= 2 for f in fam)
    held = set(fam[np.asarray(t["held_out"], bool)])
    assert held == {str(k) for k in range(13, 24)}


def test_gate_k0():
    assert ho.gate_k0(np.array([1.0, 1.5, 2.0, 2.5]))["passes"]  # rms 1.84, all matched
    assert not ho.gate_k0(np.array([2.5, 2.5]))["passes"]  # rms 2.5 > 2.0
    # 1 of 5 unmatched (NaN): fraction 0.8 passes; 2 of 5 fails
    assert ho.gate_k0(np.array([0.5, 0.5, 0.5, 0.5, np.nan]))["passes"]
    assert not ho.gate_k0(np.array([0.5, 0.5, 0.5, 4.0, np.nan]))["passes"]
    assert not ho.gate_k0(np.array([]))["passes"]


def test_frame_check_is_per_axis_and():
    assert ho.frame_check([0.4, -0.4, 0.1], [0.0, 0.2, -0.1])["passes"]
    assert not ho.frame_check([0.6, 0.7, 0.6], [0.0, 0.0, 0.0])["passes"]
    assert not ho.frame_check([0.0, 0.0, 0.0], [-0.6, -0.7, 0.0])["passes"]
    assert not ho.frame_check([np.nan], [np.nan])["passes"]


def test_readings_need_every_usable_model_to_miss():
    d = {"a": np.array([1.0, 4.0, 4.0, np.nan]), "b": np.array([4.0, 1.0, 4.0, 5.0])}
    r = ho.readings(d, ["a", "b"])
    assert list(r) == ["no_flag", "no_flag", "lensing_residual_input", "lensing_residual_input"]
    assert list(ho.readings(d, ["a"])) == [
        "no_flag",
        "lensing_residual_input",
        "lensing_residual_input",
        "lensing_residual_input",
    ]
    assert list(ho.readings(d, [])) == ["no_flag"] * 4

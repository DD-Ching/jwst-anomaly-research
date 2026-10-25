"""Offline checks of the MACS2129 CLASH held-out inputs and rule 5 label (PLAN.md)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import macs2129_clash_heldout as ho  # noqa: E402


def test_load_images_matches_slice_table5_counts():
    t = ho.load_images()
    fam = np.asarray(t["system"], str)
    assert len(t) == 53 and len(set(fam)) == 17  # SLICE Table 5: 17 sources, 53 images
    assert "15" not in set(fam)
    assert len(set(fam[np.asarray(t["z_spec"], bool)])) == 13  # 13 spectroscopic sources
    assert set(fam[np.asarray(t["held_out"], bool)]) == {"16", "17", "18"}
    assert not np.any(np.asarray(t["z_spec"], bool) & np.asarray(t["held_out"], bool))


def test_label_b_renames_only_the_flag():
    r = np.array(["no_flag", "lensing_residual_input"])
    assert list(ho.label_b(r)) == ["no_flag", "redshift_residual_input"]


def test_readings_reused_need_every_usable_model_to_miss():
    d = {"a": np.array([4.0, 4.0]), "b": np.array([1.0, 4.0])}
    assert list(ho.label_b(ho.ho.readings(d, ["a", "b"]))) == [
        "no_flag",
        "redshift_residual_input",
    ]


def test_reading_b_uses_only_spec_rows_of_the_zspec_offsets():
    spec = np.array([True, False, True])
    d = {"a": np.array([4.0, 9.0, 1.0]), "b": np.array([5.0, 9.0, 4.0])}
    assert list(ho.reading_b(d, spec, ["a", "b"])) == ["redshift_residual_input", "no_flag"]
    assert list(ho.reading_b(d, spec, [])) == ["no_flag", "no_flag"]

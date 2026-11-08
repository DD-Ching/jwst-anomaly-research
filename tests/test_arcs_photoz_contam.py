"""Frozen rule functions of scripts/arcs_photoz_contam.py (written before the run)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import arcs_photoz_contam as pc  # noqa: E402


def test_spec_class():
    assert pc.spec_class(None, 0.43) is None
    assert pc.spec_class(0.44, 0.43) == "member_or_foreground"
    assert pc.spec_class(0.20, 0.43) == "member_or_foreground"
    assert pc.spec_class(0.70, 0.43) == "background"
    assert pc.member_cut(0.43) == 0.43 + 0.03 * 1.43


def test_label_order():
    assert pc.label("member_or_foreground", "off_sequence", False) == "contaminated"
    assert pc.label("background", "red_sequence_like", True) == "background_spec"
    assert pc.label(None, "red_sequence_like", None) == "contaminated"
    assert pc.label(None, "red_sequence_like", True) == "contaminated"
    assert pc.label(None, "red_sequence_like", False) == "undetermined"
    assert pc.label(None, "off_sequence", False) == "background_photometric"
    assert pc.label(None, "off_sequence", None) == "background_photometric"
    assert pc.label(None, "off_sequence", True) == "undetermined"
    assert pc.label(None, None, None) == "undetermined"


def test_outcome():
    assert pc.outcome(["contaminated", "contaminated", "undetermined"]).startswith("open part")
    assert pc.outcome(["background_spec", "background_photometric", "undetermined"]).startswith(
        "photo-z contamination not supported"
    )
    assert (
        pc.outcome(["contaminated", "background_spec", "background_spec"]) == "mixed / undetermined"
    )
    assert pc.outcome(["undetermined"] * 3) == "mixed / undetermined"


def test_gate():
    assert pc.gate(0.8, 0.3)
    assert not pc.gate(0.79, 0.1)
    assert not pc.gate(0.9, 0.31)
    assert not pc.gate(float("nan"), 0.1)


def test_red_sequence_fit_recovers_line():
    rng = np.random.default_rng(1)
    mag = rng.uniform(18, 23, 300)
    col = 1.2 - 0.03 * (mag - 20) + rng.normal(0, 0.05, 300)
    col[:20] += 1.0  # outliers
    rs = pc.fit_red_sequence(mag, col)
    assert abs(rs["slope"] + 0.03) < 0.01
    assert abs(rs["sigma"] - 0.05) < 0.01
    assert abs(pc.rs_offset(rs, 20.0, 1.2 + 0.1) - 2.0) < 0.5


def test_ned_spec_z_skips_photometric():
    rows = [
        {"z": 0.5, "zflag": "PHOT", "sep_arcsec": 0.1},
        {"z": 0.7, "zflag": "", "sep_arcsec": 0.5},
    ]
    assert pc.ned_spec_z(rows) == 0.7
    assert pc.ned_spec_z(rows[:1]) is None


def test_combine_spec():
    assert pc.combine_spec(None, None, 0.43) == (None, False)
    assert pc.combine_spec(0.44, None, 0.43) == ("member_or_foreground", False)
    assert pc.combine_spec(None, 0.9, 0.43) == ("background", False)
    assert pc.combine_spec(1.2, 0.44, 0.43) == (None, True)

"""Unit tests of the Abell 2744 UNCOVER image-list copy (PLAN.md Inputs)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import lensmodel

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_uncover as au  # noqa: E402
from a2744_uncoverref import image_list, lenstool_id  # noqa: E402


def test_lenstool_id_keeps_the_knot_level_source():
    for uid in ("1.11", "1.13", "78.21", "83.1", "83.3", "12.2"):
        assert lensmodel.image_family(lenstool_id(uid)) == au.uncover_system(uid)
    assert lenstool_id("78.22") == "78.2b"
    assert lenstool_id("83.3") == "83c"


def test_image_list_flags_spec_and_drops_none():
    raw = Table(
        {
            "ID": [1.11, 1.12, 83.1, 83.2, 90.1],
            "RA": [3.5, 3.6, 3.7, 3.8, 3.9],
            "Dec": [-30.4] * 5,
            "z": [2.6, 2.6, 4.0, 4.0, 0.0],
            "Redshift_type": ["spec", "spec", "phot", "phot", "None"],
        }
    )
    t = image_list(raw)
    assert list(t["image_id"]) == ["1.1a", "1.1b", "83a", "83b"]
    assert list(t["z_spec"]) == [True, True, False, False]
    assert np.all(np.isfinite(t["ra"]))


def test_ambiguous_ids_fail_loudly():
    with pytest.raises(ValueError):
        lenstool_id("34.1" + "0")
    raw = Table(
        {
            "ID": [1.11, 1.11],
            "RA": [3.5, 3.6],
            "Dec": [-30.4] * 2,
            "z": [2.6, 2.6],
            "Redshift_type": ["spec", "spec"],
        }
    )
    with pytest.raises(ValueError):
        image_list(raw)

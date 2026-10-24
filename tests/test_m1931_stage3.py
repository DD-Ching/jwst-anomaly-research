"""Frozen stage 3 footprint rule of scripts/m1931_stage3.py (offline)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
pytest.importorskip("photutils")
import m1931_stage3 as s3  # noqa: E402


def test_footprint_class_outside_masked_and_kept():
    d = np.ones((10, 20))
    d[5, 7] = np.nan
    assert s3.footprint_class("absent", 3.2, 4.1, d) == "absent"
    assert s3.footprint_class("absent", 7.4, 4.6, d) == s3.OFF_IMAGE  # nearest pixel NaN
    assert s3.footprint_class("recovered", -0.6, 2.0, d) == s3.OFF_IMAGE
    assert s3.footprint_class("recovered", 19.4, 9.4, d) == "recovered"
    assert s3.footprint_class("recovered", 19.6, 9.0, d) == s3.OFF_IMAGE
    assert s3.footprint_class("absent", np.nan, 1.0, d) == s3.OFF_IMAGE


def test_off_image_is_not_testable():
    import m0417_arcref as ar

    assert s3.OFF_IMAGE not in ar.TESTABLE
    g = ar.gate_c2_shared(
        ["recovered"] * 8 + [s3.OFF_IMAGE] * 5,
        [False] * 13,
        [str(i) for i in range(13)],
        [(-1, -1)] * 13,
    )
    assert g["passes"] is False  # no arc-referenced control
    assert g["all"]["n_testable"] == 8

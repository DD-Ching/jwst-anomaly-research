"""Tests of the MACS0416 CANUCS-reference script: the image-list copy (PLAN.md ASSUMPTION)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m0416_canucsref as m  # noqa: E402


def test_constants_name_the_planned_field_and_model():
    assert (m.FIELD, m.MODEL, m.OBS, m.TILE) == ("macs0416", "macs0416-canucs", "o004", "t002")
    assert m.OUT == Path("results/lens_consistency/macs0416_canucsref")


@pytest.mark.network
def test_written_images_are_the_303_model_images_with_fixed_z(tmp_path):
    tab = Table.read(m.write_images(tmp_path), format="ascii.ecsv")
    assert len(tab) == 303 and bool(np.all(tab["z_spec"]))
    assert np.all(np.isfinite(np.asarray(tab["z"], float)))

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import astrodeep_zout as az  # noqa: E402


def test_to_zout_joins_and_converts():
    pos = Table({"ID": [1, 2, 3], "RAJ2000": [10.0, 10.1, 10.2], "DEJ2000": [-1.0, -1.1, -1.2]})
    z = Table(
        {
            "ID": [3, 1, 2],
            "zbest": [2.0, 1.0, 0.05],
            "e_zbest": [0.0, 0.1, 0.1],
            "relFlag": [1.0, 1.0, 0.0],
        }
    )
    out = az.to_zout(pos, z)
    by_id = {int(r["id"]): r for r in out}
    assert set(by_id) == {1, 3}  # relFlag 0 dropped
    assert np.isclose(by_id[1]["ra"], 10.0)
    assert np.isclose(by_id[1]["z_phot"], 1.0)
    assert np.isclose(by_id[1]["z160"], 1.0 - az.SIQR_TO_SIGMA * 0.1)
    assert np.isclose(by_id[1]["z840"], 1.0 + az.SIQR_TO_SIGMA * 0.1)
    assert by_id[3]["z160"] == by_id[3]["z840"] == 2.0  # spectroscopic, SIQR 0
    assert out.meta["provenance"] == "derived"


def test_floor_at_zero():
    pos = Table({"ID": [1], "RAJ2000": [0.0], "DEJ2000": [0.0]})
    z = Table({"ID": [1], "zbest": [0.1], "e_zbest": [1.0], "relFlag": [1.0]})
    assert az.to_zout(pos, z)["z160"][0] == az.Z_FLOOR


@pytest.mark.network
def test_fetch_m0717_live():
    pos, z = az.fetch("m0717")
    assert len(pos) == len(z) == 4082
    assert len(az.to_zout(pos, z)) == 2866

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
            "f_zbest": [0.0, 0.0, 0.0],
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


def test_z160_floor():
    pos = Table({"ID": [1], "RAJ2000": [0.0], "DEJ2000": [0.0]})
    z = Table({"ID": [1], "zbest": [0.1], "e_zbest": [1.0], "relFlag": [1.0], "f_zbest": [0.0]})
    assert az.to_zout(pos, z)["z160"][0] == az.Z_FLOOR


def test_spectroscopic_rows_kept_whatever_relflag():
    pos = Table({"ID": [1, 2], "RAJ2000": [0.0, 0.1], "DEJ2000": [0.0, 0.1]})
    z = Table(
        {
            "ID": [1, 2],
            "zbest": [0.55, 1.2],
            "e_zbest": [np.nan, 0.2],
            "relFlag": [0.0, 0.0],
            "f_zbest": [1.0, 0.0],
        }
    )
    out = az.to_zout(pos, z)
    assert list(out["id"]) == [1]
    assert out["z160"][0] == out["z840"][0] == 0.55


@pytest.mark.network
def test_fetch_m0717_live():
    pos, z = az.fetch("m0717")
    assert len(pos) == len(z) == 4082
    assert (
        len(az.to_zout(pos, z)) == 2928
    )  # relFlag 1 (2,866) plus 62 spectroscopic rows with relFlag 0

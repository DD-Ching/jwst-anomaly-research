"""Frozen decision rule of scripts/rxc0018_helios.py (offline)."""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location(
    "rxc0018_helios", ROOT / "scripts" / "rxc0018_helios.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def img(d, ratio=0.6, sep=2.0):
    return {"delay_yr": d, "mu_ratio": ratio, "sep_centre_arcsec": sep}


@pytest.mark.parametrize(
    ("path_ok", "sep", "others", "expected"),
    [
        (False, 0.0, [img(3.0)], "stopped: path"),
        (True, 0.31, [img(3.0)], "stopped: no SN image"),
        (True, float("nan"), [img(3.0)], "stopped: no SN image"),
        (True, 0.1, [img(3.0)], "agrees"),
        (True, 0.1, [img(-45.0), img(3.0), img(20.0, sep=30.0)], "agrees"),
        (True, 0.1, [img(1.0, sep=5.0)], "agrees"),
        (True, 0.1, [img(10.0)], "agrees"),
        (True, 0.1, [img(10.01)], "differs"),
        (True, 0.1, [img(0.99)], "differs"),
        (True, 0.1, [img(3.0, sep=5.01)], "differs"),
        (True, 0.1, [img(0.5, sep=20.0), img(3.0)], "differs"),
        (True, 0.1, [img(-45.0)], "differs"),
        (True, 0.1, [], "undetermined"),
        (True, 0.1, [img(-0.5), img(3.0)], "conflicts"),
        (True, 0.1, [img(-1.0, ratio=-0.5), img(3.0)], "conflicts"),
        (True, 0.1, [img(-1.01), img(3.0)], "agrees"),
        (True, 0.1, [img(-0.5, ratio=0.49), img(3.0)], "agrees"),
    ],
)
def test_decide(path_ok, sep, others, expected):
    assert mod.decide(path_ok, sep, others) == expected


def test_published_centre_near_inputs():
    ra, dec = mod.published_centre()
    assert 4.63 < ra < 4.64
    assert 16.43 < dec < 16.44
    assert mod.HELIOS_RADEC[1] > dec  # Helios lies north of the published next images

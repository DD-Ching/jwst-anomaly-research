import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import arcs_rxj2129_dja as m  # noqa: E402


def test_dja_pa_north_up_east_left():
    # +x is west (PA 90 axis), +y north (PA 0); 45 deg counter-clockwise from +x is north-west
    assert np.allclose(m.dja_pa(np.radians([0.0, 90.0, 45.0, -45.0])), [90.0, 0.0, 135.0, 45.0])


def test_convention_gate():
    assert m.convention_passes(9.9, 5)
    assert not m.convention_passes(10.1, 5)
    assert not m.convention_passes(float("nan"), 0)
    assert not m.convention_passes(1.0, 0)


def test_to_pipeline_columns_and_scale():
    dja = Table(
        {
            "id": [1, 2],
            "ra": [322.41, 322.42],
            "dec": [0.08, 0.09],
            "flux_iso": [100.0, -1.0],
            "fluxerr_iso": [1.0, 1.0],
            "a_image": [4.0, 2.0],
            "b_image": [1.0, 2.0],
            "theta_image": [0.0, np.pi / 2],
            "area_iso": [30, 10],
            "mag_auto": [24.0, np.nan],
        }
    )
    t = m.to_pipeline(dja, 0.031)
    assert np.allclose(t["ellipticity"], [0.75, 0.0])
    assert np.allclose(t["semimajor_sigma"], [4.0 * 0.04 / 0.031, 2.0 * 0.04 / 0.031])
    assert np.allclose(t["sky_orientation"], [90.0, 0.0])
    assert t["isophotal_abmag"][0] == 24.0 and np.isnan(t["isophotal_abmag"][1])
    assert t["is_extended"].all()

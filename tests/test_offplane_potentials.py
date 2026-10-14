from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

from jwst_anomaly import lensmodel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location(
    "offplane_potentials", ROOT / "scripts" / "offplane_potentials.py"
)
op = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(op)


def _model():
    comps = [
        lensmodel.DPIE("halo", 0.0, 0.0, 0.2, 0.0, 10.0, 1000.0, 1000.0, 0.4),
        lensmodel.DPIE("fg", 10.0, 0.0, 0.0, 0.0, 0.1, 5.0, 120.0, 0.4),
        lensmodel.DPIE("blend", 0.0, 10.0, 0.0, 0.0, 0.1, 5.0, 80.0, 0.4),
        lensmodel.DPIE("member", -10.0, 0.0, 0.0, 0.0, 0.1, 5.0, 80.0, 0.4),
    ]
    return lensmodel.LensModel(comps, 64.0, -24.0, FlatLambdaCDM(H0=70, Om0=0.3))


def test_match_potentials_flags_offplane_and_photoz():
    model = _model()
    ra, dec = model.to_sky([10.0, 0.0, -10.0, 30.0], [0.0, 10.0, 0.0, 30.0])
    cat = Table(
        {
            "SOURCE": [1, 2, 3, 4],
            "RA": ra,
            "DEC": dec + np.array([0.1, 0.0, 0.0, 0.0]) / 3600.0,
            "Z_SPEC": [0.11, 2.5, 0.41, 0.2],
            "Z_ML": [0.07, 0.38, 0.40, 0.2],
            "REF_SPECZ": ["a", "b", "c", "d"],
        }
    )
    out = op.match_potentials(model, cat, 0.4, 0.5, 0.02)
    assert list(out["name"]) == ["halo", "fg", "blend", "member"]
    assert list(out["cat_id"]) == [-1, 1, 2, 3]  # the halo has no source within 0.5"
    assert np.isnan(out["z_spec"][0]) and not out["offplane"][0]
    assert list(out["offplane"][1:]) == [True, True, False]
    assert list(out["photoz_agrees"][1:]) == [True, False, True]
    assert abs(out["sep_arcsec"][1] - 0.1) < 1e-3

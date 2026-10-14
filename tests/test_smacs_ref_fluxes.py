"""Tests for scripts/smacs_ref_fluxes.py (synthetic stamps and a stubbed lens model; no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("smacs_ref_fluxes", _DIR / "smacs_ref_fluxes.py")
srf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(srf)

PIX = 0.03  # arcsec


def _stamper(sources):
    """``stamp(ra, dec)``: unit-error noise-free images with Gaussian blobs at (ra, dec) keys."""
    ax = np.arange(-4.5, 4.5 + 1e-9, PIX)
    xx, yy = np.meshgrid(ax, ax)

    def stamp(ra, dec):
        img = np.zeros_like(xx)
        for (sra, sdec), amp in sources.items():
            dx, dy = (sra - ra) * 3600, (sdec - dec) * 3600
            img += amp * np.exp(-((xx - dx) ** 2 + (yy - dy) ** 2) / (2 * 0.05**2))
        return img, np.full_like(xx, 0.01), xx, yy

    return stamp


def test_measure_predicts_from_the_sibling_and_needs_a_mu_factor(monkeypatch):
    # one catalogued image (|mu| 100); a predicted counter-image (|mu| 5) 20" away on empty sky
    cat = Table({"image_id": ["9.1"], "system": ["9"], "ra": [0.0], "dec": [0.0]})
    stamps = {"f200w": _stamper({(0.0, 0.0): 1.0})}
    bt = Table({"image_id": ["9.1"], "magnification": [100.0], "beta_x": [0.0], "beta_y": [0.0]})
    imgs = Table(
        {
            "ra": [0.0, 20 / 3600],
            "dec": [0.0, 0.0],
            "x": [0.0, 20.0],
            "y": [0.0, 0.0],
            "magnification": [100.0, -5.0],
        }
    )
    monkeypatch.setattr(srf.lensmodel, "backtrace_images", lambda m, c, z: bt)
    monkeypatch.setattr(srf.lensmodel, "find_images", lambda m, g, bx, by, z: imgs)
    sib = srf.sibling_flux(cat, stamps)
    f = sib["f200w"]["9.1"]["flux"]
    assert f > 0 and sib["f200w"]["9.1"]["conc"] > 0.9  # compact blob
    out = srf.measure(None, None, cat, 2.0, sib, stamps)
    near, far = out["images"]
    assert near["d_cat"] == 0.0 and "f200w" not in near  # the catalogued image itself
    _, e0 = srf.lc.aperture_snr(*stamps["f200w"](20 / 3600, 0.0), 0.0, 0.0)
    pred = f / 100 * 5 / e0
    r = far["f200w"]
    assert r["pred_snr_from"]["9.1"] == pytest.approx(pred, abs=0.06)
    assert r["need_mu_factor"] == pytest.approx(10 / pred, rel=0.02)
    assert sib["f200w"]["9.1"]["gates_ok"] is True
    assert r["pred_snr_gated"] == r["pred_snr_from"]["9.1"]
    assert r["best_snr"] == pytest.approx(0.0, abs=1e-6)
    assert out["mu_catalogued"] == {"9.1": 100.0}


def test_sibling_gates_reject_a_neighbour_peak():
    # the catalogued spot is faint; a bright neighbour 0.3" away dominates the 0.4" peak search
    cat = Table({"image_id": ["9.1"], "system": ["9"], "ra": [0.0], "dec": [0.0]})
    stamps = {"f200w": _stamper({(0.0, 0.0): 0.1, (0.0, 0.3 / 3600): 2.0})}
    s = srf.sibling_flux(cat, stamps)["f200w"]["9.1"]
    assert s["peak_over_position"] > 2 and s["gates_ok"] is False

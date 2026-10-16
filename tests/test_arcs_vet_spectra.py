"""Tests for scripts/arcs_vet_spectra.py pure functions (synthetic; no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_vet_spectra", _DIR / "arcs_vet_spectra.py")
avs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(avs)


def test_expected_fluxes_scale_with_ew():
    e = avs.expected_fluxes(25.0, 8000.0, 1.0)
    assert np.isclose(e["EW_rest_100"] / e["EW_rest_30"], 100 / 30, rtol=1e-2)
    # AB 25 at 8000 Å: f_lambda ≈ 1.70e-19; EW_obs 60 Å -> ≈ 1.0e-17
    assert np.isclose(e["EW_rest_30"], 1.70e-19 * 60, rtol=0.02)


def test_window_result_detects_and_respects_sky_flag():
    lam = np.arange(4750.0, 9350.0, 1.25)
    snr = np.zeros_like(lam)
    k = np.argmin(np.abs(lam - 8800))
    snr[k] = 7.0
    m = {
        "lam": lam,
        "snr": snr,
        "sky_snr": np.zeros_like(lam),
        "sigma_flux": np.full_like(lam, 1e-19),
    }
    r = avs.window_result(m, 8700, 8900)
    assert r["detected"] and abs(r["max_snr_lambda_A"] - 8800) < 1.5
    assert np.isclose(r["limit_5sigma_cgs_median"], 5e-19)
    m["sky_snr"] = np.where(np.abs(lam - 8800) < 5, np.inf, 0.0)
    assert not avs.window_result(m, 8700, 8900)["detected"]
    assert avs.window_result(m, 9400, 9500)["in_range"] is False

"""Rules of scripts/h0pe_mu_check.py (frozen plan results/lens_consistency/h0pe_mu/PLAN.md)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import h0pe_mu_check as h  # noqa: E402


def test_ln_sigma_symmetric_case():
    assert math.isclose(h.ln_sigma((1.0, math.e - 1, 1 - 1 / math.e)), 1.0)


def test_offset_zero_for_photometry_itself():
    assert math.isclose(h.offset_mag(h.PHOT), 0.0, abs_tol=1e-12)
    assert all(math.isclose(h.image_z(h.PHOT, i), 0.0, abs_tol=1e-12) for i in h.IMAGES)


def test_mass_sheet_recovers_injected_lambda():
    # mass sheet on the photometry: model mu = lambda^2 mu_phot, model dt = dt_phot / lambda
    lam = 1.2
    model = {i: (h.PHOT[i][0] * lam**2, 0.1, 0.1) for i in h.IMAGES}
    model |= {k: (h.PHOT[k][0] / lam, 1.0, 1.0) for k in ("dt_ab", "dt_cb")}
    r = h.h0_mst(model)
    assert math.isclose(r["lambda_mu"], lam)
    assert math.isclose(r["lambda_t"], lam)
    assert math.isclose(r["h0_mst"], h.H0_FID)
    assert r["in_range"]


def test_ln_sigma_abs_swaps_errors_for_negative_delays():
    assert math.isclose(h.ln_sigma_abs((-100.0, 5.0, 10.0)), h.ln_sigma((100.0, 10.0, 5.0)))
    assert math.isclose(h.ln_sigma_abs((100.0, 5.0, 10.0)), h.ln_sigma((100.0, 5.0, 10.0)))


def test_ratio_z_common_term_cancels_and_sign():
    # a model with the photometric ratios has ratio z 0 whatever its absolute scale
    model = {i: (h.PHOT[i][0] * 2.0, 0.1, 0.1) for i in h.IMAGES}
    assert math.isclose(h.ratio_z(model, "b"), 0.0, abs_tol=1e-12)
    # removing the shared calibration term shrinks the error, so |z| grows
    model["b"] = (h.PHOT["b"][0] * 3.0, 0.1, 0.1)
    z_rm = h.ratio_z(model, "b")
    z_q = h.ratio_z(model, "b", remove_common=False)
    assert z_rm < z_q < 0  # model over-predicts b / a, observed lower
    assert abs(z_rm) > abs(z_q)


def test_published_run_outcomes_pinned():
    res, posthoc = h.run()
    assert res["outcome"] == "N"
    assert res["rule1"] == "method-dependent"
    assert res["n_models_offset_positive"] == 7
    assert res["n_models_h0_mst_in_range"] == 0
    assert len(res["rows"]) == 7
    assert math.isclose(res["pooled_offset_mag"], 0.633, abs_tol=5e-4)
    assert math.isclose(posthoc["pooled_offset_mag_without_model6"], 0.453, abs_tol=5e-4)

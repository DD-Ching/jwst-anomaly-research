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


def test_run_outcome_and_rule1_are_strings():
    res = h.run()
    assert res["outcome"] in {"N", "F"}
    assert res["rule1"] in {"reproduces", "method-dependent"}
    assert len(res["rows"]) == 7

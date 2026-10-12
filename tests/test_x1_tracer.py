"""Offline tests for scripts/x1_tracer.py (X1 S1 checks, docs/state/units/X1.md)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import wh_farside as wf  # noqa: E402
import x1_tracer as xt  # noqa: E402

# roundoff warnings near the photon sphere; accuracy is checked against the ODE route instead
pytestmark = pytest.mark.filterwarnings("ignore::scipy.integrate.IntegrationWarning")


@pytest.mark.parametrize("b", [6.0, 8.0, 20.0])
def test_schwarzschild_deflection_matches_darwin(b):
    assert abs(xt.deflection_schw_numeric(b) / xt.deflection_schw_exact(b) - 1) < 1e-6


def test_weak_field_deflection_limit():
    # 4/b + 15 pi / (4 b^2) to second order
    b = 2000.0
    assert abs(xt.deflection_schw_exact(b) - (4 / b + 15 * math.pi / (4 * b * b))) < 1e-8


@pytest.mark.parametrize("a", [2.1, 5.0])
def test_b_crit_two_routes(a):
    assert abs(xt.b_crit_ode(a, tol=1e-6) / wf.b_crit(a) - 1) < 1e-4
    below, above = 0.999 * wf.b_crit(a), 1.001 * wf.b_crit(a)
    assert xt.turning_u(below, a) is None and xt.turning_u(above, a) is not None


@pytest.mark.parametrize("a,frac", [(3.0, 0.5), (7.0, 0.99)])
def test_far_side_quadrature_matches_ode(a, frac):
    b = frac * wf.b_crit(a)
    assert abs(xt.phi_to_emitter(b, a, 9.0, "far")[0] / xt.phi_far_ode(b, a, 9.0) - 1) < 1e-6


def test_near_side_branches_meet_at_fold():
    bf = xt._fold_b(0.0, 9.0) * (1 - 1e-10)
    pe, pa = xt.phi_to_emitter(bf, 0.0, 9.0, "near")
    assert abs(pe - pa) < 1e-4
    assert xt.phi_to_emitter(1.01 * xt._fold_b(0.0, 9.0), 0.0, 9.0, "near") == []


def test_face_on_primary_ring_near_r_plus_one():
    tr = xt.spot_track(0.0, 9.0, 179.5, "near", n_phase=12, primary_only=True)
    rad = np.hypot(tr[:, 1], tr[:, 2])
    assert abs(rad.mean() - 10.0) < 0.5


@pytest.mark.parametrize("a", [2.1, 7.0])
def test_far_side_images_inside_b_crit(a):
    tr = xt.spot_track(a, 9.0, 155.0, "far", n_phase=12)
    assert np.nanmax(np.hypot(tr[:, 1], tr[:, 2])) < wf.b_crit(a)
    assert np.all(tr[:, 3] > 0)

"""Offline tests for scripts/x1_s4b.py (independent check of the tracer's image weights)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import wh_farside as wf  # noqa: E402
import x1_s4b  # noqa: E402
import x1_tracer as xt  # noqa: E402


def _geometry(phase: float, incl: float = 155.0):
    inc = math.radians(incl)
    o = np.array([math.sin(inc), 0.0, math.cos(inc)])
    n = np.array([math.cos(phase), math.sin(phase), 0.0])
    proj = n - (n @ o) * o
    return o, n, proj / np.linalg.norm(proj), math.acos(float(n @ o))


def test_hamilton_hit_matches_tracer_angle_and_is_null():
    b = 0.9 * wf.b_crit(2.1)
    phi, ca, sa = x1_s4b.hit_hamilton(b, 2.1, 9.0, "far")
    assert abs(phi / xt.phi_to_emitter(b, 2.1, 9.0, "far")[0] - 1.0) < 1e-9
    assert ca > 0.0  # far side: the backward ray moves away from the throat
    assert abs(ca * ca + sa * sa - 1.0) < 1e-9


def test_weak_field_magnification_is_one(monkeypatch):
    # emitter far out (R = 2000 M) at gamma = 90 deg: no lensing to O(M / R), so screen area per
    # cross-section is 1
    monkeypatch.setattr(x1_s4b, "R0", 1e4)
    re = 2000.0
    o, _, e, _ = _geometry(math.pi / 2, incl=90.0)
    b = 0.999 * re
    phi, _, _ = x1_s4b.hit_hamilton(b, 0.0, re, "near", 0)
    n = math.cos(phi) * o + math.sin(phi) * e  # the emitter is where this ray lands
    m, miss, _ = x1_s4b.mag_jacobian(0.0, re, "near", 0, b * e, o, n, 1e-4)
    assert miss < 1e-9
    assert abs(m - 1.0) < 1e-3


def test_far_image_weight_matches_tracer():
    a, re = 2.1, 9.0
    o, n, e, gamma = _geometry(0.0)
    tables = xt._tables(a, re, "far", wf.b_crit(a) * (1 - 1e-9))
    images, _ = xt._images(a, re, "far", gamma, tables, 0)
    assert len(images) == 2
    for bs, sgn, dphidb, br in images:
        m, miss, _ = x1_s4b.mag_jacobian(a, re, "far", br, sgn * bs * e, o, n, 1e-4)
        assert miss < 1e-9
        assert abs(m / (xt.image_mag(re, bs, dphidb, gamma) / re**2) - 1.0) < 1e-5

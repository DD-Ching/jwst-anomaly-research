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

# far side near b_c: scipy may report round-off; accuracy is checked against the ODE route instead
FAR_WARN = pytest.mark.filterwarnings("ignore::scipy.integrate.IntegrationWarning")


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


@FAR_WARN
@pytest.mark.parametrize("a,frac", [(3.0, 0.5), (7.0, 0.99)])
def test_far_side_quadrature_matches_ode(a, frac):
    b = frac * wf.b_crit(a)
    assert abs(xt.phi_to_emitter(b, a, 9.0, "far")[0] / xt.phi_ode(b, a, 9.0, "far") - 1) < 1e-6


def test_near_side_branches_meet_at_fold():
    bf = xt._fold_b(0.0, 9.0) * (1 - 1e-10)
    pe, pa = xt.phi_to_emitter(bf, 0.0, 9.0, "near")
    assert abs(pe - pa) < 1e-4
    assert xt.phi_to_emitter(1.01 * xt._fold_b(0.0, 9.0), 0.0, 9.0, "near") == []


def test_face_on_primary_ring_near_r_plus_one():
    tr = xt.spot_track(0.0, 9.0, 179.5, "near", n_phase=12, primary_only=True)
    rad = np.hypot(tr[:, 1], tr[:, 2])
    assert abs(rad.mean() - 10.0) < 0.5


@FAR_WARN
@pytest.mark.parametrize("a", [2.1, 7.0])
def test_far_side_images_inside_b_crit(a):
    tr = xt.spot_track(a, 9.0, 155.0, "far", n_phase=12)
    assert np.nanmax(np.hypot(tr[:, 1], tr[:, 2])) < wf.b_crit(a)
    assert np.all(tr[:, 3] > 0)


@pytest.mark.parametrize("a,br", [(5.0, 0), (5.0, 1), (0.0, 1)])
def test_near_side_just_above_b_crit_matches_ode(a, br):
    b = 1.001 * wf.b_crit(a)
    q = xt.phi_to_emitter(b, a, 9.0, "near")[br]
    assert abs(q / xt.phi_ode(b, a, 9.0, "near", br) - 1) < 1e-6


def test_emitter_inside_throat_is_rejected():
    with pytest.raises(ValueError):
        xt.phi_to_emitter(3.0, 7.0, 6.0, "far")


def _t_ode(b, a, re, side, br):
    """Coordinate time to the emitter by ODE rays, on the same relative scale as time_to_emitter."""

    def rhs(_, y):
        r = y[0]
        rr = math.sqrt(r * r + a * a)
        wp = b * b * (-2 / rr**3 + 6 / rr**4) * r / rr
        return [y[1], -0.5 * wp, 1 / (1 - 2 / rr)]

    r_e = math.sqrt(re * re - a * a) * (1 if side == "near" else -1)

    def ev(_, y):
        return y[0] - r_e

    ev.terminal = side == "far" or br == 0
    ev.direction = 1.0 if (side == "near" and br == 1) else -1.0
    rr0 = math.sqrt(xt.R_START**2 + a * a)
    y0 = [xt.R_START, -math.sqrt(xt._pot(1 / rr0, b)), 0.0]
    sol = xt.integrate.solve_ivp(rhs, (0, 1e7), y0, events=ev, rtol=1e-12, atol=1e-12)
    return sol.y_events[0][0][2] + xt._t_open(b, a, 1 / rr0)


@FAR_WARN
@pytest.mark.parametrize(
    "a,side,br,b",
    [(0.0, "near", 0, 3.0), (0.0, "near", 1, 8.0), (5.0, "near", 1, 7.0), (3.0, "far", 0, 4.0)],
)
def test_light_travel_time_matches_ode(a, side, br, b):
    assert abs(xt.time_to_emitter(b, a, 9.0, side)[br] - _t_ode(b, a, 9.0, side, br)) < 1e-6


def test_primary_delay_is_roemer_at_large_radius():
    # weak field: the primary image's delay difference between the near and far points of the orbit
    # is the Roemer term 2 R sin(i) (Shapiro corrections are ~ ln(R) M, < 1 % of it here)
    re, inc = 400.0, 60.0
    tables = xt._tables(0.0, re, "near", 3.0 * re)
    t = {}
    for gamma in (math.radians(90 - inc), math.radians(90 + inc)):
        imgs, _ = xt._images(0.0, re, "near", gamma, tables)
        (img,) = [m for m in imgs if m[1] > 0]  # the primary (same-side) image
        t[gamma] = xt.time_to_emitter(img[0], 0.0, re, "near")[img[3]]
    (t_near, t_far) = t.values()
    assert abs((t_far - t_near) / (2 * re * math.sin(math.radians(inc))) - 1) < 0.01


@FAR_WARN
@pytest.mark.parametrize("a,side", [(0.0, "near"), (5.0, "far")])
def test_delayed_track_keeps_mean_flux_and_period(a, side):
    plain = xt.spot_track(a, 9.0, 155.0, side, n_phase=72)
    delayed = xt.spot_track(a, 9.0, 155.0, side, n_phase=72, delays=True)
    assert np.isfinite(delayed[:, 1:4]).all()
    assert abs(delayed[:, 3].mean() / plain[:, 3].mean() - 1) < 0.1
    # time delays move light between phases but cannot change the range much at this radius
    assert abs(delayed[:, 3].max() / plain[:, 3].max() - 1) < 0.2


def test_sub_keplerian_lowers_the_flux_contrast():
    fast = xt.spot_track(0.0, 9.0, 155.0, "near", n_phase=36, eta=1.0)
    slow = xt.spot_track(0.0, 9.0, 155.0, "near", n_phase=36, eta=0.6)
    assert slow[:, 3].max() / slow[:, 3].min() < fast[:, 3].max() / fast[:, 3].min()


def test_s2_observables_without_static_light():
    import x1_s2

    tr = xt.spot_track(0.0, 9.0, 155.0, "near", n_phase=36)
    obs = x1_s2.observables(tr)
    c = tr[:, 1:3] * x1_s2.THETA_G
    loop = np.hypot(*(c - c.mean(0)).T).mean()
    assert obs[0] == pytest.approx((loop, math.log(tr[:, 3].max() / tr[:, 3].min())))
    # static light always dilutes the flux contrast
    assert obs[:, 1].max() == pytest.approx(obs[0, 1])

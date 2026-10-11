"""Event-by-event time accounting for a traversable wormhole with one mouth carried to a star.

Goal-round-1 audit (docs/goal_round1/README.md). Nothing here says such a wormhole can be built or
kept open: that is the unsolved premise (energy conditions, docs/goal_round1/README.md section 1.4).
This only answers the bookkeeping question: *if* a short, rigid throat exists and one mouth rides
with the traveller, which Earth events, traveller times and Earth-frame times are connected, and
when does a closed causal loop appear?

Model (ASSUMPTIONS, the Morris–Thorne–Yurtsever 1988 idealisation):
- flat ambient spacetime, Earth frame S; mouth A at x = 0 for all t;
- mouth B rides the ship from x = 0 to x = D; the throat is short and its interior static, so an
  event on B at B's proper time tau is identified with the event on A at Earth time t = tau
  (both clocks start at 0 at departure); traversal time through the throat is neglected;
- the ship either moves at constant speed (instant acceleration) or with constant proper
  acceleration a to the midpoint and deceleration after it.

Units: years and light-years (c = 1).
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass

G_EARTH_LY_PER_YR2 = 1.0323  # 9.81 m/s^2 in ly / yr^2 (c = 1)


@dataclass(frozen=True)
class Trip:
    """Earth-frame duration and proper (traveller) duration of a one-way trip of length D."""

    distance_ly: float
    earth_years: float
    proper_years: float

    @property
    def shift_years(self) -> float:
        """Time shift between the mouths after arrival: Earth-frame time minus mouth B time."""
        return self.earth_years - self.proper_years


def coast_trip(distance_ly: float, beta: float) -> Trip:
    if not 0 < beta < 1:
        raise ValueError("beta must be in (0, 1)")
    t = distance_ly / beta
    return Trip(distance_ly, t, t * math.sqrt(1.0 - beta * beta))


def accel_trip(distance_ly: float, accel_ly_yr2: float = G_EARTH_LY_PER_YR2) -> Trip:
    """Constant proper acceleration to the midpoint, then constant deceleration."""
    if accel_ly_yr2 <= 0 or distance_ly <= 0:
        raise ValueError("distance and acceleration must be positive")
    a = accel_ly_yr2
    tau_half = math.acosh(1.0 + a * distance_ly / 2.0) / a
    t_half = math.sinh(a * tau_half) / a
    return Trip(distance_ly, 2.0 * t_half, 2.0 * tau_half)


def throat_to_earth(dest_earth_frame_time: float, trip: Trip) -> float:
    """Earth event reached by entering mouth B (at rest at x = D after arrival) at Earth-frame
    time t'. B's clock reads t' - shift, and A's event with the same reading is t = t' - shift."""
    if dest_earth_frame_time < trip.earth_years:
        raise ValueError("mouth B has not arrived yet")
    return dest_earth_frame_time - trip.shift_years


def earth_to_throat(earth_time: float, trip: Trip) -> float:
    """Destination event (Earth-frame time) reached by entering mouth A at Earth time t, for t at or
    after the traveller's arrival as seen through the throat (t >= proper trip time)."""
    if earth_time < trip.proper_years:
        raise ValueError("before arrival in throat-time: mouth B is still in flight")
    return earth_time + trip.shift_years


def loop_margin_years(
    shift_years: float, distance_ly: float, other_shift_years: float | None = None
):
    """Time by which the fastest closed route returns *after* it started (> 0: no causal loop).

    Routes with one carried wormhole (shift s) and ordinary light (speed 1): Earth -> light -> B ->
    throat -> Earth returns D - s later. With a second short connection of shift s2 between the
    same places (e.g. a wormhole moved slowly, s2 ~ 0): going out through the one with the larger
    shift and back through the other returns |s - s2| earlier, i.e. margin -|s - s2| < 0 whenever
    the shifts differ. A margin <= 0 means a closed causal curve."""
    margins = {"light_plus_throat": distance_ly - shift_years}
    if other_shift_years is not None:
        margins["two_throats"] = -abs(shift_years - other_shift_years)
    return margins


def timeline(trip: Trip, stay_years: float = 1.0, ship_return_beta: float | None = None) -> dict:
    """Named events with Earth-frame time t (S), traveller proper time tau and the mother's age
    increment (she stays at mouth A on Earth, so her age increment = Earth time t of the event she
    takes part in)."""
    s = trip.shift_years
    ev = {
        "depart": {"t_earth_frame": 0.0, "tau_traveller": 0.0, "mother_elapsed": 0.0},
        "arrive": {
            "t_earth_frame": trip.earth_years,
            "tau_traveller": trip.proper_years,
            # through the throat the arrival is seen by the mother at Earth time = tau
            "mother_elapsed_via_throat": trip.proper_years,
            "mother_elapsed_via_light": trip.earth_years + trip.distance_ly,
        },
    }
    # first two-way call after arrival: zero-length throat, so both clocks read tau on each side
    ev["first_call_after_arrival"] = {
        "t_earth_frame_destination": trip.earth_years,
        "earth_time_at_mother": trip.proper_years,
        "round_trip_latency_years": 0.0,
    }
    # return home through the throat after staying `stay_years` (traveller proper time)
    tau_ret = trip.proper_years + stay_years
    ev["home_via_throat"] = {
        "tau_traveller": tau_ret,
        "mother_elapsed": tau_ret,  # exits mouth A at Earth time = tau_ret
        "t_earth_frame_departure_from_destination": trip.earth_years + stay_years,
    }
    if ship_return_beta is not None:
        back = coast_trip(trip.distance_ly, ship_return_beta)
        ev["home_via_ship"] = {
            "tau_traveller": tau_ret + back.proper_years,
            "mother_elapsed": trip.earth_years + stay_years + back.earth_years,
        }
    ev["loop_margin_years"] = loop_margin_years(s, trip.distance_ly, other_shift_years=0.0)
    ev["shift_years"] = s
    ev["trip"] = asdict(trip)
    return ev


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--distance", type=float, nargs="+", default=[20.0, 100.0])
    ap.add_argument("--beta", type=float, nargs="+", default=[0.9, 0.99, 0.999])
    a = ap.parse_args(argv)
    out = []
    for d in a.distance:
        trips = [("coast", b, coast_trip(d, b)) for b in a.beta]
        trips.append(("1g", None, accel_trip(d)))
        for kind, b, tr in trips:
            ev = timeline(tr, stay_years=1.0, ship_return_beta=0.99)
            out.append({"distance_ly": d, "profile": kind, "beta": b, **ev})
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

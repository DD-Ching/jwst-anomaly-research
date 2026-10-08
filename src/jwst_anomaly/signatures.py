"""Registry of exotic-lens signatures and the survey adapters they run on (D-054).

Each signature of D-047 (W1 negative-mass radial pair with an empty umbra, W2 Ellis pair without a
deflector, W3 inverted-microlensing light curve, W5 count deficit) is one :class:`Signature` entry:
how it is predicted and injected (``exotic_sim``, ``simulated``), which screens implement it, on
which kind of survey data it runs, and where its vetting rules and limits are recorded. A survey
enters through a thin adapter that satisfies :class:`LightCurveSurvey` or
:class:`CatalogueSurvey`; screens take the adapter, never a survey-specific file layout.

Exotic physics is a hypothesis: a screen flag is an anomaly, not evidence, until every ordinary
explanation has been tested (/vet-candidate), and a null result is reported as a limit.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from astropy.table import Table

from jwst_anomaly import exotic_sim, schema

DATA_KINDS = ("catalogue", "light_curve", "image")


@dataclass(frozen=True)
class Signature:
    """One exotic signature: prediction, injection, screens, vetting and limits (D-047, D-054)."""

    code: str  # D-047 label, e.g. "W3"
    name: str
    data_kinds: tuple[str, ...]  # survey data it needs, from DATA_KINDS
    predict: Callable | None  # exotic_sim closed form (``simulated`` output)
    inject: Callable | None  # exotic_sim injector for injection-recovery
    screens: tuple[str, ...]  # implementations, "path:entry" (scripts are not importable)
    ordinary_mimics: tuple[str, ...]  # what vetting must rule out first, cheapest first
    limits_doc: str  # where its limits are recorded
    decisions: tuple[str, ...]

    def __post_init__(self) -> None:
        bad = set(self.data_kinds) - set(DATA_KINDS)
        if bad:
            raise ValueError(f"{self.code}: unknown data kinds {sorted(bad)}")


REGISTRY: dict[str, Signature] = {}


def register(sig: Signature) -> Signature:
    if sig.code in REGISTRY:
        raise ValueError(f"signature {sig.code} is already registered")
    REGISTRY[sig.code] = sig
    return sig


def get(code: str) -> Signature:
    try:
        return REGISTRY[code]
    except KeyError:
        raise KeyError(f"no signature {code!r}; registered: {sorted(REGISTRY)}") from None


def for_kind(kind: str) -> list[Signature]:
    """Signatures that can run on survey data of ``kind``."""
    if kind not in DATA_KINDS:
        raise ValueError(f"unknown data kind {kind!r}; one of {DATA_KINDS}")
    return [s for s in REGISTRY.values() if kind in s.data_kinds]


_LIMITS = "docs/exotic_limits.md"

register(
    Signature(
        code="W1",
        name="negative-mass lens: radially stretched images beside an empty umbra",
        data_kinds=("catalogue", "image"),
        predict=exotic_sim.solve_images,
        inject=exotic_sim.inject_images,
        screens=(
            "scripts/exotic_screens.py:radial",
            "scripts/inject_radial.py",
            "scripts/orphan_pairs.py",
        ),
        ordinary_mimics=(
            "diffraction spikes and stripes",
            "cluster radial arcs",
            "edge-on or irregular galaxies",
            "voids / underdensities",
        ),
        limits_doc=_LIMITS,
        decisions=("D-031", "D-047", "D-049", "D-051", "D-053"),
    )
)
register(
    Signature(
        code="W2",
        name="Ellis-wormhole image pair with no visible deflector",
        data_kinds=("catalogue", "image"),
        predict=exotic_sim.solve_images,
        inject=exotic_sim.inject_images,
        screens=("scripts/orphan_pairs.py", "scripts/inject_pairs.py"),
        ordinary_mimics=(
            "knots of one galaxy",
            "physical companions and groups",
            "faint or dark ordinary lens galaxies",
            "chance SED matches",
        ),
        limits_doc=_LIMITS,
        decisions=("D-047", "D-048", "D-051"),
    )
)
register(
    Signature(
        code="W3",
        name="inverted microlensing: flux vanishes between caustic spikes",
        data_kinds=("light_curve",),
        predict=exotic_sim.light_curve,
        inject=exotic_sim.inject_light_curve,
        screens=("scripts/dimming_screen.py",),
        ordinary_mimics=(
            "binary-lens caustic crossings",
            "blending and photometric systematics",
            "eclipsing and other variable stars",
            "parallax and xallarap",
            "persistence and saturated-star wings",
        ),
        limits_doc=_LIMITS,
        decisions=("D-047", "D-052", "D-054"),
    )
)
register(
    Signature(
        code="W5",
        name="background-count deficit inside about theta_E",
        data_kinds=("catalogue",),
        predict=exotic_sim.count_ratio,
        inject=None,
        screens=(),
        ordinary_mimics=("masks and bright-star halos", "deblending", "cosmic variance"),
        limits_doc=_LIMITS,
        decisions=("D-047",),
    )
)


# ----------------------------------------------------------------------------- survey adapters

LIGHT_CURVE_COLUMNS = ("time", "mag", "mag_err", "band")


@runtime_checkable
class LightCurveSurvey(Protocol):
    """A time-domain survey: an event or source list and one light curve per entry.

    ``events()`` has at least ``event_id, ra, dec`` (deg), plus any published fit parameters;
    ``light_curve(event_id)`` returns :func:`standard_light_curve` output. ``efficiency(t_e)``
    is the survey's published detection efficiency for an event time scale (days), or None when
    the survey publishes none (limits then need injection-recovery on the survey's cadence).
    """

    name: str

    def events(self) -> Table: ...

    def light_curve(self, event_id: str) -> Table: ...

    def efficiency(self, t_e_days: float) -> float | None: ...


@runtime_checkable
class CatalogueSurvey(Protocol):
    """An imaging survey: a source catalogue and the area it covers."""

    name: str

    def catalogue(self) -> Table: ...

    def area_deg2(self) -> float: ...


def standard_light_curve(time, mag, mag_err, band, source: str, time_system: str) -> Table:
    """One light curve in the shared layout, ``observed`` (finite rows only, time-sorted).

    ``time_system`` names the time column (e.g. "HJD - 2450000"), recorded in ``meta``.
    """
    t = np.asarray(time, float)
    m = np.asarray(mag, float)
    e = np.asarray(mag_err, float)
    if not (t.shape == m.shape == e.shape):
        raise ValueError("time, mag and mag_err must have one shape")
    b = np.broadcast_to(np.asarray(band, str), t.shape)
    ok = np.isfinite(t) & np.isfinite(m) & np.isfinite(e) & (e > 0)
    order = np.argsort(t[ok], kind="stable")
    out = Table(
        {
            "time": t[ok][order],
            "mag": m[ok][order],
            "mag_err": e[ok][order],
            "band": b[ok][order],
        }
    )
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value, source=source, time_system=time_system
    )
    return out

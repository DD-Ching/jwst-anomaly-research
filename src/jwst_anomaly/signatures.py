"""Registry of exotic-lens signatures and the survey adapters they run on (D-054).

Each signature of D-047 (W1 negative-mass radial pair with an empty umbra, W2 Ellis pair without a
deflector, W3 inverted-microlensing light curve, W5 count deficit) is one :class:`Signature` entry:
how it is predicted and injected (``exotic_sim``, ``simulated``), which screens implement it, on
which kind of survey data it runs, and where its vetting rules and limits are recorded. A new
survey enters through a thin adapter that satisfies :class:`LightCurveSurvey` or
:class:`CatalogueSurvey`, and new screens take the adapter. The JWST screens registered below
predate the layer and still read JWST level-3 catalogues directly.

Exotic physics is a hypothesis: a screen flag is an anomaly, not evidence, until every ordinary
explanation has been tested (/vet-candidate), and a null result is reported as a limit.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from functools import partial
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
    predict: Callable | None  # exotic_sim closed form, exotic lens parameters bound
    inject: Callable | None  # exotic_sim injector, exotic lens parameters bound
    screens: tuple[str, ...]  # implementations: "script [subcommand]" (scripts are not importable)
    ordinary_mimics: tuple[str, ...]  # what vetting must rule out first, cheapest first
    limits_doc: str  # where its limits are recorded ("" while none exist)
    decisions: tuple[str, ...]
    lens: dict = field(default_factory=dict)  # exotic_sim lens parameters (n, sign)

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
_NEG = {"n": 1.0, "sign": -1}  # negative mass (Kitamura+2013 n = 1, eps < 0)
_ELLIS = {"n": 2.0, "sign": 1}  # Ellis wormhole (n = 2)


def _bound(fn: Callable, lens: dict) -> Callable:
    """``fn`` with the signature's exotic lens parameters fixed (exotic_sim defaults are an
    ordinary positive point mass)."""
    return partial(fn, **lens)


register(
    Signature(
        code="W1",
        name="negative-mass lens: radially stretched images beside an empty umbra",
        data_kinds=("catalogue", "image"),
        predict=_bound(exotic_sim.solve_images, _NEG),
        inject=_bound(exotic_sim.inject_images, _NEG),
        screens=(
            "scripts/exotic_screens.py radial",
            "scripts/inject_radial.py",
            "scripts/exotic_screens.py shear",
            "scripts/inject_shear.py",
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
        lens=_NEG,
    )
)
register(
    Signature(
        code="W2",
        name="Ellis-wormhole image pair with no visible deflector",
        data_kinds=("catalogue", "image"),
        predict=_bound(exotic_sim.solve_images, _ELLIS),
        inject=_bound(exotic_sim.inject_images, _ELLIS),
        screens=(
            "scripts/exotic_screens.py fluxratio",
            "scripts/orphan_pairs.py",
            "scripts/inject_pairs.py",
            "scripts/w12_lenscats.py",  # published lens catalogues, no visible deflector
        ),
        ordinary_mimics=(
            "knots of one galaxy",
            "physical companions and groups",
            "faint or dark ordinary lens galaxies",
            "chance SED matches",
            "catalogue position errors and blended lens light",
        ),
        limits_doc=_LIMITS,
        decisions=("D-047", "D-048", "D-051", "D-056"),
        lens=_ELLIS,
    )
)
register(
    Signature(
        code="W3",
        name="inverted microlensing: flux vanishes between caustic spikes",
        data_kinds=("light_curve",),
        predict=_bound(exotic_sim.light_curve, _NEG),
        inject=_bound(exotic_sim.inject_light_curve, _NEG),
        screens=(
            "scripts/dimming_screen.py",
            "scripts/w3_microlensing.py",
            "scripts/w3_moa.py",  # MOA-II light curves before any bump cut (D-060)
        ),
        ordinary_mimics=(
            "binary-lens caustic crossings",
            "blending and photometric systematics",
            "eclipsing and other variable stars",
            "parallax and xallarap",
            "persistence and saturated-star wings",
        ),
        limits_doc=_LIMITS,
        decisions=("D-047", "D-052", "D-054", "D-057", "D-060"),
        lens=_NEG,
    )
)
register(
    Signature(
        code="W5",
        name="background-count deficit inside about theta_E",
        data_kinds=("catalogue",),
        predict=_bound(exotic_sim.count_ratio, _NEG),
        inject=None,
        screens=(),
        ordinary_mimics=("masks and bright-star halos", "deblending", "cosmic variance"),
        limits_doc="",  # no W5 screen or limit yet
        decisions=("D-047",),
        lens=_NEG,
    )
)


# ----------------------------------------------------------------------------- survey adapters

LIGHT_CURVE_COLUMNS = schema.LIGHT_CURVE_COLUMNS
LIGHT_CURVE_FLUX_COLUMNS = schema.LIGHT_CURVE_FLUX_COLUMNS


@runtime_checkable
class LightCurveSurvey(Protocol):
    """A time-domain survey: an event or source list and one light curve per entry.

    ``events()`` has at least ``event_id, ra, dec`` (deg), plus any published fit parameters;
    ``light_curve(event_id)`` returns :func:`standard_light_curve` output, or
    :func:`standard_flux_light_curve` output for difference-imaging surveys whose flux relative to
    a reference image can be negative (MOA-II, D-060). ``efficiency(t_e)``
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

    ``time_system`` names the time column (e.g. "HJD - 2450000"), recorded in ``meta`` with the
    number of rows dropped as non-finite or with non-positive errors (``n_dropped``).
    """
    t = np.asarray(time, float)
    m = np.asarray(mag, float)
    e = np.asarray(mag_err, float)
    if not (t.shape == m.shape == e.shape):
        raise ValueError("time, mag and mag_err must have one shape")
    b = np.broadcast_to(np.asarray(band, str), t.shape)
    ok = np.isfinite(t) & np.isfinite(m) & np.isfinite(e) & (e > 0)
    order = np.argsort(t[ok], kind="stable")
    cols = (t, m, e, b)
    out = Table({k: c[ok][order] for k, c in zip(LIGHT_CURVE_COLUMNS, cols, strict=True)})
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=source,
        time_system=time_system,
        n_dropped=int((~ok).sum()),  # non-finite or non-positive-error rows (cadence audit)
    )
    return out


def standard_flux_light_curve(
    time,
    flux,
    flux_err,
    band,
    source: str,
    time_system: str,
    flux_unit: str,
    flux_kind: str = "difference",
) -> Table:
    """One flux light curve in the shared layout, ``observed`` (finite rows only, time-sorted).

    For difference-imaging photometry the flux is relative to a reference image and can be negative,
    so it has no magnitude. ``flux_unit`` and ``flux_kind`` (e.g. "difference", "detrended
    difference") are recorded in ``meta`` with ``time_system`` and ``n_dropped`` (non-finite rows or
    non-positive errors).
    """
    t = np.asarray(time, float)
    f = np.asarray(flux, float)
    e = np.asarray(flux_err, float)
    if not (t.shape == f.shape == e.shape):
        raise ValueError("time, flux and flux_err must have one shape")
    b = np.broadcast_to(np.asarray(band, str), t.shape)
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0)
    order = np.argsort(t[ok], kind="stable")
    cols = (t, f, e, b)
    out = Table({k: c[ok][order] for k, c in zip(LIGHT_CURVE_FLUX_COLUMNS, cols, strict=True)})
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=source,
        time_system=time_system,
        flux_unit=flux_unit,
        flux_kind=flux_kind,
        n_dropped=int((~ok).sum()),
    )
    return out

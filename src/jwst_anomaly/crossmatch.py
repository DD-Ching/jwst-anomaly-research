"""Cross-check targets against external catalogs (SIMBAD, NED, Gaia DR3). Owner: bootstrap unit 5.

Primary backends are batched; a per-object loop is only the last-resort NED fallback
(timings and rejected alternatives: DECISIONS.md D-006):

* ``simbad``: CDS XMatch against ``simbad`` (one request); fallback SIMBAD TAP via
  ``Simbad.query_region``.
* ``ned``: NED TAP ``NEDTAP.objdir``, 100 OR'd ADQL cones per request; fallback
  ``Ned.query_region`` per target (at most 100 targets).
* ``gaia``: CDS XMatch against VizieR ``I/355/gaiadr3`` (one request); fallback ESA Gaia archive
  TAP (pyvo) on ``gaiadr3.gaia_source``, 100 cones per request.

Raw rows are normalized to one object table per service and associated with the targets
locally (``astropy.coordinates.search_around_sky``), so all backends share one separation
definition. Every request has a client timeout and is retried with exponential backoff; a
service whose backends all fail is recorded in ``meta["services_failed"]`` and the others
continue. Normalized results are cached as ECSV under ``paths.cache_dir() / "crossmatch"``
(astroquery's own HTTP cache is bypassed so ``meta["query_utc"]`` is the real query time).

``query_matches`` returns the long format (one row per target-object pair, ``MATCH_COLUMNS``);
``crossmatch`` returns the ``schema.XMATCH_COLUMNS`` summary. Catalog content is "observed".
``is_star`` and ``is_lens_related`` are derived flags whose thresholds are ASSUMPTIONS
(documented next to each constant). Gaia positions are at epoch J2016.0 and are not propagated
to the JWST epoch.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table, unique, vstack

from jwst_anomaly import paths, schema

log = logging.getLogger(__name__)

SERVICES = ("simbad", "ned", "gaia")

# Backend chain per service: the first entry is the primary, later entries are fallbacks.
DEFAULT_BACKENDS: dict[str, tuple[str, ...]] = {
    "simbad": ("xmatch", "tap"),
    "ned": ("tap", "cone"),
    "gaia": ("xmatch", "tap"),
}

# Exact catalog/service identifiers recorded in meta["catalog"] (see SOURCES.md, unit 5).
CATALOGS: dict[tuple[str, str], str] = {
    ("simbad", "xmatch"): "CDS XMatch cat2='simbad'",
    ("simbad", "tap"): "SIMBAD TAP basic + alltypes (simbad.cds.unistra.fr)",
    ("ned", "tap"): "NED TAP NEDTAP.objdir (ned.ipac.caltech.edu/tap)",
    ("ned", "cone"): "NED object search, astroquery.ipac.ned Ned.query_region",
    ("gaia", "xmatch"): "CDS XMatch cat2='vizier:I/355/gaiadr3' (Gaia DR3, epoch J2016.0)",
    ("gaia", "tap"): "ESA Gaia archive TAP gaiadr3.gaia_source (Gaia DR3, epoch J2016.0)",
}

NED_TAP_URL = "https://ned.ipac.caltech.edu/tap"
GAIA_TAP_URL = "https://gea.esac.esa.int/tap-server/tap"

# --- Classification constants -------------------------------------------------------------
# SIMBAD object types, frozen from SIMBAD TAP table ``otypedef`` on 2026-10-07 with
# ``derive_simbad_otype_sets`` (a network test re-derives them to detect drift). Definitions:
# https://simbad.cds.unistra.fr/guide/otypes.htx
# fmt: off
SIMBAD_STAR_OTYPES = frozenset({
    "*", "**", "AB*", "Ae*", "BD*", "BS*", "BY*", "Be*", "C*", "CV*", "Ce*", "EB*", "El*",
    "Em*", "Er*", "Ev*", "HB*", "HS*", "HV*", "HXB", "Ir*", "LM*", "LP*", "LXB", "MS*",
    "Ma*", "Mi*", "N*", "No*", "OH*", "Or*", "PM*", "Pe*", "Psr", "Pu*", "RC*", "RG*",
    "RR*", "RS*", "RV*", "Ro*", "S*", "SB*", "SX*", "Sy*", "TT*", "V*", "WD*", "WR*",
    "WV*", "XB*", "Y*O", "a2*", "bC*", "cC*", "dS*", "gD*", "pA*", "s*b", "s*r", "s*y",
    "sg*",
})
# fmt: on
SIMBAD_LENS_OTYPES = frozenset({"gLS", "LS?", "gLe", "Le?", "LeI", "LI?", "LeG", "LeQ", "Lev"})
# ASSUMPTION: these sub-branches of SIMBAD's "*" branch are not point-like stars.
_SIMBAD_NON_STELLAR_BRANCHES = frozenset({"PN", "SN*", "Pl", "out"})
# NED object types (https://ned.ipac.caltech.edu/help/ui/nearposn-list_objecttypes):
# "G_Lens" lensed image of a galaxy, "Q_Lens" lensed image of a QSO. NED's "*" means
# "Star or Point Source", so NED does not feed ``is_star``.
NED_LENS_TYPES = frozenset({"G_Lens", "Q_Lens"})
# ASSUMPTION: a Gaia DR3 source is an astrometric star if parallax/error >= 5 or the
# proper-motion significance (chi^2 with the pmra-pmdec correlation, square-rooted) >= 5.
GAIA_MIN_PARALLAX_SNR = 5.0
GAIA_MIN_PM_SNR = 5.0
# Services whose classifications feed is_star (the nearest of their matches decides).
STAR_SERVICES = ("simbad", "gaia")

STAR_CRITERIA = (
    "ASSUMPTION: is_star = the nearest SIMBAD/Gaia match is a star: SIMBAD main otype in "
    "SIMBAD_STAR_OTYPES (non-candidate '*' branch minus PN/SN*/Pl/out) or Gaia DR3 "
    f"parallax_over_error >= {GAIA_MIN_PARALLAX_SNR:g} or pm_snr >= {GAIA_MIN_PM_SNR:g}; "
    "NED types are not used"
)
LENS_CRITERIA = (
    "is_lens_related = any match with a SIMBAD otype (main or other) in SIMBAD_LENS_OTYPES "
    "(grv > gLS branch incl. candidates, plus Lev) or a NED type in {G_Lens, Q_Lens}"
)
BEST_MATCH_RULE = "smallest separation; ties broken by the order of `services`"

# Long format returned by query_matches (one row per target-object pair).
MATCH_COLUMNS = (
    "source_uid",
    "service",
    "match_id",
    "match_type",  # SIMBAD main otype | NED type | Gaia "astrometric_star"/"unclassified"
    "all_types",  # SIMBAD otypes joined by "|" (main + other types); "" for NED/Gaia
    "lens_types",  # lens-related codes of this match joined by "|"
    "sep_arcsec",
    "match_ra",
    "match_dec",
    "is_star",
    "is_lens_related",
    "redshift",  # NED preferred redshift
    "parallax",  # Gaia DR3 [mas]
    "parallax_error",
    "pmra",  # Gaia DR3 [mas/yr]
    "pmra_error",
    "pmdec",
    "pmdec_error",
    "pmra_pmdec_corr",
    "parallax_over_error",  # derived
    "pm_snr",  # derived
)
# Extra summary columns beyond schema.XMATCH_COLUMNS. n_<service> == -1: not queried or failed.
SUMMARY_EXTRA_COLUMNS = ("n_simbad", "n_ned", "n_gaia", "is_lens_related", "lens_types")

_STR_COLS = ("match_id", "match_type", "all_types")
_FLOAT_COLS = (
    "match_ra",
    "match_dec",
    "redshift",
    "parallax",
    "parallax_error",
    "pmra",
    "pmra_error",
    "pmdec",
    "pmdec_error",
    "pmra_pmdec_corr",
)

# Raw service column -> normalized object column, per (service, backend).
RAW_COLUMNS: dict[tuple[str, str], dict[str, str]] = {
    ("simbad", "xmatch"): {
        "match_id": "main_id",
        "match_type": "otype",
        "all_types": "other_types",
        "match_ra": "ra",
        "match_dec": "dec",
    },
    ("simbad", "tap"): {
        "match_id": "main_id",
        "match_type": "otype",
        "all_types": "alltypes.otypes",
        "match_ra": "ra",
        "match_dec": "dec",
    },
    ("ned", "tap"): {
        "match_id": "prefname",
        "match_type": "prefphytype",
        "match_ra": "ra",
        "match_dec": "dec",
        "redshift": "z",
    },
    ("ned", "cone"): {
        "match_id": "Object Name",
        "match_type": "Type",
        "match_ra": "RA",
        "match_dec": "DEC",
        "redshift": "Redshift",
    },
    ("gaia", "xmatch"): {
        "match_id": "DR3Name",
        "match_ra": "RAdeg",
        "match_dec": "DEdeg",
        "parallax": "Plx",
        "parallax_error": "e_Plx",
        "pmra": "pmRA",
        "pmra_error": "e_pmRA",
        "pmdec": "pmDE",
        "pmdec_error": "e_pmDE",
        "pmra_pmdec_corr": "pmRApmDEcor",
    },
    ("gaia", "tap"): {
        "match_id": "source_id",
        "match_ra": "ra",
        "match_dec": "dec",
        "parallax": "parallax",
        "parallax_error": "parallax_error",
        "pmra": "pmra",
        "pmra_error": "pmra_error",
        "pmdec": "pmdec",
        "pmdec_error": "pmdec_error",
        "pmra_pmdec_corr": "pmra_pmdec_corr",
    },
}
_ID_PREFIX = {("gaia", "tap"): "Gaia DR3 "}  # match the VizieR "DR3Name" form

_ADQL_CHUNK = 100  # cones per OR'd ADQL request (~2 s per NED TAP chunk, measured 2026-10-07)
# Targets per request; backends not listed send all targets in one request.
_CHUNK_SIZE = {("ned", "tap"): _ADQL_CHUNK, ("gaia", "tap"): _ADQL_CHUNK, ("ned", "cone"): 1}
_MAX_TARGETS = {("ned", "cone"): 100}  # cap for the per-target NED fallback
_TAP_MAXREC = 100_000  # explicit MAXREC for TAP queries so truncation is detectable
_CACHE_FORMAT = 1


class CrossmatchError(RuntimeError):
    """Raised when every requested service failed."""


class _PermanentError(RuntimeError):
    """A backend failure that a retry cannot fix (limits, truncation)."""


# Not retried: deterministic failures (bad input, limits, schema changes).
_NO_RETRY = (_PermanentError, ValueError, TypeError, KeyError)


@dataclass
class _Options:
    timeout_s: float
    retries: int
    backoff_s: float


@dataclass
class _ServiceResult:
    objects: Table | None
    backend: str | None = None
    query_utc: str | None = None
    from_cache: bool = False
    elapsed_s: float = 0.0
    errors: list[str] = field(default_factory=list)


# --- Public API ---------------------------------------------------------------------------


def crossmatch(
    targets: Table,
    radius_arcsec: float = 1.0,
    services: Sequence[str] = ("simbad", "ned", "gaia"),
    *,
    backends: Mapping[str, str | Sequence[str]] | None = None,
    cache: bool = True,
    cache_dir: str | Path | None = None,
    cache_max_age_days: float | None = 30.0,
    refresh: bool = False,
    timeout_s: float = 60.0,
    retries: int = 2,
    backoff_s: float = 2.0,
) -> Table:
    """Match ``targets`` (``schema.TARGET_COLUMNS``) against external services.

    Returns one summary row per target with ``schema.XMATCH_COLUMNS`` plus
    ``SUMMARY_EXTRA_COLUMNS``. ``meta["provenance"] == "observed"`` (external catalog content)
    with services and query dates recorded in ``meta``. Equivalent to
    ``summarize_matches(targets, query_matches(targets, ...))``; call those two directly to keep
    the per-match long format. Keyword arguments are documented in :func:`query_matches`.
    """
    matches = query_matches(
        targets,
        radius_arcsec,
        services,
        backends=backends,
        cache=cache,
        cache_dir=cache_dir,
        cache_max_age_days=cache_max_age_days,
        refresh=refresh,
        timeout_s=timeout_s,
        retries=retries,
        backoff_s=backoff_s,
    )
    return summarize_matches(targets, matches)


def query_matches(
    targets: Table,
    radius_arcsec: float = 1.0,
    services: Sequence[str] = SERVICES,
    *,
    backends: Mapping[str, str | Sequence[str]] | None = None,
    cache: bool = True,
    cache_dir: str | Path | None = None,
    cache_max_age_days: float | None = 30.0,
    refresh: bool = False,
    timeout_s: float = 60.0,
    retries: int = 2,
    backoff_s: float = 2.0,
) -> Table:
    """Return every external object within ``radius_arcsec`` of each target (``MATCH_COLUMNS``).

    ``backends`` overrides the backend chain of a service, e.g. ``{"gaia": "tap"}``.
    ``cache`` reads/writes normalized results under ``cache_dir`` (default
    ``paths.cache_dir() / "crossmatch"``); a backend's cache entry is used if it is younger than
    ``cache_max_age_days`` (None: no expiry), otherwise that backend is queried, and a fallback's
    cache entry is only consulted after the backends before it failed. ``refresh=True`` skips
    cache reads. Each HTTP request gets ``timeout_s`` and up to ``retries`` retries after
    ``backoff_s * 2**k`` seconds.
    Raises :class:`CrossmatchError` if every requested service fails.
    """
    _check_targets(targets)
    services = _check_services(services)
    chains = _resolve_backends(services, backends)
    if not (math.isfinite(radius_arcsec) and radius_arcsec > 0):
        raise ValueError(f"radius_arcsec must be positive and finite, got {radius_arcsec}")
    if retries < 0 or backoff_s < 0 or not timeout_s > 0:
        raise ValueError("need retries >= 0, backoff_s >= 0 and timeout_s > 0")
    radius_arcsec = float(radius_arcsec)
    ra, dec = _deg(targets["ra"]), _deg(targets["dec"])
    opts = _Options(timeout_s=float(timeout_s), retries=int(retries), backoff_s=float(backoff_s))
    cdir = Path(cache_dir) if cache_dir is not None else paths.cache_dir() / "crossmatch"

    coords = SkyCoord(ra, dec, unit="deg")
    results: dict[str, _ServiceResult] = {}
    for svc in services if len(targets) else ():
        results[svc] = _run_service(
            svc,
            chains[svc],
            coords,
            radius_arcsec,
            opts,
            cache=cache,
            refresh=refresh,
            cache_dir=cdir,
            max_age_days=cache_max_age_days,
        )
    ok = [s for s in services if s in results and results[s].objects is not None]
    failed = [s for s in services if s in results and results[s].objects is None]
    if len(targets) and not ok:
        detail = "; ".join(f"{s}: {' | '.join(results[s].errors)}" for s in failed)
        raise CrossmatchError(f"all external services failed: {detail}")

    parts = [
        _associate(targets, coords, _classify(svc, results[svc].objects), svc, radius_arcsec)
        for svc in ok
    ]
    matches = vstack(parts, metadata_conflicts="silent") if parts else _empty_matches(targets)
    matches = matches[list(MATCH_COLUMNS)]
    matches.meta.clear()
    matches.meta.update(
        _meta(targets, services, ok, failed, results, radius_arcsec, cache=cache, refresh=refresh)
    )
    return schema.validate(matches, MATCH_COLUMNS, name="crossmatch.query_matches")


def summarize_matches(targets: Table, matches: Table) -> Table:
    """Reduce ``query_matches`` output to one ``schema.XMATCH_COLUMNS`` row per target.

    ``is_known_object``: any match from a service that answered. ``best_match_*``: the
    match with the smallest separation (ties: order of ``meta["services_requested"]``).
    ``is_star``: the nearest SIMBAD/Gaia match is a star (``STAR_CRITERIA``).
    ``n_<service>`` is -1 when that service was not queried or failed; this needs
    ``matches.meta["services_requested"]`` and ``["services_ok"]`` as set by ``query_matches``.
    """
    _check_targets(targets)
    if "services_requested" not in matches.meta or "services_ok" not in matches.meta:
        raise ValueError("matches.meta lacks services_requested/services_ok (use query_matches)")
    requested = list(matches.meta["services_requested"])
    services_ok = set(matches.meta["services_ok"])
    rank = {s: i for i, s in enumerate(requested)}
    row_of = {uid: i for i, uid in enumerate(targets["source_uid"])}
    n = len(targets)

    n_matches = np.zeros(n, dtype=int)
    per_service = {s: np.zeros(n, dtype=int) for s in SERVICES}
    is_lens = np.zeros(n, dtype=bool)
    lens_types: list[set[str]] = [set() for _ in range(n)]
    # Keys are (sep, service rank, match row); the minimum wins.
    best: list[tuple[float, int, int] | None] = [None] * n
    best_star_service: list[tuple[float, int, int] | None] = [None] * n

    for j, m in enumerate(matches):
        i = row_of[m["source_uid"]]
        svc = str(m["service"])
        n_matches[i] += 1
        per_service[svc][i] += 1
        if m["is_lens_related"]:
            is_lens[i] = True
            lens_types[i].update(f"{svc}:{t}" for t in str(m["lens_types"]).split("|") if t)
        key = (float(m["sep_arcsec"]), rank.get(svc, len(rank)), j)
        if best[i] is None or key < best[i]:
            best[i] = key
        if svc in STAR_SERVICES and (best_star_service[i] is None or key < best_star_service[i]):
            best_star_service[i] = key

    def best_value(col: str, default: Any) -> list[Any]:
        return [matches[b[2]][col] if b is not None else default for b in best]

    summary = Table()
    summary["source_uid"] = targets["source_uid"].copy()
    summary["n_matches"] = n_matches
    summary["is_known_object"] = n_matches > 0
    summary["is_star"] = np.array(
        [bool(matches[b[2]]["is_star"]) if b is not None else False for b in best_star_service],
        dtype=bool,
    )
    for col in ("match_id", "match_type", "service"):
        values = [str(v) for v in best_value(col, "")]
        summary[f"best_match_{col.removeprefix('match_')}"] = np.array(values, dtype=str)
    summary["best_match_sep_arcsec"] = np.array(best_value("sep_arcsec", np.nan), dtype=float)
    for svc in SERVICES:
        summary[f"n_{svc}"] = per_service[svc] if svc in services_ok else np.full(n, -1)
    summary["is_lens_related"] = is_lens
    summary["lens_types"] = np.array([",".join(sorted(s)) for s in lens_types], dtype=str)
    summary.meta.update(matches.meta)
    return schema.validate(summary, schema.XMATCH_COLUMNS, name="crossmatch")


def derive_simbad_otype_sets(otypedef: Table) -> tuple[frozenset[str], frozenset[str]]:
    """Return ``(stellar, lens_related)`` otype codes from SIMBAD TAP table ``otypedef``.

    Stellar (ASSUMPTION): non-candidate types whose ``path`` starts in the ``*`` branch,
    excluding the PN, SN*, Pl and out sub-branches. Lens-related: every type in the
    ``grv > gLS`` branch (candidates included) plus ``Lev`` (lensing event).
    """
    stellar, lens = set(), set()
    for row in otypedef:
        otype = str(row["otype"]).strip()
        branch = [p.strip() for p in str(row["path"]).split(">")]
        if branch[0] == "*" and not int(row["is_candidate"]):
            if not _SIMBAD_NON_STELLAR_BRANCHES.intersection(branch):
                stellar.add(otype)
        if branch[0] == "grv" and ("gLS" in branch or otype == "Lev"):
            lens.add(otype)
    return frozenset(stellar), frozenset(lens)


def gaia_astrometric_snr(
    parallax: np.ndarray,
    parallax_error: np.ndarray,
    pmra: np.ndarray,
    pmra_error: np.ndarray,
    pmdec: np.ndarray,
    pmdec_error: np.ndarray,
    pmra_pmdec_corr: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(parallax_over_error, pm_snr)``; NaN where the solution has no astrometry.

    ``pm_snr = sqrt(mu^T C^-1 mu)`` with ``C`` built from the errors and the pmra-pmdec
    correlation (NaN correlation treated as 0).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        plx_snr = np.asarray(parallax, float) / np.asarray(parallax_error, float)
        x = np.asarray(pmra, float) / np.asarray(pmra_error, float)
        y = np.asarray(pmdec, float) / np.asarray(pmdec_error, float)
        rho = np.nan_to_num(np.asarray(pmra_pmdec_corr, float), nan=0.0)
        chi2 = (x**2 - 2 * rho * x * y + y**2) / (1 - rho**2)
        pm_snr = np.sqrt(chi2)
    return plx_snr, pm_snr


# --- Service orchestration ----------------------------------------------------------------


def _run_service(
    service: str,
    chain: Sequence[str],
    coords: SkyCoord,
    radius_arcsec: float,
    opts: _Options,
    *,
    cache: bool,
    refresh: bool,
    cache_dir: Path,
    max_age_days: float | None,
) -> _ServiceResult:
    """Walk the backend chain: fresh cache entry, else query; first success wins."""
    ra, dec = coords.ra.deg, coords.dec.deg
    result = _ServiceResult(objects=None)
    for backend in chain:
        path = _cache_path(cache_dir, service, backend, ra, dec, radius_arcsec)
        if cache and not refresh:
            hit = _cache_read(path, max_age_days)
            if hit is not None:
                log.info("%s: using cached %s result %s", service, backend, path.name)
                result.objects, result.backend, result.from_cache = hit, backend, True
                result.query_utc = hit.meta["query_utc"]
                return result
        query_utc = datetime.now(UTC).isoformat(timespec="seconds")
        t0 = time.perf_counter()
        try:
            objs = _query_backend(service, backend, coords, radius_arcsec, opts)
        except Exception as exc:  # noqa: BLE001 -- any service error must not stop the others
            result.errors.append(f"{backend}: {type(exc).__name__}: {exc}"[:500])
            log.warning("%s/%s failed: %s: %s", service, backend, type(exc).__name__, exc)
            continue
        finally:
            result.elapsed_s += time.perf_counter() - t0
        result.objects, result.backend, result.query_utc = objs, backend, query_utc
        if cache:
            _cache_write(path, objs, service, backend, query_utc, radius_arcsec)
        return result
    return result


def _query_backend(
    service: str, backend: str, coords: SkyCoord, radius_arcsec: float, opts: _Options
) -> Table:
    """Send the backend's requests (retrying each one) and normalize the combined rows."""
    key = (service, backend)
    limit = _MAX_TARGETS.get(key)
    if limit is not None and len(coords) > limit:
        raise _PermanentError(f"{service}/{backend} is capped at {limit} targets ({len(coords)})")
    fetch = _FETCHERS[key]
    parts = []
    for k, chunk in enumerate(_chunks(coords, _CHUNK_SIZE.get(key) or len(coords))):
        parts.append(
            _with_retries(
                lambda chunk=chunk: fetch(chunk, radius_arcsec, opts),
                retries=opts.retries,
                backoff_s=opts.backoff_s,
                what=f"{service}/{backend} request {k + 1}",
            )
        )
    nonempty = [p for p in parts if len(p)] or parts[:1]
    raw = nonempty[0] if len(nonempty) == 1 else vstack(nonempty, metadata_conflicts="silent")
    return _normalize(service, backend, raw)


def _with_retries(fn: Callable[[], Table], *, retries: int, backoff_s: float, what: str) -> Table:
    attempt = 0
    while True:
        try:
            return fn()
        except _NO_RETRY:
            raise
        except Exception as exc:  # noqa: BLE001 -- remote services raise many exception types
            if attempt >= retries:
                raise
            delay = backoff_s * 2**attempt
            attempt += 1
            log.warning(
                "%s attempt %d/%d failed (%s: %s); retrying in %.1f s",
                *(what, attempt, retries + 1, type(exc).__name__, exc, delay),
            )
            time.sleep(delay)


# --- Backends: one request each (monkeypatched in offline tests) ---------------------------


def _fetch_xmatch(cat2: str) -> Callable[[SkyCoord, float, _Options], Table]:
    def fetch(coords: SkyCoord, radius_arcsec: float, opts: _Options) -> Table:
        from astroquery.xmatch import XMatchClass

        upload = Table(
            {"xm_idx": np.arange(len(coords)), "xm_ra": coords.ra.deg, "xm_dec": coords.dec.deg}
        )
        xmatch = XMatchClass()
        xmatch.TIMEOUT = opts.timeout_s
        return xmatch.query(
            cat1=upload,
            cat2=cat2,
            max_distance=radius_arcsec * u.arcsec,
            colRA1="xm_ra",
            colDec1="xm_dec",
            cache=False,  # this module caches normalized results itself
        )

    return fetch


def _fetch_simbad_tap(coords: SkyCoord, radius_arcsec: float, opts: _Options) -> Table:
    from astroquery.simbad import SimbadClass

    simbad = SimbadClass()
    # SimbadClass(timeout=...) only bounds async execution; give its HTTP session a timeout.
    _default_timeout(simbad._session, opts.timeout_s)
    simbad.add_votable_fields("otype", "alltypes")
    try:
        # astroquery sends one OR'd-cone ADQL query for <=300 centers and a TAP upload above.
        return simbad.query_region(coords, radius=radius_arcsec * u.arcsec)
    finally:
        SimbadClass.clear_cache()  # astroquery lru-caches TAP results per service object
        simbad._session.close()


def _fetch_tap(url: str, table: str, key: tuple[str, str]) -> Callable[..., Table]:
    """Synchronous pyvo TAP query with OR'd cones (pyvo rather than astroquery.gaia's TapPlus,
    whose connections cannot be given a timeout)."""

    def fetch(coords: SkyCoord, radius_arcsec: float, opts: _Options) -> Table:
        import pyvo
        from pyvo.utils.http import create_session

        cols = ", ".join(RAW_COLUMNS[key].values())
        query = f"SELECT {cols} FROM {table} WHERE {_adql_cones(coords, radius_arcsec)}"
        with _default_timeout(create_session(), opts.timeout_s) as session:
            result = pyvo.dal.TAPService(url, session=session).run_sync(query, maxrec=_TAP_MAXREC)
        out = result.to_table()
        if len(out) >= _TAP_MAXREC or getattr(result, "query_status", None) == "OVERFLOW":
            raise _PermanentError(f"{url} result truncated at {len(out)} rows")
        return out

    return fetch


def _fetch_ned_cone(coords: SkyCoord, radius_arcsec: float, opts: _Options) -> Table:
    from astroquery import cache_conf
    from astroquery.ipac.ned import Ned

    cols = list(RAW_COLUMNS[("ned", "cone")].values())
    previous, Ned.TIMEOUT = Ned.TIMEOUT, opts.timeout_s  # astroquery reads Ned.TIMEOUT
    try:
        with cache_conf.set_temp("cache_active", False):
            parts = [Ned.query_region(c, radius=radius_arcsec * u.arcsec) for c in coords]
    finally:
        Ned.TIMEOUT = previous
    parts = [p[cols] for p in parts if len(p)]
    if not parts:
        return Table(names=cols, dtype=(str, str, float, float, float))
    return parts[0] if len(parts) == 1 else vstack(parts, metadata_conflicts="silent")


_FETCHERS: dict[tuple[str, str], Callable[[SkyCoord, float, _Options], Table]] = {
    ("simbad", "xmatch"): _fetch_xmatch("simbad"),
    ("simbad", "tap"): _fetch_simbad_tap,
    ("ned", "tap"): _fetch_tap(NED_TAP_URL, "NEDTAP.objdir", ("ned", "tap")),
    ("ned", "cone"): _fetch_ned_cone,
    ("gaia", "xmatch"): _fetch_xmatch("vizier:I/355/gaiadr3"),
    ("gaia", "tap"): _fetch_tap(GAIA_TAP_URL, "gaiadr3.gaia_source", ("gaia", "tap")),
}


def _default_timeout(session: Any, timeout_s: float) -> Any:
    """Give a ``requests.Session`` a default timeout (pyvo/astroquery TAP calls set none)."""
    request = session.request

    def request_with_timeout(method: str, url: str, **kwargs: Any) -> Any:
        kwargs.setdefault("timeout", timeout_s)
        return request(method, url, **kwargs)

    session.request = request_with_timeout
    return session


def _chunks(coords: SkyCoord, size: int) -> list[SkyCoord]:
    return [coords[k : k + size] for k in range(0, len(coords), max(size, 1))]


def _adql_cones(coords: SkyCoord, radius_arcsec: float) -> str:
    """OR'd ADQL ``CONTAINS`` cones (one per position) on columns ``ra``/``dec``."""
    r = radius_arcsec / 3600.0
    return " OR ".join(
        f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {a:.9f}, {d:.9f}, {r:.9f})) = 1"
        for a, d in zip(coords.ra.deg, coords.dec.deg, strict=True)
    )


# --- Normalization, classification, association -------------------------------------------


def _normalize(service: str, backend: str, raw: Table) -> Table:
    """Map raw service rows to the object columns (``_STR_COLS`` + ``_FLOAT_COLS``)."""
    mapping = RAW_COLUMNS[(service, backend)]
    missing = [c for c in mapping.values() if c not in raw.colnames]
    if missing:
        raise KeyError(f"{service}/{backend} response lacks columns {missing}")
    n = len(raw)
    objs = Table()
    for col in _STR_COLS:
        objs[col] = _str_values(raw[mapping[col]]) if col in mapping else np.full(n, "", str)
    for col in _FLOAT_COLS:
        objs[col] = _float_values(raw[mapping[col]]) if col in mapping else np.full(n, np.nan)
    prefix = _ID_PREFIX.get((service, backend))
    if prefix:
        objs["match_id"] = np.array([prefix + v for v in objs["match_id"]], dtype=str)
    if n:  # XMatch and overlapping cones repeat objects
        objs = unique(objs, keys=["match_id", "match_ra", "match_dec"])
    return objs


def _classify(service: str, objs: Table) -> Table:
    """Add ``is_star``, ``is_lens_related``, ``lens_types`` and the Gaia S/N columns."""
    objs = objs.copy()
    n = len(objs)
    plx_snr, pm_snr = gaia_astrometric_snr(
        *(np.asarray(objs[c]) for c in ("parallax", "parallax_error", "pmra", "pmra_error")),
        *(np.asarray(objs[c]) for c in ("pmdec", "pmdec_error", "pmra_pmdec_corr")),
    )
    is_star = np.zeros(n, dtype=bool)
    lens: list[str] = [""] * n
    if service == "simbad":
        is_star = np.isin(objs["match_type"], list(SIMBAD_STAR_OTYPES))
        for k, (main, other) in enumerate(zip(objs["match_type"], objs["all_types"], strict=True)):
            codes = {main, *str(other).split("|")}
            lens[k] = "|".join(sorted(codes & SIMBAD_LENS_OTYPES))
    elif service == "ned":
        lens = [t if t in NED_LENS_TYPES else "" for t in objs["match_type"]]
    elif service == "gaia":
        with np.errstate(invalid="ignore"):
            is_star = (plx_snr >= GAIA_MIN_PARALLAX_SNR) | (pm_snr >= GAIA_MIN_PM_SNR)
        objs["match_type"] = np.where(is_star, "astrometric_star", "unclassified")
    objs["is_star"] = is_star
    objs["lens_types"] = np.array(lens, dtype=str)
    objs["is_lens_related"] = np.array([bool(s) for s in lens], dtype=bool)
    objs["parallax_over_error"] = plx_snr
    objs["pm_snr"] = pm_snr
    return objs


def _associate(
    targets: Table, coords: SkyCoord, objs: Table, service: str, radius_arcsec: float
) -> Table:
    """Pair targets with objects within the radius (one row per pair)."""
    good = np.isfinite(objs["match_ra"]) & np.isfinite(objs["match_dec"])
    objs = objs[good]
    if len(objs) == 0:
        return _empty_matches(targets)
    obj_coords = SkyCoord(objs["match_ra"], objs["match_dec"], unit="deg")
    i_t, i_o, sep, _ = search_around_sky(coords, obj_coords, radius_arcsec * u.arcsec)
    order = np.lexsort((sep.arcsec, i_t))
    i_t, i_o, sep = i_t[order], i_o[order], sep[order]
    pairs = objs[i_o]
    pairs["source_uid"] = targets["source_uid"][i_t]
    pairs["service"] = np.full(len(pairs), service, dtype="<U6")
    pairs["sep_arcsec"] = sep.arcsec
    return pairs[list(MATCH_COLUMNS)]


def _empty_matches(targets: Table) -> Table:
    t = Table()
    t["source_uid"] = targets["source_uid"][:0]
    for col in MATCH_COLUMNS[1:]:
        if col in ("is_star", "is_lens_related"):
            t[col] = np.zeros(0, dtype=bool)
        elif col in (*_STR_COLS, "service", "lens_types"):
            t[col] = np.zeros(0, dtype="<U1")
        else:
            t[col] = np.zeros(0, dtype=float)
    return t


def _str_values(col: Any) -> np.ndarray:
    mask = np.ma.getmaskarray(col)
    values = [
        "" if m else (v.decode() if isinstance(v, bytes) else str(v)).strip()
        for v, m in zip(np.asarray(col), mask, strict=True)
    ]
    return np.array(values, dtype=str) if values else np.array([], dtype=str)


def _float_values(col: Any) -> np.ndarray:
    return np.ma.filled(np.ma.asarray(col).astype(float), np.nan)


# --- Input checks, cache, metadata --------------------------------------------------------


def _check_targets(targets: Table) -> None:
    missing = [c for c in schema.TARGET_COLUMNS if c not in targets.colnames]
    if missing:
        raise ValueError(f"crossmatch targets: missing required columns {missing}")
    ra, dec = _deg(targets["ra"]), _deg(targets["dec"])
    bad = ~(np.isfinite(ra) & np.isfinite(dec))
    if bad.any():
        raise ValueError(f"crossmatch targets: non-finite ra/dec in {int(bad.sum())} rows")
    if len(set(targets["source_uid"])) != len(targets):
        raise ValueError("crossmatch targets: source_uid values must be unique")


def _check_services(services: Sequence[str]) -> list[str]:
    if isinstance(services, str):
        services = [services]
    out = list(dict.fromkeys(s.lower() for s in services))
    unknown = [s for s in out if s not in SERVICES]
    if unknown or not out:
        raise ValueError(f"services must be a non-empty subset of {SERVICES}, got {services}")
    return out


def _resolve_backends(
    services: Sequence[str], backends: Mapping[str, str | Sequence[str]] | None
) -> dict[str, tuple[str, ...]]:
    chains = {s: DEFAULT_BACKENDS[s] for s in services}
    for svc, chain in (backends or {}).items():
        svc = svc.lower()
        chain = (chain,) if isinstance(chain, str) else tuple(chain)
        chain = tuple(b.lower() for b in chain)
        if (
            svc not in DEFAULT_BACKENDS
            or not chain
            or any((svc, b) not in _FETCHERS for b in chain)
        ):
            valid = sorted(f"{s}:{b}" for s, b in _FETCHERS)
            raise ValueError(f"unknown backend {svc}:{chain}; valid: {valid}")
        if svc in chains:
            chains[svc] = chain
    return chains


def _deg(col: Any) -> np.ndarray:
    unit = getattr(col, "unit", None)
    values = np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)
    return (values * unit).to_value(u.deg) if unit is not None else values


def _cache_path(
    cache_dir: Path, service: str, backend: str, ra: np.ndarray, dec: np.ndarray, radius: float
) -> Path:
    h = hashlib.sha256(
        json.dumps({"v": _CACHE_FORMAT, "s": service, "b": backend, "r": round(radius, 6)}).encode()
    )
    h.update(np.ascontiguousarray(np.round(np.c_[ra, dec], 7), dtype="<f8").tobytes())
    return cache_dir / f"{service}_{backend}_{h.hexdigest()[:20]}.ecsv"


def _cache_read(path: Path, max_age_days: float | None) -> Table | None:
    """Return a cached object table, or None if missing, expired or unreadable."""
    if not path.exists():
        return None
    try:
        objs = Table.read(path, format="ascii.ecsv")
        age = datetime.now(UTC) - datetime.fromisoformat(objs.meta["query_utc"])
        age_days = age.total_seconds() / 86400
        if max_age_days is not None and age_days > max_age_days:
            return None
        out = Table()
        for col in _STR_COLS:  # ECSV may read empty strings back as masked
            out[col] = _str_values(objs[col])
        for col in _FLOAT_COLS:
            out[col] = _float_values(objs[col])
        out.meta.update(objs.meta)
        return out
    except Exception as exc:  # noqa: BLE001 -- a corrupt or outdated entry is just a miss
        log.warning("ignoring unreadable cache entry %s: %s", path, exc)
        return None


def _cache_write(
    path: Path, objs: Table, service: str, backend: str, query_utc: str, radius: float
) -> None:
    out = objs.copy()
    out.meta.clear()
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=CATALOGS[(service, backend)],
        service=service,
        backend=backend,
        query_utc=query_utc,
        radius_arcsec=radius,
    )
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        out.write(tmp, format="ascii.ecsv", overwrite=True)
        os.replace(tmp, path)
    except OSError as exc:
        log.warning("could not write cache entry %s: %s", path, exc)


def _software_versions() -> dict[str, str]:
    out = {}
    for pkg in ("astroquery", "pyvo", "astropy"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = "unknown"
    return out


def _meta(
    targets: Table,
    services: Sequence[str],
    ok: Sequence[str],
    failed: Sequence[str],
    results: Mapping[str, _ServiceResult],
    radius_arcsec: float,
    *,
    cache: bool,
    refresh: bool,
) -> dict[str, Any]:
    used = {s: results[s].backend for s in ok}
    queried = ", ".join(f"{s} [{CATALOGS[(s, b)]}]" for s, b in used.items()) or "none"
    upstream = targets.meta.get("source") or "unspecified targets"
    return {
        "provenance": schema.Provenance.OBSERVED.value,
        "source": (
            f"external catalogs {queried}; cone radius {radius_arcsec:g} arcsec around "
            f"{len(targets)} targets from: {upstream}"
        ),
        "radius_arcsec": radius_arcsec,
        "services_requested": list(services),
        "services_ok": list(ok),
        "services_failed": list(failed),
        "service_errors": {s: " | ".join(r.errors) for s, r in results.items() if r.errors},
        "backend": used,
        "catalog": {s: CATALOGS[(s, b)] for s, b in used.items()},
        "query_utc": {s: results[s].query_utc for s in ok},
        "elapsed_s": {s: round(r.elapsed_s, 3) for s, r in results.items()},
        "from_cache": {s: results[s].from_cache for s in ok},
        "cache": {"enabled": cache, "refresh": refresh},
        "software": _software_versions(),
        "star_criteria": STAR_CRITERIA,
        "lens_criteria": LENS_CRITERIA,
        "best_match_rule": BEST_MATCH_RULE,
        "caveats": (
            "Gaia DR3 positions are epoch J2016.0 (no proper-motion propagation); n_matches "
            "counts catalog entries, so one object listed by several services counts once per "
            "service; an unmatched target is only 'not in these catalogs within the radius'"
        ),
    }

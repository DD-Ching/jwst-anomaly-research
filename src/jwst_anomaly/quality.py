"""Source-level quality gate applied before anomaly ranking (decision D-011).

The first real-data run showed top-ranked lists dominated by detections on low-weight or
edge regions and by sources sharper than the PSF (hot pixels, cosmic-ray residuals). This
module flags them for the whole catalog so ranking only sees sources whose measurements can
be trusted. Flags are heuristics with ASSUMPTION thresholds (documented below and in D-011);
they decide what is ranked, not what is real.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from jwst_anomaly import schema
from jwst_anomaly.cutouts import WeightMap, _as_degrees
from jwst_anomaly.features import (
    DEFAULT_MIN_SNR,
    _detected,
    _snr_ok,
    column_as_float,
    discover_bands,
)

# ASSUMPTIONS (D-011). Relative weight: median positive WHT of the image = 1.
DEFAULT_MIN_REL_WEIGHT = 0.5
# Distance to the nearest zero-weight region or image border.
DEFAULT_MIN_EDGE_ARCSEC = 1.0
# A point source has CI_50_30 = aper50/aper30 ~ 0.5/0.3 ~ 1.67 (the apertures enclose 50% and
# 30% of the PSF). Real sources are at least that extended; on program 2736 F200W the 1st/5th
# percentiles are 1.50/1.58. Well below means sharper than the PSF: a detector artifact.
DEFAULT_MAX_ARTIFACT_CI = 1.45
# More uncovered sources than this means the weight map is for the wrong image.
MAX_UNCOVERED_FRACTION = 0.9

REASONS_D011 = ("no_coverage", "low_weight", "edge", "sharper_than_psf")  # image tests
REASONS_D014 = ("low_snr", "single_band")  # detection-confirmation tests
REASONS = REASONS_D011 + REASONS_D014


def assess_sources(
    sources: Table,
    weight_map: WeightMap | None,
    *,
    ref_band: str | None = None,
    min_rel_weight: float = DEFAULT_MIN_REL_WEIGHT,
    min_edge_arcsec: float = DEFAULT_MIN_EDGE_ARCSEC,
    max_artifact_ci: float = DEFAULT_MAX_ARTIFACT_CI,
    min_snr: float = DEFAULT_MIN_SNR,
    min_detection_snr: float | None = None,
    require_multiband: bool = False,
    confirm_column: str | None = None,
) -> Table:
    """Flag each merged source (``schema.SOURCE_COLUMNS``); return ``schema.QUALITY_COLUMNS``.

    ``weight_map`` is the reference band's coarse WHT (``cutouts.sample_weight_map``); without
    it only the PSF-sharpness test runs and the weight/edge columns are NaN.
    ``quality_reason`` is a comma-separated subset of :data:`REASONS` ("" when ok):
    ``no_coverage`` is off the image or on a zero-weight cell. The sharpness test only uses
    reference-band detections with aper50 S/N >= ``min_snr`` and a positive CI, so noise is
    not called an artifact. Raises ``ValueError`` if the map's band differs from the
    reference band or the map covers almost none of the sources (wrong image).

    Detection-confirmation tests (D-014, off by default): ``low_snr`` flags sources whose best
    aper50 S/N over the bands they are detected in is below ``min_detection_snr``. The best band,
    not the reference band, counts, so red dropouts undetected in the reference band survive. With
    ``require_multiband``, ``single_band`` flags sources detected in one band only unless an
    independent detection confirms them: a finite value in ``confirm_column`` (e.g. the matched-
    photometry separation of a catalog built from a stacked detection image).
    """
    schema.validate(sources, schema.SOURCE_COLUMNS, name="sources")
    n = len(sources)
    band = (ref_band or _ref_band(sources)).upper()
    ra = _as_degrees(sources["ra"])
    dec = _as_degrees(sources["dec"])

    if weight_map is not None:
        if weight_map.band.upper() != band:
            raise ValueError(f"weight map is {weight_map.band}, sources' ref band is {band}")
        rel, edge, covered = weight_map.at(ra, dec)
        if n and np.mean(~covered) > MAX_UNCOVERED_FRACTION:
            raise ValueError(
                f"weight map {weight_map.uri} covers only {int(covered.sum())} of {n} sources"
            )
        uncovered = ~covered
        low = covered & ~(rel >= min_rel_weight)
        near_edge = covered & ~low & (edge < min_edge_arcsec)
        rel = np.where(covered, rel, 0.0)
    else:
        rel = np.full(n, np.nan)
        edge = np.full(n, np.nan)
        uncovered = low = near_edge = np.zeros(n, dtype=bool)

    ci_name = schema.band_column(band, "CI_50_30")
    if ci_name in sources.colnames:
        ci = column_as_float(sources, ci_name)
        snr_ok = _snr_ok(sources, band.lower(), "aper50", min_snr)
        if snr_ok is None:  # no error column: cannot tell noise from artifacts
            snr_ok = np.zeros(n, dtype=bool)
        with np.errstate(invalid="ignore"):
            sharp = _detected(sources, band.lower()) & snr_ok & (ci > 0) & (ci < max_artifact_ci)
    else:
        sharp = np.zeros(n, dtype=bool)

    low_snr = np.zeros(n, dtype=bool)
    snr_test = "off"
    if min_detection_snr is not None:
        passed = np.zeros(n, dtype=bool)
        unmeasured: list[str] = []
        for b in discover_bands(sources):
            ok = _snr_ok(sources, b, "aper50", min_detection_snr)
            if ok is None:  # S/N unknown in this band: a detection there is not penalised
                unmeasured.append(b)
                ok = np.ones(n, dtype=bool)
            passed |= ok & _detected(sources, b)
        if len(unmeasured) == len(discover_bands(sources)):
            snr_test = "skipped: no <band>_aper50_abmag_err columns"
        else:
            low_snr = ~passed
            snr_test = "on" + (f" (no S/N for {unmeasured}: detections pass)" if unmeasured else "")
    single = np.zeros(n, dtype=bool)
    confirmation = None
    if require_multiband:
        n_bands = column_as_float(sources, "n_bands")
        confirmed = np.zeros(n, dtype=bool)
        if confirm_column is not None:
            if confirm_column not in sources.colnames:
                raise ValueError(f"confirm_column {confirm_column!r} not in sources")
            confirmed = np.isfinite(column_as_float(sources, confirm_column))
            info = sources.meta.get("matched_photometry") or {}
            confirmation = f"{confirm_column} from {info.get('catalog', 'an external catalog')}"
        single = (n_bands <= 1) & ~confirmed

    reasons = [
        ",".join(name for name, on in zip(REASONS, flags, strict=True) if on)
        for flags in zip(uncovered, low, near_edge, sharp, low_snr, single, strict=True)
    ]
    out = Table(
        {
            "source_uid": sources["source_uid"],
            "rel_weight": rel,
            "edge_dist_arcsec": edge,
            "sharper_than_psf": sharp,
            "quality_ok": np.array([r == "" for r in reasons], dtype=bool),
            "quality_reason": reasons,
        }
    )
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=(
            f"merged sources ({n} rows); "
            + (
                f"WHT map of {weight_map.uri} ({weight_map.step}-px cells)"
                if weight_map
                else "no WHT"
            )
            + (f"; single-band confirmation: {confirmation}" if confirmation else "")
        ),
        thresholds={
            "min_rel_weight": min_rel_weight,
            "min_edge_arcsec": min_edge_arcsec,
            "max_artifact_ci": max_artifact_ci,
            "min_snr": min_snr,
            "min_detection_snr": min_detection_snr,
            "snr_test": snr_test,
            "require_multiband": require_multiband,
            "confirm_column": confirm_column,
            "provenance": schema.Provenance.ASSUMPTION.value,
        },
        ref_band=band,
    )
    return schema.validate(out, schema.QUALITY_COLUMNS, name="quality")


def summarize(quality: Table) -> dict[str, int]:
    """Counts per reason plus totals, for reports."""
    counts = {"n_sources": len(quality), "n_ok": int(np.sum(quality["quality_ok"]))}
    for reasons in quality["quality_reason"]:
        for reason in filter(None, str(reasons).split(",")):
            counts[reason] = counts.get(reason, 0) + 1
    return counts


def _ref_band(sources: Table) -> str:
    bands = {str(b) for b in np.unique(np.asarray(sources["ref_band"]).astype(str))}
    if len(bands) != 1:
        raise ValueError(f"sources have no single ref_band: {sorted(bands)}")
    return bands.pop()

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
from jwst_anomaly.cutouts import WeightMap

# ASSUMPTIONS (D-011). Relative weight: median positive WHT of the image = 1.
DEFAULT_MIN_REL_WEIGHT = 0.5
# Distance to the nearest zero-weight region or image border.
DEFAULT_MIN_EDGE_ARCSEC = 1.0
# A point source has CI_50_30 = aper50/aper30 ~ 0.5/0.3 ~ 1.67 (the apertures enclose 50% and
# 30% of the PSF). Real sources are at least that extended; on program 2736 F200W the 1st/5th
# percentiles are 1.50/1.58. Well below means sharper than the PSF: a detector artifact.
DEFAULT_MAX_ARTIFACT_CI = 1.45


def assess_sources(
    sources: Table,
    weight_map: WeightMap | None,
    *,
    ref_band: str | None = None,
    min_rel_weight: float = DEFAULT_MIN_REL_WEIGHT,
    min_edge_arcsec: float = DEFAULT_MIN_EDGE_ARCSEC,
    max_artifact_ci: float = DEFAULT_MAX_ARTIFACT_CI,
) -> Table:
    """Flag each merged source (``schema.SOURCE_COLUMNS``); return ``schema.QUALITY_COLUMNS``.

    ``weight_map`` is the reference band's coarse WHT (``cutouts.sample_weight_map``); without
    it only the PSF-sharpness test runs and weight/edge columns are NaN. ``quality_reason`` is
    a comma-separated subset of ``no_coverage, low_weight, edge, sharper_than_psf`` ("" when
    ok); ``no_coverage`` means off the image or on a zero-weight cell.
    """
    schema.validate(sources, schema.SOURCE_COLUMNS, name="sources")
    n = len(sources)
    band = (ref_band or _ref_band(sources)).lower()
    ra = np.asarray(sources["ra"], dtype=float)
    dec = np.asarray(sources["dec"], dtype=float)

    if weight_map is not None:
        rel, edge = weight_map.at(ra, dec)
        uncovered = (rel == 0) & (edge == 0)
        low = ~uncovered & ~(rel >= min_rel_weight)
        near_edge = ~uncovered & ~low & (edge < min_edge_arcsec)
    else:
        rel = np.full(n, np.nan)
        edge = np.full(n, np.nan)
        uncovered = low = near_edge = np.zeros(n, dtype=bool)

    ci_col = schema.band_column(band, "CI_50_30")
    det_col = schema.band_column(band, "detected")
    if ci_col in sources.colnames:
        ci = np.asarray(np.ma.filled(sources[ci_col], np.nan), dtype=float)
        detected = (
            np.asarray(sources[det_col], dtype=bool)
            if det_col in sources.colnames
            else np.isfinite(ci)
        )
        with np.errstate(invalid="ignore"):
            sharp = detected & np.isfinite(ci) & (ci < max_artifact_ci)
    else:
        sharp = np.zeros(n, dtype=bool)

    reasons = []
    names = ("no_coverage", "low_weight", "edge", "sharper_than_psf")
    for flags in zip(uncovered, low, near_edge, sharp, strict=True):
        reasons.append(",".join(name for name, on in zip(names, flags, strict=True) if on))
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
            + (f"WHT map of {weight_map.uri} (step {weight_map.step})" if weight_map else "no WHT")
        ),
        thresholds={
            "min_rel_weight": min_rel_weight,
            "min_edge_arcsec": min_edge_arcsec,
            "max_artifact_ci": max_artifact_ci,
            "provenance": schema.Provenance.ASSUMPTION.value,
        },
        ref_band=band.upper(),
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

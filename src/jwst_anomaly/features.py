"""Derived per-source features for anomaly ranking. Owner: bootstrap unit 3.

Design and rejected alternatives: DECISIONS.md D-004. In short:

* colors come from the ``aper50`` AB magnitudes of wavelength-adjacent bands. The pipeline sets each
  filter's aperture radii from its own encircled-energy (EE) curve, so equal-EE colors are
  PSF-independent for unresolved sources (for extended sources they are biased red, because redder
  apertures are larger);
* morphology comes from one table-level reference band, so pixel-unit features share a pixel scale;
* measurements below an S/N floor (``min_snr``) are treated as undefined: on real data the
  heavy tails of noise-dominated pipeline statistics otherwise fill the top of the ranking;
* band-level missingness is encoded explicitly (``blue_dropout``, ``red_dropout``, ``n_gaps``);
  everything else that is undefined stays NaN here. Imputation happens in :mod:`jwst_anomaly.rank`.
"""

from __future__ import annotations

import re
import warnings
from collections import Counter
from collections.abc import Callable, Sequence

import numpy as np
from astropy.table import Table

from . import schema

DEFAULT_APERTURE = "aper50"
# JWST imaging filter names: digits are the pivot wavelength in units of 10 nm (F200W = 2.0 um).
_BAND_RE = re.compile(r"^(f\d{3,4}(?:w2|w|m|n|c))_")

DEFAULT_MIN_SNR = 3.0
# 2.5 / ln(10): S/N = this / AB-magnitude error, to first order.
_MAG_ERR_TO_SNR = 2.5 / np.log(10.0)

NAN_POLICY = (
    "build_features never imputes. A feature is NaN when it is undefined for a source: band not "
    "detected, negative flux (NaN AB magnitude), a non-positive argument to log10, or the band's "
    "aperture S/N (from the AB magnitude error) below min_snr, so that pure-noise measurements are "
    "not ranked as anomalies (bands without an error column are listed in "
    "meta['min_snr_skipped_bands']). Band-level missingness is encoded explicitly by "
    "blue_dropout, red_dropout and n_gaps, which ignore min_snr. Features that are NaN for every "
    "source are dropped and listed in meta['dropped_features']. rank.score_anomalies imputes the "
    "median (robust z = 0) and reports the NaN count per source as n_missing."
)


def band_wavelength_key(band: str) -> tuple[int, str]:
    """Sort key ordering filter names by wavelength (``F090W`` < ``F150W`` < ``F1000W``)."""
    m = re.match(r"f(\d+)", band.lower())
    return (int(m.group(1)) if m else 10**9, band.lower())


def discover_bands(table: Table, aperture: str = DEFAULT_APERTURE) -> list[str]:
    """Lower-case bands that have a ``<band>_detected`` or ``<band>_<aperture>_abmag`` column.

    Returned in wavelength order.
    """
    candidates = {m.group(1) for c in table.colnames if (m := _BAND_RE.match(c.lower()))}
    present = [
        b
        for b in candidates
        if schema.band_column(b, "detected") in table.colnames
        or schema.band_column(b, f"{aperture}_abmag") in table.colnames
    ]
    return sorted(present, key=band_wavelength_key)


def column_as_float(table: Table, name: str) -> np.ndarray:
    """Column values as a new float array; masked/non-finite entries become NaN, units dropped."""
    col = table[name]
    data = np.ma.asarray(getattr(col, "value", col), dtype=float)
    out = np.array(np.ma.filled(data, np.nan), dtype=float, copy=True)  # never a view of `table`
    out[~np.isfinite(out)] = np.nan
    return out


def _detected(table: Table, band: str) -> np.ndarray:
    name = schema.band_column(band, "detected")
    if name in table.colnames:
        values = column_as_float(table, name)  # bool/int/float flag; masked or NaN -> not detected
        return np.isfinite(values) & (values != 0)
    # Fallback when the merge did not write a detection flag: any finite float value for that band
    # (integer columns are skipped because they may hold sentinels such as label = -1).
    detected = np.zeros(len(table), dtype=bool)
    for c in table.colnames:
        if c.lower().startswith(f"{band}_") and table[c].dtype.kind == "f":
            detected |= np.isfinite(column_as_float(table, c))
    return detected


def _band_values(table: Table, band: str, quantity: str) -> np.ndarray:
    name = schema.band_column(band, quantity)
    return column_as_float(table, name) if name in table.colnames else np.full(len(table), np.nan)


def _snr_ok(table: Table, band: str, aperture: str, min_snr: float | None) -> np.ndarray | None:
    """True where the band's aperture S/N >= min_snr; ``None`` if the error column is missing."""
    name = schema.band_column(band, f"{aperture}_abmag_err")
    if min_snr is None:
        return np.ones(len(table), dtype=bool)
    if name not in table.colnames:
        return None
    err = column_as_float(table, name)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = _MAG_ERR_TO_SNR / err
    return snr >= min_snr  # NaN error -> False


def _log10_positive(x: np.ndarray) -> np.ndarray:
    out = np.full(x.shape, np.nan)
    ok = x > 0
    out[ok] = np.log10(x[ok])
    return out


def _identity(x: np.ndarray) -> np.ndarray:
    return x


def build_features(
    sources: Table,
    *,
    bands: Sequence[str] | None = None,
    ref_band: str | None = None,
    aperture: str = DEFAULT_APERTURE,
    min_snr: float | None = DEFAULT_MIN_SNR,
) -> Table:
    """Compute numeric features from a merged source table (``schema.SOURCE_COLUMNS``).

    Returns ``schema.FEATURE_ID_COLUMNS`` plus feature columns, with
    ``meta["feature_spec"]`` mapping each feature to a one-line definition.
    ``meta["provenance"] == "derived"``.

    ``bands`` defaults to every band found in the column names; ``ref_band`` (the band that
    morphology and ``ref_mag`` come from) defaults to the most common ``sources["ref_band"]``;
    ``aperture`` selects the magnitude used for colors and ``ref_mag``; measurements from a band
    whose ``aperture`` S/N is below ``min_snr`` are NaN (``None`` disables this). The NaN policy is
    :data:`NAN_POLICY` (also stored in ``meta["nan_policy"]``).
    """
    schema.validate(sources, schema.SOURCE_COLUMNS, name="sources")
    if len(sources) == 0:
        raise ValueError("sources: empty table")
    available = discover_bands(sources, aperture)
    if isinstance(bands, str):
        bands = [bands]
    band_list = sorted(
        dict.fromkeys(b.lower() for b in (available if bands is None else bands)),
        key=band_wavelength_key,
    )
    if not band_list:
        raise ValueError("sources: no per-band columns found (expected e.g. f200w_detected)")
    if ref_band is None:
        counts = Counter(
            str(v).lower() for v in sources["ref_band"] if v is not np.ma.masked and str(v)
        )
        ref_band = sorted(counts, key=lambda b: (-counts[b], b))[0] if counts else ""
    ref = ref_band.lower()
    unknown = [b for b in dict.fromkeys([*band_list, ref]) if b not in available]
    if unknown:
        raise ValueError(f"sources: no columns for band(s) {unknown}; available: {available}")
    mag = f"{aperture}_abmag"
    ref_name = ref.upper()

    snr_ok: dict[str, np.ndarray] = {}
    snr_skipped: list[str] = []
    for b in dict.fromkeys([*band_list, ref]):
        ok = _snr_ok(sources, b, aperture, min_snr)
        if ok is None:
            snr_skipped.append(b)
            ok = np.ones(len(sources), dtype=bool)
        snr_ok[b] = ok
    if snr_skipped:
        warnings.warn(
            f"min_snr={min_snr} not applied to {snr_skipped}: no {aperture}_abmag_err column",
            stacklevel=2,
        )

    def measured(band: str, quantity: str) -> np.ndarray:
        return np.where(snr_ok[band], _band_values(sources, band, quantity), np.nan)

    snr_note = "; NaN if undefined"
    if min_snr is not None:
        snr_note += f" or {aperture} S/N < {min_snr}"
    columns: dict[str, np.ndarray] = {}
    spec: dict[str, str] = {}

    # Photometry: reference magnitude and adjacent-band colors.
    columns["ref_mag"] = measured(ref, mag)
    spec["ref_mag"] = f"{aperture} AB magnitude in {ref_name}{snr_note}"
    for blue, red in zip(band_list[:-1], band_list[1:], strict=True):
        name = f"color_{blue}_{red}"
        columns[name] = measured(blue, mag) - measured(red, mag)
        spec[name] = f"{aperture} AB magnitude {blue.upper()} minus {red.upper()}{snr_note}"

    # Morphology in the reference band (pipeline/photutils quantities).
    morph: list[tuple[str, str, Callable[[np.ndarray], np.ndarray], str]] = [
        ("ref_log_isophotal_area", "isophotal_area", _log10_positive, "log10 isophotal area, pix2"),
        ("ref_ellipticity", "ellipticity", _identity, "ellipticity, 1 - semiminor/semimajor"),
        ("ref_log_ci_50_30", "CI_50_30", _log10_positive, "log10 CI_50_30 (aper50/aper30 flux)"),
        ("ref_log_ci_70_50", "CI_70_50", _log10_positive, "log10 CI_70_50 (aper70/aper50 flux)"),
        ("ref_sharpness", "sharpness", _identity, "DAOFind sharpness statistic"),
        ("ref_roundness", "roundness", _identity, "DAOFind roundness statistic"),
        ("ref_log_nn_dist", "nn_dist", _log10_positive, "log10 nearest-neighbour distance (pix)"),
    ]
    for name, quantity, transform, text in morph:
        columns[name] = transform(measured(ref, quantity))
        spec[name] = f"{text} in {ref_name}{snr_note}"

    # Detection pattern: explicit band-level missingness indicators. Every undetected band is
    # counted in exactly one of the three (n_bands = n_total - blue - red - gaps is therefore not
    # a separate feature: it would count the same missing band twice).
    det = np.column_stack([_detected(sources, b) for b in band_list])
    n_det = det.sum(axis=1)
    nb = len(band_list)
    any_det = n_det > 0
    blue = np.where(any_det, np.argmax(det, axis=1), np.nan)
    red = np.where(any_det, np.argmax(det[:, ::-1], axis=1), np.nan)
    pattern = f"of {nb} bands ({', '.join(band_list)}); NaN if detected in none"
    columns["blue_dropout"] = blue
    spec["blue_dropout"] = f"undetected bands bluer than the bluest detection, {pattern}"
    columns["red_dropout"] = red
    spec["red_dropout"] = f"undetected bands redder than the reddest detection, {pattern}"
    columns["n_gaps"] = np.where(any_det, nb - n_det - blue - red, np.nan)
    spec["n_gaps"] = f"undetected bands between the bluest and reddest detection, {pattern}"

    dropped = [n for n, v in columns.items() if not np.isfinite(v).any()]
    out = Table()
    out["source_uid"] = sources["source_uid"]
    for name, values in columns.items():
        if name not in dropped:
            out[name] = values
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"build_features of: {sources.meta['source']}",
        feature_spec={n: spec[n] for n in columns if n not in dropped},
        nan_policy=NAN_POLICY,
        bands=band_list,
        ref_band=ref_name,
        aperture=aperture,
        min_snr=min_snr,
        min_snr_skipped_bands=snr_skipped,
        dropped_features=dropped,
    )
    return schema.validate(out, schema.FEATURE_ID_COLUMNS, name="features")

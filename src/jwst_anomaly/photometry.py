"""Matched-aperture photometry from an external catalog (decision D-013).

Vetting showed that colours from the JWST pipeline's encircled-energy apertures are biased with
source size: those apertures have different angular sizes in SW and LW bands, so extended
galaxies look red (CHANGELOG 2026-10-07, ``scripts/feature_size_bias.py``). The DAWN JWST
Archive (DJA, grizli) catalogs measure every band in the same circular apertures defined on one
detection image. This module fetches such a catalog (checksum-verified, cached), and joins its
aperture magnitudes to the run's merged sources by position, as ``<band>_<label>_abmag`` and
``_abmag_err`` columns that ``features.build_features(aperture=<label>)`` uses for colours.
The pipeline catalogs stay the source list; only the colours change.
"""

from __future__ import annotations

import hashlib
import os
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import paths, schema
from jwst_anomaly.features import discover_bands

AB_ZP_UJY = 23.9  # AB magnitude of 1 microjansky
MAG_ERR_PER_SNR = 2.5 / np.log(10)  # sigma_mag ~ 1.0857 / (S/N)
# SEP aperture flags (sep.h): measurements with these bits are discarded. APER_HASMASKED
# (0x20) is kept: grizli masks neighbouring segments on purpose (APERMASK=True).
SEP_APER_TRUNC, SEP_APER_ALLMASKED, SEP_APER_NONPOSITIVE = 0x10, 0x40, 0x80
BAD_FLAGS = SEP_APER_TRUNC | SEP_APER_ALLMASKED | SEP_APER_NONPOSITIVE
DEFAULT_MATCH_RADIUS_ARCSEC = 0.2


def fetch_catalog(url: str, sha256: str, cache_dir: Path | None = None) -> Path:
    """Download ``url`` once into the cache (atomic), verifying its sha256 every time."""
    cache = cache_dir or paths.cache_dir() / "external"
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / url.rstrip("/").rsplit("/", 1)[-1]
    if not target.exists():
        fd, tmp = tempfile.mkstemp(dir=cache, suffix=".part")
        os.close(fd)
        try:
            with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as out:  # noqa: S310
                while chunk := resp.read(1 << 20):
                    out.write(chunk)
            os.replace(tmp, target)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    if digest != sha256.lower():
        raise ValueError(f"{target}: sha256 {digest} does not match the configured {sha256}")
    return target


def load_dja_catalog(path: str | Path, bands: list[str], aperture: int = 1) -> Table:
    """DJA grizli ``*_phot.fits``: positions plus AB magnitudes in aperture ``aperture``.

    Returns ``id, ra, dec`` and, per band present, ``<band>_mag``/``_mag_err`` (µJy fluxes to
    AB; NaN for non-positive flux or bad SEP flags) and the aperture diameter in ``meta``.
    """
    t = Table.read(path)
    out = Table(
        {"id": t["id"], "ra": np.asarray(t["ra"], float), "dec": np.asarray(t["dec"], float)}
    )
    diam = t.meta.get(f"ASEC_{aperture}")
    for b in bands:
        f, e, fl = (
            f"{b}_flux_aper_{aperture}",
            f"{b}_fluxerr_aper_{aperture}",
            f"{b}_flag_aper_{aperture}",
        )
        if f not in t.colnames:
            continue
        unit = str(t[f].unit or "")
        if unit.lower() not in ("ujy", "µjy"):
            raise ValueError(f"{path}: {f} has unit {unit!r}, expected uJy")
        flux = np.asarray(np.ma.filled(t[f], np.nan), float)
        err = np.asarray(np.ma.filled(t[e], np.nan), float)
        flags = (
            np.asarray(np.ma.filled(t[fl], 0), int) if fl in t.colnames else np.zeros(len(t), int)
        )
        good = np.isfinite(flux) & (flux > 0) & ((flags & BAD_FLAGS) == 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[f"{b}_mag"] = np.where(good, AB_ZP_UJY - 2.5 * np.log10(flux), np.nan)
            out[f"{b}_mag_err"] = np.where(good, MAG_ERR_PER_SNR * err / flux, np.nan)
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=f"DJA catalog {Path(path).name}, aperture {aperture}",
        aperture_diameter_arcsec=float(diam) if diam is not None else None,
    )
    return out


def join_matched_photometry(
    sources: Table,
    catalog: Table,
    label: str,
    *,
    radius_arcsec: float = DEFAULT_MATCH_RADIUS_ARCSEC,
) -> Table:
    """Copy of ``sources`` with ``<band>_<label>_abmag``/``_abmag_err`` from ``catalog``.

    Each source takes its nearest catalog object within ``radius_arcsec`` (unmatched: NaN), and
    ``<label>_match_sep_arcsec`` records the separation. Only bands present in both tables are
    added. Provenance stays ``derived`` (a positional join of observed measurements).
    """
    schema.validate(sources, schema.SOURCE_COLUMNS, name="sources")
    if not label.isidentifier() or "_" in label:
        raise ValueError(f"label {label!r} must be a short alphanumeric identifier")
    src = SkyCoord(sources["ra"], sources["dec"], unit="deg")
    cat = SkyCoord(catalog["ra"], catalog["dec"], unit="deg")
    idx, sep, _ = src.match_to_catalog_sky(cat)
    sep_arcsec = sep.arcsec
    matched = sep_arcsec <= radius_arcsec
    out = sources.copy(copy_data=True)
    out[f"{label}_match_sep_arcsec"] = np.where(matched, sep_arcsec, np.nan)
    added = []
    for b in discover_bands(sources):
        if f"{b}_mag" not in catalog.colnames:
            continue
        for q, col in (("abmag", f"{b}_mag"), ("abmag_err", f"{b}_mag_err")):
            values = np.asarray(catalog[col], float)[idx]
            out[schema.band_column(b, f"{label}_{q}")] = np.where(matched, values, np.nan)
        added.append(b)
    if not added:
        raise ValueError(f"no common bands between sources and {catalog.meta.get('source')}")
    out.meta = dict(sources.meta)
    out.meta["source"] = (
        f"{sources.meta.get('source', '')} + matched-aperture photometry "
        f"({catalog.meta.get('source')}, {int(matched.sum())} of {len(out)} matched within "
        f'{radius_arcsec}")'
    )
    out.meta["matched_photometry"] = {
        "label": label,
        "bands": added,
        "n_matched": int(matched.sum()),
        "radius_arcsec": radius_arcsec,
        "aperture_diameter_arcsec": catalog.meta.get("aperture_diameter_arcsec"),
        "catalog": catalog.meta.get("source"),
    }
    return schema.validate(out, schema.SOURCE_COLUMNS, name="sources+photometry")

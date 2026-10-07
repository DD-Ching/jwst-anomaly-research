"""Matched-aperture photometry from an external catalog (decision D-013).

Vetting showed that colours from the JWST pipeline's encircled-energy apertures are biased with
source size: those apertures have different angular sizes in SW and LW bands, so extended
galaxies look red (CHANGELOG 2026-10-07, ``scripts/feature_size_bias.py``). The DAWN JWST
Archive (DJA, grizli) catalogs measure every band in the same circular apertures defined on one
detection image. This module fetches such a catalog (checksum-verified before it enters the
cache), and joins its aperture magnitudes one-to-one to the run's merged sources by position, as
``<band>_<label>_abmag`` and ``_abmag_err`` columns that ``features.build_features`` uses for
colours. The pipeline catalogs stay the source list; only the colours change. The pinned URL and
sha256 in the run config are the reproducibility record of the external file.
"""

from __future__ import annotations

import hashlib
import os
import tempfile
import urllib.request
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import paths, schema
from jwst_anomaly.acquire import MAX_DOWNLOAD_BYTES, sha256_file
from jwst_anomaly.catalog import _match_one_to_one
from jwst_anomaly.features import discover_bands

AB_ZP_UJY = 23.9  # AB magnitude of 1 microjansky
# SEP aperture flags (sep.h): measurements with these bits are discarded. APER_HASMASKED
# (0x20) is kept: grizli masks neighbouring segments on purpose (APERMASK=True).
SEP_APER_TRUNC, SEP_APER_ALLMASKED, SEP_APER_NONPOSITIVE = 0x10, 0x40, 0x80
BAD_FLAGS = SEP_APER_TRUNC | SEP_APER_ALLMASKED | SEP_APER_NONPOSITIVE
DEFAULT_MATCH_RADIUS_ARCSEC = 0.2
# Labels that would overwrite pipeline magnitude columns.
RESERVED_LABELS = {"aper30", "aper50", "aper70", "total", "isophotal", "aper"}


def match_sep_column(label: str) -> str:
    """Name of the separation column a join with ``label`` writes (NaN = unmatched)."""
    return f"{label}_match_sep_arcsec"


def fetch_catalog(
    url: str,
    sha256: str,
    cache_dir: Path | None = None,
    max_bytes: int = MAX_DOWNLOAD_BYTES,
) -> Path:
    """Download ``url`` into the cache once, verified against ``sha256`` before it is kept.

    The cached name carries the checksum prefix, so different files with the same name do not
    collide. A cached file that no longer matches is replaced by a fresh download. Files larger
    than ``max_bytes`` are refused (CLAUDE.md: >200 MB needs a stated reason, i.e. an explicit
    ``max_bytes`` in the config).
    """
    sha256 = sha256.lower()
    cache = cache_dir or paths.cache_dir() / "external"
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / f"{sha256[:12]}_{url.rstrip('/').rsplit('/', 1)[-1]}"
    if target.exists():
        if sha256_file(target) == sha256:
            return target
        target.unlink()  # stale or corrupt: fetch again
    fd, tmp = tempfile.mkstemp(dir=cache, suffix=".part")
    os.close(fd)
    try:
        digest = hashlib.sha256()
        size = 0
        with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as out:  # noqa: S310
            length = resp.headers.get("Content-Length") if hasattr(resp, "headers") else None
            if length is not None and int(length) > max_bytes:
                raise ValueError(f"{url}: {int(length)} bytes exceeds max_bytes {max_bytes}")
            while chunk := resp.read(1 << 20):
                size += len(chunk)
                if size > max_bytes:
                    raise ValueError(f"{url}: more than max_bytes {max_bytes}")
                digest.update(chunk)
                out.write(chunk)
        if digest.hexdigest() != sha256:
            raise ValueError(
                f"{url}: sha256 {digest.hexdigest()} does not match the configured {sha256} "
                "(truncated or changed download)"
            )
        os.replace(tmp, target)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return target


def load_dja_catalog(path: str | Path, bands: list[str], aperture: int = 1) -> Table:
    """DJA grizli ``*_phot.fits``: positions plus AB magnitudes in aperture ``aperture``.

    Returns ``id, ra, dec`` and, per requested band, ``<band>_mag``/``_mag_err`` (µJy fluxes to
    AB; NaN for non-positive flux or bad SEP flags). Raises ``ValueError`` if a requested band is
    missing or a flux column is not in µJy. Provenance ``derived`` (converted and filtered).
    """
    t = Table.read(path)
    missing = [b for b in bands if f"{b}_flux_aper_{aperture}" not in t.colnames]
    if missing:
        raise ValueError(f"{path}: no aperture-{aperture} photometry for bands {missing}")
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
            # Same convention as the JWST pipeline's abmag_err, so one S/N inversion fits both.
            out[f"{b}_mag_err"] = np.where(good, 2.5 * np.log10(1.0 + err / flux), np.nan)
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"DJA catalog {Path(path).name}, aperture {aperture} (uJy -> AB, SEP flags applied)",
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

    Pairs are one-to-one within ``radius_arcsec``, closest first (``catalog._match_one_to_one``),
    so two fragments never share one catalog object's photometry; unmatched rows get NaN and
    ``<label>_match_sep_arcsec`` records the separation. Every band of ``sources`` must exist in
    ``catalog``. Provenance stays ``derived`` (a positional join of derived magnitudes).
    """
    schema.validate(sources, schema.SOURCE_COLUMNS, name="sources")
    if not (label.isalnum() and label[0].isalpha()) or label.lower() in RESERVED_LABELS:
        raise ValueError(
            f"label {label!r} must be alphanumeric and not one of {sorted(RESERVED_LABELS)}"
        )
    bands = discover_bands(sources)
    missing = [b for b in bands if f"{b}_mag" not in catalog.colnames]
    if missing:
        raise ValueError(f"{catalog.meta.get('source')}: no photometry for sample bands {missing}")
    clash = [
        c
        for b in bands
        for c in (
            schema.band_column(b, f"{label}_abmag"),
            schema.band_column(b, f"{label}_abmag_err"),
        )
        if c in sources.colnames
    ]
    if clash:
        raise ValueError(f"label {label!r} would overwrite existing columns {clash[:3]}")
    src = SkyCoord(sources["ra"], sources["dec"], unit="deg")
    cat = SkyCoord(catalog["ra"], catalog["dec"], unit="deg")
    i_src, i_cat, sep, contested = _match_one_to_one(src, cat, radius_arcsec * u.arcsec)
    out = sources.copy(copy_data=True)
    sep_col = np.full(len(out), np.nan)
    sep_col[i_src] = sep
    out[match_sep_column(label)] = sep_col
    for b in bands:
        for q, col in (("abmag", f"{b}_mag"), ("abmag_err", f"{b}_mag_err")):
            values = np.full(len(out), np.nan)
            values[i_src] = np.asarray(catalog[col], float)[i_cat]
            out[schema.band_column(b, f"{label}_{q}")] = values
    out.meta = dict(sources.meta)
    out.meta["source"] = (
        f"{sources.meta.get('source', '')} + matched-aperture photometry "
        f"({catalog.meta.get('source')}, {len(i_src)} of {len(out)} matched one-to-one within "
        f'{radius_arcsec}")'
    )
    out.meta["matched_photometry"] = {
        "label": label,
        "bands": bands,
        "n_matched": len(i_src),
        "n_contested": int(contested),
        "radius_arcsec": radius_arcsec,
        "aperture_diameter_arcsec": catalog.meta.get("aperture_diameter_arcsec"),
        "catalog": catalog.meta.get("source"),
    }
    return schema.validate(out, schema.SOURCE_COLUMNS, name="sources+photometry")


def joined_columns(joined: Table, label: str) -> Table:
    """The join's own columns (``source_uid``, separation, ``<band>_<label>_*``) for saving."""
    keep = ["source_uid", match_sep_column(label)] + [
        c for c in joined.colnames if f"_{label}_abmag" in c
    ]
    out = Table(joined[keep], copy=True)
    out.meta = {
        "provenance": schema.Provenance.DERIVED.value,
        "source": joined.meta.get("source", ""),
        "matched_photometry": joined.meta.get("matched_photometry", {}),
    }
    return out

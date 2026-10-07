"""Ingest JWST pipeline source catalogs (``_cat.ecsv``). Owner: bootstrap unit 2.

``load_pipeline_catalog`` reads one level-3 catalog (one band); ``merge_bands`` positionally
merges several bands of the same field into one source table. Design and limitations are in
DECISIONS.md D-003. In short: each band is detected independently by the pipeline (no forced
photometry), so a missing band means "not detected/deblended there", and colors are approximate.
"""

from __future__ import annotations

import copy
import hashlib
import re
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Column, MaskedColumn, Table

from jwst_anomaly import schema

# Level-3 product name: <obs_id>_cat.ecsv with obs_id = jwPPPPP-<assoc>_<target>_<instr>_<optics>,
# e.g. jw02736-o001_t001_nircam_clear-f200w or jw02736-o002_t001_miri_f770w.
_CAT_NAME = re.compile(
    r"(?P<obs_id>jw\d{5}-[a-z]\d{3,4}_[a-z0-9]+_[a-z]+_(?P<optics>[a-z0-9-]+))_cat\.ecsv$",
    re.IGNORECASE,
)
_FILTER = re.compile(r"f\d{3,4}[wmn]2?", re.IGNORECASE)  # F200W, F150W2, F470N, F1000W
_BAND = re.compile(r"[A-Z0-9]+")  # band names become column-name and source_uid parts
_FIELD_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")  # safe in file names (cutouts, reports)
# Values for "not detected in this band"; photutils segment labels start at 1.
_FILLS = {"f": np.nan, "i": -1, "u": 0, "b": False, "U": "", "S": b""}
_INPUT_PROVENANCE = {schema.Provenance.OBSERVED, schema.Provenance.DERIVED}
_MIN_PAIRS_FOR_WARNINGS = 10


def load_pipeline_catalog(path: str | Path, *, band: str | None = None) -> Table:
    """Read one level-3 ``_cat.ecsv`` into a normalized single-band table.

    Keeps the pipeline columns, adds ``schema.BAND_CATALOG_COLUMNS`` (``ra``/``dec`` in deg),
    and records the band, pipeline/photutils versions and input file in ``meta``.
    ``meta["provenance"] == "observed"``.

    The band is parsed from the MAST file name; ``band`` is required when the name has no single
    filter (e.g. a NIRCam wide+narrow pair) and must agree with it otherwise. Added ``meta`` keys:
    ``band``, ``obs_id``, ``input_file``, ``input_sha256`` (of the bytes parsed),
    ``jwst_version``, ``photutils_version`` and ``pixel_scale_arcsec`` (derived: affine fit of sky
    offsets vs pixel centroids; pixel-unit columns such as ``nn_dist`` or ``semimajor_sigma``
    differ by ~2x between NIRCam SW and LW). The pipeline's own ``meta`` (versions, aperture
    parameters, ``abvega_offset``) is kept.
    """
    path = Path(path)
    raw = path.read_bytes()
    table = Table.read(raw.decode("utf-8"), format="ascii.ecsv")  # text, so hash == parsed bytes
    sky = table["sky_centroid"] if "sky_centroid" in table.colnames else None
    if "label" not in table.colnames or not isinstance(sky, SkyCoord):
        raise ValueError(f"{path.name}: not a JWST source catalog (needs label and sky_centroid)")

    obs_id, parsed_band = _parse_catalog_name(path.name)
    if band is not None:
        band = _normalize_band(band)
        if parsed_band and parsed_band != band:
            raise ValueError(f"{path.name}: file name says {parsed_band}, but band={band!r}")
    band = band or parsed_band
    if not band:
        raise ValueError(
            f"{path.name}: cannot infer a single filter from the file name; pass band="
        )

    icrs = sky.icrs
    for index, (name, values) in enumerate((("ra", icrs.ra.deg), ("dec", icrs.dec.deg)), start=1):
        col = Column(values, name=name, unit=u.deg, description=f"ICRS {name} of sky_centroid")
        if name in table.colnames:
            table.replace_column(name, col)
        else:
            table.add_column(col, index=table.colnames.index("label") + index)

    version = table.meta.get("version") or {}
    table.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=path.name,
        band=band,
        obs_id=obs_id,
        input_file=path.name,
        input_sha256=hashlib.sha256(raw).hexdigest(),
        jwst_version=str(version.get("jwst", "")),
        photutils_version=str(version.get("photutils", "")),
        pixel_scale_arcsec=_pixel_scale_arcsec(table, icrs),
    )
    return schema.validate(table, schema.BAND_CATALOG_COLUMNS, name=path.name)


def merge_bands(
    catalogs: Mapping[str, Table],
    ref_band: str,
    radius_arcsec: float = 0.1,
    *,
    include_unmatched: bool = True,
    field_id: str | None = None,
    columns: Sequence[str] | None = None,
) -> Table:
    """Positionally merge single-band catalogs ``{band: table}`` into one source table.

    Returns ``schema.SOURCE_COLUMNS`` plus per-band columns named with
    ``schema.band_column(band, quantity)`` and per-band match separations.
    ``meta["provenance"] == "derived"`` (``"simulated"`` if any input catalog is simulated).

    Matching is one-to-one per band: all pairs within ``radius_arcsec`` (astropy
    ``search_around_sky``) are accepted greedily in order of increasing separation. Bands are
    processed ``ref_band`` first, then by wavelength. Each row's ``ra``/``dec`` is the position of
    its anchoring detection, named in the per-row ``anchor_band`` column: ``ref_band`` when
    detected there, otherwise the first band that detected it. With ``include_unmatched=True``
    (default) detections without a counterpart start new rows, so the table is the union of all
    bands (e.g. F444W-only sources); with ``False`` only reference-band sources are kept.
    The ``ref_band`` column is the merge's reference band on every row.

    ``source_uid`` is ``<field_id>_<anchor band>_<anchor label>``, e.g.
    ``jw02736-o001_t001_nircam_f200w_42``; ``field_id`` defaults to the reference catalog's
    ``obs_id`` without its filter part. It is stable for a fixed set of input files
    (``meta["inputs"]`` records their sha256) and changes if the pipeline is re-run.

    Per band ``b``: ``b_detected``, ``b_sep_arcsec`` (to the row position; NaN if not detected),
    and every scalar pipeline column (or only ``label``/``ra``/``dec`` plus ``columns``). Bands
    without a detection hold NaN (float), -1 (signed int), 0 (unsigned int), False (bool) or ""
    (str); use ``b_detected``. ``n_bands`` counts detections; ``ref_nn_sep_arcsec`` is the
    distance to the nearest reference-band detection not in the row (small values flag
    split/blended detections). Rows with non-finite positions are dropped with a warning.
    ``meta["match_stats"]`` has per-band counts and, for matches to reference-band rows, the
    median separation and median offset (band minus reference). Warnings are raised for large
    offsets, low match fractions and mixed pipeline versions.
    """
    if not np.isfinite(radius_arcsec) or radius_arcsec <= 0:
        raise ValueError(f"radius_arcsec must be positive, got {radius_arcsec!r}")
    if isinstance(columns, str):
        raise TypeError("columns must be a sequence of column names, not a string")
    cats, n_dropped = _check_catalogs(catalogs)
    ref = _normalize_band(ref_band)
    if ref not in cats:
        raise ValueError(f"ref_band {ref_band!r} not in catalogs {sorted(cats)}")
    if columns is not None:
        unknown = sorted(set(columns) - {c for cat in cats.values() for c in cat.colnames})
        if unknown:
            raise ValueError(f"columns not found in any catalog: {unknown}")
    field = field_id or _default_field_id(cats[ref])
    if not _FIELD_ID.fullmatch(field):
        raise ValueError(f"field_id {field!r} must match {_FIELD_ID.pattern}")
    bands = sorted(cats, key=_band_sort_key)  # output column order
    order = [ref] + [b for b in bands if b != ref]  # matching order
    radius = radius_arcsec * u.arcsec

    ref_cat = cats[ref]
    n_ref = len(ref_cat)
    row_ra = np.asarray(ref_cat["ra"], dtype=float)
    row_dec = np.asarray(ref_cat["dec"], dtype=float)
    row_band = [ref] * n_ref
    row_label = [int(v) for v in ref_cat["label"]]
    assign = {ref: np.arange(n_ref)}  # band -> catalog row index per merged row (-1: none)
    seps = {ref: np.zeros(n_ref)}
    self_median = 0.0 if n_ref else float("nan")
    stats: dict[str, dict] = {  # same keys as the other bands; the reference matches itself
        ref: {
            "n_catalog": n_ref,
            "n_matched": n_ref,
            "n_matched_to_ref": n_ref,
            "n_unmatched": 0,
            "n_contested": 0,
            "median_sep_arcsec": self_median,
            "median_dra_arcsec": self_median,
            "median_ddec_arcsec": self_median,
        }
    }

    for band in order[1:]:
        cat = cats[band]
        cat_coord = SkyCoord(cat["ra"], cat["dec"], unit=u.deg)
        row_coord = SkyCoord(row_ra, row_dec, unit=u.deg)
        rows, idx, sep, n_contested = _match_one_to_one(row_coord, cat_coord, radius)
        a = np.full(len(row_ra), -1)
        a[rows] = idx
        s = np.full(len(row_ra), np.nan)
        s[rows] = sep
        unmatched = np.setdiff1d(np.arange(len(cat)), idx)
        stats[band] = {
            "n_catalog": len(cat),
            "n_matched": len(rows),
            "n_unmatched": len(unmatched),
            "n_contested": n_contested,
        } | _ref_match_stats(row_coord, cat_coord, rows, idx, sep, n_ref)
        if include_unmatched and len(unmatched):
            row_ra = np.concatenate([row_ra, np.asarray(cat["ra"], dtype=float)[unmatched]])
            row_dec = np.concatenate([row_dec, np.asarray(cat["dec"], dtype=float)[unmatched]])
            row_band += [band] * len(unmatched)
            row_label += [int(v) for v in cat["label"][unmatched]]
            a = np.concatenate([a, unmatched])
            s = np.concatenate([s, np.zeros(len(unmatched))])
        assign[band], seps[band] = a, s
        _warn_if_poor_match(band, ref, stats[band], n_ref, radius_arcsec)
    for band in order:
        stats[band]["n_no_position"] = n_dropped[band]

    n_rows = len(row_ra)
    for band in order:  # rows appended after a band was matched have no detection in it
        pad = n_rows - len(assign[band])
        assign[band] = np.concatenate([assign[band], np.full(pad, -1)])
        seps[band] = np.concatenate([seps[band], np.full(pad, np.nan)])

    out = Table()
    uids = [f"{field}_{b.lower()}_{lab}" for b, lab in zip(row_band, row_label, strict=True)]
    out["source_uid"] = Column(np.array(uids, dtype=str))
    out["ra"] = Column(row_ra, unit=u.deg, description="ICRS RA of the anchoring detection")
    out["dec"] = Column(row_dec, unit=u.deg, description="ICRS Dec of the anchoring detection")
    out["ref_band"] = Column(np.full(n_rows, ref), description="reference band of the merge")
    out["n_bands"] = Column(
        np.sum([assign[b] >= 0 for b in bands], axis=0).astype(int),
        description="number of bands with a detection",
    )
    out["anchor_band"] = Column(
        np.array(row_band, dtype=str), description="band whose detection gives ra/dec"
    )
    out["ref_nn_sep_arcsec"] = Column(
        _ref_nn_sep(SkyCoord(row_ra, row_dec, unit=u.deg), ref_cat, assign[ref]),
        unit=u.arcsec,
        description=f"separation to the nearest {ref} detection not merged into this row",
    )
    for band in bands:
        _add_band_columns(out, band, cats[band], assign[band], seps[band], columns)

    _warn_if_mixed_versions(cats)
    out.meta.update(
        provenance=_merged_provenance(cats),
        source="; ".join(str(cats[b].meta["source"]) for b in order),
        field_id=field,
        ref_band=ref,
        bands=bands,
        match={
            "method": "one-to-one, greedy by separation (astropy search_around_sky)",
            "radius_arcsec": float(radius_arcsec),
            "band_order": order,
            "include_unmatched": bool(include_unmatched),
            "offset_correction": "none (median offsets reported in match_stats)",
        },
        match_stats=stats,
        inputs={b: copy.deepcopy(dict(cats[b].meta)) for b in order},
    )
    return schema.validate(out, schema.SOURCE_COLUMNS, name="merge_bands")


def _normalize_band(band: str) -> str:
    name = str(band).strip().upper()
    if not _BAND.fullmatch(name):
        raise ValueError(f"band {band!r} must be letters and digits only, e.g. 'F200W'")
    return name


def _parse_catalog_name(name: str) -> tuple[str, str | None]:
    """``(obs_id, band)`` from a MAST catalog file name; ``""``/``None`` where not parseable."""
    match = _CAT_NAME.search(name)
    if not match:
        return "", None
    filters = [t for t in match["optics"].split("-") if _FILTER.fullmatch(t)]
    return match["obs_id"].lower(), filters[0].upper() if len(filters) == 1 else None


def _pixel_scale_arcsec(table: Table, sky: SkyCoord) -> float:
    """Mean pixel scale (arcsec/pix) of the image the catalog was measured on, or NaN."""
    if "xcentroid" not in table.colnames or "ycentroid" not in table.colnames:
        return float("nan")
    x = np.asarray(table["xcentroid"], dtype=float)
    y = np.asarray(table["ycentroid"], dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sky.ra.deg) & np.isfinite(sky.dec.deg)
    if ok.sum() < 3:
        return float("nan")
    dlon, dlat = sky[ok][0].spherical_offsets_to(sky[ok])
    design = np.column_stack([x[ok], y[ok], np.ones(ok.sum())])
    jac = []
    for target in (dlon.arcsec, dlat.arcsec):
        coef, _, rank, _ = np.linalg.lstsq(design, target, rcond=None)
        if rank < 3:
            return float("nan")
        jac.append(coef[:2])
    return float(np.sqrt(abs(np.linalg.det(np.array(jac)))))


def _check_catalogs(catalogs: Mapping[str, Table]) -> tuple[dict[str, Table], dict[str, int]]:
    """Validate inputs; return them keyed by normalized band with ``ra``/``dec`` as float deg.

    Rows with non-finite positions are dropped (with a warning) and counted per band. Inputs
    are not modified.
    """
    if not catalogs:
        raise ValueError("catalogs is empty")
    cats: dict[str, Table] = {}
    dropped: dict[str, int] = {}
    for band, cat in catalogs.items():
        key = _normalize_band(band)
        if key in cats:
            raise ValueError(f"band {band!r} given twice (band names are case-insensitive)")
        schema.validate(cat, schema.BAND_CATALOG_COLUMNS, name=f"{key} catalog")
        meta_band = cat.meta.get("band")
        if meta_band and _normalize_band(meta_band) != key:
            raise ValueError(f"catalog passed as {key} has meta['band'] == {meta_band!r}")
        cat = Table(cat, copy=False)  # QTable -> Table; new column container, shared data
        for name in ("ra", "dec"):
            cat.replace_column(name, Column(_degrees(cat[name]), unit=u.deg))
        ok = np.isfinite(cat["ra"]) & np.isfinite(cat["dec"])
        dropped[key] = int((~ok).sum())
        if dropped[key]:
            warnings.warn(
                f"{key}: dropping {dropped[key]} rows with non-finite ra/dec", stacklevel=3
            )
            cat = cat[ok]
        if len(np.unique(cat["label"])) != len(cat):
            raise ValueError(f"{key} catalog: labels are not unique")
        cats[key] = cat
    return cats, dropped


def _degrees(col: Column) -> np.ndarray:
    """Column values in degrees (unitless columns are taken as degrees; masked -> NaN)."""
    if isinstance(col, MaskedColumn):
        col = col.filled(np.nan)
    values = np.asarray(col, dtype=float)
    return values if col.unit is None else (values * col.unit).to_value(u.deg)


def _merged_provenance(cats: Mapping[str, Table]) -> str:
    labels = {str(cat.meta["provenance"]) for cat in cats.values()}
    if schema.Provenance.SIMULATED in labels:
        return schema.Provenance.SIMULATED.value  # synthetic sources must stay identifiable
    if not labels <= _INPUT_PROVENANCE:
        raise ValueError(f"input catalogs must be observed, derived or simulated, got {labels}")
    return schema.Provenance.DERIVED.value


def _band_sort_key(band: str) -> tuple[int, str]:
    """Order filters by nominal wavelength (F090W < F444W < F770W < F1000W); others last."""
    match = re.match(r"F(\d{3,4})", band)
    return (int(match[1]) if match else 10**6, band)


def _default_field_id(ref_cat: Table) -> str:
    obs_id = str(ref_cat.meta.get("obs_id") or "")
    if "_" not in obs_id:
        raise ValueError("reference catalog has no obs_id in meta; pass field_id=")
    return obs_id.rsplit("_", 1)[0]  # drop the optical-element part, e.g. "_clear-f200w"


def _match_one_to_one(
    rows: SkyCoord, cat: SkyCoord, radius: u.Quantity
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Greedy one-to-one pairs within ``radius``, closest first; ties broken by index.

    Returns matched row indices, catalog indices, separations (arcsec) and the number of
    contested matches (the row or the detection had another candidate within ``radius``).
    """
    empty = np.array([], dtype=int)
    if len(rows) == 0 or len(cat) == 0:
        return empty, empty, np.array([]), 0
    i_row, i_cat, d2d, _ = search_around_sky(rows, cat, radius)
    sep = d2d.arcsec
    taken_row = np.zeros(len(rows), dtype=bool)
    taken_cat = np.zeros(len(cat), dtype=bool)
    keep = []
    for k in np.lexsort((i_cat, i_row, sep)):  # primary key: separation
        if not (taken_row[i_row[k]] or taken_cat[i_cat[k]]):
            taken_row[i_row[k]] = taken_cat[i_cat[k]] = True
            keep.append(k)
    keep = np.array(keep, dtype=int)
    n_row = np.bincount(i_row, minlength=len(rows))[i_row[keep]]
    n_cat = np.bincount(i_cat, minlength=len(cat))[i_cat[keep]]
    return i_row[keep], i_cat[keep], sep[keep], int(np.sum((n_row > 1) | (n_cat > 1)))


def _ref_match_stats(
    rows: SkyCoord,
    cat: SkyCoord,
    i_row: np.ndarray,
    i_cat: np.ndarray,
    sep: np.ndarray,
    n_ref: int,
) -> dict:
    """Count, median separation and median offsets (band minus reference, arcsec) of the
    matches to reference-band rows (the first ``n_ref`` rows)."""
    to_ref = i_row < n_ref
    stats = {
        "n_matched_to_ref": int(to_ref.sum()),
        "median_sep_arcsec": float("nan"),
        "median_dra_arcsec": float("nan"),
        "median_ddec_arcsec": float("nan"),
    }
    if to_ref.any():
        dlon, dlat = rows[i_row[to_ref]].spherical_offsets_to(cat[i_cat[to_ref]])
        stats["median_sep_arcsec"] = float(np.median(sep[to_ref]))
        stats["median_dra_arcsec"] = float(np.median(dlon.arcsec))  # = dRA * cos(Dec)
        stats["median_ddec_arcsec"] = float(np.median(dlat.arcsec))
    return stats


def _warn_if_poor_match(band: str, ref: str, stats: dict, n_ref: int, radius: float) -> None:
    """Warn on a large median offset, or on a match fraction too low to measure one."""
    n_possible = min(stats["n_catalog"], n_ref)
    n_to_ref = stats["n_matched_to_ref"]
    offset = np.hypot(stats["median_dra_arcsec"], stats["median_ddec_arcsec"])
    if n_to_ref >= _MIN_PAIRS_FOR_WARNINGS and offset > 0.25 * radius:
        warnings.warn(
            f"{band}: median offset from {ref} is {offset:.3f} arcsec "
            f"(> 25% of the {radius} arcsec match radius); check the astrometry",
            stacklevel=3,
        )
    if n_possible >= _MIN_PAIRS_FOR_WARNINGS and n_to_ref < 0.2 * n_possible:
        warnings.warn(
            f"{band}: only {n_to_ref} of {n_possible} possible matches to {ref} within "
            f"{radius} arcsec; check the astrometry and that both catalogs cover the same field",
            stacklevel=3,
        )


def _warn_if_mixed_versions(cats: Mapping[str, Table]) -> None:
    versions = {
        b: (c.meta.get("jwst_version", ""), c.meta.get("photutils_version", ""))
        for b, c in cats.items()
    }
    if len(set(versions.values())) > 1:
        warnings.warn(
            f"bands come from different jwst/photutils versions {versions}; "
            "photometry may not be directly comparable",
            stacklevel=3,
        )


def _ref_nn_sep(rows: SkyCoord, ref_cat: Table, own: np.ndarray) -> np.ndarray:
    """Per row: separation (arcsec) to the nearest reference detection other than its own."""
    out = np.full(len(rows), np.nan)
    if len(rows) == 0 or len(ref_cat) == 0:
        return out
    ref_coord = SkyCoord(ref_cat["ra"], ref_cat["dec"], unit=u.deg)
    idx1, d1, _ = rows.match_to_catalog_sky(ref_coord)
    out[:] = d1.arcsec
    is_own = idx1 == own
    out[is_own] = np.nan
    if len(ref_cat) >= 2 and is_own.any():
        _, d2, _ = rows[is_own].match_to_catalog_sky(ref_coord, nthneighbor=2)
        out[is_own] = d2.arcsec
    return out


def _add_band_columns(
    out: Table,
    band: str,
    cat: Table,
    assign: np.ndarray,
    seps: np.ndarray,
    columns: Sequence[str] | None,
) -> None:
    detected = assign >= 0
    out[schema.band_column(band, "detected")] = Column(detected, description=f"detected in {band}")
    out[schema.band_column(band, "sep_arcsec")] = Column(
        seps, unit=u.arcsec, description=f"{band} detection to row position"
    )
    names = ["label", "ra", "dec"] + [
        c for c in (cat.colnames if columns is None else columns) if c not in ("label", "ra", "dec")
    ]
    for name in names:
        # Columns only: SkyCoord mixins (sky_centroid, sky_bbox_*) are covered by ra/dec.
        if name in cat.colnames and isinstance(cat[name], Column):
            values = _take(cat[name], assign)
            if values is not None:
                out[schema.band_column(band, name)] = values


def _take(col: Column, assign: np.ndarray) -> Column | None:
    """``col[assign]`` with the band's "not detected" fill where ``assign < 0``.

    Masked input values stay masked for non-float columns (floats use NaN). Returns None for
    multi-dimensional or unsupported dtypes.
    """
    kind = col.dtype.kind
    if col.ndim != 1 or kind not in _FILLS:
        return None
    fill = _FILLS[kind]
    data = np.asarray(col.filled(fill) if isinstance(col, MaskedColumn) else col)
    have = assign >= 0
    values = np.full(len(assign), fill, dtype=data.dtype)
    values[have] = data[assign[have]]
    attrs = {"unit": col.unit, "description": col.description}
    if isinstance(col, MaskedColumn) and kind != "f" and np.any(col.mask):
        mask = np.zeros(len(assign), dtype=bool)
        mask[have] = np.asarray(col.mask)[assign[have]]
        return MaskedColumn(values, mask=mask, **attrs)
    return Column(values, **attrs)

"""Image cutouts and quality flags for ranked targets. Owner: bootstrap unit 4.

Reads a level-3 ``_i2d.fits`` (local path, ``s3://``, ``https://`` or ``mast:`` URI) with
``astropy.io.fits`` and cuts each target with ``astropy.nddata.Cutout2D`` (decision D-005).
Remote files are opened through ``fsspec`` with a small read-ahead block, and each cutout
reads its rows as one contiguous byte range, so a 3" NIRCam cutout costs a few MB instead
of the 1.8 GB file.

Quality flags are heuristics derived from SCI and WHT pixels (i2d files have no DQ array).
They mark cutouts to treat with care; they are not a vetting verdict.
"""

from __future__ import annotations

import hashlib
import re
import time
import warnings
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse
from urllib.request import url2pathname

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D, NoOverlapError
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import distance_transform_edt

from jwst_anomaly import paths, schema

MAST_DOWNLOAD_URL = "https://mast.stsci.edu/api/v0.1/Download/file?uri="

QUALITY_FLAGS = {
    "ok": "none of the conditions below",
    "outside": "no pixel with data in the cutout box (target off the image or its footprint)",
    "edge": "box crosses the image array or the coverage boundary (no-data pixels on its border)",
    "nan_center": "no-data pixel within core_radius_arcsec of the target",
    "nan": "fraction of no-data pixels in the box exceeds nan_frac_max",
    "low_weight": "median core WHT below low_weight_frac x the image's typical WHT",
}

# The table is labeled "observed" (its files hold archive pixels); these columns are derived.
_COLUMN_PROVENANCE = {
    "path": schema.Provenance.OBSERVED.value,
    "band": schema.Provenance.OBSERVED.value,
    "frac_nan": schema.Provenance.DERIVED.value,
    "on_edge": schema.Provenance.DERIVED.value,
    "quality_flag": schema.Provenance.DERIVED.value,
    "x": schema.Provenance.DERIVED.value,
    "y": schema.Provenance.DERIVED.value,
    "wht_rel": schema.Provenance.DERIVED.value,
}
_OUTSIDE = {
    "path": "",
    "frac_nan": 1.0,
    "on_edge": True,
    "quality_flag": "outside",
    "wht_rel": np.nan,
}

_WHOLE_WHT_BYTES = 8 * 2**20  # read a WHT this small in one go rather than sampling rows

# Parent-header keywords copied into each cutout file (when present).
_PRIMARY_KEYS = ("TELESCOP", "INSTRUME", "DETECTOR", "FILTER", "PUPIL", "PROGRAM", "CAL_VER")
_SCI_KEYS = ("BUNIT", "PIXAR_SR", "PIXAR_A2", "PHOTMJSR")
_FILTER_RE = re.compile(r"^F\d{3}[A-Z]+\d*$")


def make_cutouts(
    image_uri: str,
    targets: Table,
    size_arcsec: float = 3.0,
    out_dir: Path | None = None,
    *,
    extensions: Sequence[str] | None = None,
    band: str | None = None,
    core_radius_arcsec: float = 0.2,
    nan_frac_max: float = 0.1,
    low_weight_frac: float = 0.5,
    weight_ref: float | None = None,
    n_weight_rows: int = 32,
    storage_options: dict[str, Any] | None = None,
    block_size: int = 2**16,
) -> Table:
    """Cut ``targets`` (``schema.TARGET_COLUMNS``) out of one level-3 ``_i2d.fits``.

    ``image_uri`` may be a local path or a cloud URI (``s3://stpubdata/...``); cloud
    reads should fetch only the needed bytes. Returns ``schema.CUTOUT_COLUMNS``.

    Each cutout is a square of ``size_arcsec`` (rounded to an odd number of pixels so the
    target sits on the central pixel) written as a small FITS file with the cutout WCS under
    ``out_dir/<image stem>/<uid>.fits`` (default ``paths.outputs_dir() / "cutouts"``; uids
    that are not safe file names get a hash suffix). It holds SCI plus ``extensions``
    (2-D float extensions of the image; default WHT when present). One row per target,
    including targets off the image (``quality_flag == "outside"``, empty ``path``, no file).

    Columns: ``frac_nan`` is the fraction of no-data pixels in the box (non-finite SCI,
    outside the array, or WHT <= 0); ``on_edge`` is the ``edge`` condition of
    ``QUALITY_FLAGS``; ``quality_flag`` joins the ``QUALITY_FLAGS`` tokens with ``,``.
    Extra columns: ``x``/``y`` (0-based target position in the image) and ``wht_rel``
    (median core WHT over the reference weight; NaN without a WHT extension). The table is
    labeled observed; ``meta["column_provenance"]`` marks the derived columns.

    ``band`` overrides the header (``FILTER``, or ``PUPIL`` for NIRCam pupil-wheel filters).
    ``weight_ref`` is the image's typical WHT; by default it is the median positive WHT of
    ``n_weight_rows`` evenly spaced rows (all of WHT if it is small). ``storage_options`` go
    to ``fsspec`` (``s3://`` defaults to anonymous access); ``block_size`` is its read-ahead
    size in bytes.
    """
    uri = _resolve_uri(image_uri)
    root = Path(out_dir) if out_dir is not None else paths.outputs_dir() / "cutouts"
    dest = root / _image_stem(uri)
    ra = _as_degrees(targets["ra"])
    dec = _as_degrees(targets["dec"])
    uids = [str(s) for s in targets["source_uid"]]
    file_names = _file_names(uids)

    t0 = time.perf_counter()
    rows: list[dict[str, Any]] = []
    with _open_image(uri, storage_options, block_size) as (hdul, fileobj):
        primary = hdul[0].header
        # Look HDUs up by name: iterating the HDUList would read every remote header.
        sci = _get_hdu(hdul, "SCI")
        if sci is None:
            raise ValueError(f"{image_uri}: no SCI extension")
        sci_section = _image_section(sci, sci.shape, "SCI", image_uri)
        wcs = _read_wcs(sci.header, image_uri)
        band = band or _band_from_header(primary, sci.header, image_uri)
        scale_arcsec = proj_plane_pixel_scales(wcs) * 3600.0  # (x, y) arcsec / pixel
        shape = (_odd(size_arcsec / scale_arcsec[1]), _odd(size_arcsec / scale_arcsec[0]))
        core_px = max(1.0, core_radius_arcsec / float(np.mean(scale_arcsec)))
        wht_hdu = _get_hdu(hdul, "WHT")
        wht_section = (
            _image_section(wht_hdu, sci.shape, "WHT", image_uri) if wht_hdu is not None else None
        )
        if extensions is None:
            exts = ["WHT"] if wht_section is not None else []
        else:
            exts = [e.upper() for e in extensions if e.upper() != "SCI"]
        sections = {"SCI": sci_section}
        for ext in dict.fromkeys(exts):
            hdu = _get_hdu(hdul, ext)
            if hdu is None:
                raise ValueError(f"{image_uri}: requested extension {ext} not found")
            sections[ext] = _image_section(hdu, sci.shape, ext, image_uri)
        ref = None if weight_ref is None else float(weight_ref)  # computed when first needed

        valid_pos = np.isfinite(ra) & np.isfinite(dec)
        xs = np.full(len(uids), np.nan)
        ys = np.full(len(uids), np.nan)
        if valid_pos.any():
            coords = SkyCoord(ra[valid_pos], dec[valid_pos], unit="deg", frame="icrs")
            xs[valid_pos], ys[valid_pos] = wcs.world_to_pixel(coords)

        for i, uid in enumerate(uids):
            xy = (float(xs[i]), float(ys[i]))
            row = {"source_uid": uid, "band": band, "x": xy[0], "y": xy[1]}
            cut = _cut(sci_section, xy, shape, wcs)
            if cut is None:
                rows.append(row | _OUTSIDE)
                continue
            wht = _cut(wht_section, xy, shape).data if wht_section is not None else None
            nodata = ~np.isfinite(cut.data)
            if wht is not None:
                nodata |= ~(np.isfinite(wht) & (wht > 0))
            if nodata.all():
                rows.append(row | _OUTSIDE)
                continue
            if ref is None:
                ref = _reference_weight(wht_section, n_weight_rows)
            row |= _quality(cut, nodata, wht, ref, core_px, nan_frac_max, low_weight_frac)
            data = {"SCI": cut.data}
            for ext, section in sections.items():
                if ext != "SCI":
                    data[ext] = wht if ext == "WHT" else _cut(section, xy, shape).data
            path = dest / f"{file_names[i]}.fits"
            _write_cutout(path, data, cut.wcs, primary, sci.header, row, uri, ra[i], dec[i])
            row["path"] = path.as_posix()
            rows.append(row)

        cache = getattr(fileobj, "cache", None)
        fetch = {
            "bytes": getattr(cache, "total_requested_bytes", None),
            "requests": getattr(cache, "miss_count", None),
            "seconds": round(time.perf_counter() - t0, 3),
        }
        image_shape = list(sci.shape)

    out = Table(
        rows=rows or None,
        names=(*schema.CUTOUT_COLUMNS, "x", "y", "wht_rel"),
        dtype=(str, str, str, float, bool, str, float, float, float),
    )
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=str(image_uri),
        column_provenance=dict(_COLUMN_PROVENANCE),
        size_arcsec=float(size_arcsec),
        size_pix=list(shape),
        pixel_scale_arcsec=[float(s) for s in scale_arcsec],
        image_shape=image_shape,
        weight_ref=float("nan") if ref is None else ref,
        quality_flags=dict(QUALITY_FLAGS),
        fetch=fetch,
    )
    return schema.validate(out, schema.CUTOUT_COLUMNS, name="cutouts")


@dataclass(frozen=True)
class WeightMap:
    """Coarse WHT map of one level-3 image, built by :func:`sample_weight_map` (D-011).

    Cell ``(i, j)`` summarises the block of pixels ``[i*step, (i+1)*step) x [j*step, (j+1)*step)``
    by the median WHT along one full-resolution row through it (row ``i*step + step//2``).
    It gates whole catalogs cheaply; cutouts measure exact per-pixel quality for the top k.
    """

    values: np.ndarray  # median WHT per cell; NaN where unreadable
    step: int  # cell size in pixels
    shape: tuple[int, int]  # full-resolution (ny, nx)
    wcs: WCS
    band: str
    reference_weight: float  # median of the positive cells
    pixel_scale_arcsec: float
    uri: str
    valid: np.ndarray  # cells with positive, finite weight
    # (2, ny, nx): (row, col) of each cell's nearest invalid cell; -1 or n means off-image.
    nearest_invalid: np.ndarray

    def at(self, ra_deg: Any, dec_deg: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(rel_weight, edge_dist_arcsec, covered)`` at sky positions.

        ``covered`` is False off the image and on zero-weight cells (where rel_weight and the
        distance are 0). The distance runs from the source's own (fractional) position to the
        nearest invalid cell's border, so it is not quantised to whole cells.
        """
        ra = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        dec = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        x, y = self.wcs.world_to_pixel_values(ra, dec)
        x, y = np.asarray(x, float), np.asarray(y, float)
        ny, nx = self.values.shape
        inside = (
            np.isfinite(x)
            & np.isfinite(y)
            & (x >= -0.5)
            & (y >= -0.5)
            & (x <= self.shape[1] - 0.5)
            & (y <= self.shape[0] - 0.5)
        )
        # Fractional cell coordinates: cell centres sit at j*step + (step-1)/2.
        centre = (self.step - 1) / 2
        fx = np.where(inside, (x - centre) / self.step, 0.0)
        fy = np.where(inside, (y - centre) / self.step, 0.0)
        ix = np.clip(np.rint(fx), 0, nx - 1).astype(int)
        iy = np.clip(np.rint(fy), 0, ny - 1).astype(int)
        covered = inside & self.valid[iy, ix]
        rel = np.where(covered, self.values[iy, ix] / self.reference_weight, 0.0)
        jy, jx = self.nearest_invalid[0][iy, ix], self.nearest_invalid[1][iy, ix]
        cells = np.clip(np.hypot(fx - jx, fy - jy) - 0.5, 0.0, None)
        edge = np.where(covered, cells * self.step * self.pixel_scale_arcsec, 0.0)
        return rel, edge, covered


def sample_weight_map(
    image_uri: str,
    step: int | None = None,
    *,
    grid_arcsec: float = 1.0,
    storage_options: dict[str, Any] | None = None,
    block_size: int = 2**16,
) -> WeightMap:
    """Coarse WHT map of a level-3 image from one full row per cell (D-011).

    The cell size is ``step`` pixels, or ``grid_arcsec`` converted with the image's pixel
    scale, so SW, LW and MIRI maps have the same angular resolution. Remote rows are fetched
    as concurrent byte ranges (``fs.cat_ranges``); a 1.8 GB NIRCam mosaic at 1" costs ~6 MB.
    WHT arrays up to ``_WHOLE_WHT_BYTES`` are read whole. Raises ``ValueError`` when no cell
    has positive weight (wrong product or empty image).
    """
    uri = _resolve_uri(str(image_uri))
    with _open_image(uri, storage_options, block_size) as (hdul, fileobj):
        sci = _get_hdu(hdul, "SCI")
        wht_hdu = _get_hdu(hdul, "WHT")
        if sci is None or wht_hdu is None:
            raise ValueError(f"{image_uri}: needs SCI and WHT extensions")
        shape = (int(sci.header["NAXIS2"]), int(sci.header["NAXIS1"]))
        wcs = _read_wcs(sci.header, uri)
        band = _band_from_header(hdul[0].header, sci.header, uri)
        scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial)) * 3600.0)
        if step is None:
            if not grid_arcsec > 0:
                raise ValueError("grid_arcsec must be > 0")
            step = max(1, int(round(grid_arcsec / scale)))
        if step < 1:
            raise ValueError("step must be >= 1")
        wht = _image_section(wht_hdu, shape, "WHT", uri)
        rows = list(range(step // 2, shape[0], step))
        itemsize = abs(int(wht_hdu.header.get("BITPIX", -32))) // 8
        if shape[0] * shape[1] * itemsize <= _WHOLE_WHT_BYTES:
            full = np.asarray(wht[:, :], dtype=float)
            row_data = [full[r] for r in rows]
        else:
            row_data = _read_rows(wht_hdu, wht, rows, shape[1], fileobj)
    values = np.stack([_block_median(r, step) for r in row_data])
    valid = np.isfinite(values) & (values > 0)
    if not valid.any():
        raise ValueError(f"{image_uri}: WHT has no positive weight")
    reference = float(np.median(values[valid]))
    # Nearest invalid cell for every cell; padding makes the image border count as invalid.
    _, idx = distance_transform_edt(np.pad(valid, 1, constant_values=False), return_indices=True)
    nearest = idx[:, 1:-1, 1:-1] - 1
    return WeightMap(
        values, step, shape, wcs, band, reference, scale, str(image_uri), valid, nearest
    )


def _block_median(row: np.ndarray, step: int) -> np.ndarray:
    """Median of each ``step``-wide block of a row (NaN-aware; the last block may be short)."""
    n = -(-row.size // step)
    padded = np.full(n * step, np.nan)
    padded[: row.size] = row
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN blocks -> NaN
        return np.nanmedian(padded.reshape(n, step), axis=1)


def _read_rows(
    hdu: Any, section: Any, rows: list[int], width: int, fileobj: Any
) -> list[np.ndarray]:
    """Full rows of an uncompressed image HDU: concurrent byte ranges when remote."""
    fs = getattr(fileobj, "fs", None)
    bitpix = int(hdu.header.get("BITPIX", -32))
    if fs is not None and bitpix in (-32, -64) and not isinstance(hdu, fits.CompImageHDU):
        try:
            start = int(hdu.fileinfo()["datLoc"])
            rowbytes = width * abs(bitpix) // 8
            starts = [start + r * rowbytes for r in rows]
            blobs = fs.cat_ranges(
                [fileobj.path] * len(rows), starts, [s + rowbytes for s in starts]
            )
            dtype = ">f4" if bitpix == -32 else ">f8"
            return [np.frombuffer(b, dtype=dtype).astype(float) for b in blobs]
        except Exception:  # any surprise: fall back to one sequential read per row
            pass
    return [np.asarray(section[r, :], dtype=float) for r in rows]


class _RowStripSection(fits.Section):
    """``Section`` that reads full-width rows as one byte range, then slices columns.

    ``ImageHDU.section[y0:y1, x0:x1]`` issues one read per row; over HTTP/S3 each read can
    become a request. Rows are contiguous in the file, so ``[y0:y1, :]`` is a single read.
    Subclassing ``Section`` keeps ``Cutout2D``'s lazy-section support.
    """

    def __getitem__(self, key: Any) -> np.ndarray:
        if (
            isinstance(key, tuple)
            and len(key) == 2
            and all(isinstance(k, slice) for k in key)
            and key[1] != slice(None)
        ):
            return super().__getitem__((key[0], slice(None)))[:, key[1]]
        return super().__getitem__(key)


def _get_hdu(hdul: fits.HDUList, name: str) -> Any:
    try:
        return hdul[name]
    except KeyError:
        return None


def _image_section(hdu: Any, shape: tuple[int, ...], name: str, uri: str) -> Any:
    """Lazy section of a 2-D float image extension with SCI's shape (checked from the header)."""
    if (
        not isinstance(hdu, fits.ImageHDU)
        or hdu.header.get("NAXIS") != 2
        or tuple(hdu.shape) != tuple(shape)
        or hdu.header.get("BITPIX", 0) > 0
    ):
        raise ValueError(f"{uri}: extension {name} must be a 2-D float image shaped like SCI")
    if isinstance(hdu, fits.CompImageHDU):  # tile-compressed: its own section is efficient
        return hdu.section
    return _RowStripSection(hdu)


@contextmanager
def _open_image(
    uri: str, storage_options: dict[str, Any] | None, block_size: int
) -> Iterator[tuple[fits.HDUList, Any]]:
    """Open lazily: memory-mapped for local files, ``fsspec`` read-ahead for remote ones."""
    if not _is_remote(uri):
        with fits.open(uri, memmap=True, lazy_load_hdus=True) as hdul:
            yield hdul, None
        return
    import fsspec  # optional dependency: the ``cloud`` extra

    opts = dict(storage_options or {})
    if uri.startswith("s3://"):
        opts.setdefault("anon", True)
    fs, path = fsspec.core.url_to_fs(uri, **opts)
    with fs.open(path, "rb", block_size=block_size, cache_type="readahead") as fileobj:
        with fits.open(fileobj, lazy_load_hdus=True) as hdul:
            yield hdul, fileobj


def _is_remote(uri: str) -> bool:
    scheme = urlparse(uri).scheme
    return len(scheme) > 1 and scheme != "file"  # one letter = Windows drive


def _resolve_uri(uri: str) -> str:
    """``mast:`` URI -> MAST download URL (serves byte ranges); ``file://`` -> local path."""
    uri = str(uri)
    if uri.startswith("mast:"):
        return MAST_DOWNLOAD_URL + uri
    if uri.startswith("file://"):  # fits.open would copy a file:// URL into the astropy cache
        return url2pathname(urlparse(uri).path)
    return uri


def _image_stem(uri: str) -> str:
    name = uri
    if _is_remote(uri):
        parsed = urlparse(uri)
        name = parse_qs(parsed.query).get("uri", [parsed.path])[0]
    name = name.replace("\\", "/").rstrip("/").rsplit("/", 1)[-1]
    return _safe_name(re.sub(r"\.fits(\.gz)?$", "", name, flags=re.IGNORECASE))


def _safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]", "_", text) or "_"


def _file_names(uids: list[str]) -> list[str]:
    """One file stem per uid, unique even on case-insensitive file systems.

    Uids that are not already safe names, or that collide, get a short hash of the uid;
    repeated uids get a counter.
    """
    used: set[str] = set()
    names = []
    for uid in uids:
        name = _safe_name(uid)
        if name != uid or name.lower() in used:
            name = f"{name}-{hashlib.sha256(uid.encode()).hexdigest()[:8]}"
        base, k = name, 1
        while name.lower() in used:
            k += 1
            name = f"{base}-{k}"
        used.add(name.lower())
        names.append(name)
    return names


def _odd(n: float) -> int:
    k = max(1, int(round(float(n))))
    return k if k % 2 else k + 1


def _as_degrees(col: Any) -> np.ndarray:
    """Coordinates in degrees (unitless means degrees); masked entries become NaN."""
    if isinstance(col, u.Quantity):
        return np.asarray(col.to_value(u.deg), dtype=float)
    values = np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)
    unit = getattr(col, "unit", None)
    return (values * u.Unit(unit)).to_value(u.deg) if unit is not None else values


def _read_wcs(header: fits.Header, uri: str) -> WCS:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)
    if not wcs.has_celestial or wcs.naxis != 2:
        raise ValueError(f"{uri}: SCI header has no 2-D celestial FITS WCS")
    return wcs


def _band_from_header(primary: fits.Header, sci: fits.Header, uri: str) -> str:
    """``FILTER``, unless ``PUPIL`` holds a filter (NIRCam F162M, F470N, ...)."""
    pupil = str(primary.get("PUPIL", sci.get("PUPIL", ""))).strip().upper()
    if _FILTER_RE.match(pupil):
        return pupil
    filt = str(primary.get("FILTER", sci.get("FILTER", ""))).strip().upper()
    if not filt:
        raise ValueError(f"{uri}: no FILTER keyword; pass band=")
    return filt


def _reference_weight(wht: Any | None, n_rows: int) -> float:
    """Median positive WHT over ``n_rows`` evenly spaced full rows (one read each).

    Rows sit at bin centres, so the zero-weight padding at the array's top/bottom is skipped.
    Arrays up to ``_WHOLE_WHT_BYTES`` are read whole instead (one request, exact median).
    NaN without a WHT extension or positive weights.
    """
    if wht is None:
        return float("nan")
    ny = wht.shape[0]
    n = max(1, n_rows)
    if np.prod(wht.shape) * np.dtype(wht.dtype).itemsize <= _WHOLE_WHT_BYTES:
        sample = np.asarray(wht[:, :], dtype=float).ravel()
    else:
        rows = np.unique(((np.arange(n) + 0.5) * ny / n).astype(int))
        sample = np.concatenate([np.asarray(wht[int(r), :], dtype=float).ravel() for r in rows])
    good = sample[np.isfinite(sample) & (sample > 0)]
    return float(np.median(good)) if good.size else float("nan")


def _cut(
    section: Any, xy: tuple[float, float], shape: tuple[int, int], wcs: WCS | None = None
) -> Cutout2D | None:
    if not np.all(np.isfinite(xy)):
        return None
    try:
        return Cutout2D(section, xy, shape, wcs=wcs, mode="partial", fill_value=np.nan, copy=True)
    except NoOverlapError:
        return None


def _quality(
    cut: Cutout2D,
    nodata: np.ndarray,
    wht: np.ndarray | None,
    weight_ref: float,
    core_px: float,
    nan_frac_max: float,
    low_weight_frac: float,
) -> dict[str, Any]:
    ny, nx = nodata.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    cx, cy = cut.input_position_cutout
    core = (xx - cx) ** 2 + (yy - cy) ** 2 <= core_px**2
    border = np.ones_like(nodata)
    border[1:-1, 1:-1] = False
    on_edge = bool(nodata[border].any())
    frac_nan = float(nodata.mean())
    wht_rel = np.nan
    if wht is not None and np.isfinite(weight_ref) and weight_ref > 0:
        wht_rel = float(np.median(np.where(nodata[core], 0.0, wht[core]))) / weight_ref
    tokens = [
        name
        for name, hit in (
            ("edge", on_edge),
            ("nan_center", bool(nodata[core].any())),
            ("nan", frac_nan > nan_frac_max),
            ("low_weight", bool(wht_rel < low_weight_frac)),
        )
        if hit
    ]
    return {
        "frac_nan": frac_nan,
        "on_edge": on_edge,
        "quality_flag": ",".join(tokens) or "ok",
        "wht_rel": wht_rel,
    }


def _write_cutout(
    path: Path,
    data: dict[str, np.ndarray],
    wcs: WCS,
    primary: fits.Header,
    sci_header: fits.Header,
    row: dict[str, Any],
    uri: str,
    ra: float,
    dec: float,
) -> None:
    head = fits.Header()
    for key in _PRIMARY_KEYS:
        if key in primary:
            head[key] = primary[key]
    # Variable-length strings get no comment (avoids truncation warnings) and are made ASCII.
    head["SRCUID"] = _ascii(row["source_uid"])
    head["TARG_RA"] = (float(ra), "[deg] target right ascension")
    head["TARG_DEC"] = (float(dec), "[deg] target declination")
    head["BAND"] = _ascii(row["band"])
    head["PROVENAN"] = (schema.Provenance.OBSERVED.value, "pixels copied from the parent image")
    head["QUALFLAG"] = row["quality_flag"]
    head["ORIGURI"] = _ascii(uri)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs_header = wcs.to_header()
    hdus: list[Any] = [fits.PrimaryHDU(header=head)]
    for ext, arr in data.items():
        header = wcs_header.copy()
        if ext == "SCI":
            for key in _SCI_KEYS:
                if key in sci_header:
                    header[key] = sci_header[key]
        hdus.append(fits.ImageHDU(np.asarray(arr, dtype=np.float32), header=header, name=ext))
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(path, overwrite=True)


def _ascii(text: str) -> str:
    """FITS header strings must be printable ASCII; escape anything else."""
    escaped = str(text).encode("ascii", "backslashreplace").decode("ascii")
    return "".join(c if " " <= c <= "~" else repr(c)[1:-1] for c in escaped)

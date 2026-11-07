"""Cutout images: asinh-stretched PNG panels and contact sheets. Owner: bootstrap unit 4.

Reuses ``astropy.visualization`` for the stretch and matplotlib's object-oriented API
(``Figure``, Agg) so nothing touches pyplot's global state. By default every panel gets its
own display interval (``scale="panel"``), so panel brightness is not comparable between
panels; ``scale="band"`` or ``"all"`` shares one interval within a band or across the sheet.
Panels are shown in stored pixel orientation unless ``north_up=True``, which resamples each
cutout from its own WCS (:func:`north_up`). These images are for looking at sources, not for
measuring them.
"""

from __future__ import annotations

import math
import re
import textwrap
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ImageNormalize, simple_norm
from astropy.wcs import WCS, FITSFixedWarning
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy import ndimage

NO_DATA_COLOR = "#3a5f9f"  # NaN / no-coverage pixels; distinct from black sky and white sources
FLAG_COLOR = "#c8780a"  # frame of panels whose quality_flag is not "ok" (flag text is the label)
MARK_COLOR = "#e8b04a"  # ticks marking the target position
INK = "#222222"
MUTED = "#6b6b6b"
FOOTNOTE = (
    "asinh stretch (a={a}) over {scope} {lo}-{hi} percentile range; {orient}; blue = no data; "
    "ticks mark the target; orange frame = quality flag set"
)
SCALES = ("panel", "band", "all")
# ASSUMPTION: a resampled pixel is shown only if (almost) all of its bilinear weight lands on
# valid input; the margin absorbs float round-off of an exact 0/90/180/270 degree rotation.
_MIN_COVER = 0.99
_SCOPE = {"panel": "each panel's", "band": "each band's shared", "all": "one shared"}
ORIENT_NORTH_UP = "north up, east left"
ORIENT_STORED = "stored pixel orientation (not north-up)"
# Chosen by eye on program-2736 NIRCam/MIRI cutouts: cores stay resolved while PSF wings,
# faint neighbours and sky noise remain visible (see DECISIONS.md D-005).
DEFAULT_PERCENTILES = (1.0, 99.9)
DEFAULT_ASINH_A = 0.01


def normalize(
    data: np.ndarray,
    *,
    percentiles: tuple[float, float] = DEFAULT_PERCENTILES,
    asinh_a: float = DEFAULT_ASINH_A,
) -> ImageNormalize | None:
    """Asinh ``simple_norm`` over the finite pixels of ``data``; None if none are finite.

    ``clip`` stays False so NaN pixels stay masked and render in ``NO_DATA_COLOR``.
    """
    data = np.asarray(data, dtype=float)
    if not np.isfinite(data).any():
        return None
    lo, hi = percentiles
    return simple_norm(data, "asinh", asinh_a=asinh_a, min_percent=lo, max_percent=hi)


def north_up(data: np.ndarray, wcs: WCS | fits.Header, *, order: int = 1) -> np.ndarray:
    """Resample ``data`` about its central pixel so north is up and east is left.

    The local sky Jacobian at the central pixel comes from ``wcs`` itself (one-pixel steps
    converted to east/north offsets with ``SkyCoord.spherical_offsets_to``), so rotation,
    parity flips, skew and the meridian convergence of a mosaic WCS are all handled. Output
    has the input shape and the input pixel area; the central pixel (the cutout target) stays
    central; uncovered corners are NaN.

    "North" is ICRS north whatever the WCS frame (galactic or ecliptic WCSs are converted).
    ``order=1`` (bilinear, the default) is NaN-aware: values and coverage are interpolated
    separately and any output pixel that draws on a no-data pixel is NaN, so isolated
    no-data pixels stay visible (grown by about a pixel) instead of being filled. Nearest neighbour
    (``order=0``) is not the default because at a general angle it skips some input pixels,
    which can hide a one-pixel artefact (hot pixel, cosmic ray) from the eye.
    Raises ``ValueError`` when ``wcs`` has no 2-D celestial WCS.
    """
    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError(f"north_up: expected a 2-D image, got shape {data.shape}")
    if isinstance(wcs, fits.Header):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(wcs)
    if not wcs.has_celestial or wcs.celestial.naxis != 2:
        raise ValueError("north_up: no 2-D celestial WCS")
    wcs = wcs.celestial
    ny, nx = data.shape
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    pts = wcs.pixel_to_world([cx, cx + 1.0, cx], [cy, cy, cy + 1.0])
    pts = SkyCoord(pts).icrs
    centre = pts[0]
    jac = np.empty((2, 2))  # rows: east, north [arcsec]; columns: +1 px in x, +1 px in y
    for col, step in enumerate((pts[1], pts[2])):
        east, north = centre.spherical_offsets_to(step)
        jac[:, col] = (east.to_value("arcsec"), north.to_value("arcsec"))
    det = float(np.linalg.det(jac))
    if not np.isfinite(det) or det == 0.0:
        raise ValueError("north_up: singular WCS at the cutout centre")
    s = math.sqrt(abs(det))
    # Output step (u, v) is (-s east, +s north); A maps it to the input step (x, y).
    a = np.linalg.solve(jac, np.diag([-s, s]))
    matrix = a[::-1, ::-1]  # ndimage works in (row, col) = (y, x)
    centre_rc = np.array([cy, cx])
    finite = np.isfinite(data)

    def resample(arr: np.ndarray, fill: float) -> np.ndarray:
        return ndimage.affine_transform(
            arr,
            matrix,
            offset=centre_rc - matrix @ centre_rc,
            output_shape=data.shape,
            order=order,
            mode="grid-constant",  # pixels own +-0.5 px, so edge pixels survive rounding
            cval=fill,
        )

    values = resample(np.where(finite, data, 0.0), 0.0)
    cover = resample(finite.astype(float), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = values / cover
    out[~(cover >= _MIN_COVER)] = np.nan

    return out


def cutout_png(
    path: str | Path,
    out_png: str | Path | None = None,
    *,
    title: str | None = None,
    ext: str = "SCI",
    percentiles: tuple[float, float] = DEFAULT_PERCENTILES,
    asinh_a: float = DEFAULT_ASINH_A,
    panel_inches: float = 3.0,
    dpi: int = 150,
    north_up: bool = False,
) -> Path:
    """Render one cutout FITS (extension ``ext``) as a PNG; default name is ``path`` + .png.

    ``north_up`` resamples the panel from the extension's WCS (:func:`north_up`); the label
    says which orientation is shown.
    """
    path = Path(path)
    out = Path(out_png) if out_png is not None else path.with_suffix(".png")
    header = fits.getheader(path, 0)
    data, orient = _load_panel(path, ext, north_up)
    if title is None:
        title = f"{header.get('SRCUID', path.stem)}  {header.get('BAND', '')}"
    flag = str(header.get("QUALFLAG", ""))
    fig = Figure(figsize=(panel_inches, panel_inches + 0.5), dpi=dpi, layout="constrained")
    ax = fig.subplots()
    norm = normalize(data, percentiles=percentiles, asinh_a=asinh_a)
    lines = [title, flag] + ([orient] if north_up else [])
    _draw_panel(ax, data, lines, flag, norm, _chars(panel_inches))
    return _save(fig, out)


def contact_sheet(
    cutouts: Table,
    out_png: str | Path,
    *,
    ranks: Mapping[str, int] | None = None,
    ncols: int = 5,
    title: str | None = None,
    percentiles: tuple[float, float] = DEFAULT_PERCENTILES,
    asinh_a: float = DEFAULT_ASINH_A,
    panel_inches: float = 2.0,
    dpi: int = 120,
    north_up: bool = False,
    scale: str = "panel",
) -> Path:
    """Grid of cutouts from ``make_cutouts`` output (``schema.CUTOUT_COLUMNS``).

    One band: panels wrap into ``ncols`` columns. Several bands: one row per source and one
    column per band (ordered by the wavelength in the band name). Sources are ordered by
    ``rank`` (a ``rank`` column, else the ``ranks`` mapping), else by first appearance.
    Panels are labeled with source_uid, rank, band and quality_flag; rows without a file
    (``outside``) are drawn as empty panels carrying their flag.

    ``north_up`` resamples every panel from its SCI WCS (:func:`north_up`); a panel whose
    file has no celestial WCS is drawn as stored and labeled so. ``scale`` picks the display
    interval: ``"panel"`` (each panel its own), ``"band"`` (one per band, so brightness is
    comparable between sources in a band) or ``"all"`` (one for the whole sheet). Shared
    percentiles are taken over the pooled finite pixels, so one bright source sets the range.
    """
    if len(cutouts) == 0:
        raise ValueError("contact_sheet: no cutouts to draw")
    if scale not in SCALES:
        raise ValueError(f"contact_sheet: scale must be one of {SCALES}, got {scale!r}")
    uids = _ordered_sources(cutouts, ranks)
    rank_of = _rank_lookup(cutouts, ranks)
    bands = sorted({str(b) for b in cutouts["band"]}, key=_band_sort_key)
    cells = _pick_cells(cutouts)

    if len(bands) > 1:
        nrows, ncol = len(uids), len(bands)
        layout = [(uid, band) for uid in uids for band in bands]
    else:
        ncol = max(1, min(ncols, len(uids)))
        nrows = math.ceil(len(uids) / ncol)
        layout = [(uid, bands[0]) for uid in uids]

    fig = Figure(
        figsize=(ncol * panel_inches, nrows * (panel_inches + 0.4) + 0.7),
        dpi=dpi,
        layout="constrained",
    )
    axes = fig.subplots(nrows, ncol, squeeze=False).ravel()
    wrap = _chars(panel_inches)
    panels: dict[tuple[str, str], tuple[np.ndarray, str]] = {}
    if scale != "panel":  # pooled intervals need every panel first; "panel" reads one at a time
        for key in layout:
            row = cells.get(key)
            if row is not None and row["path"]:
                panels[key] = _load_panel(Path(str(row["path"])), "SCI", north_up)
    shared = _shared_norms(panels, scale, percentiles, asinh_a)
    for ax, (uid, band) in zip(axes, layout, strict=False):
        row = cells.get((uid, band))
        rank = rank_of.get(uid)
        head = f"{uid}" + (f"  #{rank}" if rank is not None else "")
        if row is None:
            _draw_panel(ax, None, [head, f"{band}: no cutout"], "missing", None, wrap)
            continue
        flag = str(row["quality_flag"])
        if not row["path"]:
            data, orient = None, ""
        elif (uid, band) in panels:
            data, orient = panels[(uid, band)]
        else:
            data, orient = _load_panel(Path(str(row["path"])), "SCI", north_up)
        lines = [head, f"{band}  {flag}"]
        if north_up and orient != ORIENT_NORTH_UP and data is not None:
            lines.append(orient)  # only panels that differ from the sheet's footnote
        if data is None:
            norm = None
        elif scale == "panel":
            norm = normalize(data, percentiles=percentiles, asinh_a=asinh_a)
        else:
            norm = shared.get("" if scale == "all" else band)
        _draw_panel(ax, data, lines, flag, norm, wrap)
    for ax in axes[len(layout) :]:
        ax.set_axis_off()
    if title:
        fig.suptitle(title, fontsize=10, color=INK)
    note = FOOTNOTE.format(
        a=asinh_a,
        scope=_SCOPE[scale],
        lo=percentiles[0],
        hi=percentiles[1],
        orient=ORIENT_NORTH_UP if north_up else ORIENT_STORED,
    )
    fig.supxlabel(textwrap.fill(note, _chars(ncol * panel_inches)), fontsize=7, color=MUTED)
    return _save(fig, Path(out_png))


def _load_panel(path: Path, ext: str, north: bool) -> tuple[np.ndarray, str]:
    """(data, orientation label) of one cutout extension, resampled north-up when asked."""
    data, header = fits.getdata(path, extname=ext, header=True)
    if not north:
        return data, ORIENT_STORED
    try:
        return north_up(data, header), ORIENT_NORTH_UP
    except ValueError:
        return data, "no WCS: " + ORIENT_STORED


def _shared_norms(
    panels: Mapping[tuple[str, str], tuple[np.ndarray, str]],
    scale: str,
    percentiles: tuple[float, float],
    asinh_a: float,
) -> dict[str, ImageNormalize]:
    """One norm per band (``scale="band"``) or one under key "" (``"all"``) over pooled pixels."""
    if scale == "panel":
        return {}
    groups: dict[str, list[np.ndarray]] = {}
    for (_, band), (data, _) in panels.items():
        pix = np.asarray(data, dtype=float).ravel()
        groups.setdefault("" if scale == "all" else band, []).append(pix[np.isfinite(pix)])
    out = {}
    for key, parts in groups.items():
        norm = normalize(np.concatenate(parts), percentiles=percentiles, asinh_a=asinh_a)
        if norm is not None:
            out[key] = norm
    return out


def _draw_panel(
    ax: Axes,
    data: np.ndarray | None,
    lines: Sequence[str],
    flag: str,
    norm: ImageNormalize | None,
    wrap: int,
) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    # Flags are comma-joined without spaces; add them so long flag lists can wrap.
    label = "\n".join(textwrap.fill(line.replace(",", ", "), wrap) for line in lines)
    ax.set_title(label, fontsize=7, color=INK, pad=3)
    if data is None or norm is None:
        ax.set_facecolor("#e6e6e6")
        ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize=8, color=MUTED)
    else:
        cmap = matplotlib.colormaps["gray"].with_extremes(bad=NO_DATA_COLOR)
        ax.imshow(data, origin="lower", cmap=cmap, norm=norm, interpolation="nearest")
        _mark_center(ax, data.shape)
    flagged = flag not in ("ok", "")
    for spine in ax.spines.values():
        spine.set_color(FLAG_COLOR if flagged else "#999999")
        spine.set_linewidth(2.0 if flagged else 0.6)


def _chars(inches: float) -> int:
    """Characters of 7 pt text that fit in ``inches`` (about 15 per inch), at least 20."""
    return max(20, int(inches * 15))


def _mark_center(ax: Axes, shape: tuple[int, ...]) -> None:
    """Short ticks pointing at the central pixel, leaving the source itself unobscured."""
    ny, nx = shape[:2]
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    inner, outer = 0.18 * min(nx, ny), 0.32 * min(nx, ny)
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ax.plot(
            [cx + dx * inner, cx + dx * outer],
            [cy + dy * inner, cy + dy * outer],
            color=MARK_COLOR,
            linewidth=1.0,
        )
    ax.set_xlim(-0.5, nx - 0.5)
    ax.set_ylim(-0.5, ny - 0.5)


def _pick_cells(cutouts: Table) -> dict[tuple[str, str], Any]:
    """(uid, band) -> row. With several images per band the first row with a file wins."""
    cells: dict[tuple[str, str], Any] = {}
    for row in cutouts:
        key = (str(row["source_uid"]), str(row["band"]))
        if key not in cells or (row["path"] and not cells[key]["path"]):
            cells[key] = row
    return cells


def _rank_lookup(cutouts: Table, ranks: Mapping[str, int] | None) -> dict[str, int]:
    """uid -> rank from a ``rank`` column, else ``ranks``; masked or NaN ranks are skipped."""
    if "rank" in cutouts.colnames:
        pairs = zip(cutouts["source_uid"], cutouts["rank"], strict=True)
    else:
        pairs = (ranks or {}).items()
    out = {}
    for uid, rank in pairs:
        if not np.ma.is_masked(rank) and np.isfinite(float(rank)):
            out.setdefault(str(uid), int(rank))
    return out


def _ordered_sources(cutouts: Table, ranks: Mapping[str, int] | None) -> list[str]:
    first_seen = list(dict.fromkeys(str(u) for u in cutouts["source_uid"]))
    rank_of = _rank_lookup(cutouts, ranks)
    # sorted() is stable, so unranked sources keep their first-appearance order.
    return sorted(first_seen, key=lambda uid: rank_of.get(uid, math.inf))


def _band_sort_key(band: str) -> tuple[float, str]:
    match = re.search(r"(\d+)", band)
    return (float(match.group(1)) if match else math.inf, band)


def _save(fig: Figure, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    return out

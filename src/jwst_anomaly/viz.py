"""Cutout images: asinh-stretched PNG panels and contact sheets. Owner: bootstrap unit 4.

Reuses ``astropy.visualization`` for the stretch and matplotlib's object-oriented API
(``Figure``, Agg) so nothing touches pyplot's global state. Every panel gets its own
display interval, so panel brightness is not comparable between panels; these images are
for looking at sources, not for measuring them.
"""

from __future__ import annotations

import math
import re
import textwrap
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ImageNormalize, simple_norm
from matplotlib.axes import Axes
from matplotlib.figure import Figure

NO_DATA_COLOR = "#3a5f9f"  # NaN / no-coverage pixels; distinct from black sky and white sources
FLAG_COLOR = "#c8780a"  # frame of panels whose quality_flag is not "ok" (flag text is the label)
MARK_COLOR = "#e8b04a"  # ticks marking the target position
INK = "#222222"
MUTED = "#6b6b6b"
FOOTNOTE = (
    "asinh stretch (a={a}) over each panel's {lo}-{hi} percentile range; blue = no data; "
    "ticks mark the target; orange frame = quality flag set"
)
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
) -> Path:
    """Render one cutout FITS (extension ``ext``) as a PNG; default name is ``path`` + .png."""
    path = Path(path)
    out = Path(out_png) if out_png is not None else path.with_suffix(".png")
    header = fits.getheader(path, 0)
    data = fits.getdata(path, extname=ext)
    if title is None:
        title = f"{header.get('SRCUID', path.stem)}  {header.get('BAND', '')}"
    flag = str(header.get("QUALFLAG", ""))
    fig = Figure(figsize=(panel_inches, panel_inches + 0.5), dpi=dpi, layout="constrained")
    ax = fig.subplots()
    norm = normalize(data, percentiles=percentiles, asinh_a=asinh_a)
    _draw_panel(ax, data, [title, flag], flag, norm, _chars(panel_inches))
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
) -> Path:
    """Grid of cutouts from ``make_cutouts`` output (``schema.CUTOUT_COLUMNS``).

    One band: panels wrap into ``ncols`` columns. Several bands: one row per source and one
    column per band (ordered by the wavelength in the band name). Sources are ordered by
    ``rank`` (a ``rank`` column, else the ``ranks`` mapping), else by first appearance.
    Panels are labeled with source_uid, rank, band and quality_flag; rows without a file
    (``outside``) are drawn as empty panels carrying their flag.
    """
    if len(cutouts) == 0:
        raise ValueError("contact_sheet: no cutouts to draw")
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
    for ax, (uid, band) in zip(axes, layout, strict=False):
        row = cells.get((uid, band))
        rank = rank_of.get(uid)
        head = f"{uid}" + (f"  #{rank}" if rank is not None else "")
        if row is None:
            _draw_panel(ax, None, [head, f"{band}: no cutout"], "missing", None, wrap)
            continue
        flag = str(row["quality_flag"])
        data = fits.getdata(str(row["path"]), extname="SCI") if row["path"] else None
        norm = (
            normalize(data, percentiles=percentiles, asinh_a=asinh_a) if data is not None else None
        )
        _draw_panel(ax, data, [head, f"{band}  {flag}"], flag, norm, wrap)
    for ax in axes[len(layout) :]:
        ax.set_axis_off()
    if title:
        fig.suptitle(title, fontsize=10, color=INK)
    note = FOOTNOTE.format(a=asinh_a, lo=percentiles[0], hi=percentiles[1])
    fig.supxlabel(textwrap.fill(note, _chars(ncol * panel_inches)), fontsize=7, color=MUTED)
    return _save(fig, Path(out_png))


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

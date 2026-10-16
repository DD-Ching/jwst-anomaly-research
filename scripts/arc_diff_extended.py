"""Local-annulus test for ``arc_diff.py`` hits on smooth extended light (Abell 2744, 2026-10-16).

The peer-population cut (``arc_diff_peers.py``) compares a hit's S with S at compact epoch-1
peaks of similar brightness. Whole-chain injections (CHANGELOG 2026-10-16) showed that it vetoes
clean point brightenings on smooth galaxy outskirts and ICL: the peer spread of their brightness
bin is set by compact sources. This module splits hits by their epoch-1 position:

- **compact**: epoch 1 has a local maximum within ``--peak-arcsec`` of the hit whose value is at
  least ``--min-contrast`` × (1.4826 MAD) above the median of epoch 1 in the annulus
  ``--ann-in``–``--ann-out`` (the mosaics are not background-subtracted, so a ratio to the
  annulus median fails; first version, 2026-10-16);
  the peer cut stays;
- **extended**: otherwise. ``z_ann`` = (S_hit − median S in the annulus) / (1.4826 MAD of S in
  the annulus), where S_hit is the signed extreme of S within ``--peak-arcsec`` of the hit; S
  pixels within ``--ann-in`` of other hits are masked so neighbours do not set the scale. The cut
  is |z_ann| ≥ ``--min-z`` (needs ≥ ``--min-ann-px`` finite annulus pixels; fewer means no test,
  i.e. pass).

Every number is ``derived``; the thresholds are ASSUMPTIONs (calibrated on injections, see
CHANGELOG).

Usage: ``arc_diff_extended.py --maps M.npz --hits T.ecsv [--id-col hit] --out ext.ecsv``. The hits
table needs ``ra``, ``dec`` and an id column.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy import ndimage

DEFAULTS = {
    "peak_arcsec": 0.15,
    "ann_in": 0.3,
    "ann_out": 0.9,
    "min_contrast": 5.0,
    "min_z": 3.0,
    "min_ann_px": 30,
}


def _disk(r_px: float) -> np.ndarray:
    n = int(np.ceil(r_px))
    yy, xx = np.mgrid[-n : n + 1, -n : n + 1]
    return yy**2 + xx**2 <= r_px**2


def measure(
    s: np.ndarray,
    e1: np.ndarray,
    yx: np.ndarray,
    scale_arcsec: float,
    peak_arcsec: float = DEFAULTS["peak_arcsec"],
    ann_in: float = DEFAULTS["ann_in"],
    ann_out: float = DEFAULTS["ann_out"],
    min_contrast: float = DEFAULTS["min_contrast"],
    min_ann_px: int = DEFAULTS["min_ann_px"],
) -> list[dict]:
    """Compactness and annulus z for each (y, x) pixel position in ``yx`` (shape (n, 2))."""
    rp, ri, ro = (v / scale_arcsec for v in (peak_arcsec, ann_in, ann_out))
    # Mask S around every hit (the annulus scale must not include neighbouring hits).
    hit_mask = np.zeros(s.shape, bool)
    yi = np.clip(np.rint(yx[:, 0]).astype(int), 0, s.shape[0] - 1) if len(yx) else []
    xi = np.clip(np.rint(yx[:, 1]).astype(int), 0, s.shape[1] - 1) if len(yx) else []
    hit_mask[yi, xi] = True
    if len(yx):
        hit_mask = ndimage.binary_dilation(hit_mask, structure=_disk(ri))
    e1f = np.where(np.isfinite(e1), e1, -np.inf)
    lmax = e1f == ndimage.maximum_filter(e1f, size=3, mode="nearest")
    n = int(np.ceil(ro))
    yy, xx = np.mgrid[-n : n + 1, -n : n + 1]
    rr = np.hypot(yy, xx)
    ann = (rr >= ri) & (rr <= ro)
    core = rr <= rp
    out = []
    for y, x in np.rint(yx).astype(int):
        r = {"compact": False, "contrast": np.nan, "s_hit": np.nan, "z_ann": np.nan, "n_ann": 0}
        if not (0 <= y < s.shape[0] and 0 <= x < s.shape[1]):
            out.append(r)
            continue
        y0, y1, x0, x1 = y - n, y + n + 1, x - n, x + n + 1
        sl = (slice(max(y0, 0), min(y1, s.shape[0])), slice(max(x0, 0), min(x1, s.shape[1])))
        ks = (slice(sl[0].start - y0, sl[0].stop - y0), slice(sl[1].start - x0, sl[1].stop - x0))
        sw, ew, mw, lw = s[sl], e1[sl], hit_mask[sl], lmax[sl]
        a_, c_ = ann[ks], core[ks]
        e_ann = ew[a_ & np.isfinite(ew)]
        e_core = ew[c_ & lw & np.isfinite(ew)]
        if e_ann.size and e_core.size:
            med = float(np.median(e_ann))
            mad = 1.4826 * float(np.median(np.abs(e_ann - med)))
            pk = float(e_core.max())
            # Peak height above the local level in units of the annulus scatter: the mosaics are
            # not background-subtracted, so a ratio to the annulus median fails on faint sources.
            r["contrast"] = (pk - med) / mad if mad > 0 else np.inf
            r["compact"] = bool(pk > med and r["contrast"] >= min_contrast)
        sc = sw[c_ & np.isfinite(sw)]
        if sc.size:
            r["s_hit"] = float(sc.flat[int(np.argmax(np.abs(sc)))])
        sa = sw[a_ & ~mw & np.isfinite(sw)]
        r["n_ann"] = int(sa.size)
        if sa.size >= min_ann_px:
            med = float(np.median(sa))
            sig = 1.4826 * float(np.median(np.abs(sa - med)))
            if sig > 0:
                r["z_ann"] = (r["s_hit"] - med) / sig
        out.append(r)
    return out


def ext_pass(row: dict, min_z: float = DEFAULTS["min_z"]) -> bool:
    """Annulus cut for an extended hit (no test when the annulus has too few pixels; a hit
    without a finite S near it fails)."""
    if not np.isfinite(row["s_hit"]):
        return False
    z = row["z_ann"]
    return bool(not np.isfinite(z) or abs(z) >= min_z)


def run(maps: Path, tab: Table, band: str, id_col: str, **kw) -> Table:
    with np.load(maps) as m:
        s = m[f"{band}_s"]
        e1 = m[f"{band}_e1img"]
        w = WCS(fits.Header.fromstring(str(m[f"{band}_wcs"])))
    scale = float(proj_plane_pixel_scales(w)[0] * 3600.0)
    x, y = w.all_world2pix(np.asarray(tab["ra"], float), np.asarray(tab["dec"], float), 0)
    rows = measure(s, e1, np.column_stack([y, x]), scale, **kw)
    out = Table(rows=rows, names=list(rows[0].keys())) if rows else Table()
    if rows:
        out.add_column(tab[id_col], index=0, name=id_col)
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--hits", type=Path, required=True)
    ap.add_argument("--id-col", default="hit")
    ap.add_argument("--band", default="F356W")
    for k, v in DEFAULTS.items():
        if k != "min_z":
            ap.add_argument("--" + k.replace("_", "-"), type=type(v), default=v)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    kw = {k: getattr(a, k) for k in DEFAULTS if k != "min_z"}
    out = run(a.maps, Table.read(a.hits), a.band, a.id_col, **kw)
    out.meta = {"provenance": "derived", "thresholds": f"ASSUMPTION: {kw}"}
    out.write(a.out, overwrite=True)


if __name__ == "__main__":
    main()

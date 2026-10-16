"""Peer-population test for ``arc_diff.py`` hits: S of a hit against S at compact epoch-1 peaks of
similar brightness in the same map (Abell 2744 UNCOVER vs ALT tiling, 2026-10-16).

Compact sources show a systematic difference between programs (PSF mismatch), so a hit's S is
compared with the S of its brightness peers. For each hit: ``log_peak`` = log10 of the largest
epoch-1 pixel within ``--peak-arcsec``; peers = local maxima of epoch 1 (``--peak-box`` px window,
positive, ≥ ``--edge-px`` from the border and from pixels without data) whose log10 value lies
within ``--dex`` of it, excluding peers within ``--self-arcsec`` of the hit; each peer's S is the
signed extreme of S within ``--peak-arcsec`` of it (peers whose window touches a pixel without S
are dropped); ``z_peers`` = (S_hit − median) / MAD (unscaled, as in ``passes.ecsv`` of
2026-10-16). The cut used in the tiling: |z_peers| ≥ ``--min-z`` or fewer than ``--min-peers``
peers (no peers = no test).

The 2026-10-16 tiling computed its peers in an uncommitted step. This reconstruction does **not**
reproduce it (``--compare`` on t−1+0 and t+0+0: cut agreement 0.35–0.54, median |Δz| 1.2–2.6, 2–23×
more peers; failed approach, 2026-10-16), so ``z_stored`` is the one to use: the hit's S against
the stored peer median and MAD of the ``--ref-passes`` rows of the same tile nearest in
``log_peak`` (median of the ``--ref-k`` nearest). Every number is ``derived``; the thresholds are
ASSUMPTIONs.

Usage: ``arc_diff_peers.py --maps M.npz --hits hits.ecsv --tile t-1+0 [--ref-passes passes.ecsv]
[--compare passes.ecsv] --out peers.ecsv``.
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


def window(a: np.ndarray, y: int, x: int, r: int) -> np.ndarray:
    """(2r + 1)² window clipped to the array (empty when the centre lies outside it)."""
    if not (0 <= y < a.shape[0] and 0 <= x < a.shape[1]):
        return a[:0, :0]
    return a[max(y - r, 0) : y + r + 1, max(x - r, 0) : x + r + 1]


def signed_extreme(win: np.ndarray) -> float:
    """The finite value of largest |value| (nan if none)."""
    ok = np.isfinite(win)
    if not ok.any():
        return float("nan")
    return float(win.flat[int(np.argmax(np.where(ok, np.abs(win), -1.0)))])


def peaks(e1: np.ndarray, box: int, edge_px: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Local maxima of the epoch-1 image: (y, x, log10 value)."""
    edge_px = max(int(edge_px), 1)
    img = np.where(np.isfinite(e1), e1, -np.inf)
    mx = ndimage.maximum_filter(img, size=box, mode="constant", cval=-np.inf)
    bad = ~np.isfinite(e1)
    near_bad = ndimage.binary_dilation(bad, iterations=edge_px) if bad.any() else bad
    ok = (img == mx) & (img > 0) & ~near_bad
    ok[:edge_px] = ok[-edge_px:] = False
    ok[:, :edge_px] = ok[:, -edge_px:] = False
    y, x = np.nonzero(ok)
    return y, x, np.log10(e1[y, x])


def peer_table(
    s: np.ndarray,
    e1: np.ndarray,
    w: WCS,
    hits: Table,
    *,
    peak_arcsec: float = 0.06,
    peak_box: int = 3,
    edge_px: int = 5,
    dex: float = 0.1,
    self_arcsec: float = 0.3,
    min_z: float = 3.0,
    min_peers: int = 1,
) -> Table:
    scale = float(proj_plane_pixel_scales(w)[0] * 3600.0)
    py, px, plog = peaks(e1, peak_box, edge_px)
    rr = max(int(np.ceil(peak_arcsec / scale)), 0)
    # Each peer's S: the signed extreme of S within peak_arcsec of the peak. Unlike the hit's own
    # S (which skips pixels without S), a peer window touching such a pixel is dropped (nan).
    gap = ndimage.maximum_filter((~np.isfinite(s)).astype(np.uint8), size=2 * rr + 1) > 0
    sa = np.nan_to_num(s)
    hi = ndimage.maximum_filter(sa, size=2 * rr + 1)
    lo = ndimage.minimum_filter(sa, size=2 * rr + 1)
    ext = np.where(gap, np.nan, np.where(np.abs(hi) >= np.abs(lo), hi, lo))
    ps = ext[py, px]
    rows = []
    for h in hits:
        x, y = (float(v) for v in w.all_world2pix(h["ra"], h["dec"], 0))
        xi, yi = int(round(x)), int(round(y))
        win = window(e1, yi, xi, rr)
        lp = float(np.log10(np.nanmax(win))) if np.nanmax(win, initial=-np.inf) > 0 else np.nan
        s_hit = signed_extreme(window(s, yi, xi, rr))
        sel = (np.abs(plog - lp) <= dex) & (np.hypot(px - x, py - y) * scale > self_arcsec)
        vals = ps[sel]
        vals = vals[np.isfinite(vals)]
        n = len(vals)
        if n:
            med = float(np.median(vals))
            mad = float(np.median(np.abs(vals - med)))
            z = (s_hit - med) / mad if mad > 0 else np.nan
        else:
            med = mad = z = np.nan
        keep = n < min_peers or not np.isfinite(z) or abs(z) >= min_z
        rows.append(
            {
                "hit": h["hit_id"],
                "ra": float(h["ra"]),
                "dec": float(h["dec"]),
                "s_hit": s_hit,
                "log_peak": lp,
                "n_peers": n,
                "peer_med_s": med,
                "peer_mad_s": mad,
                "z_peers": z,
                "peer_pass": bool(keep),
            }
        )
    t = Table(rows=rows) if rows else Table()
    t.meta = {"provenance": "derived", "thresholds": "ASSUMPTION: see arc_diff_peers.py"}
    return t


def stored_z(t: Table, ref: Table, tile: str, k: int = 3, min_z: float = 3.0) -> Table:
    """Add ``z_stored`` / ``stored_pass``: S against the stored peer statistics of ``ref``.

    A tile without usable stored rows is an error (it would pass every hit silently); a hit without
    ``log_peak`` (no positive epoch-1 pixel) gets nan and passes (no test), as in the tiling.
    """
    ref = ref[(ref["tile"] == tile) & np.isfinite(ref["peer_mad_s"]) & (ref["peer_mad_s"] > 0)]
    if not len(ref):
        raise ValueError(f"no stored peer statistics for tile {tile!r}")
    z = np.full(len(t), np.nan)
    if len(t):
        for i, r in enumerate(t):
            if not np.isfinite(r["log_peak"]):
                continue
            j = np.argsort(np.abs(ref["log_peak"] - r["log_peak"]))[:k]
            z[i] = (r["s_hit"] - np.median(ref["peer_med_s"][j])) / np.median(ref["peer_mad_s"][j])
    t["z_stored"] = z
    t["stored_pass"] = ~np.isfinite(z) | (np.abs(z) >= min_z)
    return t


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--hits", type=Path, required=True)
    ap.add_argument("--band", default="F356W")
    ap.add_argument("--peak-arcsec", type=float, default=0.06)
    ap.add_argument("--peak-box", type=int, default=3)
    ap.add_argument("--edge-px", type=int, default=5)
    ap.add_argument("--dex", type=float, default=0.1)
    ap.add_argument("--self-arcsec", type=float, default=0.3)
    ap.add_argument("--min-z", type=float, default=3.0)
    ap.add_argument("--min-peers", type=int, default=1)
    ap.add_argument("--all-hits", action="store_true", help="include near_bright_core hits")
    ap.add_argument("--compare", type=Path, help="passes.ecsv with stored z_peers")
    ap.add_argument("--ref-passes", type=Path, help="passes.ecsv: peer stats for z_stored")
    ap.add_argument("--ref-k", type=int, default=3)
    ap.add_argument("--tile", help="tile name in --compare / --ref-passes")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    with np.load(a.maps) as m:
        s = m[f"{a.band}_s"]
        e1 = m[f"{a.band}_e1img"]
        w = WCS(fits.Header.fromstring(str(m[f"{a.band}_wcs"])))
    hits = Table.read(a.hits)
    if not a.all_hits and "near_bright_core" in hits.colnames:
        hits = hits[~np.asarray(hits["near_bright_core"], bool)]
    t = peer_table(
        s,
        e1,
        w,
        hits,
        peak_arcsec=a.peak_arcsec,
        peak_box=a.peak_box,
        edge_px=a.edge_px,
        dex=a.dex,
        self_arcsec=a.self_arcsec,
        min_z=a.min_z,
        min_peers=a.min_peers,
    )
    if a.ref_passes:
        if not a.tile:
            ap.error("--ref-passes needs --tile")
        t = stored_z(t, Table.read(a.ref_passes), a.tile, a.ref_k, a.min_z)
    if a.compare:
        ref = Table.read(a.compare)
        ref = ref[ref["tile"] == a.tile] if a.tile else ref
        out = []
        for r in ref:
            d = np.hypot((t["ra"] - r["ra"]) * np.cos(np.radians(r["dec"])), t["dec"] - r["dec"])
            if len(d) and d.min() * 3600 < 0.1:
                j = int(np.argmin(d))
                out.append((r["z_peers"], t["z_peers"][j], r["n_peers"], t["n_peers"][j]))
        if out:
            o = np.array(out, float)
            ok = np.isfinite(o[:, 0]) & np.isfinite(o[:, 1])
            agree = np.mean((np.abs(o[ok, 0]) >= a.min_z) == (np.abs(o[ok, 1]) >= a.min_z))
            print(
                f"compare: {len(o)} matched, cut agreement {agree:.3f}, "
                f"median |dz| {np.median(np.abs(o[ok, 0] - o[ok, 1])):.2f}, "
                f"median n_peers ratio {np.median(o[:, 3] / np.maximum(o[:, 2], 1)):.2f}"
            )
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        t.write(a.out, overwrite=True)


if __name__ == "__main__":
    main()

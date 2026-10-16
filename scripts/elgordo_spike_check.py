"""Spike-rotation check for El Gordo ``arc_diff`` review entries (vetting; pilot, 2026-10-16).

For each review entry, the direction from the box's brightest point source (the spike source) to
the entry is compared with the source's spike directions in the mosaic pixel frame. If the entry
lies on a spike, the significance map S is sampled across the spike at the entry. A spike that
rotated by the epochs' PA_V3 difference ``dpa`` moves sideways by d·tan(dpa) at distance d, which
gives opposite-sign S on the two sides of the spike line.

Rules (ASSUMPTIONs, fixed before this script was run on the entries):
- ``on_spike``: the angle from the source to the entry is within 2° of a spike direction, and the
  distance is ≤ 6″.
- ``rotation_dipole``: on that spike, the S values sampled ±1.5 px across the spike have opposite
  signs and each |S| ≥ 2, in at least one band.

Output: ``spike.json`` (derived). The spike source is the brightest F444W pixel, located in
both bands (implementation fix after the first run, which picked a different source per band).
The spike directions are measured from the epoch-1 image itself: the angles of the 8 brightest
rays of the source between 1″ and 3″ (azimuthal profile maxima).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def ray_angles(img: np.ndarray, x0: float, y0: float, r_px: tuple[float, float], n: int = 8):
    """Angles (deg, pixel frame) of the ``n`` strongest azimuthal maxima between r_px radii."""
    yy, xx = np.indices(img.shape, dtype=float)
    r = np.hypot(xx - x0, yy - y0)
    ang = np.degrees(np.arctan2(yy - y0, xx - x0)) % 360
    m = (r >= r_px[0]) & (r <= r_px[1]) & np.isfinite(img)
    bins = np.arange(0, 361, 1.0)
    prof = np.array(
        [
            np.nanmedian(img[m & (ang >= a) & (ang < a + 1)])
            if np.any(m & (ang >= a) & (ang < a + 1))
            else -np.inf
            for a in bins[:-1]
        ]
    )
    out = []
    order = np.argsort(prof)[::-1]
    for i in order:
        a = bins[i] + 0.5
        if all(min(abs(a - b), 360 - abs(a - b)) > 10 for b in out):
            out.append(float(a))
        if len(out) == n:
            break
    return sorted(out)


def classify(angle_to_entry: float, dist_arcsec: float, spikes: list[float]) -> tuple[bool, float]:
    d = min(min(abs(angle_to_entry - s), 360 - abs(angle_to_entry - s)) for s in spikes)
    return bool(d <= 2.0 and dist_arcsec <= 6.0), float(d)


def across(S: np.ndarray, x: float, y: float, ang_deg: float, off_px: float = 1.5):
    nx, ny = -np.sin(np.radians(ang_deg)), np.cos(np.radians(ang_deg))
    vals = []
    for s in (-off_px, off_px):
        xi, yi = int(round(x + s * nx)), int(round(y + s * ny))
        vals.append(float(S[yi, xi]))
    return vals


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--maps", type=Path, required=True)
    p.add_argument("--review", type=Path, required=True)
    p.add_argument("--dpa", type=float, required=True, help="PA_V3 difference (deg), headers")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    m = np.load(a.maps, allow_pickle=True)
    entries = json.loads(a.review.read_text())["review"]
    out = {"dpa_deg": a.dpa, "bands": {}, "entries": {}}
    # one spike source for both bands: the brightest F444W pixel (fix: per-band maxima picked
    # different sources, the F150W one being the h006 star)
    w4 = WCS(fits.Header.fromstring(str(m["F444W_wcs"])))
    j4, i4 = np.unravel_index(np.nanargmax(m["F444W_e1img"]), m["F444W_e1img"].shape)
    src_sky = w4.pixel_to_world(i4, j4)
    for b in ("F150W", "F444W"):
        w = WCS(fits.Header.fromstring(str(m[b + "_wcs"])))
        img, S = m[b + "_e1img"], m[b + "_s"]
        xs, ys = w.world_to_pixel(src_sky)
        r = 3
        xi, yi = int(round(float(xs))), int(round(float(ys)))
        sub = img[yi - r : yi + r + 1, xi - r : xi + r + 1]
        dj, di = np.unravel_index(np.nanargmax(sub), sub.shape)
        y0, x0 = yi - r + dj, xi - r + di  # local peak of the same source in this band
        px = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600
        spikes = ray_angles(img, float(x0), float(y0), (1.0 / px, 3.0 / px))
        src = w.pixel_to_world(x0, y0)
        out["bands"][b] = {
            "source_radec": [round(src.ra.deg, 6), round(src.dec.deg, 6)],
            "spike_angles_pixel_deg": [round(s, 1) for s in spikes],
        }
        for e in entries:
            x, y = w.world_to_pixel_values(e["ra"], e["dec"])
            ang = float(np.degrees(np.arctan2(y - y0, x - x0)) % 360)
            dist = float(np.hypot(x - x0, y - y0) * px)
            on, dang = classify(ang, dist, spikes)
            r = out["entries"].setdefault(e["id"], {})
            r[b] = {
                "angle_from_source_deg": round(ang, 1),
                "dist_arcsec": round(dist, 2),
                "nearest_spike_offset_deg": round(dang, 1),
                "on_spike": on,
                "expected_shift_arcsec": round(dist * np.tan(np.radians(abs(a.dpa))), 3),
            }
            if on:
                v = across(S, float(x), float(y), ang)
                r[b]["S_across_spike"] = [round(t, 2) for t in v]
                r[b]["rotation_dipole"] = bool(v[0] * v[1] < 0 and min(map(abs, v)) >= 2)
    for r in out["entries"].values():
        r["on_spike_any_band"] = any(r[b]["on_spike"] for b in ("F150W", "F444W"))
        r["rotation_dipole_any_band"] = any(
            r[b].get("rotation_dipole", False) for b in ("F150W", "F444W")
        )
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

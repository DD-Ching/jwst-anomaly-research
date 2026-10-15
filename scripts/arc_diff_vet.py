"""Vet one ``arc_diff.py`` hit across several epoch pairs: forced fluxes and arc coherence.

For each map file written by ``arc_diff.py --maps`` (one per epoch pair) and each band, this
reports at a sky position:

- the forced aperture flux of the difference (nJy; ``d`` is epoch 2 minus epoch 1 after matching),
  its blank-sky noise (random apertures where epoch 1 is below its median) and its **arc noise**
  (robust spread of the same aperture on pixels whose epoch-1 aperture flux is 0.7-1.5 × the
  hit's, at least ``--exclude-arcsec`` away), with the empirical two-sided tail fraction;
- the matched-filter significance S at the hit and at the **neighbouring knots**: local maxima of
  epoch 1 above the 90th percentile within ``--neighbour-arcsec`` (closer than ``--min-sep-arcsec``
  excluded, since they are the hit itself).

A point transient changes one PSF-sized spot. A hit is **arc-coherent** in a band when at least
``--coherent-frac`` of its neighbouring knots have S of the same sign with
|S| ≥ ``--coherent-ratio`` × |S_hit|: the residual then belongs to the arc segment (a matching or
background residual of the extended light), not to a point source. The verdict for a pair uses
the band with the largest |S_hit|. Every number is ``derived``; thresholds are ASSUMPTIONs.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage

# Aperture radius (arcsec) per band; ≈ 1.5-2.5 × PSF FWHM (ASSUMPTION).
AP_RADIUS = {
    "F090W": 0.08,
    "F115W": 0.08,
    "F150W": 0.08,
    "F200W": 0.08,
    "F277W": 0.12,
    "F356W": 0.12,
    "F410M": 0.15,
    "F444W": 0.15,
}
SR_PER_ARCSEC2 = (np.pi / 180 / 3600) ** 2
NJY_PER_MJY = 1e15


def _wcs(text: str) -> WCS:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(fits.Header.fromstring(text))


def disk_sum(img: np.ndarray, r_px: float) -> np.ndarray:
    """Sum of ``img`` in a disk of radius ``r_px`` around every pixel (NaN counts as 0)."""
    n = int(np.ceil(r_px))
    yy, xx = np.mgrid[-n : n + 1, -n : n + 1]
    ker = (np.hypot(yy, xx) <= r_px).astype(float)
    return ndimage.convolve(np.nan_to_num(img), ker, mode="constant")


def neighbour_knots(
    e1: np.ndarray, x0: float, y0: float, r_max: float, r_min: float, size: int
) -> list[tuple[int, int, float]]:
    """Local maxima of ``e1`` above the 90th percentile within ``r_max`` px, beyond ``r_min`` px."""
    yy, xx = np.mgrid[: e1.shape[0], : e1.shape[1]]
    dist = np.hypot(xx - x0, yy - y0)
    region = dist <= r_max
    finite = np.where(np.isfinite(e1), e1, -np.inf)
    peaks = finite == ndimage.maximum_filter(finite, size=size)
    peaks &= region & (dist >= r_min) & (finite > np.nanpercentile(e1[region], 90))
    return [(int(y), int(x), float(dist[y, x])) for y, x in zip(*np.nonzero(peaks), strict=True)]


def coherence(s_hit: float, s_knots: list[float], ratio: float) -> float:
    """Fraction of knots with S of the hit's sign and |S| ≥ ratio × |S_hit| (NaN if no knots)."""
    if not s_knots:
        return float("nan")
    k = np.asarray(s_knots)
    return float(np.mean((np.sign(k) == np.sign(s_hit)) & (np.abs(k) >= ratio * abs(s_hit))))


def vet_band(maps, band: str, ra: float, dec: float, a, rng: np.random.Generator) -> dict | None:
    """One band's measurements; None when the position is off the map or its S window is all NaN."""
    w = _wcs(str(maps[f"{band}_wcs"]))
    x0, y0 = (float(v) for v in w.all_world2pix(ra, dec, 0))
    scale = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600
    d, s, e1 = maps[f"{band}_d"], maps[f"{band}_s"], maps[f"{band}_e1img"]
    yi, xi = int(round(y0)), int(round(x0))
    r = max(1, int(round(a.peak_arcsec / scale)))
    if not (r <= yi < d.shape[0] - r and r <= xi < d.shape[1] - r):
        return None
    win = s[yi - r : yi + r + 1, xi - r : xi + r + 1]
    if not np.isfinite(win).any():
        return None
    s_hit = float(win.flat[np.nanargmax(np.abs(win))])
    r_ap = AP_RADIUS[band] / scale
    to_njy = scale**2 * SR_PER_ARCSEC2 * NJY_PER_MJY
    fd, fe = disk_sum(d, r_ap) * to_njy, disk_sum(e1, r_ap) * to_njy
    f0, e0 = float(fd[yi, xi]), float(fe[yi, xi])
    yy, xx = np.mgrid[: d.shape[0], : d.shape[1]]
    pad = int(np.ceil(r_ap)) + 3
    inner = (xx >= pad) & (yy >= pad) & (xx < d.shape[1] - pad) & (yy < d.shape[0] - pad)
    inner &= np.isfinite(d) & np.isfinite(e1)
    far = np.hypot(xx - x0, yy - y0) * scale >= a.exclude_arcsec

    def robust(mask: np.ndarray) -> tuple[float, float, int]:
        idx = np.argwhere(mask)
        if len(idx) < 50:
            return float("nan"), float("nan"), len(idx)
        idx = idx[rng.choice(len(idx), min(3000, len(idx)), replace=False)]
        v = fd[idx[:, 0], idx[:, 1]]
        sig = 1.4826 * float(np.median(np.abs(v - np.median(v))))
        return sig, float(np.mean(np.abs(v) >= abs(f0))), int(mask.sum())

    sky_sig, _, n_sky = robust(inner & far & (e1 <= np.nanmedian(e1)))
    # Arc noise needs light at the hit; with e0 ≤ 0 (sky) it is undefined and reported as NaN.
    similar = (fe > 0.7 * e0) & (fe < 1.5 * e0) if e0 > 0 else np.zeros_like(far)
    arc_sig, arc_p, n_arc = robust(inner & far & similar)
    knots = neighbour_knots(
        e1,
        x0,
        y0,
        a.neighbour_arcsec / scale,
        a.min_sep_arcsec / scale,
        2 * int(round(r_ap)) + 1,
    )
    s_knots = [float(s[y, x]) for y, x, _ in knots]
    return {
        "flux_njy": f0,
        "e1_aperture_njy": e0,
        "sky_sigma_njy": sky_sig,
        "n_sky_pixels": n_sky,
        "arc_sigma_njy": arc_sig,
        "arc_tail_fraction": arc_p,
        "n_arc_pixels": n_arc,
        "s_hit": s_hit,
        "knots": [
            {"dist_arcsec": round(dist * scale, 3), "s": round(sv, 2)}
            for (_, _, dist), sv in sorted(zip(knots, s_knots, strict=True), key=lambda t: t[0][2])
        ],
        "coherent_fraction": coherence(s_hit, s_knots, a.coherent_ratio),
    }


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", nargs=2, action="append", required=True, metavar=("PAIR", "NPZ"))
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--peak-arcsec", type=float, default=0.06, help="S search radius at the hit")
    ap.add_argument("--neighbour-arcsec", type=float, default=1.0)
    ap.add_argument("--min-sep-arcsec", type=float, default=0.25)
    ap.add_argument("--exclude-arcsec", type=float, default=0.5)
    ap.add_argument("--coherent-ratio", type=float, default=0.5)
    ap.add_argument("--coherent-frac", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)
    pairs = {}
    for pair, path in a.maps:
        with np.load(path) as m:
            bands = sorted({k.rsplit("_", 1)[0] for k in m.files if k.endswith("_wcs")})
            res = {b: vet_band(m, b, a.ra, a.dec, a, rng) for b in bands}
        res = {b: r for b, r in res.items() if r is not None}
        if not res:
            print(pair.ljust(10), "position not covered")
            continue
        top = max(res, key=lambda b: abs(res[b]["s_hit"]))
        cf = res[top]["coherent_fraction"]
        pairs[pair] = {
            "bands": res,
            "top_band": top,
            "arc_coherent": bool(np.isfinite(cf) and cf >= a.coherent_frac),
        }
        print(
            pair.ljust(10),
            " | ".join(
                f"{b} {r['flux_njy']:+5.1f}±{r['arc_sigma_njy']:.1f} S{r['s_hit']:+.1f}"
                f" coh {r['coherent_fraction']:.2f}"
                for b, r in res.items()
            ),
            "| arc-coherent" if pairs[pair]["arc_coherent"] else "",
        )
    out = {
        "provenance": "derived",
        "position": {"ra": a.ra, "dec": a.dec},
        "thresholds_assumption": {
            k: getattr(a, k)
            for k in (
                "peak_arcsec",
                "neighbour_arcsec",
                "min_sep_arcsec",
                "exclude_arcsec",
                "coherent_ratio",
                "coherent_frac",
            )
        }
        | {"aperture_radius_arcsec": AP_RADIUS},
        "pairs": pairs,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

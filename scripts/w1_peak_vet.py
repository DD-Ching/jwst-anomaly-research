"""Streak / bright-star test of one W1 shear-screen peak (frozen rule: results/w1_shear_more/PLAN.md
addendum A).

    python scripts/w1_peak_vet.py --field macs1423 --summary <inject_shear summary.json> \
        --i2d <F200W _i2d.fits URI> --out results/w1_shear_more/vet_macs1423.json

Recomputes S at the screen maximum without (V1) used sources within ``STREAK_TOL`` of stripe
columns of the F200W cutout and (V2) used sources near the brightest catalogue source whose axis
points at it, and applies :func:`outcome`. Shapes and pixels ``observed``; S, p ``derived``;
thresholds are ASSUMPTIONs.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from jwst_anomaly import cutouts, lensmodel

sys.path.insert(0, str(Path(__file__).resolve().parent))
import inject_shear as ish  # noqa: E402

ir, es = ish.ir, ish.es

CUTOUT_ARCSEC = 30.0  # cutout used to find stripe columns (ASSUMPTION)
STRIPE_NMAD = 5.0  # a column is a stripe if |median - cutout median| > this x MAD (ASSUMPTION)
STREAK_TOL = 1.0  # arcsec from a stripe column (ASSUMPTION)
STAR_SEARCH = 10.0  # arcsec from the peak to look for the bright star (ASSUMPTION)
STAR_RADIUS = 15.0  # arcsec from the star for V2 (ASSUMPTION)
STAR_PA_TOL = 15.0  # deg between the source axis and the direction to the star (ASSUMPTION)
P_RULE = 0.05  # the D-053 p_rot threshold


def stripe_columns(img: np.ndarray, nmad: float = STRIPE_NMAD) -> np.ndarray:
    """Indices of columns whose median deviates from the cutout median by more than ``nmad``
    times the MAD of the column medians (scaled to σ). Off-mosaic fill (zero or NaN) is ignored;
    with a zero MAD no column is flagged."""
    img = np.where(img == 0.0, np.nan, np.asarray(img, float))
    with np.errstate(all="ignore"):
        med = np.nanmedian(img, axis=0)
    c = np.nanmedian(img)
    ok = np.isfinite(med)
    mad = 1.4826 * np.nanmedian(np.abs(med[ok] - np.nanmedian(med[ok]))) if ok.any() else 0.0
    if not mad > 0:
        return np.zeros(0, int)
    return np.flatnonzero(ok & (np.abs(med - c) > nmad * mad))


def axis_points_at(pa_deg, pa_to_star_deg, tol: float = STAR_PA_TOL) -> np.ndarray:
    """True where an axis (deg, mod 180) is within ``tol`` of the direction to the star."""
    return lensmodel.axis_offset_deg(pa_deg, pa_to_star_deg) <= tol


def outcome(s_clean: float, p_clean: float, b_max: float) -> str:
    """Frozen rule (addendum A): explained if S without V1 ∪ V2 is <= the B-mode max or its
    p_rot >= 0.05."""
    if s_clean <= b_max or p_clean >= P_RULE:
        return "explained by streak/spike contamination (screen-level)"
    return "not explained: to /vet-candidate"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--field", required=True, choices=ish.FIELD_NAMES)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--i2d", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--cache",
        type=Path,
        default=None,
        help="inject_shear.py output dir of this field (reuses its photo-z and Gaia files)",
    )
    args = ap.parse_args(argv)
    full = json.loads(args.summary.read_text())
    if not {"field", "screen", "assumptions"} <= set(full):
        raise SystemExit(f"{args.summary} is not an inject_shear.py per-field summary")
    if full["field"] != args.field:
        raise SystemExit(f"summary is for {full['field']}, not {args.field}")
    summ = full["screen"]
    ra_pk, dec_pk = summ["s_max_radec"]
    spec = ish.field_spec(args.field)
    with tempfile.TemporaryDirectory() as tmp:
        cache = args.cache or Path(tmp)
        cache.mkdir(parents=True, exist_ok=True)
        model, shapes, extra, _, _ = ir.load_field(args.field, cache, spec)
        # the settings the summary's screen ran with (legacy defaults for older summaries)
        screen = {**ish.LEGACY_SCREEN_OPTIONS, **full["assumptions"]["screen"]}
        sargs = SimpleNamespace(**{**ish.SHEAR_DEFAULTS, **screen})
        sargs.psf_chi = complex(*screen["psf_chi"])
        res = es.shear_screen(model, shapes, sargs, ir.psf_sigma_px(shapes), extra)
        tg = Table({"source_uid": ["peak"], "ra": [ra_pk], "dec": [dec_pk]})
        cut = cutouts.make_cutouts(
            args.i2d,
            tg,
            size_arcsec=CUTOUT_ARCSEC,
            out_dir=Path(tmp),
            storage_options={"anon": True},
        )
        with fits.open(cut["path"][0]) as h:
            img = h["SCI"].data.astype(float)
            w = WCS(h["SCI"].header)
    use = res["use"]
    ra = np.asarray(shapes["ra"], float)
    dec = np.asarray(shapes["dec"], float)
    # V1: stripe columns of the cutout
    cols = stripe_columns(img)
    px_scale = proj_plane_pixel_scales(w)[0] * 3600.0
    sx, sy = w.world_to_pixel_values(ra, dec)
    inside = (sy > -0.5) & (sy < img.shape[0] - 0.5) & (sx > -0.5) & (sx < img.shape[1] - 0.5)
    if len(cols):
        dist = np.min(np.abs(sx[:, None] - cols[None, :]), axis=1) * px_scale
    else:
        dist = np.full(len(ra), np.inf)
    v1 = use & inside & (dist <= STREAK_TOL)
    # V2: the brightest catalogue source near the peak and sources pointing at it
    sep_pk = np.hypot((ra - ra_pk) * np.cos(np.deg2rad(dec_pk)), dec - dec_pk) * 3600.0
    mag = np.asarray(shapes["mag"], float)
    near = np.flatnonzero((sep_pk <= STAR_SEARCH) & np.isfinite(mag))
    v2 = np.zeros(len(ra), bool)
    star = None
    if len(near):  # no catalogue source near the peak: V2 removes nothing
        star = int(near[np.argmin(mag[near])])
        src = SkyCoord(ra, dec, unit="deg")
        sep_star = src.separation(src[star]).arcsec
        pa_to = np.mod(src.position_angle(src[star]).deg, 180.0)
        v2 = (
            use
            & (sep_star <= STAR_RADIUS)
            & (np.arange(len(ra)) != star)
            & axis_points_at(np.asarray(shapes["pa_obs"], float), pa_to)
        )
    # S at the peak centre, with and without the vetoed rows
    x, y = res["xy"]
    cx, cy = model.to_frame(np.array([ra_pk]), np.array([dec_pk]))
    e_used = res["e"][use]
    null = res["null_max"]

    def s_at(drop: np.ndarray) -> tuple[float, float, int]:
        keep = ~drop[use]
        apm = _aperture_at(x[keep], y[keep], cx, cy, sargs)
        s, _ = apm.snr(e_used[keep])
        return float(s[0]), float(np.mean(null >= s[0])), int(keep.sum())

    none = np.zeros(len(ra), bool)
    rows = {}
    for name, drop in (("all", none), ("no_v1", v1), ("no_v2", v2), ("no_v1_v2", v1 | v2)):
        s, p, n = s_at(drop)
        rows[name] = {"s": s, "p_rot": p, "n_used_field": n}
    b_max = float(summ["b_max"])
    out = {
        "field": args.field,
        "plan": "results/w1_shear_more/PLAN.md addendum A",
        "peak_radec": [ra_pk, dec_pk],
        "b_max": b_max,
        "n_null": int(len(null)),
        "stripe_columns_px": [int(c) for c in cols],
        "pixel_scale_arcsec": px_scale,
        "n_v1": int(v1.sum()),
        "n_v2": int(v2.sum()),
        "star": None
        if star is None
        else {"ra": float(ra[star]), "dec": float(dec[star]), "mag": float(mag[star])},
        "star_sep_from_peak_arcsec": None if star is None else float(sep_pk[star]),
        "s_at_peak": rows,
        "outcome": outcome(rows["no_v1_v2"]["s"], rows["no_v1_v2"]["p_rot"], b_max),
        "provenance": {"shapes": "observed", "s_p": "derived", "thresholds": "assumption"},
    }
    args.out.write_text(json.dumps(es.json_safe(out), indent=1) + "\n")
    print(json.dumps(out["s_at_peak"]), out["outcome"])
    return 0


def _aperture_at(x, y, cx, cy, sargs):
    """The screen's aperture-mass statistic at the given centre(s)."""
    return es.ApertureMass(
        x, y, cx, cy, sargs.aperture_arcsec, sargs.filter, sargs.r_min_arcsec, sargs.min_sources
    )


if __name__ == "__main__":
    raise SystemExit(main())

"""Forced photometry under a straight stripe (diffraction spike) crossing the aperture.

``hit_lightcurve.py`` subtracts the median of an annulus. A stripe that crosses the aperture but
covers only part of the annulus leaves a positive residual, which looks like a brightening in the
one visit whose roll puts a spike there. This script models such a stripe and removes it:

1. read a box of each covering mosaic (S3 byte ranges, as ``hit_lightcurve.py``); pixel offsets
   are converted to sky offsets (east, north; arcsec) with the box WCS, so one sky geometry is
   used for every mosaic;
2. in one reference mosaic (``--ref-obs``, the one with the stripe), find the stripe position
   angle: the angle whose perpendicular-distance profile of the pixels in ``FIT_ANNULUS`` has
   the largest variance (a Radon-like search over ``ANGLE_STEP_DEG`` steps);
3. in **every** mosaic, subtract a 1-D profile B(d) across that sky angle (d = perpendicular
   distance from the line through the position), estimated as the median of the pixels whose
   along-stripe distance lies in ``ALONG_RANGE`` (so the aperture and its surroundings do not
   enter their own model);
4. measure the same aperture/annulus photometry as ``hit_lightcurve.py`` on the corrected box,
   with the random-aperture error of the corrected box.

Static light (the host, the arc) is subtracted the same way in every epoch, so epochs are
compared **within** a band, as in ``hit_lightcurve.py``. ``--inject`` puts Gaussian point sources
of known flux on the stripe in the reference mosaic, measured with and without the source, and
compares the recovered flux with the same Gaussian measured alone. Because the profile excludes
|along| < ``ALONG_RANGE[0]``, a source never enters its own model and this ratio is 1 by
construction: it checks only that, not the stripe removal (found after the h061 run; the frozen
rule kept the gate). ``stripe_control`` (the same positions without injection) and
``line_profile.ecsv`` are the post hoc checks of the removal. Only the reference mosaic's stripe
angle is removed: other mosaics' own spikes at other angles stay in their fluxes.
Mosaic reads are not retried; an unreadable mosaic is skipped and reported on stderr.
Every number is ``derived``; the radii, ranges and steps are ASSUMPTIONs.

    python scripts/stripe_phot.py --pos h061 110.82204720838848 -73.44918269902084 \\
        --ref-obs jw04043-o001_t001_nircam_clear-f444w --inject 40 \\
        --out results/transients/smacs_arcdiff/h061_spike
"""

from __future__ import annotations

import argparse
import json
import sys
import zlib
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arc_diff import read_box  # noqa: E402
from hit_lightcurve import (  # noqa: E402
    ANNULUS,
    NJY_PER_MJYSR_SR,
    aperture_photometry,
    covering_products,
    empirical_error,
)

BOX_HALF_ARCSEC = 3.0  # ASSUMPTION: box half-size read per mosaic
FIT_ANNULUS = (0.6, 2.5)  # ASSUMPTION: radii (arcsec) whose pixels define the stripe angle
ALONG_RANGE = (0.8, 2.5)  # ASSUMPTION: |along-stripe distance| (arcsec) of the profile pixels
ANGLE_STEP_DEG = 0.5  # ASSUMPTION: angle grid of the stripe search
GAUSS_FWHM = {"F444W": 0.145, "F115W": 0.040, "F090W": 0.034, "F150W": 0.050}  # ASSUMPTION


def sky_offsets(wcs, shape: tuple[int, int], ra: float, dec: float) -> tuple[np.ndarray, ...]:
    """East and north offsets (arcsec) of every pixel of a box from (ra, dec)."""
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    r, d = wcs.all_pix2world(xx.ravel(), yy.ravel(), 0)
    e = ((r - ra + 180.0) % 360.0 - 180.0) * np.cos(np.radians(dec)) * 3600.0
    n = (d - dec) * 3600.0
    return e.reshape(shape), n.reshape(shape)


def perp_along(e: np.ndarray, n: np.ndarray, pa_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Perpendicular and along-line distances for a line through the origin at sky PA (E of N)."""
    t = np.radians(pa_deg)
    ue, un = np.sin(t), np.cos(t)  # unit vector along the line
    return e * un - n * ue, e * ue + n * un


def stripe_angle(
    sci: np.ndarray, e: np.ndarray, n: np.ndarray, scale: float
) -> tuple[float, float]:
    """PA (deg, 0-180) of the strongest straight stripe through the origin, and its score.

    The score is the count-weighted variance of the binned median profile across the line.
    """
    r = np.hypot(e, n)
    m = (r >= FIT_ANNULUS[0]) & (r <= FIT_ANNULUS[1]) & np.isfinite(sci)
    v, ee, nn = sci[m], e[m], n[m]
    best = (0.0, -np.inf)
    for pa in np.arange(0.0, 180.0, ANGLE_STEP_DEG):
        d, _ = perp_along(ee, nn, pa)
        b = np.round(d / scale).astype(int)
        b -= b.min()
        cnt = np.bincount(b)
        ok = cnt >= 5
        med = np.array([np.median(v[b == i]) if ok[i] else np.nan for i in range(len(cnt))])
        w = cnt[ok]
        mu = np.sum(w * med[ok]) / np.sum(w)
        score = float(np.sum(w * (med[ok] - mu) ** 2) / np.sum(w))
        if score > best[1]:
            best = (float(pa), score)
    return best


def stripe_model(
    sci: np.ndarray, e: np.ndarray, n: np.ndarray, pa: float, scale: float
) -> np.ndarray:
    """B(d): median across-stripe profile from pixels with |along| in ``ALONG_RANGE``."""
    d, a = perp_along(e, n, pa)
    b = np.round(d / (0.5 * scale)).astype(int)  # half-pixel bins across the stripe
    sel = (np.abs(a) >= ALONG_RANGE[0]) & (np.abs(a) <= ALONG_RANGE[1]) & np.isfinite(sci)
    model = np.full(sci.shape, np.nan)
    for i in np.unique(b):
        mi = b == i
        s = sel & mi
        if s.sum() >= 5:
            model[mi] = np.median(sci[s])
    return model


def gaussian(e: np.ndarray, n: np.ndarray, e0: float, n0: float, fwhm: float) -> np.ndarray:
    """Unit-sum-per-arcsec² Gaussian (multiply by flux × pixel area in arcsec² to get pixels)."""
    s = fwhm / 2.3548
    return np.exp(-((e - e0) ** 2 + (n - n0) ** 2) / (2 * s * s)) / (2 * np.pi * s * s)


def phot(sci, e, n, scale, ra, dec, wcs, e0=0.0, n0=0.0, r_as=0.15):
    """``hit_lightcurve`` aperture photometry at sky offset (e0, n0) (box pixel units)."""
    ra1 = ra + e0 / 3600.0 / np.cos(np.radians(dec))
    x, y = (float(v) for v in wcs.all_world2pix(ra1, dec + n0 / 3600.0, 0))
    return aperture_photometry(sci, x, y, r_as / scale, (ANNULUS[0] / scale, ANNULUS[1] / scale))


def corrected(box: dict, ra: float, dec: float, pa: float) -> tuple[np.ndarray, ...]:
    e, n = sky_offsets(box["wcs"], box["sci"].shape, ra, dec)
    model = stripe_model(box["sci"], e, n, pa, box["scale"])
    return box["sci"] - np.nan_to_num(model), e, n, model


def sheet(boxes: dict, ra: float, dec: float, pa: float, band: str, path: Path) -> None:
    """Raw / stripe model / corrected stamps of every mosaic of ``band`` (rows = dates)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sel = sorted(
        (v for v in boxes.values() if str(v[0]["band"]) == band), key=lambda v: v[0]["mjd"]
    )
    fig, axs = plt.subplots(len(sel), 3, figsize=(7.5, 2.5 * len(sel)), squeeze=False)
    for row, (p, box) in zip(axs, sel, strict=True):
        corr, _, _, model = corrected(box, ra, dec, pa)
        x, y = (float(v) for v in box["wcs"].all_world2pix(ra, dec, 0))
        lo, hi = np.nanpercentile(box["sci"], [5, 99.5])
        for ax, img, t in zip(
            row, (box["sci"], model, corr), ("raw", "stripe model", "corrected"), strict=True
        ):
            ax.imshow(img, origin="lower", cmap="gray", vmin=lo, vmax=hi)
            ax.plot(x, y, "r+", ms=8)
            ax.set_title(f"{p['program']} MJD {p['mjd']:.0f} {t}", fontsize=7)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.suptitle(f"{band}: stripe PA {pa:.1f} deg E of N; red + = position", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def line_profile(boxes: dict, ra: float, dec: float, pa: float, band: str, r_as: float) -> Table:
    """Post hoc: raw aperture flux (nJy) along the stripe line through the position, per mosaic.

    A stripe gives a ridge along the whole line; a source at the position gives a local peak
    above it. Steps of ``0.1`` arcsec over ±2.5 arcsec (ASSUMPTIONs).
    """
    rows = []
    t = np.radians(pa)
    for p, box in boxes.values():
        if str(p["band"]) != band:
            continue
        to_njy = box["pixar_sr"] * NJY_PER_MJYSR_SR
        e, n = sky_offsets(box["wcs"], box["sci"].shape, ra, dec)
        for along in np.round(np.arange(-2.5, 2.51, 0.1), 2):
            e0, n0 = along * np.sin(t), along * np.cos(t)
            f = phot(box["sci"], e, n, box["scale"], ra, dec, box["wcs"], e0, n0, r_as)
            rows.append(
                (str(p["obs_id"]), str(p["program"]), float(p["mjd"]), float(along), f * to_njy)
            )
    return Table(rows=rows, names=("obs_id", "program", "mjd", "along_arcsec", "flux_njy"))


def excess_z(
    ref_flux: float, ref_err: float, other_flux: list[float], other_err: list[float]
) -> float:
    """(ref − inverse-variance mean of the others) / combined error."""
    w = 1.0 / np.asarray(other_err, float) ** 2
    mean = float(np.sum(w * np.asarray(other_flux, float)) / np.sum(w))
    return float((ref_flux - mean) / np.hypot(ref_err, 1.0 / np.sqrt(np.sum(w))))


def classify_excess(
    ref_flux: float,
    ref_err: float,
    other_flux: list[float],
    other_err: list[float],
    inj_ratio: float,
) -> str:
    """Frozen rule (PLAN.md addendum 1, before the run) for the reference epoch after correction.

    ``invalid`` unless the injection check passes (median recovered / alone in [0.7, 1.3]);
    otherwise the excess is the reference flux minus the inverse-variance mean of the other epochs
    of the same band, z = excess / sqrt(ref_err² + err_mean²): ``explained`` if z < 2,
    ``open`` if z ≥ 3, ``undetermined`` in between. All cut values are ASSUMPTIONs.
    """
    if not (0.7 <= inj_ratio <= 1.3) or not other_flux:
        return "invalid"
    z = excess_z(ref_flux, ref_err, other_flux, other_err)
    if not np.isfinite(z):  # a NaN flux or error (edge, empty profile) gives no verdict
        return "invalid"
    if z < 2.0:
        return "explained"
    return "open" if z >= 3.0 else "undetermined"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pos", nargs=3, required=True, metavar=("NAME", "RA", "DEC"))
    ap.add_argument("--ref-obs", required=True, help="obs_id of the mosaic with the stripe")
    ap.add_argument("--inject", type=int, default=0, help="injections in the reference mosaic")
    ap.add_argument("--inject-njy", type=float, default=35.0)
    ap.add_argument("--radius-arcsec", type=float, default=0.15)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    name, ra, dec = a.pos[0], float(a.pos[1]), float(a.pos[2])
    a.out.mkdir(parents=True, exist_ok=True)

    prods = covering_products([(name, ra, dec)])
    boxes = {}
    for p in prods:
        try:
            boxes[str(p["obs_id"])] = (p, read_box(str(p["uri"]), ra, dec, BOX_HALF_ARCSEC))
        except Exception as exc:  # outside the mosaic or unreadable: recorded, not fatal
            print(f"skip {p['obs_id']}: {str(exc)[:100]}", file=sys.stderr)
    if a.ref_obs not in boxes:
        raise SystemExit(f"reference mosaic {a.ref_obs} not covering or unreadable; no result")
    ref_p, ref = boxes[a.ref_obs]
    e, n = sky_offsets(ref["wcs"], ref["sci"].shape, ra, dec)
    pa, score = stripe_angle(ref["sci"], e, n, ref["scale"])

    rows = []
    for oid, (p, box) in sorted(boxes.items(), key=lambda kv: (kv[1][0]["band"], kv[1][0]["mjd"])):
        to_njy = box["pixar_sr"] * NJY_PER_MJYSR_SR
        rng = np.random.default_rng(zlib.crc32(oid.encode()))
        ann = (ANNULUS[0] / box["scale"], ANNULUS[1] / box["scale"])
        r_px = a.radius_arcsec / box["scale"]
        corr, ee, nn, _ = corrected(box, ra, dec, pa)
        f0 = phot(box["sci"], ee, nn, box["scale"], ra, dec, box["wcs"], r_as=a.radius_arcsec)
        f1 = phot(corr, ee, nn, box["scale"], ra, dec, box["wcs"], r_as=a.radius_arcsec)
        rows.append(
            {
                "obs_id": oid,
                "band": str(p["band"]),
                "program": str(p["program"]),
                "mjd": float(p["mjd"]),
                "flux_njy": f0 * to_njy,
                "err_njy": empirical_error(box["sci"], r_px, ann, rng) * to_njy,
                "flux_corr_njy": f1 * to_njy,
                "err_corr_njy": empirical_error(
                    corr, r_px, ann, np.random.default_rng(zlib.crc32(oid.encode()))
                )
                * to_njy,
            }
        )
    Table(rows=rows).write(a.out / "stripe_phot.ecsv", overwrite=True)

    inj = {}
    if a.inject:
        band = str(ref_p["band"])
        fwhm = GAUSS_FWHM.get(band, 0.1)
        to_njy = ref["pixar_sr"] * NJY_PER_MJYSR_SR
        pix_as2 = ref["scale"] ** 2
        rng = np.random.default_rng(1)
        rec, alone, ctrl_raw, ctrl_corr = [], [], [], []
        for _ in range(a.inject):
            # on the stripe, 0.8-2.0 arcsec from the position along it, either side
            along = rng.uniform(0.8, 2.0) * rng.choice([-1, 1])
            t = np.radians(pa)
            e0, n0 = along * np.sin(t), along * np.cos(t)
            g = gaussian(e, n, e0, n0, fwhm) * pix_as2 * a.inject_njy / to_njy
            # the stripe model is re-estimated around the injection position, as for a real hit
            ra1 = ra + e0 / 3600.0 / np.cos(np.radians(dec))
            dec1 = dec + n0 / 3600.0
            fs = []
            for img in (ref["sci"], ref["sci"] + g):
                c, e1, n1, _ = corrected({**ref, "sci": img}, ra1, dec1, pa)
                fs.append(
                    phot(c, e1, n1, ref["scale"], ra1, dec1, ref["wcs"], r_as=a.radius_arcsec)
                )
            rec.append((fs[1] - fs[0]) * to_njy)
            e1, n1 = sky_offsets(ref["wcs"], ref["sci"].shape, ra1, dec1)
            alone.append(
                phot(g, e1, n1, ref["scale"], ra1, dec1, ref["wcs"], r_as=a.radius_arcsec) * to_njy
            )
            # post hoc control: the same position without injection, raw and corrected (the
            # corrected value is fs[0]); its profile pixels can include the tested position itself
            ctrl_raw.append(
                phot(ref["sci"], e1, n1, ref["scale"], ra1, dec1, ref["wcs"], r_as=a.radius_arcsec)
                * to_njy
            )
            ctrl_corr.append(fs[0] * to_njy)
        rec, alone = np.array(rec), np.array(alone)
        inj = {
            "n": int(a.inject),
            "flux_njy": a.inject_njy,
            "gauss_fwhm_arcsec": fwhm,
            "aperture_flux_alone_median_njy": float(np.nanmedian(alone)),
            "recovered_median_njy": float(np.nanmedian(rec)),
            "recovered_robust_std_njy": float(
                1.4826 * np.nanmedian(np.abs(rec - np.nanmedian(rec)))
            ),
            "recovered_over_alone_median": float(np.nanmedian(rec / alone)),
            "note": "the profile excludes |along| < 0.8 arcsec, so the source never enters its own"
            " model and this ratio tests only that; stripe removal is tested by stripe_control",
        }
        stripe_control = {
            "what": f"post hoc: the {a.inject} on-stripe positions, no injection, raw and corr",
            "raw_median_njy": float(np.nanmedian(ctrl_raw)),
            "raw_robust_std_njy": float(
                1.4826 * np.nanmedian(np.abs(ctrl_raw - np.nanmedian(ctrl_raw)))
            ),
            "corr_median_njy": float(np.nanmedian(ctrl_corr)),
            "corr_robust_std_njy": float(
                1.4826 * np.nanmedian(np.abs(ctrl_corr - np.nanmedian(ctrl_corr)))
            ),
            "raw_njy": [round(float(v), 2) for v in ctrl_raw],
            "corr_njy": [round(float(v), 2) for v in ctrl_corr],
        }
        inj["stripe_control"] = stripe_control
        sheet(boxes, ra, dec, pa, str(ref_p["band"]), a.out / "sheet.png")
        line_profile(boxes, ra, dec, pa, str(ref_p["band"]), a.radius_arcsec).write(
            a.out / "line_profile.ecsv", overwrite=True
        )
    verdict, z_excess = None, None
    if inj:
        band = str(ref_p["band"])
        ref_row = next(r for r in rows if r["obs_id"] == a.ref_obs)
        oth = [r for r in rows if r["band"] == band and r["obs_id"] != a.ref_obs]
        verdict = classify_excess(
            ref_row["flux_corr_njy"],
            ref_row["err_corr_njy"],
            [r["flux_corr_njy"] for r in oth],
            [r["err_corr_njy"] for r in oth],
            inj["recovered_over_alone_median"],
        )
        z_excess = excess_z(
            ref_row["flux_corr_njy"],
            ref_row["err_corr_njy"],
            [r["flux_corr_njy"] for r in oth],
            [r["err_corr_njy"] for r in oth],
        )
    summary = {
        "verdict_classify_excess": verdict,
        "excess_z": z_excess,
        "position": [name, ra, dec],
        "ref_obs": a.ref_obs,
        "stripe_pa_deg_east_of_north": pa,
        "stripe_score": score,
        "settings": {
            "box_half_arcsec": BOX_HALF_ARCSEC,
            "fit_annulus_arcsec": FIT_ANNULUS,
            "along_range_arcsec": ALONG_RANGE,
            "angle_step_deg": ANGLE_STEP_DEG,
            "radius_arcsec": a.radius_arcsec,
            "annulus_arcsec": ANNULUS,
        },
        "injection": inj,
        "command": "python scripts/stripe_phot.py "
        + " ".join(sys.argv[1:] if argv is None else argv),
    }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    for r in rows:
        print(
            f"{r['band']:6s} {r['program']} {r['mjd']:.1f}"
            f" raw {r['flux_njy']:7.1f}±{r['err_njy']:5.1f}"
            f"  corr {r['flux_corr_njy']:7.1f}±{r['err_corr_njy']:5.1f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

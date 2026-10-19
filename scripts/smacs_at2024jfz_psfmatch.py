"""AT 2024jfz (SMACS 0723): PSF-matched NIRISS 4368 vs ERO NIRCam F150W at z_s 1.55-1.70.

Plan: ``results/transients/smacs_at2024jfz/psfmatch/PLAN.md`` (frozen before the run). The four
post hoc rows of ``../niriss4368/result.json`` that ``classify_flag`` left "undetermined" (the
flags were the wing of a static source in a plain cross-instrument difference) are re-tested on
D = E - K * T, where T is resampled onto E's grid and K is a pixel kernel fitted on the box
(``arc_diff.kernel_fit``) with the predictions and the TNS position masked. The frozen function
:func:`outcome_psfmatch` decides each row. Boxes are S3 byte-range reads; nothing is stored.
Fluxes are ``derived``; positions ``model_prediction``; thresholds ASSUMPTIONs.

    python scripts/smacs_at2024jfz_psfmatch.py --out results/transients/smacs_at2024jfz/psfmatch
"""

from __future__ import annotations

import argparse
import json
import sys
import zlib
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arc_diff import kernel_fit, read_box, resample  # noqa: E402
from hit_lightcurve import (  # noqa: E402
    ANNULUS,
    NJY_PER_MJYSR_SR,
    aperture_photometry,
    empirical_error,
)
from smacs_at2024jfz import TNS_RADEC  # noqa: E402
from smacs_at2024jfz_niriss import (  # noqa: E402
    INJ_FWHM_ARCSEC,
    MAX_NAN_FRAC,
    NSIGMA,
    POSTHOC_MU_FACTOR,
    RADIUS_ARCSEC,
    SEARCH_ARCSEC,
    grid_points,
    inject_gaussian,
)

BOX_HALF_ARCSEC = 6.5  # one box holding the four predictions and the TNS position
KERNEL_HALF = 4  # ASSUMPTION: kernel half-width in E pixels (0.26 arcsec)
MASK_ARCSEC = 0.5  # ASSUMPTION: excluded from the kernel fit around predictions and TNS
CONTROL_FRAC = 0.5  # positive control: D(TNS) >= this x F_now and >= NSIGMA sigma
POSTHOC_CORE_ARCSEC = (
    0.35  # post hoc ASSUMPTION: disc points this close to the source are its core residual
)


def outcome_psfmatch(
    disc_delta: np.ndarray,
    sig: float,
    thr: float,
    source_resid: float,
    inj_delta: float,
    control_ok: bool,
) -> str:
    """Frozen rule of PLAN.md for one row: "present", "absent" or "untestable".

    ``disc_delta``: D aperture fluxes (nJy) on the 0.1 arcsec disc grid (NaN where unmeasured);
    ``sig``: sigma of D; ``thr``: the row threshold; ``source_resid``: D at the static source
    centroid; ``inj_delta``: D at the prediction after injecting the expected counter-image
    (whole chain re-run); ``control_ok``: the global positive control at the TNS position.
    """
    if not control_ok:
        return "untestable"
    d = np.asarray(disc_delta, float)
    if d.size == 0 or float(np.mean(~np.isfinite(d))) > MAX_NAN_FRAC:
        return "untestable"
    if not (np.isfinite(source_resid) and abs(source_resid) < thr):
        return "untestable"
    if not (np.isfinite(inj_delta) and inj_delta >= thr and inj_delta >= NSIGMA * sig):
        return "untestable"
    hit = np.isfinite(d) & (d >= thr) & (d >= NSIGMA * sig)
    return "present" if hit.any() else "absent"


def control_passes(d_tns: float, sig: float, f_now: float) -> bool:
    """Positive control: the transient itself is recovered in D."""
    return bool(np.isfinite(d_tns) and d_tns >= CONTROL_FRAC * f_now and d_tns >= NSIGMA * sig)


def fit_mask(shape: tuple[int, int], xy: np.ndarray, r_px: float) -> np.ndarray:
    """True where the kernel fit may use a pixel: outside ``r_px`` of every (x, y) in ``xy``."""
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    keep = np.ones(shape, bool)
    for x, y in np.atleast_2d(xy):
        keep &= (xx - x) ** 2 + (yy - y) ** 2 > r_px**2
    return keep


def psf_matched_diff(e: dict, t_sci: np.ndarray, t_var: np.ndarray, mask: np.ndarray):
    """D = E - K * T (T already on E's grid) and the kernel K; NaN where either is undefined."""
    var = e["err"] ** 2 + t_var
    k, chi2 = kernel_fit(t_sci, e["sci"], var, KERNEL_HALF, mask=mask)
    filled = np.where(np.isfinite(t_sci), t_sci, 0.0)
    conv = ndimage.convolve(filled, k, mode="constant")
    bad = ndimage.maximum_filter((~np.isfinite(t_sci)).astype(np.uint8), size=2 * KERNEL_HALF + 1)
    d = e["sci"] - conv
    d[(bad > 0) | ~np.isfinite(e["sci"])] = np.nan
    return d, k, chi2


def phot(img: np.ndarray, wcs, radec: np.ndarray, scale: float, to_njy: float) -> np.ndarray:
    """Aperture fluxes (nJy) at ``radec`` with #360's aperture and annulus."""
    r_px, ann = RADIUS_ARCSEC / scale, (ANNULUS[0] / scale, ANNULUS[1] / scale)
    xy = np.column_stack(wcs.all_world2pix(radec[:, 0], radec[:, 1], 0))
    return np.array([aperture_photometry(img, x, y, r_px, ann) for x, y in xy]) * to_njy


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    prev = json.loads((args.out.parent / "niriss4368" / "result.json").read_text())
    rows = [r for r in prev["posthoc"]["rows"] if r["outcome"] == "undetermined"]
    if [r["z_s"] for r in rows] != [1.55, 1.6, 1.65, 1.7]:
        raise RuntimeError(f"unexpected undetermined rows: {[r['z_s'] for r in rows]}")
    f_now = float(prev["posthoc"]["f_now_njy"])
    src = np.array([rows[0]["static_check"]["source_radec_t"]])
    pred = np.array([[r["ra"], r["dec"]] for r in rows])
    centre = pred.mean(axis=0)
    uris = prev["uris"]
    e = read_box(uris["E"], centre[0], centre[1], BOX_HALF_ARCSEC)
    t = read_box(uris["T"], centre[0], centre[1], BOX_HALF_ARCSEC)
    shape = e["sci"].shape
    # T to E's grid; surface brightness (MJy/sr) is conserved by interpolation
    t_sci = resample(t["sci"], t["wcs"], e["wcs"], shape)
    t_var = resample(t["err"] ** 2, t["wcs"], e["wcs"], shape)
    # resampling a finer grid onto a coarser one keeps the per-pixel variance of the finer one;
    # scale it by the pixel-area ratio (independent-pixel ASSUMPTION; enters only the fit weights)
    t_var *= (t["scale"] / e["scale"]) ** 2
    scale = e["scale"]
    to_njy = e["pixar_sr"] * NJY_PER_MJYSR_SR
    masked = np.vstack([pred, np.array([TNS_RADEC])])
    mxy = np.column_stack(e["wcs"].all_world2pix(masked[:, 0], masked[:, 1], 0))
    mask = fit_mask(shape, mxy, MASK_ARCSEC / scale)
    d, k, chi2 = psf_matched_diff(e, t_sci, t_var, mask)
    rng = np.random.default_rng(zlib.crc32(b"psfmatch"))
    r_px, ann = RADIUS_ARCSEC / scale, (ANNULUS[0] / scale, ANNULUS[1] / scale)
    sig = float(empirical_error(d, r_px, ann, rng) * to_njy)
    d_tns = float(phot(d, e["wcs"], np.array([TNS_RADEC]), scale, to_njy)[0])
    ctrl = control_passes(d_tns, sig, f_now)
    d_src = float(phot(d, e["wcs"], src, scale, to_njy)[0])
    out_rows = []
    for r in rows:
        thr = float(r["threshold_njy"])
        pt = np.array([[r["ra"], r["dec"]]])
        disc = grid_points(np.repeat(pt, 2, axis=0), radius=SEARCH_ARCSEC)
        dd = phot(d, e["wcs"], disc, scale, to_njy)
        x, y = (float(v) for v in e["wcs"].all_world2pix(r["ra"], r["dec"], 0))
        e_inj = dict(e)
        e_inj["sci"] = inject_gaussian(
            e["sci"], x, y, thr / POSTHOC_MU_FACTOR / to_njy, INJ_FWHM_ARCSEC / scale
        )
        d_inj, _, _ = psf_matched_diff(e_inj, t_sci, t_var, mask)
        inj = float(phot(d_inj, e["wcs"], pt, scale, to_njy)[0])
        verdict = outcome_psfmatch(dd, sig, thr, d_src, inj, ctrl)
        i = int(np.nanargmax(dd)) if np.isfinite(dd).any() else None
        # post hoc (not in PLAN.md): the disc without the static source's core residual
        cosd = np.cos(np.radians(r["dec"]))
        to_src = np.hypot((disc[:, 0] - src[0, 0]) * cosd, disc[:, 1] - src[0, 1]) * 3600.0
        outer = dd[to_src > POSTHOC_CORE_ARCSEC]
        d_pred = float(phot(d, e["wcs"], pt, scale, to_njy)[0])
        out_rows.append(
            {
                "z_s": r["z_s"],
                "ra": r["ra"],
                "dec": r["dec"],
                "threshold_njy": thr,
                "n_disc": int(len(disc)),
                "nan_frac": float(np.mean(~np.isfinite(dd))),
                "delta_max_njy": None if i is None else round(float(dd[i]), 1),
                "delta_max_snr": None if i is None else round(float(dd[i] / sig), 2),
                "delta_max_radec": None if i is None else np.round(disc[i], 7).tolist(),
                "delta_min_njy": round(float(np.nanmin(dd)), 1),
                "inj_delta_njy": round(inj, 1),
                "plain_delta_max_njy": r["delta_max_njy"],
                "outcome": verdict,
                "posthoc": {
                    "d_at_prediction_njy": round(d_pred, 1),
                    "n_disc_outside_core": int(outer.size),
                    "delta_max_outside_core_njy": (
                        round(float(np.nanmax(outer)), 1) if np.isfinite(outer).any() else None
                    ),
                    "outcome_outside_core": outcome_psfmatch(outer, sig, thr, 0.0, inj, ctrl),
                },
            }
        )
    stamps(e, t_sci, d, centre, pred, src, TNS_RADEC, to_njy, args.out / "stamps_psfmatch.png")
    out = {
        "provenance": {
            "positions": "model_prediction",
            "fluxes": "derived",
            "thresholds": "assumption",
        },
        "plan": "PLAN.md",
        "uris": uris,
        "box_centre_radec": np.round(centre, 7).tolist(),
        "box_shape": list(shape),
        "pixel_scale_arcsec": round(scale, 5),
        "kernel_half_px": KERNEL_HALF,
        "kernel_sum": round(float(k.sum()), 4),
        "kernel_fit_chi2": round(chi2, 3),
        "mask_arcsec": MASK_ARCSEC,
        "sigma_d_njy": round(sig, 2),
        "control": {
            "d_tns_njy": round(d_tns, 1),
            "f_now_njy": f_now,
            "passes": ctrl,
        },
        "static_source": {
            "radec": src[0].tolist(),
            "d_njy": round(d_src, 1),
        },
        "posthoc_note": "not in PLAN.md: disc points > "
        f"{POSTHOC_CORE_ARCSEC} arcsec from the static source, with the subtraction gate "
        "waived (source_resid set to 0); does not change the frozen outcome",
        "rows": out_rows,
        "window_outcome": "absent"
        if all(r["outcome"] == "absent" for r in out_rows)
        else "not absent",
    }
    (args.out / "result.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


def stamps(e, t_sci, d, centre, pred, src, tns, to_njy, path: Path) -> None:
    """E, K-free T on E's grid, and D around the predictions; predictions, source, TNS marked."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    w = e["wcs"]
    x, y = (float(v) for v in w.all_world2pix(centre[0], centre[1], 0))
    h = int(np.ceil(1.6 / e["scale"]))
    xi, yi = int(round(x)), int(round(y))
    y0, x0 = max(yi - h, 0), max(xi - h, 0)
    sl = (slice(y0, yi + h + 1), slice(x0, xi + h + 1))
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    panels = [("E (NIRISS 4368)", e["sci"]), ("T on E grid", t_sci), ("D = E - K*T", d)]
    tx, ty = w.all_world2pix(tns[0], tns[1], 0)
    th = h
    ts = (
        slice(max(int(ty) - th, 0), int(ty) + th + 1),
        slice(max(int(tx) - th, 0), int(tx) + th + 1),
    )
    for ax, (name, img) in zip(axes[:3], panels, strict=True):
        sub = img[sl] * to_njy
        lo, hi = np.nanpercentile(sub, [5, 99.5]) if name != "D = E - K*T" else (-30, 30)
        ax.imshow(sub, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        px, py = w.all_world2pix(pred[:, 0], pred[:, 1], 0)
        ax.plot(px - x0, py - y0, "r+", ms=8)
        sx, sy = w.all_world2pix(src[:, 0], src[:, 1], 0)
        ax.plot(sx - x0, sy - y0, "co", mfc="none", ms=10)
        ax.set_title(f"{name}, 3.2 arcsec, nJy/pix", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    sub = d[ts] * to_njy
    axes[3].imshow(sub, origin="lower", cmap="gray", vmin=-30, vmax=30)
    axes[3].plot(tx - ts[1].start, ty - ts[0].start, "y+", ms=10)
    axes[3].set_title("D at the TNS position (control)", fontsize=8)
    axes[3].set_xticks([])
    axes[3].set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

"""g09725 (Abell 2744 flank, 2561 o002 -> o006) difference photometry in six NIRCam bands.

Frozen plan: ``results/transients/a2744_2561_flank/sw/PLAN.md``. Reuses the S3 byte-range
``Mosaic``, ``resample`` and ``register_offset`` of ``sneos_ratio``. The rules are functions
(``measure``, ``classify``, ``outcome``); ``run`` reads the stamps and writes ``summary.json`` and
the contact sheet. Fluxes are observed (difference images); ratios are derived.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sneos_ratio import Mosaic, register_offset, resample  # noqa: E402

RA, DEC = 3.6117000138645383, -30.457647789409304  # g09725 ra_cat / dec_cat, ../hits.ecsv
EARLY, LATE = "jw02561-o002_t001", "jw02561-o006_t007"
SW, LW = ("f115w", "f150w", "f200w"), ("f277w", "f356w", "f444w")
URI = "s3://stpubdata/jwst/public/jw02561/L3/t/o{o}/{stem}_nircam_clear-{band}_i2d.fits"
STAMP_HALF_ARCSEC = 3.0
MAX_SHIFT_ARCSEC = 0.03  # ASSUMPTION (rule 1)
MASK_ARCSEC = 0.4  # rule 1
APER_ARCSEC = 0.20  # ASSUMPTION (rule 2)
ANNULUS_ARCSEC = (0.35, 0.60)  # rule 2
RING_ARCSEC, RING_N = 1.6, 16  # rule 3
MIN_SNR = 5.0  # ASSUMPTION (rule 4)
OUT = "results/transients/a2744_2561_flank/sw"


def aperture_flux(diff: np.ndarray, err: np.ndarray, p: tuple[float, float], pix: float):
    """Rule 2: background-subtracted aperture sum and its propagated error (image units)."""
    r = APER_ARCSEC / pix
    rin, rout = (a / pix for a in ANNULUS_ARCSEC)
    ann = CircularAnnulus(p, rin, rout).to_mask(method="center").to_image(diff.shape) > 0
    bkg = float(np.nanmedian(diff[ann & np.isfinite(diff)]))
    bad = ~np.isfinite(diff) | ~np.isfinite(err)
    t = aperture_photometry(
        np.where(bad, 0, diff - bkg), CircularAperture(p, r), error=np.where(bad, 0, err)
    )
    n_bad = int(np.sum(CircularAperture(p, r).to_mask(method="center").to_image(diff.shape) * bad))
    return float(t["aperture_sum"][0]), float(t["aperture_sum_err"][0]), n_bad


def measure(diff: np.ndarray, err: np.ndarray, p: tuple[float, float], pix: float) -> dict:
    """Rules 2-3: flux at ``p``, σ = max(σ_prop, 1.4826 MAD of the ring apertures)."""
    f, s_prop, n_bad = aperture_flux(diff, err, p, pix)
    ang = 2 * np.pi * np.arange(RING_N) / RING_N
    rr = RING_ARCSEC / pix
    ring = [
        aperture_flux(diff, err, (p[0] + rr * np.cos(a), p[1] + rr * np.sin(a)), pix)[0]
        for a in ang
    ]
    ring = np.asarray(ring)
    s_emp = float(1.4826 * np.median(np.abs(ring - np.median(ring))))
    return {
        "flux": f,
        "sigma_prop": s_prop,
        "sigma_emp": s_emp,
        "sigma": max(s_prop, s_emp),
        "n_bad_pix": n_bad,
    }


def classify(flux: float, sigma: float, registered: bool) -> str:
    """Rule 4 (per band)."""
    if not registered or not np.isfinite(flux) or not np.isfinite(sigma) or sigma <= 0:
        return "unregistered"
    snr = flux / sigma
    return "brighter" if snr >= MIN_SNR else "fainter" if snr <= -MIN_SNR else "none"


def outcome(classes: dict[str, str]) -> str:
    """Rules 5-6 on {band: class}."""
    if sum(classes[b] == "brighter" for b in LW) < 2:
        return "U"
    sw = [classes[b] for b in SW]
    if "brighter" in sw:
        return "SW_detected"
    if all(c == "none" for c in sw):
        return "SW_none"
    if "fainter" in sw:
        return "SW_fainter"
    return "U"


def uri(stem: str, band: str) -> str:
    return URI.format(o=stem[9:12], stem=stem, band=band)


def run(out: str = OUT) -> dict:
    import json

    rows, sheets = {}, {}
    for band in SW + LW:
        e = Mosaic(uri(EARLY, band))
        sci_e, err_e, w = e.section(RA, DEC, STAMP_HALF_ARCSEC)
        pix = e.pix
        pixar = float(e.hdul["SCI"].header["PIXAR_SR"])
        e.close()
        m = Mosaic(uri(LATE, band))
        sci_l, err_l, wl = m.section(RA, DEC, STAMP_HALF_ARCSEC + 1.0)
        m.close()
        res, rerr = (resample(x, wl, w, sci_e.shape) for x in (sci_l, err_l))
        p = tuple(float(v) for v in w.world_to_pixel_values(RA, DEC))
        yy, xx = np.indices(sci_e.shape)
        mask = np.hypot(xx - p[0], yy - p[1]) <= MASK_ARCSEC / pix
        dx, dy = register_offset(sci_e, res, mask)
        shift = float(np.hypot(dx, dy) * pix)
        diff, derr = res - sci_e, np.hypot(err_e, rerr)
        r = measure(diff, derr, p, pix)
        to_ujy = pixar * 1e12
        row = {
            "pix_arcsec": pix,
            "shift_arcsec": shift,
            "registered": shift <= MAX_SHIFT_ARCSEC,
            "df_ujy": r["flux"] * to_ujy,
            "sigma_ujy": r["sigma"] * to_ujy,
            "sigma_prop_ujy": r["sigma_prop"] * to_ujy,
            "sigma_emp_ujy": r["sigma_emp"] * to_ujy,
            "n_bad_pix": r["n_bad_pix"],
        }
        row["snr"] = row["df_ujy"] / row["sigma_ujy"] if row["sigma_ujy"] > 0 else float("nan")
        row["class"] = classify(row["df_ujy"], row["sigma_ujy"], row["registered"])
        row["limit5_ujy"] = MIN_SNR * row["sigma_ujy"]
        rows[band] = row
        sheets[band] = (sci_e, res, diff, p, pix)
    ref = rows["f356w"]["df_ujy"]
    summary = {
        "plan": f"{out}/PLAN.md",
        "position": {"ra": RA, "dec": DEC, "id": "g09725"},
        "epochs": {"early": EARLY, "late": LATE},
        "bands": rows,
        "ratio_to_f356w": {b: rows[b]["df_ujy"] / ref if ref else float("nan") for b in rows},
        "outcome": outcome({b: rows[b]["class"] for b in rows}),
        "provenance": {"fluxes": "observed (difference images)", "ratios": "derived"},
    }
    Path(out).mkdir(parents=True, exist_ok=True)
    Path(out, "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    contact_sheet(sheets, Path(out, "contact_sheet.png"))
    return summary


def contact_sheet(sheets: dict, path) -> None:
    """o002, o006 (resampled) and the difference per band; aperture circled."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(sheets), 3, figsize=(6, 2 * len(sheets)))
    for row, (band, (e, late, d, p, pix)) in zip(axes, sheets.items(), strict=True):
        v = np.nanpercentile(np.abs(d), 99.5)
        for ax, im, lab in zip(row, (e, late, d), ("o002", "o006", "o006 - o002"), strict=True):
            lo, hi = (-v, v) if lab.startswith("o006 -") else np.nanpercentile(im, (1, 99.5))
            ax.imshow(im, origin="lower", cmap="gray", vmin=lo, vmax=hi)
            ax.add_patch(plt.Circle(p, APER_ARCSEC / pix, fill=False, ec="r", lw=0.8))
            ax.set_title(f"{band} {lab}", fontsize=7)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=70)
    plt.close(fig)


if __name__ == "__main__":
    import json as _json

    s = run()
    print(s["outcome"])
    for b, r in s["bands"].items():
        print(
            b,
            _json.dumps({k: r[k] for k in ("shift_arcsec", "df_ujy", "sigma_ujy", "snr", "class")}),
        )

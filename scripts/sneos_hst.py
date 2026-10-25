"""SN Eos 101.2 / 101.1 ratio from HST WFC3/IR F110W 2024 (PID 17504) minus 2012 (CLASH)
difference images.

Plan: ``results/transients/sneos_hst/PLAN.md`` (frozen before any pixel was read). Reuses the
resampling, registration, injection, combination and decision functions
of ``scripts/sneos_ratio.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from scipy.ndimage import shift as nd_shift

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sneos_ratio as sr  # noqa: E402

APER_ARCSEC = 0.25  # ASSUMPTION (PLAN.md rule 3)
ANNULUS_ARCSEC = (0.40, 0.60)  # PLAN.md rule 3
MAX_PRESHIFT_ARCSEC = 0.3  # ASSUMPTION (PLAN.md rule 2)
N_RANDOM = 40  # PLAN.md rule 3
RANDOM_MIN_SEP_ARCSEC = 0.6  # PLAN.md rule 3
R_S1 = (1.805, 0.039)  # Coulter et al. 2026 Table S1 F110W 2024 (sneos_ratio/literature.json)
R_JWST = (1.382, 0.083)  # #541 (sneos_ratio/summary.json)

SCIENCE = ("ifaa53010", "ifaa52010", "ifaa50010", "ifaa49010", "ifaa51010", "ifaa54010")
TEMPLATE = "ibt2a2020"
TEMPLATE_CHECK = "ibt2a9030"
OUT = "results/transients/sneos_hst"


def compare(ln_r: float, sig_ln: float) -> str:
    """PLAN.md Q2: which of Table S1 (1.805) and JWST (1.382) R_HST is consistent with within 3σ
    (σ = the two σ_ln in quadrature): 'S1', 'JWST', 'both' or 'neither'."""
    if not (np.isfinite(ln_r) and np.isfinite(sig_ln)):
        return "undetermined"
    ok = [abs(ln_r - np.log(r)) <= 3.0 * np.hypot(sig_ln, s) for r, s in (R_S1, R_JWST)]
    return {(True, True): "both", (True, False): "S1", (False, True): "JWST"}.get(
        tuple(ok), "neither"
    )


def aperture_fluxes(img: np.ndarray, p1, p2, pix: float) -> tuple[float, float]:
    """PLAN.md rule 3: background-subtracted sums at p1 (101.1), p2 (101.2); the annulus
    excludes the other aperture."""
    r = APER_ARCSEC / pix
    rin, rout = (a / pix for a in ANNULUS_ARCSEC)
    yy, xx = np.indices(img.shape)
    out = []
    for p, q in ((p1, p2), (p2, p1)):
        other = np.hypot(xx - q[0], yy - q[1]) <= r
        ann = CircularAnnulus(p, rin, rout).to_mask(method="center").to_image(img.shape) > 0
        sel = ann & ~other & np.isfinite(img)
        bkg = float(np.median(img[sel]))
        bad = ~np.isfinite(img)
        t = aperture_photometry(img - bkg, CircularAperture(p, r), mask=bad)
        out.append(float(t["aperture_sum"][0]))
    return out[0], out[1]


def random_sigma(img: np.ndarray, p1, p2, pix: float, seed: int = 0) -> float:
    """PLAN.md rule 3: 1.4826 × MAD of background-subtracted aperture sums at N_RANDOM positions
    ≥ 0.6″ from both SN images, with the whole annulus inside the stamp."""
    rng = np.random.default_rng(seed)
    edge = ANNULUS_ARCSEC[1] / pix + 1
    ny, nx = img.shape
    sums = []
    while len(sums) < N_RANDOM:
        x, y = rng.uniform(edge, nx - 1 - edge), rng.uniform(edge, ny - 1 - edge)
        if min(np.hypot(x - p[0], y - p[1]) for p in (p1, p2)) * pix < RANDOM_MIN_SEP_ARCSEC:
            continue
        sums.append(_single(img, (x, y), pix))
    s = np.asarray(sums)
    return float(1.4826 * np.median(np.abs(s - np.median(s))))


def _single(img: np.ndarray, p, pix: float) -> float:
    r = APER_ARCSEC / pix
    rin, rout = (a / pix for a in ANNULUS_ARCSEC)
    ann = CircularAnnulus(p, rin, rout).to_mask(method="center").to_image(img.shape) > 0
    bkg = float(np.median(img[ann & np.isfinite(img)]))
    t = aperture_photometry(img - bkg, CircularAperture(p, r), mask=~np.isfinite(img))
    return float(t["aperture_sum"][0])


def ratio_with_error(f1, f2, s1, s2) -> dict:
    ratio = f2 / f1
    return {
        "f1": f1,
        "f2": f2,
        "sf1": s1,
        "sf2": s2,
        "snr1": f1 / s1,
        "snr2": f2 / s2,
        "ratio": ratio,
        "sratio": abs(ratio) * float(np.hypot(s1 / f1, s2 / f2)),
        "detected": bool(f1 / s1 >= sr.MIN_SNR and f2 / s2 >= sr.MIN_SNR),
    }


def injection_trials(diff, psf, p1, p2, f1, pix) -> list[dict]:
    """PLAN.md rule 5 (as sneos_ratio.injection_trials, with this plan's apertures)."""
    mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
    sep = (p2[0] - p1[0], p2[1] - p1[1])
    r = sr.INJ_RADIUS_ARCSEC / pix
    out = []
    for ang in (0.0, 120.0, 240.0):
        q1 = (
            mid[0] + r * np.cos(np.radians(ang)) - sep[0] / 2,
            mid[1] + r * np.sin(np.radians(ang)) - sep[1] / 2,
        )
        q2 = (q1[0] + sep[0], q1[1] + sep[1])
        b1, b2 = aperture_fluxes(diff, q1, q2, pix)
        for ratio in sr.INJ_RATIOS:
            g1, g2 = aperture_fluxes(sr.inject_pair(diff, psf, q1, sep, f1, ratio), q1, q2, pix)
            meas = (g2 - b2) / (g1 - b1)
            out.append(
                {
                    "angle_deg": ang,
                    "injected": ratio,
                    "measured": float(meas),
                    "recovered": sr.recovered(meas, ratio),
                }
            )
    return out


# --- run driver -------------------------------------------------------------------------------


class Drz(sr.Mosaic):
    """SCI section and header keywords of one HST ``_drz.fits`` by S3 byte range (no ERR)."""

    def __init__(self, obs: str):
        super().__init__(f"s3://stpubdata/hst/public/{obs[:4]}/{obs}/{obs}_drz.fits")
        h0, h1 = self.hdul[0].header, self.hdul["SCI"].header
        self.photflam = float(h1.get("PHOTFLAM", h0.get("PHOTFLAM")))
        self.exptime = float(h0["EXPTIME"])
        self.date = str(h0.get("DATE-OBS"))

    def section(self, ra, dec, half_arcsec):
        x, y = (int(round(float(v))) for v in self.wcs.world_to_pixel_values(ra, dec))
        h = int(np.ceil(half_arcsec / self.pix))
        ny, nx = self.hdul["SCI"].header["NAXIS2"], self.hdul["SCI"].header["NAXIS1"]
        if not (h <= x < nx - h and h <= y < ny - h):
            raise ValueError(f"section at ({ra}, {dec}) runs off the mosaic")
        sl = (slice(y - h, y + h + 1), slice(x - h, x + h + 1))
        sci = np.asarray(self.hdul["SCI"].section[sl], float) * self.photflam  # rule 1
        return sci, None, self.wcs[sl]


def difference(sci, w, tmpl: Drz, ra, dec, sn_mask, pix) -> tuple[np.ndarray, dict]:
    """PLAN.md rule 2."""
    t_sci, _, tw = tmpl.section(ra, dec, sr.STAMP_HALF_ARCSEC + 1.0)
    res = sr.resample(t_sci, tw, w, sci.shape)
    dx, dy = sr.register_offset(sci, res, sn_mask)
    info = {"offset0_arcsec": float(np.hypot(dx, dy) * pix), "preshift_px": [dx, dy]}
    if info["offset0_arcsec"] <= MAX_PRESHIFT_ARCSEC:
        # deviation: edges take the nearest value; a NaN edge (cval=nan) biased the re-measured
        # offset by up to 0.7 px
        res = nd_shift(np.nan_to_num(res), (-dy, -dx), order=1, mode="nearest")
        dx2, dy2 = sr.register_offset(sci, res, sn_mask)
        info["shift_arcsec"] = float(np.hypot(dx2, dy2) * pix)
    else:
        info["shift_arcsec"] = info["offset0_arcsec"]
    info["registered"] = info["shift_arcsec"] <= sr.MAX_SHIFT_ARCSEC
    return sci - res, info


def run(out: str = OUT) -> dict:
    import json

    from astropy.table import Table

    t = Table.read(sr.IMAGES, format="ascii.ecsv")
    t = t[np.isin(np.asarray(t["image_id"], str), ["101.1", "101.2"])]
    t.sort("image_id")
    (ra1, ra2), (de1, de2) = list(map(float, t["ra"])), list(map(float, t["dec"]))
    ram, dem = (ra1 + ra2) / 2, (de1 + de2) / 2
    star = sr.gaia_psf_star(ram, dem)
    tmpl, tchk = Drz(TEMPLATE), Drz(TEMPLATE_CHECK)
    rows, sheets = [], {}
    for obs in SCIENCE:
        m = Drz(obs)
        sci, _, w = m.section(ram, dem, sr.STAMP_HALF_ARCSEC)
        pix = m.pix
        try:
            ps, _, _ = m.section(star["ra"], star["dec"], sr.PSF_HALF_ARCSEC)
            c = ps.shape[0] // 2
            border = np.concatenate([ps[0], ps[-1], ps[:, 0], ps[:, -1]])
            psf = (
                np.nan_to_num(ps - np.nanmedian(border))
                if np.all(np.isfinite(ps[c - 2 : c + 3, c - 2 : c + 3]))
                else None
            )
        except ValueError:
            psf = None
        p1 = tuple(float(v) for v in w.world_to_pixel_values(ra1, de1))
        p2 = tuple(float(v) for v in w.world_to_pixel_values(ra2, de2))
        yy, xx = np.indices(sci.shape)
        sn = (np.hypot(xx - p1[0], yy - p1[1]) <= sr.SN_MASK_ARCSEC / pix) | (
            np.hypot(xx - p2[0], yy - p2[1]) <= sr.SN_MASK_ARCSEC / pix
        )
        diff, info = difference(sci, w, tmpl, ram, dem, sn, pix)
        f1, f2 = aperture_fluxes(diff, p1, p2, pix)
        s_rand = random_sigma(diff, p1, p2, pix)
        pois = [np.sqrt(max(f, 0.0) / (m.exptime * m.photflam)) * m.photflam for f in (f1, f2)]
        ph = ratio_with_error(
            f1, f2, float(np.hypot(s_rand, pois[0])), float(np.hypot(s_rand, pois[1]))
        )
        row = {
            "science": obs,
            "date": m.date,
            "exptime": m.exptime,
            "pix_arcsec": pix,
            "sigma_random": s_rand,
            **info,
            **ph,
        }
        if obs == SCIENCE[0]:
            d2, i2 = difference(sci, w, tchk, ram, dem, sn, pix)
            g1, g2 = aperture_fluxes(d2, p1, p2, pix)
            row["template_check"] = {
                "template": TEMPLATE_CHECK,
                "f1": g1,
                "f2": g2,
                "ratio": g2 / g1,
                **i2,
            }
        if psf is not None and psf.sum() > 0 and info["registered"]:
            f_inj = f1 if f1 > 0 else sr.MIN_SNR * ph["sf1"]
            tr = injection_trials(diff, psf, p1, p2, f_inj, pix)
            row["injections"] = tr
            row["gate_I"] = sum(x["recovered"] for x in tr) >= 5
        else:
            row["gate_I"] = False
            row["gate_I_note"] = "not run (no PSF in mosaic or not registered)"
        row["used"] = bool(
            info["registered"] and ph["detected"] and row["gate_I"] and ph["ratio"] > 0
        )
        rows.append(row)
        sheets[obs] = (sci, diff, p1, p2)
        m.close()
    tmpl.close()
    tchk.close()
    used = [r for r in rows if r["used"]]
    ln_r = np.log([r["ratio"] for r in used]) if used else np.array([])
    sig = np.array([r["sratio"] / r["ratio"] for r in used]) if used else np.array([])
    mean, s = sr.combine(ln_r, sig)
    summary = {
        "plan": f"{OUT}/PLAN.md",
        "r_model": sr.R_MODEL,
        "r_table_s1": R_S1,
        "r_jwst": R_JWST,
        "template": TEMPLATE,
        "psf_star": star,
        "rows": rows,
        "n_used": len(used),
        "ln_r_obs": mean,
        "sigma_ln": s,
        "r_obs": float(np.exp(mean)) if np.isfinite(mean) else None,
        "q1_decision": sr.decide(mean, s),
        "q2_compare": compare(mean, s),
        "provenance": {
            "fluxes": "observed (difference images)",
            "r_model": "model_prediction",
            "r_table_s1": "derived (published photometry)",
        },
    }
    Path(out).mkdir(parents=True, exist_ok=True)
    Path(out, "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    contact_sheet(sheets, Path(out, "contact_sheet.png"))
    return summary


def contact_sheet(sheets: dict, path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, len(sheets), figsize=(2.2 * len(sheets), 4.6))
    for k, (obs, (sci, diff, p1, p2)) in enumerate(sheets.items()):
        for j, (img, lab) in enumerate(((sci, "2024"), (diff, "2024 − 2012"))):
            a = ax[j, k]
            d0 = diff - np.nanmedian(diff)  # display only: remove the sky pedestal
            v = np.nanpercentile(np.abs(d0), 99.5)
            shown = img - np.nanmedian(img)
            a.imshow(shown, origin="lower", cmap="gray", vmin=-v, vmax=v * (4 if j == 0 else 1))
            for p in (p1, p2):
                a.add_patch(plt.Circle(p, APER_ARCSEC / 0.128, fill=False, color="r", lw=0.6))
            a.set_title(f"{obs} {lab}", fontsize=6)
            a.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    s = run()
    print({k: s[k] for k in ("n_used", "r_obs", "sigma_ln", "q1_decision", "q2_compare")})

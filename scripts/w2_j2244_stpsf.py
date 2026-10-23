"""STScI official PSF capability check (results/w2_j2244_stpsf/PLAN.md, frozen at 908ecca).

The WFC3-team default STDPSF model (hst_notebooks WFC3 PSF notebook, section 2.2.1), read with
photutils `GriddedPSFModel.read(format="stdpsf")` and fitted with `PSFPhotometry` on single FLT
exposures, for three development cases. Only the data interface, settings and result handling are
written here; the model and the fitter are used as they are. Run in an isolated environment with
photutils (offline, inputs from scripts/net_ledger.py downloads):

    python scripts/w2_j2244_stpsf.py DATA_DIR OUT_DIR
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling.models import Sersic2D
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from photutils.psf import GriddedPSFModel, PSFPhotometry, SourceGrouper
from scipy.signal import fftconvolve

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w2_j2244_stpsf_rules as rules  # noqa: E402

ZP = 25.93113467346076  # F160W AB, GO 17455 (results/w2_j2244_hst/runs/inj.json)
PIX = 0.1283  # arcsec per FLT pixel, nominal (r_eff conversion and rms in px)
FULL = 1014  # WFC3/IR science frame (STDPSF grid at 0, 507, 1014)
FIT_SHAPE = (11, 11)
BAD_DQ = 1 | 2 | 4 | 8 | 16 | 32 | 256 | 512
ANNULI = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 8), (8, 12)]
HALF = 15  # star stamp half-size (px)
OFFSET = (-3, 7)  # case 2: B relative to A's peak pixel (dx, dy)
APER_R = 2.0  # extended-light proxy aperture at P (px)
INJ = {"n": 4.0, "q": 0.8, "r_eff": 0.2, "pa_deg": 30.0, "mags": (20.0, 21.0, 22.0), "frac": 0.4}
EXPOSURES = {
    "case1": ["if8a06cuq", "if8a06cvq", "if8a06cwq", "if8a06cxq"],
    "A": ["if8a4jf4q", "if8a4jf5q", "if8a4jf6q", "if8a4jf7q"],
    "B": ["if8a38guq", "if8a38gvq", "if8a38gwq", "if8a38gxq"],
}
GAIA = {"case1": 2755194341830999040, "A": 2834388281729993472, "B": 299642795747798272}


def load(data_dir: Path, root: str) -> dict:
    with fits.open(data_dir / f"{root}_flt.fits") as h:
        sci, err, dq = (h[e].data.astype(float) for e in ("SCI", "ERR", "DQ"))
        hdr = h["SCI"].header
        wcs = WCS(hdr, h)
        exptime = float(h[0].header["EXPTIME"])
    x0, y0 = -int(hdr["LTV1"]), -int(hdr["LTV2"])  # sub -> full: x_full = x_sub - LTV1
    return {
        "sci": sci,
        "err": err,
        "mask": (dq.astype(int) & BAD_DQ) > 0,
        "wcs": wcs,
        "x0": x0,
        "y0": y0,
        "exptime": exptime,
        "root": root,
    }


def full_frame(e: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.zeros((FULL, FULL))
    er = np.ones((FULL, FULL))
    m = np.ones((FULL, FULL), bool)
    sl = np.s_[e["y0"] : e["y0"] + e["sci"].shape[0], e["x0"] : e["x0"] + e["sci"].shape[1]]
    d[sl], er[sl], m[sl] = e["sci"], e["err"], e["mask"]
    return d, er, m


def peak(e: dict, ra: float, dec: float) -> tuple[int, int]:
    """Brightest unmasked pixel (3x3 box sum) within 5 px of the Gaia position (sub-array frame)."""
    x, y = e["wcs"].world_to_pixel(SkyCoord(ra, dec, unit="deg"))
    xi, yi = int(round(float(x))), int(round(float(y)))
    d = np.where(e["mask"], 0.0, e["sci"])
    best, bxy = -np.inf, (xi, yi)
    for yy in range(yi - 5, yi + 6):
        for xx in range(xi - 5, xi + 6):
            s = d[yy - 1 : yy + 2, xx - 1 : xx + 2].sum()
            if s > best:
                best, bxy = s, (xx, yy)
    return bxy


def ring(shape, xc, yc, r0, r1):
    yy, xx = np.indices(shape)
    r = np.hypot(xx - xc, yy - yc)
    return (r >= r0) & (r < r1)


def bkg(d, m, xc, yc, r0, r1) -> float:
    sel = ring(d.shape, xc, yc, r0, r1) & ~m
    return float(sigma_clipped_stats(d[sel], sigma=3.0)[1])


def fit(model, d, er, m, xy: list[tuple[float, float]]):
    flux0 = []
    for x, y in xy:
        sel = ring(d.shape, x, y, 0, 3.01) & ~m
        flux0.append(d[sel].sum() / 0.943)  # centre-PSF fraction within r <= 3 px (PLAN section 1)
    init = Table({"x": [p[0] for p in xy], "y": [p[1] for p in xy], "flux": flux0})
    phot = PSFPhotometry(model, FIT_SHAPE, grouper=SourceGrouper(min_separation=20))
    res = phot(d, error=er, mask=m, init_params=init)
    resid = phot.make_residual_image(d, psf_shape=(25, 25))
    return res, resid


def col(res, *names):
    for n in names:
        if n in res.colnames:
            return np.asarray(res[n], float)
    raise KeyError(names)


def star_stats(resid, er, m, xc, yc, flux) -> dict:
    core = ring(resid.shape, xc, yc, 0, 5.0001) & ~m
    chi = float(((resid[core] / er[core]) ** 2).sum() / max(core.sum() - 3, 1))
    z, frac = [], []
    for r0, r1 in ANNULI:
        sel = ring(resid.shape, xc, yc, r0, r1) & ~m
        z.append(float(resid[sel].sum() / np.sqrt((er[sel] ** 2).sum())))
        frac.append(float(resid[sel].sum() / flux))
    return {"redchi_r5": chi, "z_annuli": z, "frac_annuli": frac}


def aperture(img, er, m, xc, yc, r):
    sel = ring(img.shape, xc, yc, 0, r + 1e-4) & ~m
    return float(img[sel].sum()), float(np.sqrt((er[sel] ** 2).sum()))


def psf_at(model, x: float, y: float) -> np.ndarray:
    """Bilinear interpolation of the 3x3 STDPSF grid at (x, y) (oversampled array)."""
    gx = np.unique(model.grid_xypos[:, 0])
    gy = np.unique(model.grid_xypos[:, 1])
    ix = int(np.clip(np.searchsorted(gx, x) - 1, 0, len(gx) - 2))
    iy = int(np.clip(np.searchsorted(gy, y) - 1, 0, len(gy) - 2))
    tx = (x - gx[ix]) / (gx[ix + 1] - gx[ix])
    ty = (y - gy[iy]) / (gy[iy + 1] - gy[iy])
    idx = {tuple(p): k for k, p in enumerate(model.grid_xypos.tolist())}
    a = model.data
    return (
        (1 - tx) * (1 - ty) * a[idx[(gx[ix], gy[iy])]]
        + tx * (1 - ty) * a[idx[(gx[ix + 1], gy[iy])]]
        + (1 - tx) * ty * a[idx[(gx[ix], gy[iy + 1])]]
        + tx * ty * a[idx[(gx[ix + 1], gy[iy + 1])]]
    )


def sersic_image(model, px: float, py: float, mag: float, shape) -> np.ndarray:
    """Noiseless Sersic (PLAN case 3) in e-/s on the detector grid: the profile on a 0.25 px lattice
    aligned with pixel centres, convolved with the oversampled ePSF at P, sampled at pixel centres.
    The flux is normalised within +-12 px of P."""
    os_ = 4
    cx, cy = int(round(px)), int(round(py))
    n = 12 * os_
    off = np.arange(-n, n + 1) / os_
    xx, yy = np.meshgrid(cx + off, cy + off)
    s = Sersic2D(
        amplitude=1.0,
        r_eff=INJ["r_eff"] / PIX,
        n=INJ["n"],
        x_0=px,
        y_0=py,
        ellip=1 - INJ["q"],
        theta=np.deg2rad(INJ["pa_deg"]),
    )(xx, yy)
    s *= 10 ** (-0.4 * (mag - ZP)) / s.sum()
    conv = fftconvolve(s, psf_at(model, px, py), mode="same")
    pix = conv[::os_, ::os_]  # lattice points at whole-pixel offsets from (cx, cy)
    out = np.zeros(shape)
    h = pix.shape[0] // 2
    out[cy - h : cy + h + 1, cx - h : cx + h + 1] = pix
    return out


def case1(model, data_dir: Path, gaia: dict) -> dict:
    rows, sky = [], []
    for k, root in enumerate(EXPOSURES["case1"]):
        e = load(data_dir, root)
        xs, ys = peak(e, *gaia["case1"])
        d, er, m = full_frame(e)
        xc, yc = xs + e["x0"], ys + e["y0"]
        d = d - bkg(d, m, xc, yc, 10, 15)
        res, resid = fit(model, d, er, m, [(xc, yc)])
        xf, yf = col(res, "x_fit")[0], col(res, "y_fit")[0]
        fl, fe = col(res, "flux_fit")[0], col(res, "flux_err")[0]
        st = star_stats(resid, er, m, xf, yf, fl)
        ap = aperture(d, er, m, xf, yf, 6.0)[0]
        c = e["wcs"].pixel_to_world(xf - e["x0"], yf - e["y0"])
        sky.append(c)
        rows.append(
            {
                "root": root,
                "x_full": xf,
                "y_full": yf,
                "flux": fl,
                "flux_err": fe,
                "aper6_over_fit": ap / fl,
                **st,
            }
        )
        if k == 0:
            stamp = (d, resid, m, xf, yf)
    ra0 = np.mean([c.ra.deg for c in sky])
    dec0 = np.mean([c.dec.deg for c in sky])
    ref = SkyCoord(ra0, dec0, unit="deg")
    sep_mas = [float(ref.separation(c).arcsec * 1000) for c in sky]
    rms_mas = float(np.sqrt(np.mean(np.square(sep_mas))))
    fl = np.array([r["flux"] for r in rows])
    out = {
        "exposures": rows,
        "sky_offsets_mas": sep_mas,
        "sky_rms_mas": rms_mas,
        "sky_rms_px": rms_mas / (PIX * 1000),
        "flux_frac_std": float(fl.std(ddof=1) / np.median(fl)),
        "flux_frac_err_median": float(np.median([r["flux_err"] for r in rows]) / np.median(fl)),
    }
    return out, stamp


def pair_frames(model, data_dir: Path, gaia: dict, k: int) -> dict:
    a = load(data_dir, EXPOSURES["A"][k])
    b = load(data_dir, EXPOSURES["B"][k])
    xa, ya = peak(a, *gaia["A"])
    xb, yb = peak(b, *gaia["B"])
    d, er, m = full_frame(a)
    xa_f, ya_f = xa + a["x0"], ya + a["y0"]
    # star B: background-subtracted stamp from its own frame, added at the fixed integer offset
    bd, ber, bm = full_frame(b)
    xb_f, yb_f = xb + b["x0"], yb + b["y0"]
    bd = bd - bkg(bd, bm, xb_f, yb_f, 10, 15)
    tx, ty = xa_f + OFFSET[0], ya_f + OFFSET[1]
    src = np.s_[yb_f - HALF : yb_f + HALF + 1, xb_f - HALF : xb_f + HALF + 1]
    dst = np.s_[ty - HALF : ty + HALF + 1, tx - HALF : tx + HALF + 1]
    d[dst] += bd[src]
    er[dst] = np.sqrt(er[dst] ** 2 + ber[src] ** 2)
    m[dst] |= bm[src]
    mx, my = (xa_f + tx) / 2, (ya_f + ty) / 2
    d = d - bkg(d, m, mx, my, 12, 18)
    # P from initial peaks and r <= 3 px aperture fluxes (fainter + 0.4 x (brighter - fainter))
    fa = aperture(d, er, m, xa_f, ya_f, 3.0)[0]
    fb = aperture(d, er, m, tx, ty, 3.0)[0]
    (fx, fy), (bx_, by_) = ((xa_f, ya_f), (tx, ty)) if fa <= fb else ((tx, ty), (xa_f, ya_f))
    px, py = fx + INJ["frac"] * (bx_ - fx), fy + INJ["frac"] * (by_ - fy)
    dpsf = psf_at(model, xb_f, yb_f) - psf_at(model, tx, ty)
    return {
        "d": d,
        "er": er,
        "m": m,
        "xy": [(xa_f, ya_f), (tx, ty)],
        "P": (px, py),
        "exptime": a["exptime"],
        "roots": (a["root"], b["root"]),
        "psf_change_B": float(np.abs(dpsf).max() / psf_at(model, tx, ty).max()),
    }


def case23(model, data_dir: Path, gaia: dict) -> tuple[dict, dict, tuple]:
    rows, inj_rows = [], []
    for k in range(4):
        f = pair_frames(model, data_dir, gaia, k)
        res, resid = fit(model, f["d"], f["er"], f["m"], f["xy"])
        xf, yf, fl = col(res, "x_fit"), col(res, "y_fit"), col(res, "flux_fit")
        stars = [star_stats(resid, f["er"], f["m"], xf[i], yf[i], fl[i]) for i in range(2)]
        s0, n0 = aperture(resid, f["er"], f["m"], *f["P"], APER_R)
        rows.append(
            {
                "roots": f["roots"],
                "x_full": xf.tolist(),
                "y_full": yf.tolist(),
                "flux": fl.tolist(),
                "P": list(f["P"]),
                "S": s0,
                "N": n0,
                "psf_change_B": f["psf_change_B"],
                "stars": stars,
            }
        )
        if k == 0:
            stamp = (f["d"], resid, f["m"], *f["P"])
        for j, mag in enumerate(INJ["mags"]):
            gal = sersic_image(model, *f["P"], mag, f["d"].shape)
            rng = np.random.default_rng(1000 + 10 * k + j)
            noisy = rng.poisson(np.clip(gal, 0, None) * f["exptime"]) / f["exptime"]
            d2 = f["d"] + noisy
            er2 = np.sqrt(f["er"] ** 2 + np.clip(gal, 0, None) / f["exptime"])
            _, resid2 = fit(model, d2, er2, f["m"], f["xy"])
            s1, n1 = aperture(resid2, er2, f["m"], *f["P"], APER_R)
            fap = aperture(gal, er2, f["m"], *f["P"], APER_R)[0]
            inj_rows.append(
                {
                    "k": k,
                    "mag": mag,
                    "S": s1,
                    "N": n1,
                    "F_ap": fap,
                    "S_noinj": s0,
                    "F_total_stamp": float(gal.sum()),
                }
            )
            if k == 0 and mag == 20.0:
                stamp_inj = (d2, resid2, f["m"], *f["P"])
    z_p = rules.combine([r["S"] for r in rows], [r["N"] for r in rows])
    inj = {}
    for mag in INJ["mags"]:
        rr = [r for r in inj_rows if r["mag"] == mag]
        z = rules.combine([r["S"] for r in rr], [r["N"] for r in rr])
        rec = (sum(r["S"] for r in rr) - sum(r["S_noinj"] for r in rr)) / sum(r["F_ap"] for r in rr)
        inj[str(mag)] = {"z_P": z, "recovered_fraction": rec, "exposures": rr}
    return {"exposures": rows, "z_P": z_p}, inj, (stamp, stamp_inj)


def png(path: Path, stamps: list, titles: list) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(len(stamps), 2, figsize=(5.2, 2.5 * len(stamps)))
    for i, (st, t) in enumerate(zip(stamps, titles, strict=True)):
        d, resid, m, xc, yc = st
        sl = np.s_[int(yc) - 14 : int(yc) + 15, int(xc) - 14 : int(xc) + 15]
        dd = np.where(m[sl], np.nan, d[sl])
        rr = np.where(m[sl], np.nan, resid[sl])
        ax[i, 0].imshow(np.arcsinh(dd / max(np.nanstd(dd), 1e-3)), origin="lower", cmap="gray")
        v = np.nanpercentile(np.abs(rr), 99)
        ax[i, 1].imshow(rr, origin="lower", cmap="RdBu_r", vmin=-v, vmax=v)
        ax[i, 0].set_title(f"{t}: data", fontsize=8)
        ax[i, 1].set_title(f"{t}: data - PSF model (+-{v:.2f} e/s)", fontsize=8)
        for a in ax[i]:
            a.set_xticks([])
            a.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=80)


def main(data_dir: str, out_dir: str) -> None:
    data_dir, out = Path(data_dir), Path(out_dir)
    model = GriddedPSFModel.read(data_dir / "STDPSF_WFC3IR_F160W.fits", format="stdpsf")
    rows = json.load(open(data_dir / "gaia3.json"))
    names = [c["name"] for c in rows["metadata"]]
    pos = {
        r[names.index("source_id")]: (r[names.index("ra")], r[names.index("dec")])
        for r in rows["data"]
    }
    gaia = {k: pos[v] for k, v in GAIA.items()}
    c1, st1 = case1(model, data_dir, gaia)
    c2, c3, (st2, st3) = case23(model, data_dir, gaia)
    h = rules.hookup(c1["sky_rms_px"])
    s_ok = rules.star_ok(
        [r["redchi_r5"] for r in c1["exposures"]], [r["z_annuli"] for r in c1["exposures"]]
    )
    p_ok = rules.pair_ok(c2["z_P"])
    i20 = rules.injection_ok(c3["20.0"]["z_P"], c3["20.0"]["recovered_fraction"])
    worth = rules.worth_validation(h, s_ok, p_ok, i20, c3["21.0"]["z_P"])
    result = {
        "label": "derived (development cases on #474 controls; official default STDPSF model "
        "fitted with photutils PSFPhotometry on FLT exposures; injections simulated)",
        "settings": {
            "fit_shape": FIT_SHAPE,
            "bad_dq": BAD_DQ,
            "annuli_px": ANNULI,
            "aperture_r_px": APER_R,
            "offset_B_px": OFFSET,
            "injection": INJ,
            "zp_ab": ZP,
        },
        "case1": c1,
        "case2": c2,
        "case3": c3,
        "decisions": {"H": h, "S": s_ok, "P": p_ok, "I_20": i20, "worth_validation": worth},
    }
    json.dump(result, open(out / "result.json", "w"), indent=1)
    open(out / "result.json", "a").write("\n")
    png(
        out / "cases.png",
        [st1, st2, st3],
        ["case 1 star", "case 2 pair", "case 3 pair + 20 mag Sersic"],
    )
    print(json.dumps(result["decisions"]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

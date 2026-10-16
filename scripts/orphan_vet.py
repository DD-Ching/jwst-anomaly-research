"""Per-pair vetting of orphan image pairs: bridge, surface brightness, colours, deflector light.

Implements the checks of ``results/orphan_pairs/abell2744/vetting/PLAN.md`` (written before the
run) for pairs of ``orphan_pairs.py``'s ``top_orphans.ecsv``. Small boxes of the level-3 mosaics
are read by S3 byte range (``arc_diff.read_box``); nothing is stored. Every number is ``derived``;
every threshold is an ASSUMPTION fixed in the plan.

- E1 bridge: mean of a 0.1"-wide strip along the joining segment (0.15" kept clear of each
  member) in the summed LW stack, against the robust std of the same strip at random positions.
- E2 surface brightness: aperture flux (0.15") over the PSF-deconvolved second-moment area; its
  error from re-measuring with random sky patches added (Monte Carlo). Members whose deconvolved
  size is below ``RESOLVED_FRAC`` x the PSF FWHM are unresolved: E2 undetermined.
- E3 colours: chi^2 of a constant flux ratio a/b over the bands.
- E4 deflector light: 0.2" aperture at the brightest strip point; 3-sigma limit.
- E6 detection: bands with S/N >= 3 per member.

    python scripts/orphan_vet.py --pairs results/orphan_pairs/abell2744/top_orphans.ecsv \\
        --ranks 1 2 3 5 6 7 9 10 11 13 14 --i2d jw02561-o001_t003 \\
        --out results/orphan_pairs/abell2744/vetting
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arc_diff import read_box  # noqa: E402

BANDS = ("f115w", "f150w", "f200w", "f277w", "f356w", "f410m", "f444w")
LW = ("f277w", "f356w", "f444w")
# PSF FWHM (arcsec) per band, JWST user documentation (NIRCam PSF page; ASSUMPTION: point values)
PSF_FWHM = {
    "f115w": 0.040,
    "f150w": 0.050,
    "f200w": 0.066,
    "f277w": 0.092,
    "f356w": 0.116,
    "f410m": 0.137,
    "f444w": 0.145,
}
HALF_BOX = 5.0  # arcsec
R_AP = 0.15  # member aperture radius (plan)
R_DEFL = 0.2  # deflector aperture radius (plan)
STRIP_HALF_WIDTH = 0.05  # 0.1" wide strip (plan)
CLEAR = 0.15  # strip kept this far from each member (plan)
R_MOM = 0.3  # second-moment radius (plan)
N_RANDOM = 40  # random strip placements (plan)
N_MC = 50  # Monte Carlo sky patches for E2 errors
RESOLVED_FRAC = 0.5  # deconvolved FWHM below this x PSF FWHM: unresolved (plan)
STACK_SCALE_TOL = 1e-4  # relative pixel-scale agreement for stacking LW bands (ASSUMPTION)
STACK_PIX_TOL = 0.1  # member pixel positions must agree to this (pixels) across stacked bands
SKY_PATCH_MAX = 4.0  # an MC sky patch is source-free when |pixel| < this x sky std (ASSUMPTION)
S3 = "s3://stpubdata/jwst/public/jw{p}/L3/t/o{o}/{prefix}_nircam_clear-{b}_i2d.fits"


def uri(prefix: str, band: str) -> str:
    prog, obs = prefix[2:7], prefix.split("-o")[1][:3]
    return S3.format(p=prog, o=obs, prefix=prefix, b=band)


def sky_stats(sci: np.ndarray) -> tuple[float, float]:
    """Sigma-clipped median and robust std of the box (the mosaics carry a sky pedestal)."""
    v = sci[np.isfinite(sci)]
    for _ in range(5):
        med = np.median(v)
        s = 1.4826 * np.median(np.abs(v - med))
        v = v[np.abs(v - med) < 3 * s]
    return float(np.median(v)), float(1.4826 * np.median(np.abs(v - np.median(v))))


def ap_sum(img: np.ndarray, x: float, y: float, r: float) -> float:
    yy, xx = np.indices(img.shape)
    m = (xx - x) ** 2 + (yy - y) ** 2 <= r * r
    v = img[m]
    return float(np.nansum(v)) if np.isfinite(v).mean() > 0.9 else float("nan")


def moments(img: np.ndarray, x: float, y: float, r: float) -> float:
    """Gaussian-weighted (sigma = r/2) second-moment size sigma^2 (pixels^2), iterated centre."""
    yy, xx = np.indices(img.shape)
    sw = r / 2.0
    for _ in range(5):
        w = np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sw * sw))
        f = np.where(np.isfinite(img), img, 0.0) * w
        f = np.where(((xx - x) ** 2 + (yy - y) ** 2) <= r * r, f, 0.0)
        tot = f.sum()
        if tot <= 0:
            return float("nan")
        x, y = (f * xx).sum() / tot, (f * yy).sum() / tot
    qxx = (f * (xx - x) ** 2).sum() / tot
    qyy = (f * (yy - y) ** 2).sum() / tot
    # weighted moments of a Gaussian of width s with weight width sw: s_w^2 = s^2 sw^2/(s^2+sw^2)
    sw2 = 0.5 * (qxx + qyy)
    if sw2 >= sw * sw:
        return float("inf")
    return float(sw2 * sw * sw / (sw * sw - sw2))


def strip_mask(shape, p0, p1, half_w: float, clear: float) -> np.ndarray:
    yy, xx = np.indices(shape)
    d = np.array(p1) - np.array(p0)
    L = np.hypot(*d)
    u = d / L
    rx, ry = xx - p0[0], yy - p0[1]
    t = rx * u[0] + ry * u[1]
    n = np.abs(-rx * u[1] + ry * u[0])
    return (t >= clear) & (t <= L - clear) & (n <= half_w)


def measure_pair(row, prefix: str, rng: np.random.Generator) -> dict:
    ra = 0.5 * (row["ra_a"] + row["ra_b"])
    dec = 0.5 * (row["dec_a"] + row["dec_b"])
    boxes = {}
    with ThreadPoolExecutor(7) as pool:
        for b, box in zip(
            BANDS,
            pool.map(lambda b: read_box(uri(prefix, b), ra, dec, HALF_BOX), BANDS),
            strict=False,
        ):
            boxes[b] = box
    out = {"rank": int(row["rank"]), "bands": {}}
    lw_stack, lw_bands = None, []
    for b in BANDS:
        box = boxes[b]
        sci, w, sc = box["sci"], box["wcs"], box["scale"]
        med, sd = sky_stats(sci)
        img = sci - med
        xa, ya = (float(v) for v in w.all_world2pix(row["ra_a"], row["dec_a"], 0))
        xb, yb = (float(v) for v in w.all_world2pix(row["ra_b"], row["dec_b"], 0))
        r = R_AP / sc
        fa, fb = ap_sum(img, xa, ya, r), ap_sum(img, xb, yb, r)
        # empirical aperture error: random apertures in the box
        ny, nx = img.shape
        rs = []
        for _ in range(200):
            x, y = rng.uniform(r + 1, nx - r - 1), rng.uniform(r + 1, ny - r - 1)
            if min(np.hypot(x - xa, y - ya), np.hypot(x - xb, y - yb)) > 3 * r:
                rs.append(ap_sum(img, x, y, r))
        rs = np.array([v for v in rs if np.isfinite(v)])
        err = float(1.4826 * np.median(np.abs(rs - np.median(rs))))
        out["bands"][b] = {
            "fa": fa,
            "fb": fb,
            "err": err,
            "snr_a": fa / err,
            "snr_b": fb / err,
            "scale": sc,
        }
        if b in LW:
            if lw_stack is None:
                lw_stack, lw_scale = img.copy(), sc
                lw_xy = (xa, ya, xb, yb)
                lw_bands.append(b)
            elif (
                img.shape == lw_stack.shape
                and abs(sc / lw_scale - 1) < STACK_SCALE_TOL
                and max(abs(u - v) for u, v in zip((xa, ya, xb, yb), lw_xy, strict=True))
                < STACK_PIX_TOL
            ):
                # same grid: both members fall on the same pixel coordinates (fix after review:
                # the first run compared scales to 1e-6 and stacked F277W only)
                lw_stack = lw_stack + img
                lw_bands.append(b)
    out["lw_stack_bands"] = lw_bands
    xa, ya, xb, yb = lw_xy
    sc = lw_scale
    img = lw_stack
    # E1 bridge
    m = strip_mask(img.shape, (xa, ya), (xb, yb), STRIP_HALF_WIDTH / sc, CLEAR / sc)
    bridge = float(np.nanmean(img[m])) if m.sum() else float("nan")
    ny, nx = img.shape
    vals = []
    L = np.hypot(xb - xa, yb - ya)
    for _ in range(N_RANDOM * 5):
        ang = rng.uniform(0, np.pi)
        x0, y0 = rng.uniform(L, nx - L), rng.uniform(L, ny - L)
        p1 = (x0 + L * np.cos(ang), y0 + L * np.sin(ang))
        mm = strip_mask(img.shape, (x0, y0), p1, STRIP_HALF_WIDTH / sc, CLEAR / sc)
        if mm.sum() and min(np.hypot(x0 - xa, y0 - ya), np.hypot(x0 - xb, y0 - yb)) > L + 0.5 / sc:
            vals.append(float(np.nanmean(img[mm])))
        if len(vals) >= N_RANDOM:
            break
    vals = np.array(vals)
    bmed = float(np.median(vals))
    bsd = float(1.4826 * np.median(np.abs(vals - bmed)))
    out["E1"] = {
        "strip_mean": bridge,
        "random_median": bmed,
        "random_std": bsd,
        "z": (bridge - bmed) / bsd,
        "n_random": len(vals),
    }
    # E4 deflector light: brightest strip point, 0.2" aperture
    yy, xx = np.nonzero(m)
    if len(xx):
        k = int(np.nanargmax(img[yy, xx]))
        xd, yd = float(xx[k]), float(yy[k])
    else:
        xd, yd = 0.5 * (xa + xb), 0.5 * (ya + yb)
    rd = R_DEFL / sc
    fd = ap_sum(img, xd, yd, rd)
    rs = []
    for _ in range(200):
        x, y = rng.uniform(rd + 1, nx - rd - 1), rng.uniform(rd + 1, ny - rd - 1)
        if min(np.hypot(x - xa, y - ya), np.hypot(x - xb, y - yb)) > 3 * rd:
            rs.append(ap_sum(img, x, y, rd))
    rs = np.array([v for v in rs if np.isfinite(v)])
    ed = float(1.4826 * np.median(np.abs(rs - np.median(rs))))
    # upper bound on light at the brightest strip point (it includes member wings, so it is
    # conservative): flux + 3 sigma, as the mean per stacked band, in AB mag
    pix_njy = (sc / 206264.806) ** 2 * 1e15
    ub_njy = (max(fd, 0.0) + 3 * ed) / len(lw_bands) * pix_njy
    out["E4"] = {
        "flux_stack": fd,
        "err_stack": ed,
        "x": xd,
        "y": yd,
        "upper_njy": ub_njy,
        "mag_limit_3sig": float(31.4 - 2.5 * np.log10(ub_njy)),
    }
    # E2 sizes on the LW stack (PSF: mean LW FWHM in quadrature)
    psf_fwhm = float(np.mean([PSF_FWHM[b] for b in out["lw_stack_bands"]]))
    psf_s2 = (psf_fwhm / 2.3548 / sc) ** 2
    rm = R_MOM / sc
    res = {}
    for tag, (x, y) in (("a", (xa, ya)), ("b", (xb, yb))):
        s2 = moments(img, x, y, rm)
        f = ap_sum(img, x, y, R_AP / sc)
        mc_s2, mc_f = [], []
        sky_sd = sky_stats(img)[1]
        tries = 0
        while len(mc_s2) < N_MC and tries < 50 * N_MC:
            tries += 1
            ox, oy = rng.uniform(rm + 4, nx - rm - 4), rng.uniform(rm + 4, ny - rm - 4)
            if min(np.hypot(ox - xa, oy - ya), np.hypot(ox - xb, oy - yb)) < 2 * rm:
                continue
            pert = _shifted_patch(img, ox, oy, x, y, rm)
            # source-free sky only (fix after the first run: patches with galaxies inflated the
            # errors); the stack's correlated noise is kept
            if not np.any(pert) or np.nanmax(np.abs(pert)) > SKY_PATCH_MAX * sky_sd:
                continue
            mc_s2.append(moments(img + pert, x, y, rm))
            mc_f.append(ap_sum(img + pert, x, y, R_AP / sc))
        mc_s2 = np.array([v for v in mc_s2 if np.isfinite(v)])
        dec2 = s2 - psf_s2
        fwhm_dec = 2.3548 * np.sqrt(max(dec2, 0.0)) * sc
        res[tag] = {
            "flux": f,
            "flux_sd": float(np.std(mc_f)),
            "sigma2_px": s2,
            "sigma2_sd": float(np.std(mc_s2)) if len(mc_s2) > 3 else float("nan"),
            "psf_sigma2_px": psf_s2,
            "fwhm_deconv_arcsec": fwhm_dec,
            "resolved": bool(np.isfinite(s2) and fwhm_dec >= RESOLVED_FRAC * psf_fwhm),
        }
    out["E2"] = res
    if res["a"]["resolved"] and res["b"]["resolved"]:
        sba = res["a"]["flux"] / (res["a"]["sigma2_px"] - psf_s2)
        sbb = res["b"]["flux"] / (res["b"]["sigma2_px"] - psf_s2)
        # relative errors in quadrature (flux and deconvolved area)
        rel = np.sqrt(
            sum(
                (res[t]["flux_sd"] / res[t]["flux"]) ** 2
                + (res[t]["sigma2_sd"] / (res[t]["sigma2_px"] - psf_s2)) ** 2
                for t in "ab"
            )
        )
        ratio = sba / sbb
        out["E2"]["sb_ratio"] = float(ratio)
        out["E2"]["sb_ratio_rel_err"] = float(rel)
        out["E2"]["sb_z"] = float(np.log(ratio) / rel)
    # E3 colours: chi^2 of a constant ratio over bands with both fluxes finite
    rb = [(v["fa"], v["fb"], v["err"]) for v in out["bands"].values()]
    rb = [t for t in rb if np.isfinite(t[0]) and np.isfinite(t[1]) and t[1] > 3 * t[2]]
    if len(rb) >= 3:
        fa, fb, e = (np.array(x) for x in zip(*rb, strict=False))
        r = fa / fb
        er = np.abs(r) * np.hypot(e / np.abs(fa).clip(1e-9), e / fb)
        wmean = np.sum(r / er**2) / np.sum(1 / er**2)
        chi2 = float(np.sum(((r - wmean) / er) ** 2))
        from scipy.stats import chi2 as c2

        out["E3"] = {
            "n_bands": len(r),
            "ratio_mean": float(wmean),
            "chi2": chi2,
            "p": float(c2.sf(chi2, len(r) - 1)),
        }
    else:
        out["E3"] = {"n_bands": len(rb), "p": None}
    out["E6_detected_bands"] = {
        "a": int(sum(v["snr_a"] >= 3 for v in out["bands"].values())),
        "b": int(sum(v["snr_b"] >= 3 for v in out["bands"].values())),
    }
    out["_stamp"] = (img, (xa, ya, xb, yb, xd, yd), sc)
    return out


def _shifted_patch(img, ox, oy, x, y, rm):
    """A sky patch from (ox, oy) placed at (x, y) (integer shift), zero elsewhere."""
    out = np.zeros_like(img)
    h = int(np.ceil(rm)) + 2
    xi, yi, oxi, oyi = int(round(x)), int(round(y)), int(round(ox)), int(round(oy))
    src = img[oyi - h : oyi + h + 1, oxi - h : oxi + h + 1]
    dst = out[yi - h : yi + h + 1, xi - h : xi + h + 1]
    if src.shape == dst.shape:
        out[yi - h : yi + h + 1, xi - h : xi + h + 1] = np.nan_to_num(src)
    return out


def sheet(results: list[dict], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(2.0 * n, 2.3), squeeze=False)
    for ax, r in zip(axes[0], results, strict=True):
        img, (xa, ya, xb, yb, xd, yd), sc = r["_stamp"]
        cx, cy = 0.5 * (xa + xb), 0.5 * (ya + yb)
        h = int(1.6 / sc)
        sub = img[int(cy) - h : int(cy) + h, int(cx) - h : int(cx) + h]
        lo, hi = np.nanpercentile(sub, [5, 99.7])
        ax.imshow(sub, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        off = (int(cx) - h, int(cy) - h)
        ax.plot([xa - off[0]], [ya - off[1]], "o", mfc="none", mec="r", ms=10)
        ax.plot([xb - off[0]], [yb - off[1]], "o", mfc="none", mec="c", ms=10)
        ax.plot([xd - off[0]], [yd - off[1]], "+", color="y", ms=8)
        ax.set_title(f'#{r["rank"]} LW stack 3.2"', fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# PLAN addendum: CATS v4.1 pairs <= 4" apart with known z (fixed before measuring)
CONTROL_PAIRS = (
    ("39.1", "39.2"),
    ("37.1", "37.2"),
    ("33.1", "33.2"),
    ("4.4", "4.5"),
    ("62.1", "62.2"),
    ("40.1", "40.2"),
    ("8.1", "8.2"),
    ("3.1", "3.2"),
    ("61.1", "61.2"),
    ("10.1", "10.2"),
    ("134.1", "134.2"),
)


def e2_valid(sb_z: list[float]) -> dict:
    """PLAN addendum: E2 is usable only if the known equal-SB controls behave (all conditions)."""
    z = np.asarray([v for v in sb_z if v is not None and np.isfinite(v)], float)
    n = len(z)
    if n < 3:
        return {"valid": False, "n": n, "reason": "fewer than 3 measured controls"}
    med = float(np.median(z))
    rsd = float(1.4826 * np.median(np.abs(z - med)))
    tail = float(np.mean(np.abs(z) > 3))
    ok = abs(med) <= 1.0 and rsd <= 1.5 and tail <= 0.10
    return {"valid": bool(ok), "n": n, "median": med, "robust_std": rsd, "frac_abs_gt3": tail}


def classify(m: dict, e2_ok: bool, artifact: bool, dz_spec_ok: bool | None, sis_mag: float) -> dict:
    """PLAN addendum classification. ``m``: one measurements.json entry; ``dz_spec_ok`` False when
    two public spectroscopic redshifts differ by > 0.01 (1 + z), None when unavailable;
    ``sis_mag``: SIS-scaled deflector magnitude. Returns the class and the conditions that fired."""
    e1 = m["E1"]["z"] >= 3
    sbz = m["E2"].get("sb_z")
    e2_mis = bool(e2_ok and sbz is not None and abs(sbz) > 3)
    e3p = m["E3"].get("p")
    e3_mis = e3p is not None and e3p < 0.01
    e5_mis = dz_spec_ok is False
    fired = [
        k
        for k, v in (
            ("E1_bridge", e1),
            ("E2_sb_mismatch", e2_mis),
            ("E3_colour_mismatch", e3_mis),
            ("E5_specz_mismatch", e5_mis),
            ("E6_artifact", artifact),
        )
        if v
    ]
    if fired:
        return {"class": "ordinary", "why": fired}
    res = m["E2"]["a"]["resolved"] and m["E2"]["b"]["resolved"]
    det = min(m["E6_detected_bands"].values()) >= 3
    lacking = [
        k
        for k, v in (("E2_not_valid", not e2_ok), ("unresolved", not res), ("detection", not det))
        if v
    ]
    if lacking:
        return {"class": "insufficient", "why": lacking}
    defl_ok = m["E4"].get("mag_limit_3sig", -np.inf) >= sis_mag + 3
    if sbz is not None and abs(sbz) <= 3 and e3p is not None and e3p >= 0.01 and defl_ok:
        return {"class": "anomaly", "why": ["E2_equal_sb", "E3_equal_colour", "E4_dark"]}
    return {"class": "insufficient", "why": ["fallback"]}


VIZIER_SPEC = "J/ApJ/982/51"  # Price+25 UNCOVER spectroscopic redshifts
VIZIER_SPS = "J/ApJS/270/12"  # Wang+24 SPS redshifts (z50, zsp)
VIZIER_DR3 = "J/ApJS/270/7"  # Weaver+24 UNCOVER DR3 photometry (fluxes in 10 nJy)
CONE_ARCSEC = 1.5
SAMEZ_DZ = 0.15  # same-z neighbour: |z50 - z_pair| <= this x (1 + z_pair) (plan addendum)
ARTIFACT_RANKS = (1,)  # #269 RESULT.md sheet pass: #1 artefact-affected (visual flag)


def _cone(cat: str, ra: float, dec: float) -> list[dict]:
    """Rows of a VizieR catalogue within CONE_ARCSEC (an empty list when none)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    c = SkyCoord(ra * u.deg, dec * u.deg)
    import time

    for attempt in range(4):  # VizieR resets connections now and then
        try:
            res = Vizier(columns=["**"], row_limit=-1).query_region(
                c, radius=CONE_ARCSEC * u.arcsec, catalog=[cat]
            )
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))
    rows = []
    for tb in res:
        for row in tb:
            rc = SkyCoord(float(row["RAJ2000"]) * u.deg, float(row["DEJ2000"]) * u.deg)
            d = {"sep_arcsec": float(rc.separation(c).arcsec), "ra": float(rc.ra.deg)}
            d["dec"] = float(rc.dec.deg)
            for k in ("zsp", "z50", "z16", "z84", "zph", "Seq", "f356w", "e_f356w"):
                if k in tb.colnames and not np.ma.is_masked(row[k]):
                    d[k] = float(row[k])
            rows.append(d)
    return rows


def classify_run(out: Path, pairs: Path, ranks: list[int]) -> dict:
    """Plan addendum: E2 validity from the controls, E5/E6 cone searches, E4 vs SIS, classes.

    Reads ``measurements.json`` and ``controls.json`` in ``out``; writes ``e2_validation.json``,
    ``vizier.json`` and ``classification.json`` there (every query result stored, empty or not)."""
    meas = {m["rank"]: m for m in json.loads((out / "measurements.json").read_text())}
    ctrl = json.loads((out / "controls.json").read_text())
    sbz = [
        c["E2"].get("sb_z") for c in ctrl if c["E2"]["a"]["resolved"] and c["E2"]["b"]["resolved"]
    ]
    val = e2_valid(sbz)
    val["control_min_detected_bands"] = [min(c["E6_detected_bands"].values()) for c in ctrl]
    val["control_e1_n_random"] = [c["E1"]["n_random"] for c in ctrl]
    val["control_e3_n_bands"] = [c["E3"]["n_bands"] for c in ctrl]
    val["sb_z"] = sbz
    (out / "e2_validation.json").write_text(json.dumps(val, indent=1, default=float) + "\n")
    t = Table.read(pairs)
    viz, cls = {}, {}
    for r in t:
        k = int(r["rank"])
        if k not in ranks:
            continue
        rec = {}
        for tag in "ab":
            ra, dec = float(r[f"ra_{tag}"]), float(r[f"dec_{tag}"])
            rec[tag] = {c: _cone(c, ra, dec) for c in (VIZIER_SPEC, VIZIER_SPS, VIZIER_DR3)}
        mra, mdec = 0.5 * (r["ra_a"] + r["ra_b"]), 0.5 * (r["dec_a"] + r["dec_b"])
        rec["midpoint"] = {VIZIER_SPS: _cone(VIZIER_SPS, float(mra), float(mdec))}
        viz[k] = rec
        zp = 0.5 * (float(r["z_a"]) + float(r["z_b"]))
        nb = 0
        for row in rec["midpoint"][VIZIER_SPS]:
            da = 3600 * np.hypot(
                (row["ra"] - r["ra_a"]) * np.cos(np.radians(dec)), row["dec"] - r["dec_a"]
            )
            db = 3600 * np.hypot(
                (row["ra"] - r["ra_b"]) * np.cos(np.radians(dec)), row["dec"] - r["dec_b"]
            )
            z = row.get("z50")
            if da > 0.2 and db > 0.2 and z is not None and abs(z - zp) <= SAMEZ_DZ * (1 + zp):
                nb += 1
        spec = [
            min(rec[tag][VIZIER_SPEC], key=lambda x: x["sep_arcsec"]).get("zsp")
            if rec[tag][VIZIER_SPEC]
            else None
            for tag in "ab"
        ]
        if spec[0] is not None and spec[1] is not None:
            dz_ok = abs(spec[0] - spec[1]) <= 0.01 * (1 + 0.5 * (spec[0] + spec[1]))
        else:
            dz_ok = None
        sis = float(r["sis_sigma"])
        sis_mag = float(19.0 - 10.0 * np.log10(sis / 180.0))  # D-048 F160W FJ scaling
        m = meas[k]
        c = classify(m, val["valid"], k in ARTIFACT_RANKS, dz_ok, sis_mag)
        cls[k] = {
            **c,
            "E1_z": m["E1"]["z"],
            "E1_n_random": m["E1"]["n_random"],
            "E2_sb_z": m["E2"].get("sb_z"),
            "E3_p": m["E3"].get("p"),
            "E4_mag_upper": m["E4"]["mag_limit_3sig"],
            "sis_sigma": sis,
            "sis_mag_F160W_scaling": sis_mag,
            "E5_specz": spec,
            "E6_samez_neighbours": nb,
            "E6_detected_bands": m["E6_detected_bands"],
            "lw_stack_bands": m["lw_stack_bands"],
            "artifact_flag_visual_269": k in ARTIFACT_RANKS,
        }
    (out / "vizier.json").write_text(json.dumps(viz, indent=1) + "\n")
    summary = {"e2_valid": val["valid"], "classification": cls}
    (out / "classification.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pairs", type=Path)
    ap.add_argument("--ranks", type=int, nargs="+")
    ap.add_argument(
        "--controls",
        help="measure CONTROL_PAIRS of this lens_consistency image list (e.g. abell2744-cats)",
    )
    ap.add_argument(
        "--classify",
        action="store_true",
        help="classify from measurements.json and controls.json in --out (no image reads)",
    )
    ap.add_argument("--i2d", help="level-3 product prefix, e.g. jw02561-o001_t003")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.classify:
        summary = classify_run(args.out, args.pairs, args.ranks)
        for k, v in summary["classification"].items():
            print(k, v["class"], v["why"])
        return 0
    if args.controls:
        import lens_consistency as lc

        il, _ = lc.image_list(args.controls, {}, None)
        lc.shift_images(args.controls, il)
        pos = {
            str(i): (float(r), float(d))
            for i, r, d in zip(il["image_id"], il["ra"], il["dec"], strict=True)
        }
        rows = [
            {
                "rank": k + 1,
                "ra_a": pos[a][0],
                "dec_a": pos[a][1],
                "ra_b": pos[b][0],
                "dec_b": pos[b][1],
                "pair": f"{a}/{b}",
            }
            for k, (a, b) in enumerate(CONTROL_PAIRS)
        ]
    else:
        t = Table.read(args.pairs)
        rows = [r for r in t if int(r["rank"]) in args.ranks]
    rng = np.random.default_rng(args.seed)
    results = [measure_pair(r, args.i2d, rng) for r in rows]
    args.out.mkdir(parents=True, exist_ok=True)
    for r, row in zip(results, rows, strict=True):
        if isinstance(row, dict) and "pair" in row:
            r["pair"] = row["pair"]
    name = "controls" if args.controls else "measurements"
    sheet(results, args.out / f"stamps_lw{'_controls' if args.controls else ''}.png")
    clean = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    (args.out / f"{name}.json").write_text(json.dumps(clean, indent=1, default=float) + "\n")
    for r in clean:
        e2 = r["E2"]
        print(
            r["rank"],
            "E1 z {:.1f}".format(r["E1"]["z"]),
            "E2 res",
            e2["a"]["resolved"],
            e2["b"]["resolved"],
            "fwhm {:.3f} {:.3f}".format(
                e2["a"]["fwhm_deconv_arcsec"], e2["b"]["fwhm_deconv_arcsec"]
            ),
            "sb_z",
            round(e2.get("sb_z", float("nan")), 2),
            "E3 p",
            r["E3"].get("p"),
            "det",
            r["E6_detected_bands"],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

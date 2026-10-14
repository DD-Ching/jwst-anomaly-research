"""Observational parity test: does image B look like A or like A's mirror image, as the model says?

For two images A, B of one source, a lens model predicts the local map between them: with
``J = d beta / d theta`` at each image, ``J_A dtheta_A = J_B dtheta_B``, so offsets around B map
to offsets around A by ``M_pred = J_A^-1 J_B``. ``sign(det M_pred) = sign(mu_A mu_B)`` is the
relative parity: negative means B is a mirror image of A.

The measurement fits B's stamp with A's stamp resampled through a free affine map
``M = R(t1) diag(s1, h s2) R(t2)`` plus a shift, a flux scale and a constant, separately in the two
handedness families ``h = +1`` and ``h = -1``. ``dchi2 = chi2(wrong family) - chi2(predicted
family)``: positive means the data prefer the model's parity. Both stamps come from the same
band, so their PSFs match to first order (ASSUMPTION; the affine resampling also distorts A's PSF).

Per-pair power (simulated): B is replaced by each family's best-fit model plus noise, and the fit
is repeated. A pair is ``testable`` only when the predicted-family synthetic gives
``dchi2 > DCHI2_MIN`` and the wrong-family synthetic gives ``dchi2 < -DCHI2_MIN`` (a symmetric
or unresolved source gives about 0 for both). Classes: ``consistent`` / ``parity_residual``
(``dchi2`` beyond +-DCHI2_MIN) / ``undetermined``; ``near_critical`` when either image has
``|mu| > MAX_ABS_MU`` (its parity is model-dependent; ASSUMPTION). Image positions are observed,
``mu`` and ``M_pred`` model_prediction, fitted maps and ``dchi2`` derived, synthetics simulated.

    OMP_NUM_THREADS=1 python scripts/parity_test.py --model abell2744-bergamini23 \
        --image s3://stpubdata/jwst/public/jw02561/L3/t/o001/jw02561-o001_t003_nircam_clear-f200w_i2d.fits
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy import ndimage, optimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel, paths  # noqa: E402

FIT_RADIUS_ARCSEC = 0.6  # pixels of B used in the fit (ASSUMPTION)
TEMPLATE_HALF_ARCSEC = 1.6  # A's stamp (B uses its central FIT_RADIUS_ARCSEC); byte-range cost
N_READERS = 8  # concurrent S3 byte-range readers
HALO_FILTER_ARCSEC = 1.0  # median filter subtracted from stamps (ASSUMPTION; also trims long arcs)
SCALE_RANGE = (0.2, 5.0)  # singular values of M
MAX_SHIFT_ARCSEC = 0.3  # catalogue positions vs light peaks
DCHI2_MIN = 100.0  # family preference (ASSUMPTION; symmetric-source null spread ~ +-50, ~1200 px)
MAX_ABS_MU = 10.0  # beyond it the parity is model-dependent (ASSUMPTION)
MIN_PEAK_SNR = 10.0  # the template's peak pixel S/N
N_ANGLES = 24  # start grid per rotation angle (0..180 deg)
N_REFINE = 8  # grid starts refined by least squares
N_MIRROR_AXES = 72  # reflection axes (2.5 deg steps) tried on the other family's best map
SYN_SMOOTH_PIX = 2.0  # Gaussian sigma (pixels) of the template copy the synthetics use
SEG_SNR = 2.0  # segmentation threshold on the 1-px-smoothed S/N map (ASSUMPTION)
SEG_GROW_PIX = 2  # dilation of each segment (wings below the threshold)
PSF_SIGMA_ARCSEC = 0.028  # F200W PSF core, FWHM 0.066" / 2.355 (ASSUMPTION; other bands: scale it)
MIN_SIZE_RATIO = 2.0  # major-axis sigma / PSF sigma; below it PSF-dominated (ASSUMPTION)


def _rot(t: float) -> np.ndarray:
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s], [s, c]])


def affine(p: np.ndarray, h: int) -> tuple[np.ndarray, np.ndarray]:
    """``M`` and the shift from ``p = (t1, t2, ln s1, ln s2, dx, dy)`` in family ``h``."""
    m = _rot(p[0]) @ np.diag([np.exp(p[2]), h * np.exp(p[3])]) @ _rot(p[1])
    return m, np.array(p[4:6])


def decompose(m: np.ndarray) -> tuple[float, float, float, float]:
    """``(t1, t2, ln s1, ln s2)`` with ``m = R(t1) diag(s1, sign(det) s2) R(t2)``."""
    u, s, vt = np.linalg.svd(m)
    if np.linalg.det(u) < 0:  # move the reflection into the second singular value
        u[:, 1] *= -1
    if np.linalg.det(vt) < 0:
        vt[1, :] *= -1
    return (
        float(np.arctan2(u[1, 0], u[0, 0])),
        float(np.arctan2(vt[1, 0], vt[0, 0])),
        float(np.log(s[0])),
        float(np.log(s[1])),
    )


class Template:
    """Bilinear sampler of a stamp at model-frame offsets (arcsec)."""

    def __init__(self, sci, err, xx, yy):
        self.sci = np.nan_to_num(np.asarray(sci, float))
        self.var = np.nan_to_num(np.asarray(err, float) ** 2, nan=np.inf)
        ny, nx = self.sci.shape
        jj, ii = np.mgrid[0:ny, 0:nx]
        a = np.column_stack([np.ones(ii.size), ii.ravel(), jj.ravel()])
        cx = np.linalg.lstsq(a, np.ravel(xx), rcond=None)[0]
        cy = np.linalg.lstsq(a, np.ravel(yy), rcond=None)[0]
        self.o = np.array([cx[0], cy[0]])
        self.inv = np.linalg.inv(np.array([[cx[1], cx[2]], [cy[1], cy[2]]]))
        self.shape = self.sci.shape

    def sample(self, x, y):
        ij = self.inv @ np.vstack([np.ravel(x) - self.o[0], np.ravel(y) - self.o[1]])
        coords = ij[::-1]  # (row = j, col = i)
        inside = (
            (coords[0] >= 0)
            & (coords[0] <= self.shape[0] - 1)
            & (coords[1] >= 0)
            & (coords[1] <= self.shape[1] - 1)
        )
        val = ndimage.map_coordinates(self.sci, coords, order=1, mode="nearest")
        var = ndimage.map_coordinates(self.var, coords, order=1, mode="nearest")
        return val, var, inside


def _model(tmpl: Template, p, h, bx, by, b, bvar):
    m, sh = affine(p, h)
    d = m @ np.vstack([bx, by])
    val, var, inside = tmpl.sample(d[0] + sh[0], d[1] + sh[1])
    w = 1.0 / bvar
    a = np.column_stack([val, np.ones_like(val)]) * np.sqrt(w)[:, None]
    f, c = np.linalg.lstsq(a, b * np.sqrt(w), rcond=None)[0]
    sig2 = bvar + f * f * var
    pred = f * val + c
    return pred, sig2, inside, f


def _resid(p, tmpl, h, bx, by, b, bvar):
    pred, sig2, inside, _ = _model(tmpl, p, h, bx, by, b, bvar)
    r = (b - pred) / np.sqrt(sig2)
    r[~inside] = 10.0  # mapped outside A's stamp: penalised
    return r


def mirror_starts(p: np.ndarray, h: int) -> list[np.ndarray]:
    """Parameter vectors in family ``-h`` for the map of ``p`` composed with a reflection about
    each of ``N_MIRROR_AXES`` axes through B's centre. A source mirror-symmetric about some axis
    is fitted equally well by one of them, so the cross-family start is within half a step of a
    tie; a coarse axis grid on a sharp, bright core leaves a large spurious preference."""
    m, sh = affine(p, h)
    out = []
    for phi in np.arange(N_MIRROR_AXES) * np.pi / N_MIRROR_AXES:
        refl = _rot(phi) @ np.diag([1.0, -1.0]) @ _rot(-phi)
        out.append(np.array([*decompose(m @ refl), *sh]))
    return out


def fit_family(
    tmpl: Template, h: int, bx, by, b, bvar, m_pred: np.ndarray, seeds=(), exact=()
) -> dict:
    """Best affine fit in handedness family ``h`` from a grid of rotation starts at the predicted
    scales, plus each seed parameter vector (e.g. the other family's best fit) with its second
    rotation scanned over the same grid (its mirror axis), plus each ``exact`` start as given."""
    t1, t2, l1, l2 = decompose(m_pred)
    lo_s, hi_s = np.log(SCALE_RANGE[0]), np.log(SCALE_RANGE[1])
    l1, l2 = (float(np.clip(v, lo_s + 1e-3, hi_s - 1e-3)) for v in (l1, l2))
    lo = [-np.inf, -np.inf, lo_s, lo_s, -MAX_SHIFT_ARCSEC, -MAX_SHIFT_ARCSEC]
    hi = [np.inf, np.inf, hi_s, hi_s, MAX_SHIFT_ARCSEC, MAX_SHIFT_ARCSEC]
    angles = np.arange(N_ANGLES) * np.pi / N_ANGLES
    starts = []
    for a1 in angles:  # cheap grid over both rotations, then refine the best few
        for a2 in angles:
            p0 = np.array([t1 + a1, t2 + a2, l1, l2, 0.0, 0.0])
            r = _resid(p0, tmpl, h, bx, by, b, bvar)
            starts.append((float(r @ r), p0))
    for sp in seeds:
        for a in angles:
            p0 = np.array(sp, float).copy()
            p0[1] += a
            p0[2:4] = np.clip(p0[2:4], lo_s + 1e-3, hi_s - 1e-3)
            p0[4:6] = np.clip(p0[4:6], -MAX_SHIFT_ARCSEC + 1e-6, MAX_SHIFT_ARCSEC - 1e-6)
            r = _resid(p0, tmpl, h, bx, by, b, bvar)
            starts.append((float(r @ r), p0))
    for p0 in exact:
        p0 = np.array(p0, float).copy()
        p0[2:4] = np.clip(p0[2:4], lo_s + 1e-3, hi_s - 1e-3)
        p0[4:6] = np.clip(p0[4:6], -MAX_SHIFT_ARCSEC + 1e-6, MAX_SHIFT_ARCSEC - 1e-6)
        r = _resid(p0, tmpl, h, bx, by, b, bvar)
        starts.append((float(r @ r), p0))
    starts.sort(key=lambda t: t[0])
    best = None
    for _, p0 in starts[:N_REFINE]:
        res = optimize.least_squares(
            _resid, p0, bounds=(lo, hi), args=(tmpl, h, bx, by, b, bvar), max_nfev=80
        )
        chi2 = float(2 * res.cost)
        if best is None or chi2 < best["chi2"]:
            best = {"chi2": chi2, "p": res.x}
    m, sh = affine(best["p"], h)
    best["m"] = m
    best["shift"] = sh
    return best


def compare(tmpl: Template, bx, by, b, bvar, m_pred: np.ndarray, seeds=()) -> dict:
    """Fit both families, then refit each seeded by the other's best fit (and ``seeds``), so a
    local minimum in one family cannot fake a preference; ``dchi2 = chi2(wrong) - chi2(pred)``."""
    h_pred = 1 if np.linalg.det(m_pred) > 0 else -1
    m_alt = m_pred @ np.diag([1.0, -1.0])  # start the wrong family from a mirrored prediction
    fp = fit_family(tmpl, h_pred, bx, by, b, bvar, m_pred, seeds)
    fw = fit_family(tmpl, -h_pred, bx, by, b, bvar, m_alt, seeds)
    fp2 = fit_family(
        tmpl, h_pred, bx, by, b, bvar, m_pred, [fw["p"]], mirror_starts(fw["p"], -h_pred)
    )
    fw2 = fit_family(
        tmpl, -h_pred, bx, by, b, bvar, m_alt, [fp["p"]], mirror_starts(fp["p"], h_pred)
    )
    fp = fp2 if fp2["chi2"] < fp["chi2"] else fp
    fw = fw2 if fw2["chi2"] < fw["chi2"] else fw
    return {"h_pred": h_pred, "pred": fp, "wrong": fw, "dchi2": fw["chi2"] - fp["chi2"]}


def pair_test(tmpl: Template, bx, by, b, bvar, m_pred: np.ndarray, seed: int = 0) -> dict:
    """Real comparison plus the two synthetic power checks."""
    rng = np.random.default_rng(seed)
    real = compare(tmpl, bx, by, b, bvar, m_pred)
    # The synthetics are drawn from a smoothed copy of A: drawn from A itself, they would carry A's
    # own noise pattern, which only the generating geometry reproduces (spurious power on a
    # symmetric source). Smoothing blurs real structure too, so the power estimate is conservative.
    smooth = Template.__new__(Template)
    smooth.__dict__.update(tmpl.__dict__)
    smooth.sci = ndimage.gaussian_filter(tmpl.sci, SYN_SMOOTH_PIX)
    syn = {}
    for key in ("pred", "wrong"):
        fit = real[key]
        h = real["h_pred"] if key == "pred" else -real["h_pred"]
        model, sig2, _, _ = _model(smooth, fit["p"], h, bx, by, b, bvar)
        sb = model + rng.normal(0.0, np.sqrt(sig2))
        seeds = [real["pred"]["p"], real["wrong"]["p"]]  # includes the generating parameters
        syn[key] = compare(tmpl, bx, by, sb, bvar, m_pred, seeds)["dchi2"]
    testable = syn["pred"] > DCHI2_MIN and syn["wrong"] < -DCHI2_MIN
    d = real["dchi2"]
    if not testable:
        cls = "untestable"
    elif d > DCHI2_MIN:
        cls = "consistent"
    elif d < -DCHI2_MIN:
        cls = "parity_residual"
    else:
        cls = "undetermined"
    return {
        "class": cls,
        "dchi2": d,
        "dchi2_syn_pred": syn["pred"],
        "dchi2_syn_wrong": syn["wrong"],
        "h_pred": real["h_pred"],
        "n_pix": int(b.size),
        "chi2_pred": real["pred"]["chi2"],
        "chi2_wrong": real["wrong"]["chi2"],
        "m_fit_pred": real["pred"]["m"].round(3).tolist(),
        "m_pred": np.asarray(m_pred).round(3).tolist(),
    }


def _fit_pixels(sci, err, xx, yy, other=None):
    sel = (np.hypot(xx, yy) <= FIT_RADIUS_ARCSEC) & np.isfinite(sci) & np.isfinite(err) & (err > 0)
    if other is not None:  # pixels of neighbouring sources are not B's light
        sel &= ~other
    return xx[sel], yy[sel], sci[sel], err[sel] ** 2


def segment(sci, err, xx, yy) -> tuple[np.ndarray, np.ndarray]:
    """``(own, other)`` masks: the connected S/N > ``SEG_SNR`` component holding the brightest
    pixel within 0.15" of the stamp centre (the source), and every other component, each grown by
    ``SEG_GROW_PIX``; ``other`` excludes ``own``. Derived from the image alone (no catalogue
    segmentation map, which would need the full mosaic)."""
    ok = np.isfinite(sci) & np.isfinite(err) & (err > 0)
    snr = np.where(ok, np.nan_to_num(sci) / np.where(ok, err, 1.0), 0.0)
    det = ndimage.gaussian_filter(snr, 1.0) > SEG_SNR
    lab, _ = ndimage.label(det)
    core = (np.hypot(xx, yy) <= 0.15) & ok
    own = np.zeros(sci.shape, bool)
    if core.any():
        cand = np.where(core & (lab > 0), snr, -np.inf)
        if np.isfinite(cand).any():
            own = lab == lab.flat[int(np.argmax(cand))]
    grow = np.ones((2 * SEG_GROW_PIX + 1,) * 2, bool)
    own = ndimage.binary_dilation(own, grow)
    other = ndimage.binary_dilation((lab > 0) & ~own, grow) & ~own
    return own, other


def size_ratio(sci, xx, yy, own) -> float:
    """Major-axis sigma of the source's light (flux-weighted second moments inside ``own``,
    positive pixels only) over ``PSF_SIGMA_ARCSEC``; ~1 for a point source."""
    w = np.where(own & np.isfinite(sci), np.clip(np.nan_to_num(sci), 0, None), 0.0)
    if w.sum() <= 0:
        return np.nan
    mx, my = (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    dx, dy = xx - mx, yy - my
    cov = (
        np.array(
            [[(w * dx * dx).sum(), (w * dx * dy).sum()], [(w * dx * dy).sum(), (w * dy * dy).sum()]]
        )
        / w.sum()
    )
    return float(np.sqrt(np.linalg.eigvalsh(cov)[-1]) / PSF_SIGMA_ARCSEC)


def read_stamps(
    image_uri: str, images, want: list[int], halo_filter: float = 0.0, mask: bool = True
) -> dict:
    """``{k: (sci, err, xx, yy, peak_snr, own, other)}`` for images ``want`` (``own``/``other``
    from ``segment`` when ``mask``, else None), by ``N_READERS`` threads, each with its own
    byte-range file handle (I/O-bound; one handle is not thread-safe)."""

    def work(chunk):
        stamp, close = lc.image_stamper(image_uri, TEMPLATE_HALF_ARCSEC)
        out = {}
        try:
            for k in chunk:
                sci, err, xx, yy = stamp(float(images["ra"][k]), float(images["dec"][k]))
                if not np.isfinite(sci).any():
                    continue
                if halo_filter > 0:  # ICL / member gradients would drive the fit
                    sci = lc._median_subtract(sci, xx, yy, halo_filter)
                c = (np.hypot(xx, yy) <= 0.15) & np.isfinite(sci) & (err > 0)
                snr = float(np.max(sci[c] / err[c])) if c.any() else np.nan
                own, other = segment(sci, err, xx, yy) if mask else (None, None)
                out[k] = (sci, err, xx, yy, snr, own, other)
        finally:
            close()
        return out

    chunks = [want[i::N_READERS] for i in range(N_READERS) if want[i::N_READERS]]
    stamps: dict = {}
    with ThreadPoolExecutor(N_READERS) as ex:
        for part in ex.map(work, chunks):
            stamps.update(part)
    return stamps


def _sxy(st):
    return st[0], st[2], st[3], st[5]


def run(
    model_name: str, image_uri: str, halo_filter: float = HALO_FILTER_ARCSEC, mask: bool = True
) -> list[dict]:
    model, files, par = lc.load_model(model_name)
    images, z_m_limit = lc.image_list(model_name, files, par)
    lc.apply_frame_offset(model_name, model, images)
    z = lensmodel.image_redshifts(images, z_m_limit)
    ok = np.isfinite(z)
    images, z = images[ok], z[ok]
    x, y = model.to_frame(np.asarray(images["ra"], float), np.asarray(images["dec"], float))
    lm = model.lens_map(x, y, z)
    jac = np.stack([[lm["a11"], lm["a12"]], [lm["a21"], lm["a22"]]]).transpose(2, 0, 1)
    mu = 1.0 / np.linalg.det(jac)
    systems = np.asarray(images["system"]).astype(str)
    usable = np.abs(mu) <= MAX_ABS_MU
    # read only images in systems with >= 2 usable images (the rest are near_critical)
    want = [
        k for k in range(len(images)) if usable[k] and (usable & (systems == systems[k])).sum() > 1
    ]
    stamps = read_stamps(image_uri, images, want, halo_filter, mask)
    rows: list[dict] = []
    jobs: list[tuple[int, tuple]] = []
    for s in sorted(set(systems)):
        members = list(np.flatnonzero(systems == s))
        if len(members) < 2:
            continue
        idx = [k for k in members if k in stamps and np.isfinite(stamps[k][4])]
        ka = max(idx, key=lambda k: stamps[k][4]) if idx else None  # template: highest peak S/N
        for kb in members:
            if kb == ka:
                continue
            row = {
                "system": s,
                "image_a": str(images["image_id"][ka]) if ka is not None else None,
                "image_b": str(images["image_id"][kb]),
                "z": float(z[kb]),
                "mu_a": float(mu[ka]) if ka is not None else None,
                "mu_b": float(mu[kb]),
                "snr_a": stamps[ka][4] if ka is not None else None,
                "snr_b": stamps[kb][4] if kb in stamps else None,
            }
            if ka is None or kb not in stamps:
                row["class"] = "near_critical" if not usable[kb] else "off_image"
            elif stamps[ka][4] < MIN_PEAK_SNR:
                row["class"] = "faint"
            elif mask and not (row.setdefault("size_ratio_a", size_ratio(*_sxy(stamps[ka]))) > 0):
                row["class"] = "no_segment"  # nothing detected at A's centre
            elif mask and row["size_ratio_a"] < MIN_SIZE_RATIO:
                row["class"] = "unresolved"  # the unlensed PSF would set the handedness
            else:
                sci_a, err_a, xa, ya, _, _, other_a = stamps[ka]
                if mask:  # neighbours of A must not be mapped into B's disc
                    sci_a = np.where(other_a, 0.0, sci_a)
                tmpl = Template(sci_a, err_a, xa, ya)
                bx, by, b, bvar = _fit_pixels(*stamps[kb][:4], stamps[kb][6])
                m_pred = np.linalg.inv(jac[ka]) @ jac[kb]
                jobs.append((len(rows), (tmpl, bx, by, b, bvar, m_pred)))
            rows.append(row)
    # CPU-bound fits: one process per core (run with OMP_NUM_THREADS=1; scripts/CLAUDE.md topology)
    with ProcessPoolExecutor(os.cpu_count() or 1) as ex:
        for (i, _), res in zip(jobs, ex.map(_pair_job, [a for _, a in jobs]), strict=True):
            rows[i].update(res)
    for row in rows:
        print(json.dumps({k: v for k, v in row.items() if not k.startswith("m_")}), flush=True)
    return rows


def _pair_job(args) -> dict:
    return pair_test(*args)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", required=True, choices=sorted(lc.MODELS))
    ap.add_argument("--image", required=True, help="i2d path or s3:// URI (one band)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument(
        "--halo-filter-arcsec",
        type=float,
        default=HALO_FILTER_ARCSEC,
        help="median-filter size subtracted from every stamp (0: none)",
    )
    ap.add_argument(
        "--no-mask",
        action="store_true",
        help="no segmentation mask and no resolution gate (the first, unvalidated version)",
    )
    args = ap.parse_args(argv)
    rows = run(args.model, args.image, args.halo_filter_arcsec, not args.no_mask)
    out = args.out or paths.repo_root() / "results" / "parity" / f"{args.model}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    summary = {
        "model": args.model,
        "image": args.image,
        "thresholds": {
            "dchi2_min": DCHI2_MIN,
            "max_abs_mu": MAX_ABS_MU,
            "min_peak_snr": MIN_PEAK_SNR,
            "fit_radius_arcsec": FIT_RADIUS_ARCSEC,
            "halo_filter_arcsec": args.halo_filter_arcsec,
            "mask": not args.no_mask,
            "seg_snr": SEG_SNR,
            "min_size_ratio": MIN_SIZE_RATIO,
            "psf_sigma_arcsec": PSF_SIGMA_ARCSEC,
        },
        "counts": counts,
        "rows": rows,
    }
    out.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(counts))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

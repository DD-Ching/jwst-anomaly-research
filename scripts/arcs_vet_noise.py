"""Measurement-noise injections for ``arcs`` anti rows (vetting; plan in the output dir's PLAN.md).

For each strong-shear background row of a stored ``arcs`` table, an elliptical Gaussian with the
row's catalogue flux, ``semimajor_px`` (σ along the major axis) and ``ellipticity`` is injected
with its major axis along the z_s = 2 prediction (``aligned``) at empty positions in an annulus
around the row, on the same image (S3 byte-range cutout). The axis is re-measured with the
``arcs_bgsub_shapes.py`` method (median-filter subtraction, Gaussian-weighted moments around the
brightest pixel). P_anti is the fraction of valid injections measured ``anti``. A bright copy of
each injection (``--bright-mag``) is the method's known-case gate. ASSUMPTIONs: the Gaussian
profile, the annulus, N, the emptiness rule and the thresholds (``arcs.json``).

    python scripts/arcs_vet_noise.py --model abells1063-cats \\
        --image s3://stpubdata/jwst/public/jw03293/L3/t/o001/jw03293-o001_t001_nircam_clear-f200w_i2d.fits

Output: ``results/arcs/<model>/vetting/noise_injections.ecsv`` and ``noise_injections.json``
(derived; predictions model_prediction).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy.ndimage import median_filter

from jwst_anomaly import cutouts, paths

_spec = importlib.util.spec_from_file_location(
    "arcs_bgsub_shapes", Path(__file__).resolve().parent / "arcs_bgsub_shapes.py"
)
abs_ = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(abs_)


def gate_ok(frac_aligned_bright, min_frac: float = 0.95) -> bool:
    """Method validity: every row's bright injections come back aligned in ≥ ``min_frac``."""
    f = np.asarray(frac_aligned_bright, float)
    return bool(f.size > 0 and np.all(np.isfinite(f)) and np.all(f >= min_frac))


def noise_class(p_anti: float, p_min: float = 0.05) -> str:
    """Per target: ``noise-explicable`` if P_anti ≥ ``p_min``, else ``not noise-explicable``."""
    if not np.isfinite(p_anti):
        return "undetermined"
    return "noise-explicable" if p_anti >= p_min else "not noise-explicable"


def excess_p(p_rows, n_obs: int) -> float:
    """Poisson-binomial P(≥ n_obs successes) for independent rows with probabilities ``p_rows``."""
    dist = np.array([1.0])
    for p in np.asarray(p_rows, float):
        dist = np.convolve(dist, [1.0 - p, p])
    return float(dist[n_obs:].sum()) if n_obs < dist.size else 0.0


def excess_class(p: float, p_min: float = 0.05) -> str:
    return "excess explained by measurement noise" if p >= p_min else "excess remains"


def gaussian_stamp(shape, x0, y0, total, sig_major, ellip, theta_rad) -> np.ndarray:
    """Elliptical Gaussian with integral ``total``; ``theta_rad`` is the pixel-frame major axis."""
    yy, xx = np.indices(shape, dtype=float)
    dx, dy = xx - x0, yy - y0
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    u, v = c * dx + s * dy, -s * dx + c * dy
    sa, sb = sig_major, sig_major * (1.0 - ellip)
    return total / (2 * np.pi * sa * sb) * np.exp(-0.5 * ((u / sa) ** 2 + (v / sb) ** 2))


def measure(box: np.ndarray, filt_px: int, peak_px: int, weight_px: float):
    """Pixel-frame axis (deg) and ellipticity at the box centre, as ``arcs_bgsub_shapes``."""
    sub = box - median_filter(np.nan_to_num(box), size=filt_px)
    c = box.shape[0] // 2
    r = peak_px
    win = sub[c - r : c + r + 1, c - r : c + r + 1]
    j, i = np.unravel_index(np.nanargmax(win), win.shape)
    return abs_.weighted_axis(sub, c - r + i, c - r + j, weight_px)


def _row(job):
    path, row, a = job
    with fits.open(path) as h:
        sci = np.asarray(h["SCI"].data, float)
        wcs = WCS(h["SCI"].header)
    px = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))) * 3600.0
    filt = int(round(a["filter_arcsec"] / px)) | 1
    half = filt + 2 * int(np.ceil(3 * a["weight_px"]))
    omega = (px / 206264.806) ** 2
    to_img = 1.0 / (1e6 * omega)  # Jy -> MJy/sr summed over pixels
    xt, yt = wcs.world_to_pixel_values(row["ra"], row["dec"])
    rng = np.random.default_rng(a["seed"] + int(row["label"]))
    out = {"label": int(row["label"]), "n_valid": 0, "anti": 0, "aligned_bright": 0, "drawn": 0}
    offs, offs_b = [], []
    while out["n_valid"] < a["n_inject"] and out["drawn"] < 20 * a["n_inject"]:
        out["drawn"] += 1
        r = np.sqrt(rng.uniform(a["r_in"] ** 2, a["r_out"] ** 2)) / px
        phi = rng.uniform(0, 2 * np.pi)
        x, y = int(round(xt + r * np.cos(phi))), int(round(yt + r * np.sin(phi)))
        if x - half < 0 or y - half < 0 or x + half >= sci.shape[1] or y + half >= sci.shape[0]:
            continue
        box = sci[y - half : y + half + 1, x - half : x + half + 1]
        if not np.all(np.isfinite(box)):
            continue
        sub0 = box - median_filter(box, size=filt)
        sig = 1.4826 * np.median(np.abs(sub0 - np.median(sub0)))
        k = 4
        if np.max(np.abs(sub0[half - k : half + k + 1, half - k : half + k + 1])) >= 3 * sig:
            continue
        # pixel-frame direction of the predicted (aligned) sky PA at this position
        here = wcs.pixel_to_world(x, y)
        step = here.directional_offset_by(row["pa_pred_z2"] * u.deg, px * u.arcsec)
        x2, y2 = wcs.world_to_pixel(step)
        theta = np.arctan2(float(y2) - y, float(x2) - x)
        res = []
        for m in (row["mag"], a["bright_mag"]):
            tot = 3631.0 * 10 ** (-0.4 * m) * to_img
            ell = row["ellipticity"] if a["ellip"] is None else a["ellip"]
            stamp = gaussian_stamp(box.shape, half, half, tot, row["semimajor_px"], ell, theta)
            t_pix, e = measure(box + stamp, filt, a["peak_px"], a["weight_px"])
            res.append((t_pix, e))
        if not all(np.isfinite(v[0]) for v in res):
            continue
        out["n_valid"] += 1
        for (t_pix, _), lst in zip(res, (offs, offs_b), strict=True):
            lst.append(abs_.axis_offset(t_pix, np.degrees(theta) % 180.0))
    offs, offs_b = np.asarray(offs, float), np.asarray(offs_b, float)
    out["anti"] = int(np.sum(offs >= a["anti_deg"]))
    out["aligned"] = int(np.sum(offs <= a["aligned_deg"]))
    out["aligned_bright"] = int(np.sum(offs_b <= a["aligned_deg"]))
    out["median_offset"] = float(np.median(offs)) if offs.size else float("nan")
    return out


def run(args) -> dict:
    d = args.arcs_dir / args.model
    t = Table.read(d / "arcs.ecsv")
    t = t[np.asarray(t["strong_shear"], bool) & (t["population"] == "background")]
    stored = json.loads((d / "arcs.json").read_text())["assumptions"]
    bg = Table.read(d / "bgsub_shapes.ecsv")
    bg_class = {int(r["label"]): str(r["class_bgsub"]) for r in bg}
    a = {
        k: getattr(args, k)
        for k in (
            "filter_arcsec",
            "weight_px",
            "peak_px",
            "n_inject",
            "r_in",
            "r_out",
            "bright_mag",
            "seed",
            "ellip",
        )
    } | {"aligned_deg": stored["aligned_deg"], "anti_deg": stored["anti_deg"]}
    tg = Table({"source_uid": [str(v) for v in t["label"]], "ra": t["ra"], "dec": t["dec"]})
    with tempfile.TemporaryDirectory() as tmp:
        cut = cutouts.make_cutouts(
            args.image, tg, size_arcsec=2 * args.r_out + 3.0, out_dir=Path(tmp)
        )
        cols = ("label", "ra", "dec", "mag", "semimajor_px", "ellipticity", "pa_pred_z2")
        jobs = [
            (p, {k: float(r[k]) for k in cols}, a) for p, r in zip(cut["path"], t, strict=True) if p
        ]
        with ProcessPoolExecutor(args.procs) as ex:
            rows = list(ex.map(_row, jobs))
    out = Table(rows=[[r[k] for k in _COLS] for r in rows], names=_COLS)
    out["mag"] = [float(t["mag"][t["label"] == r][0]) for r in out["label"]]
    out["p_anti"] = out["anti"] / np.maximum(out["n_valid"], 1)
    out["frac_aligned_bright"] = out["aligned_bright"] / np.maximum(out["n_valid"], 1)
    out["class_bgsub"] = [bg_class.get(int(v), "missing") for v in out["label"]]
    out.meta.update(
        provenance="derived", source=f"{args.model} arcs.ecsv + {args.image}", assumptions=a
    )
    vd = d / "vetting"
    vd.mkdir(exist_ok=True)
    out.write(vd / f"noise_injections{args.tag}.ecsv", overwrite=True)
    gate = gate_ok(out["frac_aligned_bright"])
    n_obs = int(np.sum(out["class_bgsub"] == "anti"))
    p_ex = excess_p(out["p_anti"], n_obs)
    targets = {}
    for lab in args.targets:
        sel = out["label"] == lab
        p = float(out["p_anti"][sel][0]) if sel.any() else float("nan")
        targets[str(lab)] = {
            "p_anti": round(p, 4),
            "n_valid": int(out["n_valid"][sel][0]) if sel.any() else 0,
            "class": noise_class(p) if gate else "invalid (gate failed)",
        }
    summary = {
        "model": args.model,
        "image": args.image,
        "label": (
            "vetting test (derived); injected sources are simulations; predictions model_prediction"
        ),
        "assumptions": a,
        "gate": {
            "ok": gate,
            "min_frac_aligned_bright": round(float(np.min(out["frac_aligned_bright"])), 4),
            "rows": len(out),
        },
        "targets": targets,
        "sample": {
            "n_obs_anti_bgsub": n_obs,
            "expected_anti_measurement_noise": round(float(np.sum(out["p_anti"])), 3),
            "p_ge_n_obs": round(p_ex, 5),
            "class": excess_class(p_ex) if gate else "invalid (gate failed)",
        },
    }
    (vd / f"noise_injections{args.tag}.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


_COLS = ("label", "n_valid", "drawn", "anti", "aligned", "aligned_bright", "median_offset")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", required=True)
    p.add_argument("--image", required=True)
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--targets", type=int, nargs="+", default=[1316, 2184])
    p.add_argument("--filter-arcsec", type=float, default=1.0, help="ASSUMPTION")
    p.add_argument("--weight-px", type=float, default=3.0, help="ASSUMPTION")
    p.add_argument("--peak-px", type=int, default=3, help="ASSUMPTION")
    p.add_argument("--n-inject", type=int, default=60, help="ASSUMPTION")
    p.add_argument("--r-in", type=float, default=1.5, help="annulus inner radius, arcsec")
    p.add_argument("--r-out", type=float, default=6.0, help="annulus outer radius, arcsec")
    p.add_argument("--bright-mag", type=float, default=24.0, help="gate magnitude")
    p.add_argument("--seed", type=int, default=20261016)
    p.add_argument("--procs", type=int, default=4)
    p.add_argument(
        "--ellip", type=float, default=None, help="injected ellipticity (default: row's)"
    )
    p.add_argument("--tag", default="", help="output file suffix")
    args = p.parse_args(argv)
    print(json.dumps(run(args), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""RX J0152 family 4: revised forced photometry at the predicted third image (D-034 addendum bg).

Frozen plan: results/lens_consistency/rxj0152_fam4_photometry_v2/PLAN.md. Reads 50" cutouts of the
ten VENUS (GO 6882) NIRCam mosaics through MAST byte ranges (not stored) and writes ``result.json``.
Detection, segmentation and WCS helpers come from ``rxj0152_fam4_photometry.py`` (addendum bf).

    python scripts/rxj0152_fam4_photometry_v2.py [--cutout-dir DIR]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import tempfile
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from scipy import stats

_V1 = Path(__file__).resolve().with_name("rxj0152_fam4_photometry.py")
_spec = importlib.util.spec_from_file_location("rxj0152_fam4_photometry", _V1)
v1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v1)

OUT = v1.OUT.parent / "rxj0152_fam4_photometry_v2"
BANDS = (
    "f090w",
    "f115w",
    "f150w",
    "f200w",
    "f210m",
    "f277w",
    "f300m",
    "f356w",
    "f410m",
    "f444w",
)
REF = "f277w"
PREDICTED = v1.PREDICTED  # model_prediction
# derived: |mu 4.1| / mu_3 per model (posthoc.json of addendum be)
RATIO41 = {"lenstool": 7.165, "glafic": 5.325, "zitrin": 7.482}
REGION_RADIUS = 3.0  # ASSUMPTION: gate-K match radius
FLUX_FACTOR = 3.0  # ASSUMPTION
COLOUR_TOL = 0.3  # ASSUMPTION, dex
DETECTABLE_SIGMA = 3.0  # ASSUMPTION
MIN_BANDS = 3
EXCLUDE_RADIUS = 0.3
APER_RADIUS = 0.3
N_APER = 300
SEED = 20261025
OFFSET_ARCSEC = 1.5
OFFSET_ANGLES = tuple(range(0, 360, 40))
RECOVER_RADIUS = 1.0
V_MIN_FRACTION = 0.8  # ASSUMPTION
P_CHANCE = 0.05  # ASSUMPTION


def flux_window(f41_ref: float) -> tuple[float, float]:
    """Rule 6 F277W window for the third image."""
    lo = f41_ref / max(RATIO41.values()) / FLUX_FACTOR
    hi = FLUX_FACTOR * f41_ref / min(RATIO41.values())
    return lo, hi


def colour_match(f: dict[str, float], ref: dict[str, float], sigma: dict[str, float]) -> bool:
    """Rule 6 colour test over the bands where the scaled reference SED is >= 3 sigma_ap."""
    if not (np.isfinite(f[REF]) and f[REF] > 0 and np.isfinite(ref[REF]) and ref[REF] > 0):
        return False
    used = 0
    for b in BANDS:
        if b == REF or not np.isfinite(ref[b]) or ref[b] <= 0:
            continue
        if ref[b] * f[REF] / ref[REF] < DETECTABLE_SIGMA * sigma[b]:
            continue
        used += 1
        if not np.isfinite(f[b]) or f[b] <= 0:
            return False
        d = np.log10(f[b] / f[REF]) - np.log10(ref[b] / ref[REF])
        if abs(d) > COLOUR_TOL:
            return False
    return used >= MIN_BANDS


def chance_p(n_field: int, a_region: float, a_valid: float, n_obs: int) -> tuple[float, float]:
    """Gate X: expected chance count in the region and Poisson P(N >= n_obs)."""
    lam = n_field * a_region / (a_valid - a_region)
    if n_obs == 0:
        return lam, 1.0
    if lam == 0:  # no field rate measured: uninformative, never "significant" (post-run fix)
        return lam, float("nan")
    return lam, float(stats.poisson.sf(n_obs - 1, lam))


def v_pass(at_prediction: list[bool], all_recovered: list[bool]) -> bool:
    return all(at_prediction) and sum(all_recovered) >= V_MIN_FRACTION * len(all_recovered)


def known_case_pass(n_accepted: int) -> bool:
    return n_accepted >= 6


def outcome(v_ok: bool, c_ok: bool, d_ok: bool, n_obs: int, p: float) -> str:
    """Frozen outcome table."""
    if not (v_ok and c_ok):
        return "O3_undetermined"
    if n_obs > 0:
        return "O1_counterpart" if p <= P_CHANCE else "O1c_chance_level"
    return "O2_missing" if d_ok else "O3_undetermined"


def offset_positions(centre: tuple[float, float]) -> list[tuple[float, float]]:
    c = SkyCoord(*centre, unit="deg")
    out = []
    for a in OFFSET_ANGLES:
        p = c.directional_offset_by(a * u.deg, OFFSET_ARCSEC * u.arcsec)
        out.append((float(p.ra.deg), float(p.dec.deg)))
    return out


# --- image handling -------------------------------------------------------------------------


def fetch(cutout_dir: Path) -> dict[str, tuple[np.ndarray, WCS]]:
    from jwst_anomaly.cutouts import make_cutouts

    targets = Table({"source_uid": ["fam4"], "ra": [v1.CENTRE[0]], "dec": [v1.CENTRE[1]]})
    out = {}
    for band in BANDS:
        tab = make_cutouts(
            v1.URI.format(band=band), targets, size_arcsec=v1.SIZE_ARCSEC, out_dir=cutout_dir
        )
        with fits.open(tab["path"][0]) as h:
            out[band] = (np.array(h["SCI"].data, float), WCS(h["SCI"].header))
    return out


def forced(images: dict, positions: list[tuple[float, float]]) -> list[dict[str, float]]:
    from photutils.aperture import SkyCircularAperture, aperture_photometry

    sky = SkyCoord([p[0] for p in positions], [p[1] for p in positions], unit="deg")
    ap = SkyCircularAperture(sky, r=APER_RADIUS * u.arcsec)
    res = [{} for _ in positions]
    for b in BANDS:
        sci, wcs = images[b]
        tab = aperture_photometry(np.nan_to_num(sci), ap, wcs=wcs, method="exact")
        for i, v in enumerate(tab["aperture_sum"]):
            res[i][b] = float(v)
    return res


def aperture_noise(images: dict, segmaps: dict) -> dict[str, float]:
    from photutils.aperture import CircularAperture, aperture_photometry

    out = {}
    for b in BANDS:
        sci, wcs = images[b]
        seg = segmaps[b]
        r = APER_RADIUS / wcs.proj_plane_pixel_scales()[0].to("arcsec").value
        rng = np.random.default_rng(SEED)
        cand = rng.uniform([r, r], [sci.shape[1] - r, sci.shape[0] - r], size=(20 * N_APER, 2))
        keep = seg[cand[:, 1].astype(int), cand[:, 0].astype(int)] == 0
        cand = cand[keep]
        tab = aperture_photometry(sci, CircularAperture(cand, r), method="exact")
        sums = np.asarray(tab["aperture_sum"], float)
        sums = sums[np.isfinite(sums)][:N_APER]
        out[b] = float(sigma_clipped_stats(sums, sigma=3.0)[2])
    return out


def segmaps_for(images: dict, seg_ref: np.ndarray) -> dict[str, np.ndarray]:
    wcs_ref = images[REF][1]
    return {
        b: seg_ref
        if b == REF
        else v1.reproject_seg(seg_ref, wcs_ref, images[b][0].shape, images[b][1])
        for b in BANDS
    }


def _sep(a: tuple[float, float], b: tuple[float, float]) -> float:
    return SkyCoord(*a, unit="deg").separation(SkyCoord(*b, unit="deg")).arcsec


def in_region(p: tuple[float, float]) -> bool:
    return any(_sep(p, q) <= REGION_RADIUS for q in PREDICTED.values())


def measure(images: dict, pos: dict, sigma: dict | None = None) -> dict:
    """Rules 2-6 (and the gate X field count); segmaps and sigma are computed only without sigma."""
    sci_ref, wcs_ref = images[REF]
    seg_ref = v1.segment(sci_ref)
    cents = v1.seg_centroids(sci_ref, seg_ref, wcs_ref)
    near41 = v1.labels_near(seg_ref, wcs_ref, *pos["4.1"], EXCLUDE_RADIUS)
    near42 = v1.labels_near(seg_ref, wcs_ref, *pos["4.2"], EXCLUDE_RADIUS)
    excl = near41 | near42
    segmaps = None
    if sigma is None:
        segmaps = segmaps_for(images, seg_ref)
        sigma = aperture_noise(images, segmaps)
    labs = [lab for lab in cents if lab not in excl]
    fl = forced(images, [cents[lab] for lab in labs] + [pos["4.1"]])
    f41 = fl.pop()
    lo, hi = flux_window(f41[REF])
    region, field = [], 0
    for lab, f in zip(labs, fl, strict=True):
        in_flux = bool(lo <= f[REF] <= hi)
        col = colour_match(f, f41, sigma)
        ok = in_flux and col
        if in_region(cents[lab]):
            region.append(
                {
                    "label": lab,
                    "ra": cents[lab][0],
                    "dec": cents[lab][1],
                    "fluxes": f,
                    "in_flux_window": in_flux,
                    "colour_match": col,
                    "passes": bool(ok),
                }
            )
        elif ok:
            field += 1
    return {
        "seg_ref": seg_ref,
        "segmaps": segmaps,
        "near41": sorted(near41),
        "near42": sorted(near42),
        "sigma": sigma,
        "f41": f41,
        "window": (lo, hi),
        "region": region,
        "counterparts": [r for r in region if r["passes"]],
        "n_field": field,
    }


def inject(images: dict, segmaps: dict, labels41: list[int], pos41, target, scale: float) -> dict:
    """Rule 7: add 4.1's own segment pixels / R_m, shifted to ``target``, in every band."""
    out = {}
    for b in BANDS:
        sci, wcs = images[b]
        ys, xs = np.nonzero(np.isin(segmaps[b], labels41))
        x41, y41 = wcs.world_to_pixel(SkyCoord(*pos41, unit="deg"))
        tx, ty = wcs.world_to_pixel(SkyCoord(*target, unit="deg"))
        dx, dy = int(round(float(tx - x41))), int(round(float(ty - y41)))
        new = sci.copy()
        ny, nx = ys + dy, xs + dx
        ok = (ny >= 0) & (nx >= 0) & (ny < sci.shape[0]) & (nx < sci.shape[1])
        new[ny[ok], nx[ok]] += np.nan_to_num(sci[ys[ok], xs[ok]]) * scale
        out[b] = (new, wcs)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutout-dir", type=Path, default=None)
    args = ap.parse_args()
    # cloud disk rule: a default temporary cutout directory is deleted after the run
    cdir = args.cutout_dir or Path(tempfile.mkdtemp(prefix="fam4v2_"))
    raw = fetch(cdir)
    images = {b: (v1.background(s), w) for b, (s, w) in raw.items()}
    cat = Table.read(v1.IMAGES)
    pos = {str(r["image_id"]): (float(r["ra"]), float(r["dec"])) for r in cat}
    res = measure(images, pos)
    summary: dict = {"plan": "PLAN.md", "bands": list(BANDS), "near41": res["near41"]}
    if not res["near41"] or set(res["near41"]) & set(res["near42"]):
        summary["outcome"] = "stopped: 4.1 not isolated"
        (OUT / "result.json").write_text(json.dumps(summary, indent=1) + "\n")
        print(summary)
        return
    sigma = res["sigma"]
    # gate V
    at_pred, every, vdetail = [], [], []
    for name, p in PREDICTED.items():
        for k, t in enumerate([p, *offset_positions(p)]):
            r2 = measure(
                inject(images, res["segmaps"], res["near41"], pos["4.1"], t, 1.0 / RATIO41[name]),
                pos,
                sigma,
            )
            hit = any(_sep((c["ra"], c["dec"]), t) < RECOVER_RADIUS for c in r2["counterparts"])
            vdetail.append({"model": name, "k": k, "ra": t[0], "dec": t[1], "recovered": hit})
            every.append(hit)
            if k == 0:
                at_pred.append(hit)
    v_ok = v_pass(at_pred, every)
    # gate C
    ids = [i for i in pos if i.split(".")[0] in ("2", "3", "5", "6")]
    fc = dict(zip(ids, forced(images, [pos[i] for i in ids]), strict=True))
    cdetail = {}
    for i in ids:
        fam = i.split(".")[0]
        sib = [j for j in ids if j != i and j.split(".")[0] == fam]
        cdetail[i] = bool(any(colour_match(fc[i], fc[j], sigma) for j in sib))
    n_acc = sum(cdetail.values())
    c_ok = known_case_pass(n_acc)
    # gate D
    lo, hi = res["window"]
    d_ok = lo >= 5 * sigma[REF]
    # gate X
    seg_ref = res["seg_ref"]
    wcs_ref = images[REF][1]
    yy, xx = np.mgrid[: seg_ref.shape[0], : seg_ref.shape[1]]
    sky = wcs_ref.pixel_to_world(xx.ravel(), yy.ravel())
    reg = np.zeros(sky.size, bool)
    for p in PREDICTED.values():
        reg |= sky.separation(SkyCoord(*p, unit="deg")).arcsec <= REGION_RADIUS
    finite = np.isfinite(images[REF][0]).ravel()
    pix_area = wcs_ref.proj_plane_pixel_area().to("arcsec2").value
    a_region, a_valid = float((reg & finite).sum() * pix_area), float(finite.sum() * pix_area)
    n_obs = len(res["counterparts"])
    lam, p = chance_p(res["n_field"], a_region, a_valid, n_obs)
    oc = outcome(v_ok, c_ok, d_ok, n_obs, p)
    for c in res["counterparts"]:
        near = min(PREDICTED, key=lambda m: _sep((c["ra"], c["dec"]), PREDICTED[m]))
        exp = res["f41"][REF] / RATIO41[near]
        c.update(
            nearest_model=near,
            sep_arcsec=_sep((c["ra"], c["dec"]), PREDICTED[near]),
            within_factor3_of_model=bool(
                exp / FLUX_FACTOR <= c["fluxes"][REF] <= FLUX_FACTOR * exp
            ),
        )
    summary.update(
        f41=res["f41"],
        sigma_ap=sigma,
        window_ref=[lo, hi],
        gate_V_at_prediction=dict(zip(PREDICTED, at_pred, strict=True)),
        gate_V_recovered=f"{sum(every)}/{len(every)}",
        gate_V_detail=vdetail,
        gate_V_pass=v_ok,
        gate_C_detail=cdetail,
        gate_C_accepted=n_acc,
        gate_C_pass=c_ok,
        gate_D_pass=d_ok,
        gate_X={
            "n_field": res["n_field"],
            "a_region": a_region,
            "a_valid": a_valid,
            "lambda": lam,
            "p_ge_obs": p,
        },
        region_segments=res["region"],
        counterparts=res["counterparts"],
        outcome=oc,
        units="0.3 arcsec aperture sums of i2d SCI pixels (MJy/sr x pixels); ratios only are used",
    )
    (OUT / "result.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    if args.cutout_dir:
        np.save(cdir / "seg_ref.npy", seg_ref)
    else:
        shutil.rmtree(cdir, ignore_errors=True)
    keys = ("gate_V_at_prediction", "gate_V_recovered", "gate_C_accepted", "gate_D_pass", "gate_X")
    print(json.dumps({k: summary[k] for k in (*keys, "outcome")}, default=float))
    print("counterparts", n_obs, "region segments", len(res["region"]))


if __name__ == "__main__":
    main()

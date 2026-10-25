"""RX J0152 family 4: forced photometry at the model-predicted third image (D-034 addendum bf).

Frozen plan: results/lens_consistency/rxj0152_fam4_photometry/PLAN.md. Reads 50" cutouts of the
VENUS (GO 6882) NIRCam mosaics through MAST byte ranges (not stored) and writes ``result.json``.

    python scripts/rxj0152_fam4_photometry.py [--cutout-dir DIR]
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS

OUT = (
    Path(__file__).resolve().parents[1] / "results" / "lens_consistency" / "rxj0152_fam4_photometry"
)
IMAGES = OUT.parent / "slice_relics_rxj0152" / "slice_images.ecsv"
BANDS = ("f115w", "f150w", "f277w", "f444w")
REF = "f277w"
URI = "mast:JWST/product/jw06882-o008_t008_nircam_clear-{band}_i2d.fits"
CENTRE = (28.18100, -13.95480)
SIZE_ARCSEC = 50.0
# model_prediction (posthoc.json of addendum be): positive-parity unmatched image of family 4
PREDICTED = {
    "lenstool": (28.17894873479369, -13.952926865904491),
    "glafic": (28.17805, -13.95268),
    "zitrin": (28.17945, -13.95343),
}
# derived: (|mu 4.1| + |mu 4.2|) / mu_3 per model (image-plane tables of addendum be)
RATIO = {"lenstool": 14.05, "glafic": 10.05, "zitrin": 15.10}
REGION_RADIUS = 3.0  # ASSUMPTION: gate-K match radius
FLUX_FACTOR = 3.0  # ASSUMPTION
COLOUR_TOL = 0.3  # ASSUMPTION, dex
PAIR_RADIUS = 0.3
APER_RADIUS = 0.3
N_APER = 300
SEED = 20261025


def flux_window(f_pair_ref: float) -> tuple[float, float]:
    """Rule 5 F277W flux window for a counterpart."""
    lo = f_pair_ref / (FLUX_FACTOR * max(RATIO.values()))
    hi = FLUX_FACTOR * f_pair_ref / min(RATIO.values())
    return lo, hi


def colour_match(fluxes: dict[str, float], ref_fluxes: dict[str, float]) -> bool:
    """Rule 5 colour test: every band within COLOUR_TOL dex of the reference SED (vs F277W)."""
    if any(not np.isfinite(fluxes[b]) or fluxes[b] <= 0 for b in BANDS):
        return False
    if any(not np.isfinite(ref_fluxes[b]) or ref_fluxes[b] <= 0 for b in BANDS):
        return False
    for b in BANDS:
        if b == REF:
            continue
        d = np.log10(fluxes[b] / fluxes[REF]) - np.log10(ref_fluxes[b] / ref_fluxes[REF])
        if abs(d) > COLOUR_TOL:
            return False
    return True


def outcome(v_pass: bool, c_pass: bool, d_pass: bool, n_counterparts: int) -> str:
    """Frozen outcome table."""
    if not (v_pass and c_pass):
        return "O3_undetermined"
    if n_counterparts > 0:
        return "O1_counterpart"
    return "O2_missing" if d_pass else "O3_undetermined"


def known_case_pass(n_accepted: int, n_total: int = 12) -> bool:
    return n_total > 0 and n_accepted >= 6


def in_region(ra: float, dec: float) -> bool:
    c = SkyCoord(ra, dec, unit="deg")
    return any(
        c.separation(SkyCoord(*p, unit="deg")).arcsec <= REGION_RADIUS for p in PREDICTED.values()
    )


# --- image handling -------------------------------------------------------------------------


def fetch(cutout_dir: Path) -> dict[str, tuple[np.ndarray, WCS]]:
    from jwst_anomaly.cutouts import make_cutouts

    targets = Table({"source_uid": ["fam4"], "ra": [CENTRE[0]], "dec": [CENTRE[1]]})
    out = {}
    for band in BANDS:
        tab = make_cutouts(
            URI.format(band=band), targets, size_arcsec=SIZE_ARCSEC, out_dir=cutout_dir
        )
        with fits.open(tab["path"][0]) as h:
            sci = np.array(h["SCI"].data, float)
            out[band] = (sci, WCS(h["SCI"].header))
    return out


def background(sci: np.ndarray) -> np.ndarray:
    _, med, _ = sigma_clipped_stats(sci[np.isfinite(sci)], sigma=3.0)
    return sci - med


def segment(sci: np.ndarray):
    from photutils.segmentation import deblend_sources, detect_sources

    _, _, std = sigma_clipped_stats(sci[np.isfinite(sci)], sigma=3.0)
    img = np.nan_to_num(sci)
    conv = convolve(img, Gaussian2DKernel(2.0 * gaussian_fwhm_to_sigma))
    seg = detect_sources(conv, 2.0 * std, npixels=10)
    seg = deblend_sources(conv, seg, npixels=10, progress_bar=False)
    return seg.data


def reproject_seg(seg: np.ndarray, wcs_ref: WCS, shape: tuple[int, int], wcs: WCS) -> np.ndarray:
    """Nearest-neighbour map of the reference segmentation onto another band's pixel grid."""
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    sky = wcs.pixel_to_world(xx.ravel(), yy.ravel())
    x, y = wcs_ref.world_to_pixel(sky)
    xi, yi = np.rint(x).astype(int), np.rint(y).astype(int)
    ok = (xi >= 0) & (yi >= 0) & (xi < seg.shape[1]) & (yi < seg.shape[0])
    out = np.zeros(xi.size, int)
    out[ok] = seg[yi[ok], xi[ok]]
    return out.reshape(shape)


def seg_fluxes(images: dict, segmaps: dict) -> dict[int, dict[str, float]]:
    labels = np.unique(segmaps[REF])
    labels = labels[labels > 0]
    res = {int(lab): {} for lab in labels}
    for b in BANDS:
        sci = np.nan_to_num(images[b][0])
        sums = np.bincount(segmaps[b].ravel(), weights=sci.ravel(), minlength=labels.max() + 1)
        for lab in labels:
            res[int(lab)][b] = float(sums[lab])
    return res


def seg_centroids(sci: np.ndarray, seg: np.ndarray, wcs: WCS) -> dict[int, tuple[float, float]]:
    from scipy import ndimage

    labels = np.unique(seg)
    labels = labels[labels > 0]
    w = np.clip(np.nan_to_num(sci), 0, None)
    cents = ndimage.center_of_mass(w, seg, labels)
    out = {}
    for lab, (cy, cx) in zip(labels, cents, strict=True):
        if not np.isfinite(cx):
            continue
        s = wcs.pixel_to_world(cx, cy)
        out[int(lab)] = (float(s.ra.deg), float(s.dec.deg))
    return out


def labels_near(seg: np.ndarray, wcs: WCS, ra: float, dec: float, radius: float) -> set[int]:
    x, y = wcs.world_to_pixel(SkyCoord(ra, dec, unit="deg"))
    scale = wcs.proj_plane_pixel_scales()[0].to("arcsec").value
    yy, xx = np.mgrid[: seg.shape[0], : seg.shape[1]]
    m = (xx - x) ** 2 + (yy - y) ** 2 <= (radius / scale) ** 2
    return {int(v) for v in np.unique(seg[m]) if v > 0}


def aperture_noise(sci: np.ndarray, seg: np.ndarray, wcs: WCS) -> float:
    scale = wcs.proj_plane_pixel_scales()[0].to("arcsec").value
    r = APER_RADIUS / scale
    rng = np.random.default_rng(SEED)
    yy, xx = np.mgrid[: sci.shape[0], : sci.shape[1]]
    sums = []
    tries = 0
    while len(sums) < N_APER and tries < 50 * N_APER:
        tries += 1
        x, y = rng.uniform(r, sci.shape[1] - r), rng.uniform(r, sci.shape[0] - r)
        if seg[int(y), int(x)] > 0:
            continue
        m = (xx - x) ** 2 + (yy - y) ** 2 <= r**2
        vals = sci[m]
        if not np.all(np.isfinite(vals)):
            continue
        sums.append(vals.sum())
    _, _, std = sigma_clipped_stats(np.array(sums), sigma=3.0)
    return float(std)


def measure(images: dict) -> dict:
    """Rules 2-5 on one set of background-subtracted images."""
    sci_ref, wcs_ref = images[REF]
    seg_ref = segment(sci_ref)
    segmaps = {
        b: seg_ref
        if b == REF
        else reproject_seg(seg_ref, wcs_ref, images[b][0].shape, images[b][1])
        for b in BANDS
    }
    fl = seg_fluxes(images, segmaps)
    cents = seg_centroids(sci_ref, seg_ref, wcs_ref)
    cat = Table.read(IMAGES)
    pos = {str(r["image_id"]): (float(r["ra"]), float(r["dec"])) for r in cat}
    pair = set()
    for iid in ("4.1", "4.2"):
        pair |= labels_near(seg_ref, wcs_ref, *pos[iid], PAIR_RADIUS)
    res = {
        "seg_ref": seg_ref,
        "fluxes": fl,
        "centroids": cents,
        "pair_labels": sorted(pair),
        "pos": pos,
    }
    if not pair:
        return res
    f_pair = {b: sum(fl[lab][b] for lab in pair) for b in BANDS}
    lo, hi = flux_window(f_pair[REF])
    found = []
    for lab, (ra, dec) in cents.items():
        if lab in pair or not in_region(ra, dec):
            continue
        f = fl[lab]
        in_flux = lo <= f[REF] <= hi
        col = colour_match(f, f_pair)
        entry = {
            "label": lab,
            "ra": ra,
            "dec": dec,
            "fluxes": f,
            "in_flux_window": in_flux,
            "colour_match": col,
        }
        res.setdefault("region_segments", []).append(entry)
        if in_flux and col:
            found.append(entry)
    res.update(f_pair=f_pair, window=(lo, hi), counterparts=found)
    return res


def inject(images: dict, seg_ref: np.ndarray, pair: list[int], target: tuple[float, float]) -> dict:
    """Rule 7: add the pair's segment pixels scaled by 1/max(R), shifted to ``target``."""
    out = {}
    wcs_ref = images[REF][1]
    for b in BANDS:
        sci, wcs = images[b]
        segb = seg_ref if b == REF else reproject_seg(seg_ref, wcs_ref, sci.shape, wcs)
        mask = np.isin(segb, pair)
        ys, xs = np.nonzero(mask)
        cx, cy = xs.mean(), ys.mean()
        tx, ty = wcs.world_to_pixel(SkyCoord(*target, unit="deg"))
        dx, dy = int(round(float(tx) - cx)), int(round(float(ty) - cy))
        new = sci.copy()
        ny, nx = ys + dy, xs + dx
        ok = (ny >= 0) & (nx >= 0) & (ny < sci.shape[0]) & (nx < sci.shape[1])
        new[ny[ok], nx[ok]] += np.nan_to_num(sci[ys[ok], xs[ok]]) / max(RATIO.values())
        out[b] = (new, wcs)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutout-dir", type=Path, default=None)
    args = ap.parse_args()
    cdir = args.cutout_dir or Path(tempfile.mkdtemp(prefix="fam4_"))
    raw = fetch(cdir)
    images = {b: (background(s), w) for b, (s, w) in raw.items()}
    res = measure(images)
    wcs_ref = images[REF][1]
    summary: dict = {"plan": "PLAN.md", "bands": list(BANDS), "pair_labels": res["pair_labels"]}
    if not res["pair_labels"]:
        summary["outcome"] = "stopped: pair not detected"
        (OUT / "result.json").write_text(json.dumps(summary, indent=1) + "\n")
        return
    # gate V
    v = {}
    for name, p in PREDICTED.items():
        r2 = measure(inject(images, res["seg_ref"], res["pair_labels"], p))
        v[name] = any(
            SkyCoord(c["ra"], c["dec"], unit="deg").separation(SkyCoord(*p, unit="deg")).arcsec
            < 1.0
            for c in r2["counterparts"]
        )
    v_pass = all(v.values())
    # gate C
    accepted, details = 0, {}
    seg_ref, fl = res["seg_ref"], res["fluxes"]
    for fam in ("2", "3", "5", "6"):
        ids = [i for i in res["pos"] if i.split(".")[0] == fam]
        labs = {i: sorted(labels_near(seg_ref, wcs_ref, *res["pos"][i], PAIR_RADIUS)) for i in ids}
        sed = {
            i: ({b: sum(fl[lab][b] for lab in labs[i]) for b in BANDS} if labs[i] else None)
            for i in ids
        }
        for i in ids:
            ok = sed[i] is not None and any(
                sed[j] is not None and set(labs[j]) != set(labs[i]) and colour_match(sed[i], sed[j])
                for j in ids
                if j != i
            )
            details[i] = {"labels": labs[i], "accepted": bool(ok)}
            accepted += int(ok)
    c_pass = known_case_pass(accepted, 12)
    # gate D
    sigma = aperture_noise(images[REF][0], seg_ref, wcs_ref)
    lo, hi = res["window"]
    d_pass = lo >= 5 * sigma
    oc = outcome(v_pass, c_pass, d_pass, len(res["counterparts"]))
    summary.update(
        f_pair_MJy_sr_sum=res["f_pair"],
        window_ref=[lo, hi],
        sigma_ap_ref=sigma,
        gate_V=v,
        gate_V_pass=v_pass,
        gate_C_accepted=accepted,
        gate_C_detail=details,
        gate_C_pass=c_pass,
        gate_D_pass=d_pass,
        region_segments=res.get("region_segments", []),
        counterparts=res["counterparts"],
        outcome=oc,
        units="sums of i2d SCI pixels (MJy/sr x pixels); ratios only are used",
    )
    (OUT / "result.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    np.save(cdir / "seg_ref.npy", seg_ref)
    print(
        json.dumps({k: summary[k] for k in ("gate_V", "gate_C_accepted", "gate_D_pass", "outcome")})
    )
    print(
        "counterparts",
        len(res["counterparts"]),
        "region segments",
        len(res.get("region_segments", [])),
    )


if __name__ == "__main__":
    main()

"""AT 2024jfz (SMACS 0723): z_s 1.35-1.40 counter-image photometry in NIRISS 4368 F150W.

Plan: ``results/transients/smacs_at2024jfz/niriss4368/PLAN.md`` (frozen before the run). The
Mahler+22 ICLv2 track of the nearby counter-image (``../posthoc.json``, z_s 1.35 and 1.40, plus
points between) is gridded within ``SEARCH_ARCSEC``; each grid point is measured in the 4368 epoch
(E) and the ERO 2736 NIRCam F150W template (T) with ``hit_lightcurve.aperture_photometry``, and
the frozen function :func:`outcome` decides "present", "absent" or "untestable" (with an injection
control in E). Boxes are S3 byte-range reads; nothing is stored. Fluxes are ``derived``; the
positions ``model_prediction``; thresholds ASSUMPTIONs.

    python scripts/smacs_at2024jfz_niriss.py --out results/transients/smacs_at2024jfz/niriss4368
"""

from __future__ import annotations

import argparse
import json
import sys
import zlib
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arc_diff import read_box  # noqa: E402
from hit_lightcurve import (  # noqa: E402
    ANNULUS,
    NJY_PER_MJYSR_SR,
    aperture_photometry,
    empirical_error,
)
from smacs_at2024jfz import F0_NJY, T0_MJD, TNS_RADEC  # noqa: E402

EPOCH_OBS = "jw04368-o001_t001_niriss_clear-f150w"
TEMPLATE_OBS = "jw02736-o001_t001_nircam_clear-f150w"
Z_TRACK = (1.35, 1.40)  # posthoc.json rows whose window covers MJD 60557.96
N_TRACK = 6  # track points: the two ends and 4 between (ASSUMPTION)
SEARCH_ARCSEC = 0.5  # ASSUMPTION: model position error
GRID_ARCSEC = 0.1  # ASSUMPTION: grid step
RADIUS_ARCSEC = 0.15  # aperture radius (hit_lightcurve default)
BOX_HALF_ARCSEC = 4.0  # box half-size; the empirical error uses the whole box
MU_RATIO = 0.84  # |mu_k / mu_0| at z_s 1.35-1.40 (posthoc.json: 0.8425, 0.8375)
F_PRED_NJY = F0_NJY * MU_RATIO  # 2166 nJy; the plan rounds it to 2170
T_AP_NJY = 0.125 * 2170.0  # ASSUMPTION: aperture fraction >= 0.5, fading or mu error <= x4
INJ_NJY = 0.25 * 2170.0  # injected total flux of the control
INJ_FWHM_ARCSEC = 0.08  # ASSUMPTION: Gaussian stand-in for the drizzled NIRISS F150W PSF
NSIGMA = 5.0
MAX_NAN_FRAC = 0.2
MIN_INJ_PASS = 5  # of N_TRACK
STATIC_SEARCH_ARCSEC = 0.6  # post hoc: static-source search radius around a flag
STATIC_AP_ARCSEC = 0.8  # post hoc: aperture holding most of a point source in both instruments
POSTHOC_MU_FACTOR = 0.5  # post hoc ASSUMPTION: model |mu_k/mu_0| may be high by up to x2
STATIC_FLAG_ARCSEC = 0.3  # post hoc ASSUMPTION: a flag this close to a T source is its wing
STATIC_OFFSET_ARCSEC = 0.05  # post hoc ASSUMPTION: same source if E and T centroids agree to this
EPOCH_MJD = 60557.96  # MJD of EPOCH_OBS as listed in ../mast_epochs.json and ../posthoc.json


def classify_flag(verdict: str, static: dict | None, threshold: float) -> str:
    """Post hoc: a "present" that is the wing of a static source becomes "absent" or "undetermined".

    A flag within ``STATIC_FLAG_ARCSEC`` of a source whose E and T centroids agree within
    ``STATIC_OFFSET_ARCSEC`` is that source's cross-instrument residual. Then a counter-image
    blended with the source is excluded ("absent") only if the source's total E - T flux change
    is below ``threshold``; otherwise "undetermined". Without a usable ``static`` the verdict stays.
    """
    if verdict != "present" or not static:
        return verdict
    if (
        static["flag_to_source_arcsec"] <= STATIC_FLAG_ARCSEC
        and static["e_minus_t_offset_arcsec"] <= STATIC_OFFSET_ARCSEC
    ):
        return "absent" if abs(static["flux_e_minus_t_njy"]) < threshold else "undetermined"
    return verdict


def posthoc_outcome(
    delta: np.ndarray, sig: float, nan_frac: float, threshold: float, inj_delta: float | None = None
) -> str:
    """Post hoc rule (not in PLAN.md) for one counter-image, with a one-point injection gate.

    ``inj_delta`` is E - T at the predicted position after injecting the expected counter-image
    (``threshold / POSTHOC_MU_FACTOR`` in aperture units); it must pass the same cuts, else
    "untestable".

    ``threshold`` is ``POSTHOC_MU_FACTOR * |mu_k/mu_0| * (F_E - F_T)(TNS)``: under a single-peaked
    light curve, a counter-image whose delay lies in [0, t_E - t_0] shows a phase at which the
    event was at least as bright as at the TNS position in E (same instrument, so no aperture
    correction).
    """
    if nan_frac > MAX_NAN_FRAC:
        return "untestable"
    if inj_delta is not None and not (
        np.isfinite(inj_delta) and inj_delta >= threshold and inj_delta >= NSIGMA * sig
    ):
        return "untestable"
    d = np.asarray(delta, float)
    hit = np.isfinite(d) & (d >= threshold) & (d >= NSIGMA * sig)
    return "present" if hit.any() else "absent"


def outcome(
    delta: np.ndarray,
    sig: float,
    nan_frac_e: float,
    nan_frac_t: float,
    inj_delta: np.ndarray,
    inj_sig: float,
    t_ap: float = T_AP_NJY,
) -> str:
    """Frozen rule of PLAN.md, combined over the grid: "present", "absent" or "untestable".

    ``delta``: E - T aperture fluxes at the grid points (nJy, NaN where unmeasured); ``sig``: their
    sigma; ``inj_delta``/``inj_sig``: the same at the track points after injection into E.
    """
    if nan_frac_e > MAX_NAN_FRAC or nan_frac_t > MAX_NAN_FRAC:
        return "untestable"
    inj = np.asarray(inj_delta, float)
    passed = np.isfinite(inj) & (inj >= t_ap) & (inj >= NSIGMA * inj_sig)
    if passed.sum() < MIN_INJ_PASS:
        return "untestable"
    d = np.asarray(delta, float)
    hit = np.isfinite(d) & (d >= t_ap) & (d >= NSIGMA * sig)
    return "present" if hit.any() else "absent"


def _peak_centroid(box: dict, ra: float, dec: float, search: float) -> tuple[float, float]:
    """Sky position of the 7 x 7 px centroid of the brightest pixel near (ra, dec).

    The box median is subtracted first; raises ValueError on an empty, all-NaN or flat window.
    """
    sci = box["sci"] - np.nanmedian(box["sci"])
    x, y = (float(v) for v in box["wcs"].all_world2pix(ra, dec, 0))
    h = max(int(np.ceil(search / box["scale"])), 1)
    x0, y0 = max(int(round(x)) - h, 3), max(int(round(y)) - h, 3)
    x1 = min(int(round(x)) + h + 1, sci.shape[1] - 3)
    y1 = min(int(round(y)) + h + 1, sci.shape[0] - 3)
    win = sci[y0:y1, x0:x1]
    if win.size == 0 or not np.isfinite(win).any():
        raise ValueError("no finite pixels in the centroid window")
    py, px = np.unravel_index(np.nanargmax(win), win.shape)
    px, py = px + x0, py + y0
    w = np.clip(np.nan_to_num(sci[py - 3 : py + 4, px - 3 : px + 4]), 0, None)
    yy, xx = np.mgrid[py - 3 : py + 4, px - 3 : px + 4]
    if w.sum() <= 0:
        raise ValueError("flat centroid window")
    cx, cy = (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    ra_c, dec_c = box["wcs"].all_pix2world(cx, cy, 0)
    return float(ra_c), float(dec_c)


def static_check(boxes: dict, ra: float, dec: float) -> dict:
    """Post hoc: is a flag the wing of a static source? Same source in T and E, and its total flux.

    The brightest T pixel within ``STATIC_SEARCH_ARCSEC`` of the flag, centroided in T and then in
    E (searched near the T centroid); fluxes in ``STATIC_AP_ARCSEC`` (annulus 1.0-1.5 arcsec).
    """
    ra_t, dec_t = _peak_centroid(boxes["T"], ra, dec, STATIC_SEARCH_ARCSEC)
    ra_e, dec_e = _peak_centroid(boxes["E"], ra_t, dec_t, 0.15)
    cosd = np.cos(np.radians(dec))
    out = {
        "source_radec_t": [round(ra_t, 7), round(dec_t, 7)],
        "e_minus_t_offset_arcsec": round(
            float(np.hypot((ra_e - ra_t) * cosd, dec_e - dec_t)) * 3600, 3
        ),
        "flag_to_source_arcsec": round(float(np.hypot((ra - ra_t) * cosd, dec - dec_t)) * 3600, 3),
    }
    for k, (rk, dk) in {"E": (ra_e, dec_e), "T": (ra_t, dec_t)}.items():
        b = boxes[k]
        x, y = (float(v) for v in b["wcs"].all_world2pix(rk, dk, 0))
        f = aperture_photometry(
            b["sci"], x, y, STATIC_AP_ARCSEC / b["scale"], (1.0 / b["scale"], 1.5 / b["scale"])
        )
        out[f"flux_{k}_njy"] = round(float(f * b["pixar_sr"] * NJY_PER_MJYSR_SR), 1)
    out["flux_e_minus_t_njy"] = round(out["flux_E_njy"] - out["flux_T_njy"], 1)
    out["flux_change_frac"] = round(out["flux_e_minus_t_njy"] / out["flux_T_njy"], 4)
    return out


def track_points(p0: tuple[float, float], p1: tuple[float, float], n: int = N_TRACK) -> np.ndarray:
    """``n`` (ra, dec) points evenly from ``p0`` to ``p1`` (linear; fine over < 1 arcsec)."""
    t = np.linspace(0.0, 1.0, n)[:, None]
    return np.asarray(p0, float) * (1 - t) + np.asarray(p1, float) * t


def grid_points(
    track: np.ndarray, radius: float = SEARCH_ARCSEC, step: float = GRID_ARCSEC
) -> np.ndarray:
    """(ra, dec) of a ``step`` grid about the track midpoint, within ``radius`` of the segment."""
    ra0, dec0 = track.mean(axis=0)
    cosd = np.cos(np.radians(dec0))

    def to_xy(p):
        return np.column_stack([(p[:, 0] - ra0) * cosd * 3600.0, (p[:, 1] - dec0) * 3600.0])

    a, b = to_xy(track[[0, -1]])
    half = np.linalg.norm(b - a) / 2 + radius
    k = int(np.ceil(half / step))
    gx, gy = np.meshgrid(np.arange(-k, k + 1) * step, np.arange(-k, k + 1) * step)
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    ab = b - a
    if ab @ ab == 0:  # a single point: a disc
        t = np.zeros(len(pts))
    else:
        t = np.clip(((pts - a) @ ab) / (ab @ ab), 0.0, 1.0)
    dist = np.linalg.norm(pts - (a + t[:, None] * ab), axis=1)
    keep = pts[dist <= radius + 1e-9]
    return np.column_stack([ra0 + keep[:, 0] / 3600.0 / cosd, dec0 + keep[:, 1] / 3600.0])


def inject_gaussian(
    sci: np.ndarray, x: float, y: float, total: float, fwhm_px: float
) -> np.ndarray:
    """Copy of ``sci`` with a circular Gaussian of integral ``total`` (pixel units) at (x, y)."""
    s = fwhm_px / 2.3548
    yy, xx = np.mgrid[0 : sci.shape[0], 0 : sci.shape[1]]
    g = np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * s * s))
    return sci + total * g / g.sum()


def cloud_uri(obs_id: str) -> str:
    """S3 URI of the level-3 ``_i2d.fits`` of ``obs_id`` (MAST path lookup)."""
    from jwst_anomaly import query

    obs = query.mast_client().query_criteria(obs_id=obs_id)
    prods = query.list_products(obs, subgroups=("I2D",), cloud_uris=True)
    uris = [
        u
        for f, u in zip(
            query.str_values(prods["productFilename"]),
            query.str_values(prods["cloud_uri"]),
            strict=True,
        )
        if f == f"{obs_id}_i2d.fits" and u
    ]
    if not uris:
        raise RuntimeError(f"no cloud URI for {obs_id}")
    return uris[0]


def measure(box: dict, radec: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    """Aperture fluxes (nJy) at ``radec``, the box's empirical sigma (nJy) and pixel positions."""
    scale = box["scale"]
    to_njy = box["pixar_sr"] * NJY_PER_MJYSR_SR
    r_px, ann = RADIUS_ARCSEC / scale, (ANNULUS[0] / scale, ANNULUS[1] / scale)
    xy = np.column_stack(box["wcs"].all_world2pix(radec[:, 0], radec[:, 1], 0))
    f = np.array([aperture_photometry(box["sci"], x, y, r_px, ann) for x, y in xy]) * to_njy
    rng = np.random.default_rng(zlib.crc32(box.get("name", "").encode()))
    err = empirical_error(box["sci"], r_px, ann, rng) * to_njy
    return f, float(err), xy


def stamps(boxes: dict, centre: tuple[float, float], track: np.ndarray, path: Path) -> None:
    """E and T stamps side by side (pixel grids differ, so no E - T image); track marked."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(boxes), figsize=(4.2 * len(boxes), 4.2))
    for ax, (name, box) in zip(np.atleast_1d(axes), boxes.items(), strict=True):
        sci = box["sci"] * box["pixar_sr"] * NJY_PER_MJYSR_SR
        x, y = box["wcs"].all_world2pix(centre[0], centre[1], 0)
        h = int(np.ceil(2.0 / box["scale"]))
        xi, yi = int(round(float(x))), int(round(float(y)))
        sub = sci[max(yi - h, 0) : yi + h + 1, max(xi - h, 0) : xi + h + 1]
        lo, hi = np.nanpercentile(sub, [5, 99.5])
        ax.imshow(sub, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        tx, ty = box["wcs"].all_world2pix(track[:, 0], track[:, 1], 0)
        ax.plot(tx - max(xi - h, 0), ty - max(yi - h, 0), "r.", ms=3)
        r = SEARCH_ARCSEC / box["scale"]
        mx, my = np.mean(tx) - max(xi - h, 0), np.mean(ty) - max(yi - h, 0)
        ax.add_patch(plt.Circle((mx, my), r, fill=False, color="c", lw=0.6))
        ax.set_title(f"{name} (4 arcsec stamp, nJy/pix)", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    post = json.loads((args.out.parent / "posthoc.json").read_text())
    ends = []
    for z in Z_TRACK:
        rows = [r for r in post["rows"] if abs(r["z_s"] - z) < 1e-9 and r["testable_epochs_mjd"]]
        if len(rows) != 1:
            raise RuntimeError(f"expected one testable posthoc row at z_s {z}")
        ends.append((rows[0]["ra"], rows[0]["dec"]))
    track = track_points(*ends)
    grid = grid_points(track)
    centre = tuple(track.mean(axis=0))
    uris = {"E": cloud_uri(EPOCH_OBS), "T": cloud_uri(TEMPLATE_OBS)}
    boxes, res = {}, {}
    for k, uri in uris.items():
        b = read_box(uri, centre[0], centre[1], BOX_HALF_ARCSEC)
        b["name"] = k
        boxes[k] = b
        f, err, _ = measure(b, grid)
        ft, _, _ = measure(b, track)
        res[k] = {"grid": f, "track": ft, "err": err}
    sig = float(np.hypot(res["E"]["err"], res["T"]["err"]))
    delta = res["E"]["grid"] - res["T"]["grid"]
    # injection control: one source at a time into E, measured at the same track point
    e = boxes["E"]
    to_njy = e["pixar_sr"] * NJY_PER_MJYSR_SR
    txy = np.column_stack(e["wcs"].all_world2pix(track[:, 0], track[:, 1], 0))
    r_px = RADIUS_ARCSEC / e["scale"]
    ann = (ANNULUS[0] / e["scale"], ANNULUS[1] / e["scale"])
    inj_e = []
    for x, y in txy:
        sci = inject_gaussian(e["sci"], x, y, INJ_NJY / to_njy, INJ_FWHM_ARCSEC / e["scale"])
        inj_e.append(aperture_photometry(sci, x, y, r_px, ann) * to_njy)
    inj_delta = np.array(inj_e) - res["T"]["track"]
    nan_e = float(np.mean(~np.isfinite(res["E"]["grid"])))
    nan_t = float(np.mean(~np.isfinite(res["T"]["grid"])))
    verdict = outcome(delta, sig, nan_e, nan_t, inj_delta, sig)
    # TNS position itself (not part of the rule)
    tns = np.array([TNS_RADEC])
    tb = {k: read_box(u, *TNS_RADEC, BOX_HALF_ARCSEC) for k, u in uris.items()}
    tns_f = {k: measure(b | {"name": k}, tns)[:2] for k, b in tb.items()}
    stamps(boxes, centre, track, args.out / "stamps_track.png")
    stamps(tb, TNS_RADEC, tns, args.out / "stamps_tns.png")
    # post hoc (not in PLAN.md): every posthoc.json counter-image whose delay is in [0, t_E - t_0]
    from concurrent.futures import ThreadPoolExecutor

    if not np.any(np.isclose(post["epochs_mjd"], EPOCH_MJD, atol=0.01)):
        raise RuntimeError(f"EPOCH_MJD {EPOCH_MJD} is not an epoch of posthoc.json")
    t_e = EPOCH_MJD
    span = t_e - T0_MJD
    f_now = float(tns_f["E"][0][0] - tns_f["T"][0][0])
    rows = [r for r in post["rows"] if 0.0 <= r["delay_days"] <= span]

    def one(r: dict) -> dict:
        thr = POSTHOC_MU_FACTOR * abs(r["mu_ratio"]) * f_now
        row = {
            "z_s": r["z_s"],
            "sep_from_tns_arcsec": r["sep_from_tns_arcsec"],
            "mu_ratio": r["mu_ratio"],
            "delay_days": r["delay_days"],
            "phase_in_e_days": round(span - r["delay_days"], 1),
            "threshold_njy": round(thr, 1),
        }
        pt = np.array([[r["ra"], r["dec"]]])
        disc = grid_points(np.repeat(pt, 2, axis=0))
        m, bx = {}, {}
        try:
            for k, u in uris.items():
                bx[k] = read_box(u, r["ra"], r["dec"], BOX_HALF_ARCSEC)
                m[k] = measure(bx[k] | {"name": k}, disc)[:2]
        except ValueError as exc:  # outside a mosaic
            return row | {"outcome": "untestable", "error": str(exc)[:120]}
        d = m["E"][0] - m["T"][0]
        sg = float(np.hypot(m["E"][1], m["T"][1]))
        nanf = float(max(np.mean(~np.isfinite(m["E"][0])), np.mean(~np.isfinite(m["T"][0]))))
        # injection of the expected counter-image (aperture flux thr / factor) at the prediction
        e = bx["E"]
        to = e["pixar_sr"] * NJY_PER_MJYSR_SR
        x, y = (float(v) for v in e["wcs"].all_world2pix(r["ra"], r["dec"], 0))
        sci = inject_gaussian(
            e["sci"], x, y, thr / POSTHOC_MU_FACTOR / to, INJ_FWHM_ARCSEC / e["scale"]
        )
        f_inj = aperture_photometry(
            sci,
            x,
            y,
            RADIUS_ARCSEC / e["scale"],
            (ANNULUS[0] / e["scale"], ANNULUS[1] / e["scale"]),
        )
        inj_d = float(f_inj * to - measure(bx["T"] | {"name": "T"}, pt)[0][0])
        verdict_ph = posthoc_outcome(d, sg, nanf, thr, inj_d)
        row |= {"sigma_njy": round(sg, 2), "nan_frac": nanf, "inj_delta_njy": round(inj_d, 1)}
        if np.isfinite(d).any():
            i = int(np.nanargmax(d))
            row |= {
                "delta_max_njy": round(float(d[i]), 1),
                "delta_max_snr": round(float(d[i] / sg), 2),
            }
        row["outcome_raw"] = verdict_ph
        sc = None
        if verdict_ph == "present":
            row["delta_max_radec"] = np.round(disc[i], 7).tolist()
            try:
                sc = static_check(bx, *disc[i])
            except ValueError as exc:
                sc = None
                row["static_check_error"] = str(exc)[:120]
            if sc:
                cosd = np.cos(np.radians(r["dec"]))
                ra_s, dec_s = sc["source_radec_t"]
                sc["predicted_to_source_arcsec"] = round(
                    float(np.hypot((r["ra"] - ra_s) * cosd, r["dec"] - dec_s)) * 3600, 3
                )
                row["static_check"] = sc
        row["outcome"] = classify_flag(verdict_ph, sc, thr)
        return row | {"ra": r["ra"], "dec": r["dec"]}

    with ThreadPoolExecutor(8) as ex:
        ph = list(ex.map(one, rows))
    flagged = np.array([[p["ra"], p["dec"]] for p in ph if p.get("outcome_raw") == "present"])
    if len(flagged):
        fc = tuple(flagged.mean(axis=0))
        fb = {k: read_box(u, fc[0], fc[1], BOX_HALF_ARCSEC) for k, u in uris.items()}
        stamps(fb, fc, flagged, args.out / "stamps_posthoc_flags.png")
    i_max = int(np.nanargmax(delta))
    out = {
        "provenance": {
            "positions": "model_prediction",
            "fluxes": "derived",
            "thresholds": "assumption",
        },
        "plan": "PLAN.md",
        "uris": uris,
        "pixel_scale_arcsec": {k: round(b["scale"], 5) for k, b in boxes.items()},
        "track_radec": np.round(track, 7).tolist(),
        "n_grid": int(len(grid)),
        "sigma_njy": {
            "E": round(res["E"]["err"], 2),
            "T": round(res["T"]["err"], 2),
            "delta": round(sig, 2),
        },
        "t_ap_njy": T_AP_NJY,
        "inj_njy": INJ_NJY,
        "nan_frac": {"E": nan_e, "T": nan_t},
        "track_flux_njy": {k: np.round(res[k]["track"], 1).tolist() for k in res},
        "track_delta_njy": np.round(res["E"]["track"] - res["T"]["track"], 1).tolist(),
        "inj_delta_njy": np.round(inj_delta, 1).tolist(),
        "inj_passed": int(np.sum((inj_delta >= T_AP_NJY) & (inj_delta >= NSIGMA * sig))),
        "grid_delta_max_njy": round(float(delta[i_max]), 1),
        "grid_delta_max_radec": np.round(grid[i_max], 7).tolist(),
        "grid_delta_max_snr": round(float(delta[i_max] / sig), 2),
        "grid_delta_min_njy": round(float(np.nanmin(delta)), 1),
        "n_grid_over_t_ap": int(np.nansum(delta >= T_AP_NJY)),
        "n_grid_over_5sigma": int(np.nansum(delta >= NSIGMA * sig)),
        "outcome": verdict,
        "tns_position": {
            k: {"flux_njy": round(float(v[0][0]), 1), "err_njy": round(v[1], 1)}
            for k, v in tns_f.items()
        },
        "posthoc": {
            "note": "not in PLAN.md; single-peaked light curve ASSUMPTION; outcome_raw from "
            "posthoc_outcome (one-point injection gate), outcome after classify_flag; threshold = "
            f"{POSTHOC_MU_FACTOR} x |mu_ratio| x (F_E - F_T at the TNS position)",
            "t_e_mjd": t_e,
            "t_e_minus_t0_days": round(span, 2),
            "f_now_njy": round(f_now, 1),
            "rows": ph,
        },
    }
    (args.out / "result.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

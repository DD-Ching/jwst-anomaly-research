"""Vet Abell 2744 system 700.1 (UNCOVER 78): colours of the mirrored a/b knots, perturber check.

Plan (frozen before the run): ``results/lens_consistency/a2744_s700_vet/PLAN.md``.

* Colours: high-pass forced photometry (r = 0.15") of UNCOVER 78.x1 / 78.x2 in seven NIRCam bands of
  ``jw02561-o001_t003`` (S3 byte ranges), noise from random apertures, and a χ² of the colour
  difference with the flux ratio free. Validated first on catalogued UNCOVER multiple-image pairs.
* Perturber: DJA v7.2 eazy zout sources within 5" of the 78.2 a/b midpoint, with a crude F200W
  aperture magnitude from the mosaic.

Usage: ``python scripts/a2744_s700_vet.py [--out DIR]``, then ``--sheet-only`` (sheet.png) and
``--model-check`` (post hoc Bergamini+23 components); both read the written summary.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a2744_followups as af  # noqa: E402
import a2744_uncover as au  # noqa: E402
import lens_consistency as lc  # noqa: E402

from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

BANDS = ("f115w", "f150w", "f200w", "f277w", "f356w", "f410m", "f444w")
REF_BAND = "f200w"
HALF = 2.0  # stamp half-width, arcsec
R_AP = 0.15  # knot aperture radius, arcsec (ASSUMPTION)
RECENTRE = 0.1  # re-centring radius on the F200W high-pass peak, arcsec (ASSUMPTION)
N_RANDOM = 200
RANDOM_RING = (0.6, 1.8)  # random-aperture centres, arcsec from the target (ASSUMPTION)
P_DIFFER = 0.01  # a pair whose colour χ² has p below this "differs" (ASSUMPTION)
SNR_MIN = 10.0  # known pairs need this F200W S/N in both images (ASSUMPTION)
MAX_KNOWN = 20
KNOWN_GATE = 0.8  # fraction of known pairs that must read "consistent" (ASSUMPTION)
TARGET = "78"
TARGET_PAIRS = (("78.11", "78.12"), ("78.21", "78.22"), ("78.31", "78.32"))
PERTURBER_R = 5.0  # arcsec around the 78.2 a/b midpoint (ASSUMPTION)
Z_CLUSTER = au.Z_LENS
DZ_MEMBER = 0.1  # |z - z_cluster| for a member (ASSUMPTION)
MAG_PERTURBER = 22.0  # F200W AB brighter than this (ASSUMPTION)
R_TOTAL = 0.5  # crude total-magnitude aperture, arcsec (ASSUMPTION)
R_BG = 1.5  # background = stamp median beyond this radius, arcsec (ASSUMPTION)
SEED = 20261020
DEFAULT_OUT = Path("results/lens_consistency/a2744_s700_vet")


def aperture_noise(hp: np.ndarray, xx: np.ndarray, yy: np.ndarray, rng) -> float:
    """1.4826 × MAD of r = R_AP aperture sums at random centres 0.6–1.8" from the target."""
    r = np.sqrt(rng.uniform(RANDOM_RING[0] ** 2, RANDOM_RING[1] ** 2, N_RANDOM))
    phi = rng.uniform(0.0, 2.0 * np.pi, N_RANDOM)
    sums = []
    for cx, cy in zip(r * np.cos(phi), r * np.sin(phi), strict=True):
        ap = np.hypot(xx - cx, yy - cy) <= R_AP
        sums.append(np.nansum(hp[ap]))
    s = np.asarray(sums)
    return float(1.4826 * np.median(np.abs(s - np.median(s))))


def recentre(hp: np.ndarray, xx: np.ndarray, yy: np.ndarray) -> tuple[float, float]:
    """Offset (West, North; arcsec) of the high-pass peak within RECENTRE of the target."""
    inner = np.where((np.hypot(xx, yy) <= RECENTRE) & np.isfinite(hp), hp, -np.inf)
    if not np.isfinite(inner).any():
        return 0.0, 0.0
    k = int(np.argmax(inner))
    return float(xx.flat[k]), float(yy.flat[k])


def pair_chi2(fa, sa, fb, sb) -> dict:
    """χ² of f_a = r f_b over bands (r free, error on both sides); p from χ² with n - 1 dof."""
    fa, sa, fb, sb = (np.asarray(v, float) for v in (fa, sa, fb, sb))
    ok = np.isfinite(fa) & np.isfinite(fb) & (sa > 0) & (sb > 0)
    fa, sa, fb, sb = fa[ok], sa[ok], fb[ok], sb[ok]
    n = int(ok.sum())
    if n < 2 or not np.any(fb > 0) or not np.any(fa > 0):
        return {"n_bands": n, "ratio": None, "chi2": None, "p": None, "ratio_at_grid_edge": None}

    def chi2(r):
        r = np.asarray(r, float)[..., None]
        return np.sum((fa - r * fb) ** 2 / (sa**2 + r**2 * sb**2), axis=-1)

    # coarse log grid over 1e-3..1e3, then a bounded refinement around the best grid point
    grid = np.logspace(-3, 3, 601)
    k = int(np.argmin(chi2(grid)))
    lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
    res = optimize.minimize_scalar(chi2, bounds=(lo, hi), method="bounded")
    c = float(res.fun)
    return {
        "n_bands": n,
        "ratio": float(res.x),
        "chi2": c,
        "p": float(stats.chi2.sf(c, n - 1)),
        "ratio_at_grid_edge": k in (0, len(grid) - 1),  # the true minimum may lie outside
    }


def differs(p: float | None, at_edge: bool | None = False) -> bool | None:
    """None (not measurable) without a p value or with a ratio fit at the grid edge (post hoc)."""
    return None if p is None or at_edge else bool(p < P_DIFFER)


def classify_perturber(member: bool | None, mag: float) -> bool | None:
    """True: member brighter than MAG_PERTURBER; None ("not assessed", post hoc): membership
    unknown (no usable z) or a member without a measurable magnitude; False otherwise, including
    any source too faint to be a perturber whatever its z (post hoc, like the None state)."""
    if np.isfinite(mag) and mag >= MAG_PERTURBER:
        return False
    if member is None:
        return None
    if not member:
        return False
    if not np.isfinite(mag):
        return None
    return True  # a finite magnitude brighter than MAG_PERTURBER (fainter ones returned above)


def tri_any(flags: list[bool | None]) -> bool | None:
    """True if any flag is True; None if none is but some is None; else False."""
    if any(f is True for f in flags):
        return True
    return None if any(f is None for f in flags) else False


def decide(gate_passed: bool, any_differs: bool | None, perturber: bool | None) -> str:
    """Decision function of PLAN.md. Post hoc additions (after /code-review): the
    ``any_differs is None`` branch (a target pair not measurable is never "consistent") and the
    ``perturber is None`` branch and tail (a nearby source not assessed is never "no perturber")."""
    if not gate_passed:
        if perturber is None:
            tail = "member perturber not assessed"
        elif perturber:
            tail = "member perturber present"
        else:
            tail = "no member perturber"
        return f"undetermined (colour test not validated; {tail})"
    if any_differs:
        return "one image (H1)"
    if any_differs is None:
        return "undetermined (target colours not measurable)"
    if perturber:
        return "ordinary perturber possible (H2)"
    if perturber is None:
        return "undetermined (member perturber not assessed)"
    return "undetermined (H3)"


def measure(positions: dict[str, tuple[float, float]]) -> dict[str, dict[str, dict]]:
    """Per band and image: flux, noise and S/N at the re-centred position (S3 range reads)."""
    # re-centre on F200W first; the same offset is used in every band
    shift: dict[str, tuple[float, float]] = {}
    stamp, close = lc.image_stamper(af.I2D.format(REF_BAND), HALF)
    try:
        for iid, (ra, dec) in positions.items():
            sci, _, xx, yy = stamp(ra, dec)
            shift[iid] = recentre(au._highpass(sci, xx, yy), xx, yy)
    finally:
        close()

    def one_band(band: str) -> tuple[str, dict]:
        rng = np.random.default_rng(SEED)
        stamp, close = lc.image_stamper(af.I2D.format(band), HALF)
        res = {}
        try:
            for iid, (ra, dec) in positions.items():
                sci, _, xx, yy = stamp(ra, dec)
                hp = au._highpass(sci, xx, yy)
                cx, cy = shift[iid]
                ap = np.hypot(xx - cx, yy - cy) <= R_AP
                f = float(np.nansum(hp[ap])) if np.isfinite(hp[ap]).any() else np.nan
                s = aperture_noise(hp, xx, yy, rng)
                res[iid] = {"flux": f, "sigma": s, "snr": f / s if s > 0 else np.nan}
        finally:
            close()
        return band, res

    with ThreadPoolExecutor(len(BANDS)) as pool:
        out = dict(pool.map(one_band, BANDS))
    for band in out:
        for iid in out[band]:
            out[band][iid]["shift_arcsec"] = [round(v, 3) for v in shift[iid]]
    return out


def pair_result(phot: dict, a: str, b: str) -> dict:
    fa = [phot[band][a]["flux"] for band in BANDS]
    fb = [phot[band][b]["flux"] for band in BANDS]
    sa = [phot[band][a]["sigma"] for band in BANDS]
    sb = [phot[band][b]["sigma"] for band in BANDS]
    r = pair_chi2(fa, sa, fb, sb)
    return {"a": a, "b": b, **r, "differs": differs(r["p"], r["ratio_at_grid_edge"])}


def known_pairs(images: Table) -> list[tuple[str, str]]:
    """First two images of every spec-z UNCOVER system except the target's source."""
    raw = images
    pairs = []
    seen = set()
    for row in raw:
        iid = repr(float(row["ID"]))
        sys_ = au.uncover_system(iid)
        if sys_.split(".")[0] == TARGET or str(row["Redshift_type"]).strip() != "spec":
            continue
        if sys_ in seen:
            continue
        members = [repr(float(i)) for i in raw["ID"] if au.uncover_system(repr(float(i))) == sys_]
        seen.add(sys_)
        if len(members) >= 2:
            pairs.append((members[0], members[1]))
    return pairs


def total_mag(ra: float, dec: float, stamp) -> float:
    """Crude F200W AB magnitude: r = 0.5" aperture minus the stamp median beyond 1.5" (MJy/sr)."""
    sci, _, xx, yy = stamp(ra, dec)
    rr = np.hypot(xx, yy)
    bg = np.nanmedian(sci[rr > R_BG])
    # pixel solid angle; hypot, because the mosaic may be rotated relative to north
    pix_sr = (float(np.hypot(xx[0, 1] - xx[0, 0], yy[0, 1] - yy[0, 0])) / 206264.806) ** 2
    fjy = float(np.nansum(sci[rr <= R_TOTAL] - bg)) * 1e6 * pix_sr
    return float(-2.5 * np.log10(fjy / 3631.0)) if fjy > 0 else np.nan


def _r(v, nd: int = 3):
    return None if v is None or not np.isfinite(v) else round(float(v), nd)


def perturbers(zout: Table, ra0: float, dec0: float) -> list[dict]:
    """DJA zout sources within PERTURBER_R of (ra0, dec0) with z, member flag and F200W mag."""
    d = au._sep_arcsec(zout["ra"], zout["dec"], ra0, dec0)
    near = np.flatnonzero(d <= PERTURBER_R)
    stamp, close = lc.image_stamper(af.I2D.format(REF_BAND), 2.5)
    rows = []
    try:
        for k in near[np.argsort(d[near])]:
            zs = float(zout["z_spec"][k])
            z = zs if zs > 0 else float(zout["z_phot"][k])
            mag = total_mag(float(zout["ra"][k]), float(zout["dec"][k]), stamp)
            # failed eazy fit (z < 0) or non-finite z: membership unknown (None), not "non-member"
            usable = bool(np.isfinite(z) and z >= 0)
            member = bool(abs(z - Z_CLUSTER) <= DZ_MEMBER) if usable else None
            zp = float(zout["z_phot"][k])
            eazy_ok = bool(np.isfinite(zp) and zp >= 0)  # the interval is eazy's, even for spec-z
            row = {"dja_id": int(zout["id"][k]), "sep_arcsec": _r(d[k], 2)}
            # eazy's failure sentinels are never stored as a z or an interval
            row["z"], row["z_type"] = (
                (_r(z, 3), "spec" if zs > 0 else "phot") if usable else (None, None)
            )
            row["z160_z840_eazy"] = (
                [_r(float(zout["z160"][k]), 3), _r(float(zout["z840"][k]), 3)] if eazy_ok else None
            )
            rows.append(
                row
                | {
                    "f200w_mag_crude": _r(mag, 2),
                    "member": member,
                    "perturber": classify_perturber(member, mag),
                }
            )
    finally:
        close()
    return rows


def sheet(pos: dict, summary: dict, zout: Table, path: Path) -> None:
    """F115W / F200W / F444W high-pass and raw stamps (6" wide) around the 78.2 a/b midpoint with
    the 78 knots (circles) and the DJA sources within PERTURBER_R (squares, members red)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ra0 = float(np.mean([pos["78.21"][0], pos["78.22"][0]]))
    dec0 = float(np.mean([pos["78.21"][1], pos["78.22"][1]]))
    cos = np.cos(np.deg2rad(dec0))
    near = {p["dja_id"]: p for p in summary["perturbers"]}
    sel = np.isin(np.asarray(zout["id"]), list(near))
    bands = ("f115w", "f200w", "f444w")
    fig, axes = plt.subplots(2, len(bands), figsize=(3.4 * len(bands), 7.0))
    for j, band in enumerate(bands):
        stamp, close = lc.image_stamper(af.I2D.format(band), 3.5)
        try:
            sci, _, xx, yy = stamp(ra0, dec0)
        finally:
            close()
        for i, img in enumerate((sci, au._highpass(sci, xx, yy))):
            ax = axes[i, j]
            med = np.nanmedian(img)
            sig = 1.4826 * np.nanmedian(np.abs(img - med)) or 1.0
            ax.pcolormesh(xx, yy, np.arcsinh((img - med) / sig), cmap="gray_r", vmin=-1, vmax=5)
            for iid, (ra, dec) in pos.items():
                if iid.split(".")[0] != TARGET:
                    continue
                dx, dy = -(ra - ra0) * cos * 3600.0, (dec - dec0) * 3600.0
                if max(abs(dx), abs(dy)) < 3.0:
                    ax.scatter(dx, dy, s=50, facecolors="none", edgecolors="tab:blue")
                    if i == 0:
                        ax.annotate(iid, (dx, dy), fontsize=6, color="tab:blue")
            for r in zout[sel]:
                dx = -(float(r["ra"]) - ra0) * cos * 3600.0
                dy = (float(r["dec"]) - dec0) * 3600.0
                m = near[int(r["id"])]["member"]
                col = {True: "tab:red", False: "tab:orange"}.get(m, "tab:purple")
                ax.scatter(dx, dy, marker="s", s=70, facecolors="none", edgecolors=col)
                if i == 0:
                    zr = near[int(r["id"])]["z"]
                    ax.annotate("z?" if zr is None else f"z{zr}", (dx, dy), fontsize=6, color=col)
            ax.set_xlim(-3, 3)
            ax.set_ylim(-3, 3)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{band.upper()} {'raw' if i == 0 else 'high-pass'}", fontsize=9)
    fig.suptitle(
        "UNCOVER 78 (700.1) a/b: knots (o); DJA sources (□): red member, orange non-member,\n"
        'purple z unknown; 6", West = +x'
    )
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def lenstool_potentials(par_text: str) -> dict[str, dict]:
    """Lenstool ``potentiel`` blocks: centre (model frame, arcsec), v_disp, cut_radius_kpc."""
    import re

    out = {}
    for blk in re.split(r"\npotentiel\s+", par_text)[1:]:
        name = blk.split()[0]
        body = blk.split("\nend")[0]

        def g(key: str, body: str = body) -> float | None:
            m = re.search(rf"\b{key}\s+(\S+)", body)
            return float(m.group(1)) if m else None

        out[name] = {k: g(k) for k in ("x_centre", "y_centre", "v_disp", "cut_radius_kpc")}
    return out


def model_components(pert: list[dict], zout: Table, radius: float = 1.0) -> list[dict]:
    """Post hoc: Bergamini+23 potentials within ``radius`` of each member perturber's position."""
    model, files, _ = lc.load_model("abell2744-bergamini23")
    pots = lenstool_potentials(Path(files["best.par"]).read_text())
    names = list(pots)
    px = np.array([pots[n]["x_centre"] for n in names], float)
    py = np.array([pots[n]["y_centre"] for n in names], float)
    rows = []
    for p in pert:
        if p["perturber"] is not True:
            continue
        r = zout[np.asarray(zout["id"]) == p["dja_id"]][0]
        x, y = model.to_frame(np.array([float(r["ra"])]), np.array([float(r["dec"])]))
        d = np.hypot(px - x[0], py - y[0])
        comps = [
            {"potential": names[k], "offset_arcsec": _r(d[k], 2), **pots[names[k]]}
            for k in np.argsort(d)
            if d[k] <= radius
        ]
        rows.append({"dja_id": p["dja_id"], "components": comps})
    return rows


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--sheet-only", action="store_true", help="redraw sheet.png from summary.json")
    ap.add_argument(
        "--model-check",
        action="store_true",
        help="post hoc: add the Bergamini+23 components at each perturber to summary.json",
    )
    args = ap.parse_args(argv)
    raw = Table.read(fetch_catalog(*au.UNCOVER_V2["images"]))
    pos = {repr(float(r["ID"])): (float(r["RA"]), float(r["Dec"])) for r in raw}
    if args.model_check:
        summary = json.loads((args.out / "summary.json").read_text())
        zout = Table.read(fetch_catalog(*af.ZOUT))
        summary["post_hoc_bergamini23_components"] = model_components(summary["perturbers"], zout)
        (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
        print(json.dumps(summary["post_hoc_bergamini23_components"], indent=1))
        return summary
    if args.sheet_only:
        summary = json.loads((args.out / "summary.json").read_text())
        sheet(pos, summary, Table.read(fetch_catalog(*af.ZOUT)), args.out / "sheet.png")
        return summary
    kp = known_pairs(raw)
    ids = sorted({i for p in kp for i in p} | {i for p in TARGET_PAIRS for i in p})
    phot = measure({i: pos[i] for i in ids})

    def snr(i: str) -> float:
        return float(phot[REF_BAND][i]["snr"])

    # NaN S/N (no data) never passes: min() would let it through
    bright = [p for p in kp if all(snr(i) >= SNR_MIN for i in p)][:MAX_KNOWN]
    known = [pair_result(phot, a, b) for a, b in bright]
    judged = [k for k in known if k["differs"] is not None]
    frac = sum(not k["differs"] for k in judged) / len(judged) if judged else 0.0
    gate = bool(judged and frac >= KNOWN_GATE)
    target = [pair_result(phot, a, b) for a, b in TARGET_PAIRS]
    any_diff = tri_any([t["differs"] for t in target])
    ra_m = np.mean([pos["78.21"][0], pos["78.22"][0]])
    dec_m = np.mean([pos["78.21"][1], pos["78.22"][1]])
    pert = perturbers(Table.read(fetch_catalog(*af.ZOUT)), float(ra_m), float(dec_m))
    has_pert = tri_any([p["perturber"] for p in pert])
    out = {
        "plan": "PLAN.md",
        "bands": list(BANDS),
        "params": {
            "r_ap": R_AP,
            "recentre": RECENTRE,
            "n_random": N_RANDOM,
            "random_ring": list(RANDOM_RING),
            "p_differ": P_DIFFER,
            "snr_min": SNR_MIN,
            "known_gate": KNOWN_GATE,
            "perturber_r": PERTURBER_R,
            "dz_member": DZ_MEMBER,
            "mag_perturber": MAG_PERTURBER,
        },
        "n_known_candidates": len(kp),  # knot sub-systems "S.K", not independent sources
        "n_tested_sources": len({a.split(".")[0] for a, _ in bright}),
        "known_pairs": known,
        "known_consistent_fraction": _r(frac),
        "colour_gate_passed": gate,
        "target_pairs": target,
        "any_target_pair_differs": any_diff,
        "perturbers": pert,
        "perturber_found": has_pert,
        "not_assessed": [p["dja_id"] for p in pert if p["perturber"] is None],
        "outcome": decide(gate, any_diff, has_pert),
        "photometry": {
            band: {
                i: {k: (_r(v, 5) if k != "shift_arcsec" else v) for k, v in d.items()}
                for i, d in phot[band].items()
                if i.split(".")[0] == TARGET
            }
            for band in BANDS
        },
        "provenance": "derived (high-pass aperture photometry of observed mosaics; z: catalogue)",
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "photometry"}, indent=1))
    return out


if __name__ == "__main__":
    main()

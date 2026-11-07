"""Image-level two-epoch centroids of SMACS ``f200w_1816`` / ``f200w_384`` (2736 o001 -> 6882 o057).

Plan: ``results/smacs_jwst3/epochs_1816_384_img/PLAN.md`` (frozen in the claim commit). Cutouts are
S3 byte-range reads of the level-3 ``_i2d.fits`` (``jwst_anomaly.cutouts``); centroids are
``photutils.centroids.centroid_2dg`` on a small box; the frame tie and the null come from peers
(compact sources of similar magnitude and size). ``decide`` encodes the plan's outcome rules.

    python scripts/smacs_epochs_img.py <scratch dir> --cat-dir <dir with the CATs> \\
        --out results/smacs_jwst3/epochs_1816_384_img
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from photutils.centroids import centroid_2dg
from scipy.ndimage import fourier_shift

from jwst_anomaly import catalog, cutouts

sys.path.insert(0, str(Path(__file__).resolve().parent))
from epoch_compare import mutual_matches, offsets_mas  # noqa: E402
from smacs_epochs_1816_384 import BASELINE_YR, E1, E2, TARGETS  # noqa: E402
from smacs_epochs_1816_384 import _angle_deg as angle_deg  # noqa: E402

BANDS = ("f150w", "f444w")
S3 = {"e1": "jw02736/L3/t/o001/", "e2": "jw06882/L3/t/o057/"}
NAMES = {"e1": E1, "e2": E2}
# Frozen plan settings (ASSUMPTIONs).
SIZE_ARCSEC = 1.5
PEAK_RADIUS_ARCSEC = 0.15
BOX = 7
PEER_RADIUS_ARCSEC = 60.0
PEER_DMAG = 1.0
PEER_SIZE_FRAC = 0.5
PEER_ISOLATION_ARCSEC = 1.0
MATCH_ARCSEC = 0.3
MIN_PEERS = 10
PCTL = 95.0
ANGLE_DEG = 45.0
INJECT_MAS = 20.0
INJECT_TOL_MAS = 5.0
ZOOM_ARCSEC = 0.6  # contact-sheet panel width (display only)


def centroid_sky(data: np.ndarray, wcs: WCS) -> SkyCoord | None:
    """2-D Gaussian centroid of the brightest pixel's box near the cutout centre, or None."""
    data = np.asarray(data, dtype=float)
    if not np.isfinite(data).any():
        return None
    _, med, _ = sigma_clipped_stats(data[np.isfinite(data)])
    img = np.where(np.isfinite(data), data - med, 0.0)
    ny, nx = img.shape
    cy, cx = (ny - 1) / 2, (nx - 1) / 2
    scale = float(np.mean([q.to_value(u.arcsec) for q in wcs.proj_plane_pixel_scales()]))
    rad = PEAK_RADIUS_ARCSEC / scale
    yy, xx = np.mgrid[:ny, :nx]
    near = np.hypot(yy - cy, xx - cx) <= rad
    k = np.argmax(np.where(near, img, -np.inf))
    py, px = np.unravel_index(k, img.shape)
    h = BOX // 2
    if py - h < 0 or px - h < 0 or py + h >= ny or px + h >= nx:
        return None
    sub = img[py - h : py + h + 1, px - h : px + h + 1]
    if sub[h, h] <= 0:
        return None
    try:
        x, y = centroid_2dg(sub)
    except Exception:  # noqa: BLE001 - a failed fit is "no centroid"
        return None
    if not (np.isfinite(x) and np.isfinite(y)) or not (0 <= x < BOX and 0 <= y < BOX):
        return None
    return wcs.pixel_to_world(px - h + x, py - h + y)


def shifted(data: np.ndarray, wcs: WCS, sky: SkyCoord, dmas: tuple[float, float]) -> np.ndarray:
    """``data`` Fourier-shifted so a source at ``sky`` moves by ``dmas`` (dRA cos Dec, dDec) mas."""
    new = sky.spherical_offsets_by(dmas[0] * u.mas, dmas[1] * u.mas)
    x0, y0 = wcs.world_to_pixel(sky)
    x1, y1 = wcs.world_to_pixel(new)
    filled = np.where(np.isfinite(data), data, 0.0)
    out = fourier_shift(np.fft.fft2(filled), (float(y1 - y0), float(x1 - x0)))
    return np.real(np.fft.ifft2(out))


def band_stats(target_off: np.ndarray | None, peer_off: np.ndarray) -> dict[str, Any] | None:
    """Residuals vs the peer tie (target: peer median; each peer: leave-one-out median)."""
    n = len(peer_off)
    if target_off is None or n < MIN_PEERS:
        return None
    tie = np.median(peer_off, axis=0)
    t_res = target_off - tie
    p_res = np.array(
        [np.hypot(*(peer_off[i] - np.median(np.delete(peer_off, i, 0), axis=0))) for i in range(n)]
    )
    return {
        "n_peers": n,
        "tie_mas": tie.tolist(),
        "target_residual_mas": t_res.tolist(),
        "target_residual_total_mas": float(np.hypot(*t_res)),
        "peer_residual_median_mas": float(np.median(p_res)),
        "peer_residual_p95_mas": float(np.percentile(p_res, PCTL)),
        "n_peers_with_larger_residual": int(np.sum(p_res > np.hypot(*t_res))),
        "label": "derived",
    }


def decide(per_band: dict[str, dict[str, Any] | None]) -> str:
    """Plan outcome rules (frozen); ``per_band[b]`` is ``band_stats`` output or None."""
    res = [per_band.get(b) for b in BANDS]
    if any(r is None for r in res):
        return "undetermined"
    tot = [r["target_residual_total_mas"] for r in res]
    vec = [np.array(r["target_residual_mas"]) for r in res]
    if all(t > r["peer_residual_p95_mas"] for t, r in zip(tot, res, strict=True)) and (
        angle_deg(*vec) <= ANGLE_DEG
    ):
        return "moves: star"
    if all(t <= r["peer_residual_median_mas"] for t, r in zip(tot, res, strict=True)):
        return "no motion detected"
    return "undetermined"


def select_peers(cat1: Table, target: SkyCoord, k_target: int, cat2: Table) -> list[int]:
    """Epoch-1 F150W rows that qualify as peers (plan 'Peers'), with a mutual epoch-2 match."""
    c1 = SkyCoord(cat1["ra"], cat1["dec"], unit="deg")
    c2 = SkyCoord(cat2["ra"], cat2["dec"], unit="deg")
    mag = np.asarray(cat1["aper50_abmag"], dtype=float)
    size = np.asarray(cat1["semimajor_sigma"], dtype=float)
    near = target.separation(c1).arcsec <= PEER_RADIUS_ARCSEC
    ok = (
        near
        & (np.abs(mag - mag[k_target]) <= PEER_DMAG)
        & (np.abs(size - size[k_target]) <= PEER_SIZE_FRAC * size[k_target])
    )
    _, d2, _ = c1.match_to_catalog_sky(c1, nthneighbor=2)
    ok &= d2.arcsec > PEER_ISOLATION_ARCSEC
    ok[k_target] = False
    i1, _ = mutual_matches(c1, c2, MATCH_ARCSEC)
    cand = np.flatnonzero(ok)
    return [int(i) for i in cand[np.isin(cand, i1)]]


def _read(path: str) -> tuple[np.ndarray, WCS]:
    with fits.open(path) as h:
        return h["SCI"].data.astype(float), WCS(h["SCI"].header)


def sheet(path: Path, files: dict, cents: dict, names: list[str]) -> None:
    """Target cutouts per band and epoch, both epochs' centroids marked (ZOOM_ARCSEC wide)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(len(names), 4, figsize=(10, 2.6 * len(names)), squeeze=False)
    for r, name in enumerate(names):
        for c, (b, ep) in enumerate((x, e) for x in BANDS for e in ("e1", "e2")):
            ax = axs[r, c]
            data, wcs = _read(files[(b, ep)][name])
            lo, hi = np.nanpercentile(data, [5, 99.5])
            ax.imshow(data, origin="lower", cmap="gray", vmin=lo, vmax=hi)
            for e, mk in (("e1", "c"), ("e2", "r")):
                sky = cents[(b, e)].get(name)
                if sky is not None:
                    x, y = wcs.world_to_pixel(sky)
                    ax.plot(x, y, "+", color=mk, ms=10)
            scale = np.mean([q.to_value(u.arcsec) for q in wcs.proj_plane_pixel_scales()])
            half = ZOOM_ARCSEC / 2 / scale
            cy, cx = (np.array(data.shape) - 1) / 2
            ax.set_xlim(cx - half, cx + half)
            ax.set_ylim(cy - half, cy + half)
            ax.set_title(f"{name} {b} {ep} (+ cyan e1, red e2)", fontsize=7)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("work", type=Path, help="scratch directory for cutout FITS files")
    ap.add_argument("--cat-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    prev = json.loads((args.out.parent / "epochs_1816_384" / "summary.json").read_text())["targets"]
    cat1 = catalog.load_pipeline_catalog(args.cat_dir / E1.format(b="f150w"))
    cat2 = catalog.load_pipeline_catalog(args.cat_dir / E2.format(b="f150w"))
    c1 = SkyCoord(cat1["ra"], cat1["dec"], unit="deg")
    # Positions to cut: each target (nearest epoch-1 F150W row) and its peers.
    uids, ras, decs, groups, target_sep = [], [], [], {}, {}
    for name in TARGETS:
        t = SkyCoord(prev[name]["ra"], prev[name]["dec"], unit="deg")
        sep = t.separation(c1).arcsec
        k = int(np.argmin(sep))
        if sep[k] > MATCH_ARCSEC:
            raise SystemExit(f"{name}: no epoch-1 F150W row within {MATCH_ARCSEC} arcsec")
        target_sep[name] = float(sep[k])
        peers = select_peers(cat1, c1[k], k, cat2)
        groups[name] = [f"{name}"] + [f"{name}_p{i}" for i in peers]
        for uid, i in zip(groups[name], [k, *peers], strict=True):
            uids.append(uid)
            ras.append(float(cat1["ra"][i]))
            decs.append(float(cat1["dec"][i]))
    targets = Table({"source_uid": uids, "ra": ras, "dec": decs})
    pos = dict(zip(uids, SkyCoord(ras, decs, unit="deg"), strict=True))
    cents: dict[tuple[str, str], dict[str, SkyCoord | None]] = {}
    inject: dict[str, list[list[float]]] = {b: [] for b in BANDS}
    peer_uids = {m for members in groups.values() for m in members[1:]}

    def fetch(key: tuple[str, str]) -> dict[str, str]:
        b, ep = key
        i2d = NAMES[ep].format(b=b).replace("_cat.ecsv", "_i2d.fits")
        tab = cutouts.make_cutouts(
            f"s3://stpubdata/jwst/public/{S3[ep]}{i2d}",
            targets,
            size_arcsec=SIZE_ARCSEC,
            out_dir=args.work,
            storage_options={"anon": True},
        )
        return dict(zip(tab["source_uid"], tab["path"], strict=True))

    keys = [(b, ep) for b in BANDS for ep in ("e1", "e2")]
    with ThreadPoolExecutor(max_workers=len(keys)) as pool:  # four independent S3 reads
        allfiles = dict(zip(keys, pool.map(fetch, keys), strict=True))
    for key, fl in allfiles.items():
        cents[key] = {uid: (centroid_sky(*_read(p)) if p else None) for uid, p in fl.items()}
    for b in BANDS:
        # Synthetic validation on peers: shift the epoch-2 cutout, re-measure.
        for uid, p in allfiles[(b, "e2")].items():
            if uid not in peer_uids or not p or cents[(b, "e2")][uid] is None:
                continue
            data, wcs = _read(p)
            base = cents[(b, "e2")][uid]
            got = []
            for d in ((INJECT_MAS, 0.0), (0.0, INJECT_MAS)):
                c = centroid_sky(shifted(data, wcs, base, d), wcs)
                got.extend(offsets_mas(base, c)[0].tolist() if c is not None else [np.nan] * 2)
            inject[b].append(got)
    out: dict[str, Any] = {
        "plan": "results/smacs_jwst3/epochs_1816_384_img/PLAN.md",
        "assumptions": {
            "size_arcsec": SIZE_ARCSEC,
            "box_px": BOX,
            "peer_dmag": PEER_DMAG,
            "peer_size_frac": PEER_SIZE_FRAC,
            "peer_radius_arcsec": PEER_RADIUS_ARCSEC,
            "peer_isolation_arcsec": PEER_ISOLATION_ARCSEC,
            "min_peers": MIN_PEERS,
            "percentile": PCTL,
            "angle_deg": ANGLE_DEG,
            "inject_mas": INJECT_MAS,
            "inject_tol_mas": INJECT_TOL_MAS,
            "peak_radius_arcsec": PEAK_RADIUS_ARCSEC,
            "match_arcsec": MATCH_ARCSEC,
            "baseline_yr": BASELINE_YR,
        },
        "validation": {},
        "targets": {},
    }
    valid = True
    for b in BANDS:
        g = np.array(inject[b], dtype=float)
        ok = np.all(np.isfinite(g), axis=1) if len(g) else np.zeros(0, bool)
        med = np.median(g[ok], axis=0) if ok.any() else np.full(4, np.nan)
        want = np.array([INJECT_MAS, 0.0, 0.0, INJECT_MAS])
        passed = bool(ok.sum() >= MIN_PEERS and np.all(np.abs(med - want) <= INJECT_TOL_MAS))
        valid &= passed
        out["validation"][b] = {
            "n": int(ok.sum()),
            "median_recovered_mas": med.tolist(),
            "injected_mas": want.tolist(),
            "max_abs_error_of_median_mas": float(np.max(np.abs(med - want))),
            "passed": passed,
            "label": "simulation",
        }
    for name, members in groups.items():
        rec: dict[str, Any] = {"bands": {}}
        per_band = {}
        for b in BANDS:

            def off(uid: str, b: str = b) -> np.ndarray | None:
                a, z = cents[(b, "e1")].get(uid), cents[(b, "e2")].get(uid)
                return None if a is None or z is None else offsets_mas(a, z)[0]

            t_off = off(members[0])
            peers = [o for o in (off(m) for m in members[1:]) if o is not None]
            st = band_stats(t_off, np.array(peers).reshape(-1, 2))
            per_band[b] = st
            rec["bands"][b] = st or {
                "not_measured": f"target centroid {'ok' if t_off is not None else 'failed'}; "
                f"{len(peers)} peers with centroids ({MIN_PEERS} needed)",
                "label": "derived",
            }
            if st is not None:
                st["peer_p95_motion_mas_per_yr"] = st["peer_residual_p95_mas"] / BASELINE_YR
        if all(per_band[b] is not None for b in BANDS):
            vecs = [np.array(per_band[b]["target_residual_mas"]) for b in BANDS]
            rec["post_hoc"] = {
                "angle_between_band_residuals_deg": angle_deg(*vecs),
                "label": "derived (post hoc; not part of the frozen rule)",
            }
        rec["n_peers_selected"] = len(members) - 1
        rec["epoch1_f150w_position"] = [pos[members[0]].ra.deg, pos[members[0]].dec.deg]
        rec["epoch1_f150w_row_sep_arcsec"] = target_sep[name]
        rec["decision"] = decide(per_band) if valid else "not evaluated: validation failed"
        out["targets"][name] = rec
    args.out.mkdir(parents=True, exist_ok=True)
    sheet(args.out / "contact_sheet.png", allfiles, cents, list(TARGETS))
    (args.out / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v["decision"] for k, v in out["targets"].items()}))
    print(json.dumps(out["validation"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Per-exposure (per-dither) photometry of a target against nearby compact knots on level-2 ``_cal``
frames (MACS0416 h044 vetting).

An epoch-to-epoch excess seen in resampled ``_i2d`` differences can come from the mosaic
(resampling, outlier rejection, a bad dither) rather than the sky. This script measures the
target and comparable knots in every ``_cal`` exposure that covers them, so the excess can be
tested exposure by exposure:

1. per exposure and source: recentred aperture flux at several radii (local annulus background);
2. each source's flux is divided by its median over all exposures of that filter (``q``);
3. per exposure, the knots' median ``q`` is the exposure's common factor (throughput, PSF,
   distortion near the target); the target's ``t = q_target / median_k q_k``;
4. per epoch, the mean ``t`` and its standard error; the epoch excess of epoch B over epoch A
   is ``mean t_B / mean t_A − 1``. Every knot is run through the same steps (leave-one-out
   normalisation), so the target's excess is ranked against the knots in robust σ.

A dither-dependent deficit (a few exposures low, the rest normal) points to processing; an excess
that is the same in every exposure and grows toward the aperture core points to the PSF core
(focus, jitter, wavefront) rather than to added flux; an excess constant with radius and well
outside the knots is what added point-source flux would give. Every number is ``derived``;
thresholds are ASSUMPTIONs.

    python scripts/cal_dither_phot.py --target h044 64.022499 -24.081711 \\
        --knots-json results/transients/macs0416_arcdiff/vet_sw/knots_o211_o212.json \\
        --cal-dir s3://stpubdata/jwst/public/jw01176/jw01176211001 \\
        --cal-dir s3://stpubdata/jwst/public/jw01176/jw01176212001 --out outputs/dither
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arc_diff_knots import finite_or_none, robust_z  # noqa: E402
from persistence_check import _open, _safe, header_info, list_cal  # noqa: E402
from transient_forced import aperture_flux, recentre  # noqa: E402

RADII_ARCSEC = (0.08, 0.16, 0.24)  # ASSUMPTION: 0.08″ is the arc_diff_vet SW aperture
SW_SCALE = 0.031  # arcsec per NIRCam SW pixel
LW_SCALE = 0.063
MIN_NJY, MAX_NJY = 150.0, 1000.0  # ASSUMPTION: knot epoch-1 flux range (arc_diff_knots + bright)
EXCLUDE_ARCSEC = 0.2  # ASSUMPTION: knots this close to a target are not references
RECENTRE_PX = 2  # centroid box half-size; recentring absorbs per-exposure WCS offsets
MAX_SHIFT_PX = 2.0  # ASSUMPTION: a larger centroid shift followed a neighbour; keep the WCS pixel
OFFSET_BOX_PX = 3  # centroid box half-size for the per-exposure knot offset
MAX_OFFSET_PX = 6  # ASSUMPTION: level-2 WCS errors stay below ~0.2″ (SW) per exposure


def epoch_of(name: str) -> str:
    """``01176o211`` from ``jw01176211001_..._cal.fits`` (program, then observation number)."""
    m = re.match(r"jw(\d{5})(\d{3})\d{3}_", name.rsplit("/", 1)[-1])
    return f"{m.group(1)}o{m.group(2)}" if m else "?"


def exposure_offset(cut: dict, refs: list[str]) -> tuple[float, float, int]:
    """Median (dx, dy) of the reference knots' centroids from their WCS pixels, and how many
    knots gave one. Level-2 ``_cal`` WCS lacks the level-3 tweakreg alignment, so every source of
    an exposure is off by about the same 1–3 px; (0, 0, n) when fewer than 3 knots centroid."""
    shifts = []
    for uid in refs:
        if uid not in cut:
            continue
        sci, lx, ly = cut[uid][0], cut[uid][2], cut[uid][3]
        cx, cy = recentre(sci, lx, ly, OFFSET_BOX_PX)
        if 0 < np.hypot(cx - lx, cy - ly) <= MAX_OFFSET_PX:
            shifts.append((cx - lx, cy - ly))
    if len(shifts) < 3:
        return 0.0, 0.0, len(shifts)
    dx, dy = np.median(np.array(shifts), axis=0)
    return float(dx), float(dy), len(shifts)


def measure(
    uri: str,
    pix: dict[str, tuple[float, float]],
    long: bool,
    radii: list[float],
    refs: list[str],
):
    """Recentred fluxes of every position on one exposure: {uid: {"x", "y", "f<i>", "e<i>"}}.

    Positions are first moved by the exposure's knot offset (``exposure_offset``), then each
    source is recentred within ``MAX_SHIFT_PX`` of that."""
    scale = LW_SCALE if long else SW_SCALE
    r_px = [r / scale for r in radii]
    half = int(np.ceil(3 * max(r_px))) + 6 + int(MAX_OFFSET_PX)
    cut = {}
    with _open(uri) as hdul, warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        ny, nx = hdul["SCI"].header["NAXIS2"], hdul["SCI"].header["NAXIS1"]
        wcs = WCS(hdul["SCI"].header, naxis=2)
        for uid, (x, y) in pix.items():
            xi, yi = int(round(x)), int(round(y))
            if not (half <= xi < nx - half and half <= yi < ny - half):
                continue
            sl = np.s_[yi - half : yi + half + 1, xi - half : xi + half + 1]
            sci = np.asarray(hdul["SCI"].section[sl], dtype=float)
            err = np.asarray(hdul["ERR"].section[sl], dtype=float)
            dq = np.asarray(hdul["DQ"].section[sl])
            sci[(dq & 1) > 0] = np.nan  # DO_NOT_USE
            cut[uid] = (sci, err, x - (xi - half), y - (yi - half), xi - half, yi - half)
    dx, dy, n_ref = exposure_offset(cut, refs)
    out = {}
    for uid, (sci, err, lx, ly, x0, y0) in cut.items():
        sx, sy = lx + dx, ly + dy
        cx, cy = recentre(sci, sx, sy, RECENTRE_PX)
        # recentre returns its input unchanged on every failure path (a converged centroid equals
        # it to float precision only by accident); a fallback position biases core apertures and
        # centroids, so such a source gets no measurement on this exposure
        if (cx, cy) == (sx, sy) or np.hypot(cx - sx, cy - sy) > MAX_SHIFT_PX:
            continue
        # sky centroid in the knots' frame: remove the exposure offset before the WCS
        ra, dec = (float(v) for v in wcs.all_pix2world(cx + x0 - dx, cy + y0 - dy, 0))
        row = {
            "x": cx + x0,
            "y": cy + y0,
            "ra": ra,
            "dec": dec,
            "offset_dx": dx,
            "offset_dy": dy,
            "n_offset_refs": n_ref,
            "shift_px": float(np.hypot(cx - sx, cy - sy)),
        }
        for i, r in enumerate(r_px):
            r_in = max(2 * r, r + 3)
            f, e = aperture_flux(sci, err, cx, cy, r, r_in, r_in + 4)
            row[f"f{i}"], row[f"e{i}"] = f, e
        out[uid] = row
    return out


def normalise(tab: Table, uids: list[str], knots: list[str], col: str) -> dict[str, np.ndarray]:
    """Per-exposure ``t`` for every uid in ``uids`` (knots: leave-one-out reference)."""
    files = sorted(set(tab["file"]))
    q = {u: np.full(len(files), np.nan) for u in uids}
    idx = {f: i for i, f in enumerate(files)}
    for r in tab:
        if r["uid"] in q:
            q[r["uid"]][idx[r["file"]]] = r[col]
    for u in uids:
        med = np.nanmedian(q[u]) if np.isfinite(q[u]).any() else np.nan
        q[u] = q[u] / med if med and np.isfinite(med) and med > 0 else np.full(len(files), np.nan)
    t = {}
    for u in uids:
        ref = np.array([q[k] for k in knots if k != u])
        with np.errstate(all="ignore"):
            n_ok = np.isfinite(ref).sum(axis=0) if len(ref) else np.zeros(len(files))
            m = np.where(n_ok >= 3, np.nanmedian(ref, axis=0) if len(ref) else np.nan, np.nan)
        t[u] = q[u] / m
    t["_files"] = np.array(files)
    return t


def epoch_excess(t: np.ndarray, epochs: np.ndarray, a: str, b: str) -> tuple[float, float]:
    """(mean t_b / mean t_a − 1, its standard error from the per-exposure scatter)."""
    ta, tb = t[(epochs == a) & np.isfinite(t)], t[(epochs == b) & np.isfinite(t)]
    if len(ta) < 2 or len(tb) < 2:
        return float("nan"), float("nan")
    ma, mb = ta.mean(), tb.mean()
    sa, sb = ta.std(ddof=1) / np.sqrt(len(ta)), tb.std(ddof=1) / np.sqrt(len(tb))
    return float(mb / ma - 1), float(mb / ma * np.hypot(sa / ma, sb / mb))


def centroid_shifts(tab: Table, targets: list[str], knots: list[str], epochs: list[str]) -> dict:
    """Epoch-to-epoch shift (mas, east and north) of each source's mean sky centroid, minus the
    knots' median shift, with the target ranked against the knots in robust σ per axis.

    A core whose light centre moves between epochs relative to its neighbours is what an
    i2d difference dipole needs (a registration or PSF-shape change, or added off-centre flux)."""
    pos: dict[str, dict[str, np.ndarray]] = {}
    for u in [*targets, *knots]:
        q = tab[tab["uid"] == u]
        if not len(q):
            continue
        cosd = np.cos(np.radians(float(np.median(q["dec"]))))
        for e in epochs:
            m = (q["epoch"] == e) & np.isfinite(q["ra"])
            if m.sum() >= 2:
                xy = np.c_[q["ra"][m] * cosd * 3.6e6, q["dec"][m] * 3.6e6]
                pos.setdefault(u, {})[e] = xy.mean(axis=0)
    out = {}
    for i, a in enumerate(epochs):
        for b in epochs[i + 1 :]:
            d = {u: pos[u][b] - pos[u][a] for u in pos if a in pos[u] and b in pos[u]}
            kd = np.array([d[k] for k in knots if k in d])
            if len(kd) < 5:
                continue
            med = np.median(kd, axis=0)
            res = {}
            for u in targets:
                if u not in d:
                    continue
                rel = d[u] - med
                z = [robust_z(float(d[u][i]), kd[:, i])[2] for i in (0, 1)]
                res[u] = {
                    "east_mas": round(float(rel[0]), 2),
                    "north_mas": round(float(rel[1]), 2),
                    "z_east": round(z[0], 2),
                    "z_north": round(z[1], 2),
                }
            sig = [robust_z(0.0, kd[:, i])[1] for i in (0, 1)]
            out[f"{a}_{b}"] = {
                "knot_robust_sigma_mas": [round(x, 2) for x in sig],
                "n_knots": len(kd),
                "targets": res,
            }
    return out


def load_knots(path: Path, targets: dict[str, tuple[float, float]], bands: list[str]):
    """Knot positions (union over ``bands``) in the flux range, away from every target."""
    d = json.loads(path.read_text())
    pos: dict[tuple[float, float], None] = {}
    for band in bands:
        for k in d["bands"].get(band, {}).get("knots", []):
            if MIN_NJY <= k["e1_njy"] <= MAX_NJY:
                pos[(round(k["ra"], 6), round(k["dec"], 6))] = None
    out = {}
    for ra, dec in pos:
        cosd = np.cos(np.radians(dec))
        near = [
            np.hypot((ra - tr) * cosd, dec - td) * 3600 < EXCLUDE_ARCSEC
            for tr, td in targets.values()
        ]
        dup = [np.hypot((ra - r) * cosd, dec - d_) * 3600 < 0.1 for r, d_ in out.values()]
        if not any(near) and not any(dup):
            out[f"k{len(out):02d}"] = (ra, dec)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--target", nargs=3, action="append", required=True, metavar=("UID", "RA", "DEC")
    )
    ap.add_argument("--knots-json", type=Path, required=True)
    ap.add_argument("--cal-dir", action="append", required=True)
    ap.add_argument("--filters", nargs="+", default=["F150W", "F200W"])
    ap.add_argument("--radii-arcsec", nargs="+", type=float, default=list(RADII_ARCSEC))
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    targets = {u: (float(r), float(d)) for u, r, d in args.target}
    knots = load_knots(args.knots_json, targets, args.filters)
    if clash := set(targets) & set(knots):
        ap.error(f"target names clash with knot ids: {sorted(clash)}")
    positions = {**targets, **knots}
    files = [f for d in args.cal_dir for f in list_cal(d)]
    failed: list[str] = []
    with ThreadPoolExecutor(args.workers) as ex:
        heads = list(ex.map(lambda f: _safe(header_info, f, positions), files))
    failed += [e for _, e in heads if e]
    info = [h for h, _ in heads if h and h["filter"] in args.filters and h["pix"]]

    def one(h):
        res, err = _safe(
            measure,
            h["file"],
            h["pix"],
            "LONG" in h["detector"].upper(),
            args.radii_arcsec,
            list(knots),
        )
        if err:
            failed.append(err)
        return h, res or {}

    with ThreadPoolExecutor(args.workers) as ex:
        measured = list(ex.map(one, info))
    rows = []
    for h, res in measured:
        name = h["file"].rsplit("/", 1)[-1]
        for uid, m in res.items():
            rows.append(
                {
                    "uid": uid,
                    "file": name,
                    "epoch": epoch_of(name),
                    "detector": h["detector"],
                    "filter": h["filter"],
                    "mjd": h["mjd"],
                    **m,
                }
            )
    if not rows:
        print(json.dumps({"error": "no measurements", "failed": failed}, indent=1))
        return 1
    tab = Table(rows=rows)
    args.out.mkdir(parents=True, exist_ok=True)
    tab.meta.update({"provenance": "derived", "radii_arcsec": args.radii_arcsec, "failed": failed})
    tab.write(args.out / "dither_phot.ecsv", overwrite=True)

    # chronological, so every pair "a_b" reads earlier → later
    start = {e: float(np.nanmedian(tab["mjd"][tab["epoch"] == e])) for e in set(tab["epoch"])}
    if not all(np.isfinite(v) for v in start.values()):
        print(json.dumps({"error": "an epoch has no finite EXPSTART", "epochs": start}))
        return 1
    epochs_all = sorted(start, key=start.get)
    summary = {
        "provenance": "derived",
        "source": "level-2 _cal exposures: " + ", ".join(args.cal_dir),
        "thresholds_assumption": {
            "radii_arcsec": args.radii_arcsec,
            "knot_njy": [MIN_NJY, MAX_NJY],
            "exclude_arcsec": EXCLUDE_ARCSEC,
            "recentre_px": RECENTRE_PX,
            "offset_box_px": OFFSET_BOX_PX,
            "max_offset_px": MAX_OFFSET_PX,
            "max_shift_px": MAX_SHIFT_PX,
        },
        "targets": {u: list(p) for u, p in targets.items()},
        "knots": {u: list(p) for u, p in knots.items()},
        "failed": failed,
        "filters": {},
    }
    knot_ids = list(knots)
    uids = list(targets) + knot_ids
    for band in args.filters:
        sub = tab[tab["filter"] == band]
        if not len(sub):
            continue
        bres: dict = {
            "n_exposures": {e: len(set(sub["file"][sub["epoch"] == e])) for e in epochs_all}
        }
        fe = dict(zip(sub["file"], sub["epoch"], strict=True))
        for i, r in enumerate(args.radii_arcsec):
            t = normalise(sub, uids, knot_ids, f"f{i}")
            ep = np.array([fe[f] for f in t["_files"]])
            rres: dict = {}
            for j, a in enumerate(epochs_all):
                for b in epochs_all[j + 1 :]:
                    ex = {u: epoch_excess(t[u], ep, a, b) for u in uids}
                    kref = np.array([ex[k][0] for k in knot_ids])
                    kmed, ksig, _ = robust_z(0, kref)
                    rres[f"{a}_{b}"] = {
                        "knot_median": kmed,
                        "knot_robust_sigma": ksig,
                        "n_knots": int(np.isfinite(kref).sum()),
                        "targets": {
                            u: {
                                "excess": ex[u][0],
                                "err": ex[u][1],
                                "z_vs_knots": robust_z(ex[u][0], kref)[2],
                                # exposures with a finite t: recentring can fail, an aperture can
                                # touch a DO_NOT_USE pixel, or < 3 knots are measured there
                                "n_measured": [
                                    int(np.isfinite(t[u][ep == a]).sum()),
                                    int(np.isfinite(t[u][ep == b]).sum()),
                                ],
                            }
                            for u in targets
                        },
                    }
            # per-exposure t of the targets, to look for a dither-dependent pattern
            rres["per_exposure"] = {
                u: [
                    {
                        "file": str(f),
                        "epoch": str(e),
                        "t": None if not np.isfinite(v) else round(float(v), 4),
                    }
                    for f, e, v in zip(t["_files"], ep, t[u], strict=True)
                ]
                for u in targets
            }
            bres[f"r{r:.2f}"] = rres
        bres["centroid_mas"] = centroid_shifts(sub, list(targets), knot_ids, epochs_all)
        summary["filters"][band] = bres
    (args.out / "summary.json").write_text(
        json.dumps(finite_or_none(summary), indent=1, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            {
                b: {k: v for k, v in d.items() if k != "per_exposure"}
                for b, d in summary["filters"].items()
            },
            indent=1,
            default=str,
        )[:6000]
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

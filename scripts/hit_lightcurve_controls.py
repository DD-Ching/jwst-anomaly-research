"""Calibrate ``hit_lightcurve.py`` epoch-pair changes of hits against nearby control sources.

Program-to-program PSF, background and registration differences change aperture fluxes even of
constant sources (the compact-core class of D-027 / the Abell 2744 grid). Two steps:

``pick``: control positions = local maxima of one mosaic box whose forced flux (same aperture
as ``hit_lightcurve.py``) lies in ``--flux-range``; prints ``--pos cNN RA DEC`` arguments.

``pick-arc``: arc-pixel controls for a hit on an arc. Compact peaks miss the on-arc systematics
(PSF and registration differences act on steep extended structure). Grid positions (step
``--step-arcsec``) whose aperture flux **and** annulus median (the extended light around the
aperture) both lie within ``--match`` times the hit's own values in the same mosaic, at least
``--min-sep-arcsec`` from the hit and from each other; prints ``--pos cNN RA DEC`` arguments.

``compare``: for each epoch pair (two ``obs_id`` of one band), z = (f2 − f1) / hypot(e1, e2) per
position; the controls' median and robust std (1.4826 MAD) calibrate each hit's z, and the
control tail fraction is the share of controls at least as far from the median. Controls are
not de-duplicated against hits (a control on a hit's own knot makes the test conservative).
Every number is ``derived``; thresholds are ASSUMPTIONs.

    python scripts/hit_lightcurve_controls.py pick --uri <i2d> --ra 3.5522 --dec -30.3753
    python scripts/hit_lightcurve_controls.py pick-arc --uri <i2d> --ra 3.552148 --dec -30.376850 \\
        --half-arcsec 16 --match 1.5 --step-arcsec 0.2 --min-sep-arcsec 0.6 --n 80
    python scripts/hit_lightcurve_controls.py compare --lc vet/ctrl301/lightcurves.ecsv \\
        --pair <obs_id_1> <obs_id_2> --out vet/ctrl301/control_z.ecsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hit_lightcurve import ANNULUS, aperture_photometry  # noqa: E402
from transient_forced import robust_std  # noqa: E402

MIN_CONTROLS = 20  # ASSUMPTION: controls needed for a stable median and robust std
CONTROL_NAME = re.compile(r"c\d+")  # the names ``pick`` prints (c00, c01, ...)


def is_control(name: str) -> bool:
    return CONTROL_NAME.fullmatch(name) is not None


def pick(
    uri: str, ra: float, dec: float, half: float, flux_range, n: int, radius: float, seed: int
):
    from arc_diff import read_box
    from scipy import ndimage

    b = read_box(uri, ra, dec, half)
    s, sc = b["sci"], b["scale"]
    k = b["pixar_sr"] * 1e15
    ys, xs = np.nonzero(ndimage.maximum_filter(s, size=9) == s)
    out = []
    for y, x in zip(ys, xs, strict=True):
        f = aperture_photometry(s, x, y, radius / sc, (ANNULUS[0] / sc, ANNULUS[1] / sc)) * k
        if flux_range[0] < f < flux_range[1]:
            r, d = b["wcs"].all_pix2world(x, y, 0)
            out.append((float(r), float(d)))
    idx = np.random.default_rng(seed).permutation(len(out))[:n]
    return [out[j] for j in idx]


def annulus_median(sci: np.ndarray, x: float, y: float, ann_px: tuple[float, float]) -> float:
    """Median of the background annulus that ``aperture_photometry`` subtracts; NaN if off."""
    ny, nx = sci.shape
    r_out = ann_px[1]
    if not (r_out <= x < nx - r_out and r_out <= y < ny - r_out):
        return float("nan")
    # the same cutout and annulus test as hit_lightcurve.aperture_photometry
    x0, y0 = max(int(np.floor(x - r_out)), 0), max(int(np.floor(y - r_out)), 0)
    x1, y1 = min(int(np.ceil(x + r_out)) + 1, nx), min(int(np.ceil(y + r_out)) + 1, ny)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    d2 = (xx - x) ** 2 + (yy - y) ** 2
    v = sci[y0:y1, x0:x1][(d2 >= ann_px[0] ** 2) & (d2 <= r_out**2)]
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else float("nan")


def pick_arc(
    uri: str,
    ra: float,
    dec: float,
    half: float,
    match: float,
    step: float,
    min_sep: float,
    n: int,
    radius: float,
    seed: int,
):
    """Grid positions matching the hit's aperture flux and annulus level within ``match``."""
    from arc_diff import read_box

    if not match > 1:
        raise SystemExit(f"--match must be > 1 (got {match})")
    b = read_box(uri, ra, dec, half)
    s, sc = b["sci"], b["scale"]
    ann = (ANNULUS[0] / sc, ANNULUS[1] / sc)
    x0, y0 = (float(v) for v in b["wcs"].all_world2pix(ra, dec, 0))
    f0 = aperture_photometry(s, x0, y0, radius / sc, ann)
    a0 = annulus_median(s, x0, y0, ann)
    if not (f0 > 0 and a0 > 0):
        raise SystemExit(f"hit has flux {f0} and annulus level {a0}: no arc to match")
    stp = step / sc
    cand = []
    for y in np.arange(ann[1], s.shape[0] - ann[1], stp):
        for x in np.arange(ann[1], s.shape[1] - ann[1], stp):
            if np.hypot(x - x0, y - y0) * sc < min_sep:
                continue
            f = aperture_photometry(s, x, y, radius / sc, ann)
            a = annulus_median(s, x, y, ann)
            if f0 / match <= f <= f0 * match and a0 / match <= a <= a0 * match:
                cand.append((x, y))
    out: list[tuple[float, float]] = []
    for j in np.random.default_rng(seed).permutation(len(cand)):
        x, y = cand[j]
        if all(np.hypot(x - u, y - v) * sc >= min_sep for u, v in out):
            out.append((x, y))
        if len(out) == n:
            break
    pts = []
    for x, y in out:
        r, d = b["wcs"].all_pix2world(x, y, 0)
        pts.append((float(r), float(d)))
    if len(pts) < max(n, MIN_CONTROLS):
        print(
            f"warning: {len(pts)} arc controls (asked {n}; compare needs {MIN_CONTROLS})",
            file=sys.stderr,
        )
    return pts


def compare(lc: Table, e1: str, e2: str) -> Table:
    d = {(r["name"], r["obs_id"]): (r["flux_njy"], r["err_njy"]) for r in lc}
    names = sorted({str(n) for n in lc["name"]})
    hits = [n for n in names if not is_control(n)]
    z = {}
    for n in names:
        if (n, e1) in d and (n, e2) in d:
            (f1, s1), (f2, s2) = d[(n, e1)], d[(n, e2)]
            if np.isfinite(f1) and np.isfinite(f2) and s1 > 0 and s2 > 0:
                z[n] = ((f2 - f1) / np.hypot(s1, s2), f2 / f1 if f1 else np.nan)
    ctrl = np.array([v[0] for k, v in z.items() if is_control(k)])
    if len(ctrl) < MIN_CONTROLS:
        raise SystemExit(f"{e1} -> {e2}: only {len(ctrl)} controls")
    med = float(np.median(ctrl))
    rsd = robust_std(ctrl, min_n=MIN_CONTROLS)
    if not rsd > 0:
        raise SystemExit(f"{e1} -> {e2}: control robust std is {rsd}")
    rows = []
    for h in hits:
        if h in z:
            zh, ratio = z[h]
            rows.append(
                (
                    h,
                    e1,
                    e2,
                    zh,
                    ratio,
                    (zh - med) / rsd,
                    float(np.mean(np.abs(ctrl - med) >= abs(zh - med))),
                    len(ctrl),
                    med,
                    rsd,
                )
            )
    # control statistics are per pair, so they are columns (table meta would be lost in vstack)
    t = Table(
        rows=rows,
        names=(
            "name",
            "obs_id_1",
            "obs_id_2",
            "z",
            "ratio",
            "z_cal",
            "ctrl_tail_frac",
            "n_controls",
            "ctrl_median",
            "ctrl_rsd",
        ),
        dtype=(str, str, str, float, float, float, float, int, float, float),
    )
    t.meta["provenance"] = "derived"
    return t


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pick")
    p.add_argument("--uri", required=True)
    p.add_argument("--ra", type=float, required=True)
    p.add_argument("--dec", type=float, required=True)
    p.add_argument("--half-arcsec", type=float, default=16.0)
    p.add_argument("--flux-range", type=float, nargs=2, default=(20.0, 150.0))
    p.add_argument("--n", type=int, default=80)
    p.add_argument("--radius-arcsec", type=float, default=0.15)
    p.add_argument("--seed", type=int, default=1)
    q = sub.add_parser("pick-arc")
    q.add_argument("--uri", required=True)
    q.add_argument("--ra", type=float, required=True)
    q.add_argument("--dec", type=float, required=True)
    q.add_argument("--half-arcsec", type=float, default=16.0)
    q.add_argument("--match", type=float, default=1.5)
    q.add_argument("--step-arcsec", type=float, default=0.2)
    q.add_argument("--min-sep-arcsec", type=float, default=0.6)
    q.add_argument("--n", type=int, default=80)
    q.add_argument("--radius-arcsec", type=float, default=0.15)
    q.add_argument("--seed", type=int, default=1)
    c = sub.add_parser("compare")
    c.add_argument("--lc", type=Path, required=True)
    c.add_argument("--pair", nargs=2, action="append", required=True, metavar=("OBS1", "OBS2"))
    c.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "pick":
        pts = pick(a.uri, a.ra, a.dec, a.half_arcsec, a.flux_range, a.n, a.radius_arcsec, a.seed)
        print(" ".join(f"--pos c{i:02d} {r:.6f} {d:.6f}" for i, (r, d) in enumerate(pts)))
        return 0
    if a.cmd == "pick-arc":
        pts = pick_arc(
            a.uri,
            a.ra,
            a.dec,
            a.half_arcsec,
            a.match,
            a.step_arcsec,
            a.min_sep_arcsec,
            a.n,
            a.radius_arcsec,
            a.seed,
        )
        print(" ".join(f"--pos c{i:02d} {r:.6f} {d:.6f}" for i, (r, d) in enumerate(pts)))
        return 0
    from astropy.table import vstack

    lc = Table.read(a.lc)
    t = vstack([compare(lc, e1, e2) for e1, e2 in a.pair])
    for col in ("z", "ratio", "z_cal", "ctrl_tail_frac", "ctrl_median", "ctrl_rsd"):
        t[col].info.format = ".3f"
    t.write(a.out, overwrite=True)
    t.pprint(max_lines=-1, max_width=-1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

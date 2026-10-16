"""Registration offset between two level-3 NIRCam mosaics from their MAST source catalogues.

Each pipeline ``*_cat.ecsv`` gives ``sky_centroid`` from its own mosaic's WCS, so matching compact
sources between two catalogues measures how the two mosaics' WCSs disagree. Controls, the split
into a calibration set A and a held-out set B, the two models and the pass rules are fixed in
``results/transients/a2744_grid_alt/reg1176/PLAN.md`` (written before the first measurement).

- M0: one shift (ΔRA cosδ, ΔDec) per mosaic pair.
- M1: affine, each offset component linear in (x, y) on the tangent plane.

Both are fitted on A with an iterative 3 × 1.4826 MAD clip and scored on B. All numbers are
``derived``; the cuts are ASSUMPTIONs listed in ``PLAN.md``.

    python scripts/reg_offsets.py \\
        --pair 1176:jw01176/L3/t/o221/jw01176-o221_t010_nircam_clear-{b}_cat.ecsv \\
        --ref 2561:jw02561/L3/t/o001/jw02561-o001_t003_nircam_clear-{b}_cat.ecsv --bands f356w \\
        --out results/transients/a2744_grid_alt/reg1176
"""

from __future__ import annotations

import argparse
import json
import zlib
from pathlib import Path

import numpy as np
from astropy.table import Table

S3 = "s3://stpubdata/jwst/public/"
MAG_RANGE = (21.0, 25.5)  # ASSUMPTION (PLAN.md): centroidable, below saturation
ISOLATION_ARCSEC = 1.0  # ASSUMPTION: no other catalogue source this close
MATCH_ARCSEC = 0.6  # ASSUMPTION: nearest-neighbour match radius
EDGE_ARCSEC = 3.0  # ASSUMPTION: distance from the catalogue footprint (convex hull) edge
CLIP = 3.0
MAD = 1.4826


def read_cat(path: str) -> Table:
    """A MAST source catalogue (S3 or local)."""
    if "://" in path or path.startswith("jw"):
        import fsspec

        uri = path if "://" in path else S3 + path
        with fsspec.open(uri, "rb", anon=True) as f:
            return Table.read(f.read().decode(), format="ascii.ecsv")
    return Table.read(path, format="ascii.ecsv")


def tangent(
    ra: np.ndarray, dec: np.ndarray, ra0: float, dec0: float
) -> tuple[np.ndarray, np.ndarray]:
    """Small-field tangent-plane offsets (arcsec) from (ra0, dec0)."""
    x = (ra - ra0) * np.cos(np.radians(dec0)) * 3600.0
    y = (dec - dec0) * 3600.0
    return x, y


def _mag(t: Table) -> np.ndarray:
    for c in ("aper_total_abmag", "isophotal_abmag", "aper30_abmag"):
        if c in t.colnames:
            return np.asarray(t[c].filled(np.nan) if hasattr(t[c], "filled") else t[c], dtype=float)
    raise KeyError("no AB magnitude column")


def _isolated(x: np.ndarray, y: np.ndarray, r: float) -> np.ndarray:
    from scipy.spatial import cKDTree

    tree = cKDTree(np.c_[x, y])
    n = tree.query_ball_point(np.c_[x, y], r, return_length=True)
    return n == 1


def _hull_distance(px: np.ndarray, py: np.ndarray, hx: np.ndarray, hy: np.ndarray) -> np.ndarray:
    """Distance of points inside the convex hull of (hx, hy) to its boundary (negative outside)."""
    from scipy.spatial import ConvexHull

    pts = np.c_[hx, hy]
    hull = ConvexHull(pts)
    # hull.equations: n·p + d ≤ 0 inside; distance to each facet = -(n·p + d)
    d = -(np.c_[px, py] @ hull.equations[:, :2].T + hull.equations[:, 2])
    return d.min(axis=1)


def controls(cat: Table, ref: Table) -> Table:
    """Matched compact, isolated, unsaturated sources inside both footprints."""
    from scipy.spatial import cKDTree

    ra0 = float(np.median(ref["sky_centroid"].ra.deg))
    dec0 = float(np.median(ref["sky_centroid"].dec.deg))
    out = {}
    for name, t in (("cat", cat), ("ref", ref)):
        sc = t["sky_centroid"]
        x, y = tangent(sc.ra.deg, sc.dec.deg, ra0, dec0)
        ok = np.isfinite(x) & np.isfinite(y)
        out[name] = {"x": x, "y": y, "ok": ok, "iso": np.zeros(len(t), bool)}
        out[name]["iso"][ok] = _isolated(x[ok], y[ok], ISOLATION_ARCSEC)
        out[name]["mag"] = _mag(t)
        out[name]["ext"] = np.asarray(t["is_extended"], dtype=bool)
    c, r = out["cat"], out["ref"]
    tree = cKDTree(np.c_[c["x"][c["ok"]], c["y"][c["ok"]]])
    ci = np.flatnonzero(c["ok"])
    d, j = tree.query(np.c_[r["x"], r["y"]], distance_upper_bound=MATCH_ARCSEC)
    m = np.isfinite(d) & r["ok"]
    ri = np.flatnonzero(m)
    cj = ci[j[m]]
    keep = (
        ~r["ext"][ri]
        & ~c["ext"][cj]
        & r["iso"][ri]
        & c["iso"][cj]
        & (r["mag"][ri] >= MAG_RANGE[0])
        & (r["mag"][ri] <= MAG_RANGE[1])
        & (c["mag"][cj] >= MAG_RANGE[0])
        & (c["mag"][cj] <= MAG_RANGE[1])
    )
    ri, cj = ri[keep], cj[keep]
    edge = np.minimum(
        _hull_distance(r["x"][ri], r["y"][ri], r["x"][r["ok"]], r["y"][r["ok"]]),
        _hull_distance(r["x"][ri], r["y"][ri], c["x"][c["ok"]], c["y"][c["ok"]]),
    )
    sel = edge >= EDGE_ARCSEC
    ri, cj = ri[sel], cj[sel]
    tab = Table()
    tab["ref_id"] = np.asarray(ref["label"])[ri] if "label" in ref.colnames else ri
    tab["ra"] = ref["sky_centroid"].ra.deg[ri]
    tab["dec"] = ref["sky_centroid"].dec.deg[ri]
    tab["x"] = r["x"][ri]
    tab["y"] = r["y"][ri]
    tab["mag_ref"] = r["mag"][ri]
    tab["dx"] = c["x"][cj] - r["x"][ri]  # cat − ref, arcsec (RA cosδ)
    tab["dy"] = c["y"][cj] - r["y"][ri]
    tab["set"] = [
        "A" if zlib.crc32(f"{a:.5f},{b:.5f}".encode()) % 2 == 0 else "B"
        for a, b in zip(tab["ra"], tab["dec"], strict=True)
    ]
    tab.meta.update(ra0=ra0, dec0=dec0)
    return tab


def robust_rms(v: np.ndarray) -> float:
    v = np.asarray(v, float)
    return float(MAD * np.median(np.abs(v - np.median(v)))) if v.size else float("nan")


def fit(tab: Table, model: str) -> dict:
    """Clipped least-squares fit of (dx, dy) on set A; returns parameters and their errors."""
    a = tab[tab["set"] == "A"]
    x, y, dx, dy = (np.asarray(a[k], float) for k in ("x", "y", "dx", "dy"))
    design = np.ones((len(a), 1)) if model == "M0" else np.c_[np.ones(len(a)), x / 60.0, y / 60.0]
    use = np.ones(len(a), bool)
    for _ in range(10):
        px, *_ = np.linalg.lstsq(design[use], dx[use], rcond=None)
        py, *_ = np.linalg.lstsq(design[use], dy[use], rcond=None)
        rx, ry = dx - design @ px, dy - design @ py
        sx, sy = robust_rms(rx[use]), robust_rms(ry[use])
        new = (np.abs(rx - np.median(rx[use])) <= CLIP * sx) & (
            np.abs(ry - np.median(ry[use])) <= CLIP * sy
        )
        if np.array_equal(new, use):
            break
        use = new
    cov = np.linalg.inv(design[use].T @ design[use])
    return {
        "model": model,
        "px": px.tolist(),
        "py": py.tolist(),
        "ex": (np.sqrt(np.diag(cov)) * sx).tolist(),
        "ey": (np.sqrt(np.diag(cov)) * sy).tolist(),
        "n_used": int(use.sum()),
        "n_A": len(a),
        "rms_A": [sx, sy],
        "note": "M1 slopes are per arcmin of tangent-plane x (RA) and y (Dec)",
    }


def score(tab: Table, f: dict) -> dict:
    """Residuals of the held-out set B after subtracting the model (clipped at 3σ for the stats)."""
    b = tab[tab["set"] == "B"]
    x, y = np.asarray(b["x"], float), np.asarray(b["y"], float)
    design = (
        np.ones((len(b), 1)) if f["model"] == "M0" else np.c_[np.ones(len(b)), x / 60.0, y / 60.0]
    )
    rx = np.asarray(b["dx"], float) - design @ np.array(f["px"])
    ry = np.asarray(b["dy"], float) - design @ np.array(f["py"])
    return {
        "n_B": len(b),
        "median": [float(np.median(rx)), float(np.median(ry))],
        "rms": [robust_rms(rx), robust_rms(ry)],
        "raw_median": [float(np.median(b["dx"])), float(np.median(b["dy"]))],
    }


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pair", required=True, help="label:catalogue path with {b} for the band")
    ap.add_argument("--ref", required=True, help="label:reference catalogue path with {b}")
    ap.add_argument("--bands", nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    lab, path = args.pair.split(":", 1)
    rlab, rpath = args.ref.split(":", 1)
    args.out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for b in args.bands:
        tab = controls(read_cat(path.format(b=b)), read_cat(rpath.format(b=b)))
        tab.write(args.out / f"controls_{lab}_vs_{rlab}_{b}.ecsv", overwrite=True)
        res = {}
        for m in ("M0", "M1"):
            f = fit(tab, m)
            res[m] = {"fit": f, "B": score(tab, f)}
        summary[b] = {"n_controls": len(tab), **res}
        print(b, json.dumps(summary[b]["M0"]["B"]), "M1 B rms", summary[b]["M1"]["B"]["rms"])
        (args.out / f"summary_{lab}_vs_{rlab}.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

"""Abell 2744 system 700.1: fine source-redshift scan in Bergamini+23.

D-036 left open whether Bergamini+23 (``abell2744-bergamini23``, z of 700.1 free, best fit 1.217)
reproduces 700.1a/b/c as three separate images at any redshift; z = 1, 1.5, 2, 3, 5 did not.
Plan (frozen first): ``results/lens_consistency/a2744_s700_zscan/PLAN.md``.

At each z: one source (``a2744_uncover.wls_source`` of the three back-traced images),
``lensmodel.find_images``, then a one-to-one assignment of the catalogued images to predicted
images with the smallest total separation. A z *reproduces* the system when every catalogued
image gets its own predicted image within ``MATCH_ARCSEC`` (ASSUMPTION, the D-037 pair criterion).

Everything computed is ``model_prediction``; the image positions are ``observed``.

Usage: python scripts/a2744_s700.py [--out results/lens_consistency/a2744_s700_zscan]
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a2744_uncover as au  # noqa: E402
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

MODEL = "abell2744-bergamini23"
SYSTEM = "700.1"
Z_GRID = np.round(np.arange(1.0, 4.0001, 0.025), 3)  # source redshifts scanned
MATCH_ARCSEC = 1.0  # each catalogued image needs its own predicted image this close (ASSUMPTION)
GRID_STEP = 0.25
SHEET_BANDS = ("f150w", "f200w", "f444w")
DEFAULT_OUT = Path("results/lens_consistency/a2744_s700_zscan")


def _r(x: float, nd: int = 3) -> float | None:
    return None if not np.isfinite(x) else round(float(x), nd)


def assign(pred: Table, cat: Table) -> tuple[list[int], np.ndarray]:
    """One-to-one assignment of catalogued images to predicted ones (minimum total separation).

    Returns the predicted index per catalogued image and the separations (arcsec); an empty list
    and NaNs when fewer images are predicted than catalogued."""
    n = len(cat)
    if len(pred) < n:
        return [], np.full(n, np.nan)
    d = np.array(
        [au._sep_arcsec(pred["ra"], pred["dec"], float(c["ra"]), float(c["dec"])) for c in cat]
    )
    best = min(itertools.permutations(range(len(pred)), n), key=lambda p: d[range(n), p].sum())
    return list(best), d[range(n), list(best)]


def nearest(pred: Table, cat: Table) -> np.ndarray:
    """Separation (arcsec) of each catalogued image to its nearest predicted image (may share)."""
    if len(pred) == 0:
        return np.full(len(cat), np.nan)
    return np.array(
        [
            float(np.min(au._sep_arcsec(pred["ra"], pred["dec"], float(c["ra"]), float(c["dec"]))))
            for c in cat
        ]
    )


def verdict(rows: list[dict]) -> dict:
    """Frozen outcome: A ("fitted-z artefact") if any z reproduces the system, else B."""
    good = [r["z"] for r in rows if r["reproduced"]]
    finite = [r for r in rows if r["rms_arcsec"] is not None]
    best = min(finite, key=lambda r: r["rms_arcsec"]) if finite else None
    return {
        "outcome": "A: fitted-z artefact" if good else "B: not reproduced at any z",
        "z_reproduced": [min(good), max(good)] if good else None,
        "n_z_reproduced": len(good),
        "best_assignment": best,
    }


def solve(model, cat: Table, grid, z: float) -> tuple[Table, np.ndarray]:
    """Predicted images of the ``wls_source`` of ``cat`` at ``z``, and the back-traced sources."""
    ox, oy = model.to_frame(np.asarray(cat["ra"], float), np.asarray(cat["dec"], float))
    m = model.lens_map(ox, oy, np.full(len(cat), z))
    jac = np.stack([np.stack([m["a11"], m["a12"]], -1), np.stack([m["a21"], m["a22"]], -1)], -2)
    beta = np.stack([np.asarray(m["beta_x"], float), np.asarray(m["beta_y"], float)], -1)
    bx, by = au.wls_source(beta, jac)
    return lensmodel.find_images(model, grid, float(bx), float(by), z), beta


def validate(model, arcs: Table, grid) -> dict:
    """Known cases: the same assignment rule on every other system with a fixed (catalogued) z
    and at least two images. Returns the reproduced fraction and the failing systems."""
    systems = np.asarray(arcs["system"]).astype(str)
    zcol = np.asarray(arcs["z"], float)
    done, failed = 0, []
    for s in sorted(set(systems) - {SYSTEM}):
        cat = arcs[systems == s]
        z = float(zcol[systems == s][0])
        if len(cat) < 2 or not z > 0:
            continue
        pred, _ = solve(model, cat, grid, z)
        idx, sep = assign(pred, cat)
        done += 1
        if not (idx and np.all(sep <= MATCH_ARCSEC)):
            failed.append(s)
    return {
        "n_systems": done,
        "n_reproduced": done - len(failed),
        "fraction": _r((done - len(failed)) / done) if done else None,
        "failed": failed,
    }


def scan(model, cat: Table, grid) -> list[dict]:
    rows = []
    for z in Z_GRID:
        pred, beta = solve(model, cat, grid, float(z))
        idx, sep = assign(pred, cat)
        near = nearest(pred, cat)
        rows.append(
            {
                "z": float(z),
                "n_images": len(pred),
                "sep_arcsec": [_r(s) for s in sep],
                "rms_arcsec": _r(float(np.sqrt(np.mean(sep**2)))) if idx else None,
                "mu": [_r(float(pred["magnification"][k]), 1) for k in idx],
                "nearest_arcsec": [_r(s) for s in near],
                "source_spread_arcsec": _r(float(np.ptp(beta, axis=0).max())),
                "reproduced": bool(idx) and bool(np.all(sep <= MATCH_ARCSEC)),
            }
        )
    return rows


def _stretch(img: np.ndarray) -> np.ndarray:
    """arcsinh of the stamp in units of its robust noise (median-subtracted, MAD sigma)."""
    med = np.nanmedian(img)
    sig = 1.4826 * np.nanmedian(np.abs(img - med))
    return np.arcsinh((img - med) / (sig if sig > 0 else 1.0))


def sheet(cat: Table, rows: list[dict], model, grid, path: Path) -> None:
    """F150W / F200W / F444W stamps around 700.1a/b and 700.1c with the catalogued positions
    (circles) and the predicted images at the best-rms z (crosses)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    finite = [r for r in rows if r["rms_arcsec"] is not None]
    if not finite:
        print("no z predicts three images: no sheet", flush=True)
        return
    best = min(finite, key=lambda r: r["rms_arcsec"])
    pred, _ = solve(model, cat, grid, best["z"])
    centres = {
        "700.1a/b": (float(np.mean(cat["ra"][:2])), float(np.mean(cat["dec"][:2])), 4.0),
        "700.1c": (float(cat["ra"][2]), float(cat["dec"][2]), 4.0),
    }
    fig, axes = plt.subplots(2, len(SHEET_BANDS), figsize=(3.2 * len(SHEET_BANDS), 6.4))
    for j, band in enumerate(SHEET_BANDS):
        stamper, close = lc.image_stamper(au.S3_I2D.format(band), 4.0)
        try:
            for i, (name, (ra0, dec0, half)) in enumerate(centres.items()):
                img, _, xx, yy = stamper(ra0, dec0)
                ax = axes[i, j]
                # the mosaic is not north-up: plot on the stamper's sky offsets (West = +x)
                ax.pcolormesh(xx, yy, _stretch(img), cmap="gray_r", vmin=-1, vmax=8)
                ax.set_xlim(-half, half)
                ax.set_ylim(-half, half)
                ax.set_aspect("equal")
                cos = np.cos(np.deg2rad(dec0))
                for t, mk, col in ((cat, "o", "tab:blue"), (pred, "x", "tab:red")):
                    dx = -(np.asarray(t["ra"], float) - ra0) * cos * 3600.0
                    dy = (np.asarray(t["dec"], float) - dec0) * 3600.0
                    inside = (np.abs(dx) < half) & (np.abs(dy) < half)
                    style = {"facecolors": "none", "edgecolors": col} if mk == "o" else {"c": col}
                    ax.scatter(dx[inside], dy[inside], marker=mk, s=60, **style)
                ax.set_title(f"{name} {band.upper()}", fontsize=9)
                ax.set_xticks([])
                ax.set_yticks([])
        finally:
            close()
    fig.suptitle(f'Abell 2744 700.1: catalogued (o), Bergamini+23 at z {best["z"]} (x); 8" stamps')
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--sheet", action="store_true", help="also write sheet.png (S3 range reads)")
    args = ap.parse_args(argv)
    model, files, par = lc.load_model(MODEL)
    arcs, _ = lc.image_list(MODEL, files, par)
    cat = arcs[np.asarray(arcs["system"]).astype(str) == SYSTEM]
    cat.sort("image_id")  # a, b, c: the sheet and the summary rely on this order
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, arcs), GRID_STEP)
    known = validate(model, arcs, grid)
    print("known cases:", known, flush=True)
    rows = scan(model, cat, grid)
    out = {
        "model": MODEL,
        "system": SYSTEM,
        "images": [str(i) for i in cat["image_id"]],
        "match_arcsec": MATCH_ARCSEC,
        "grid_step_arcsec": GRID_STEP,
        "z_grid": [float(Z_GRID[0]), float(Z_GRID[-1]), 0.025],
        "known_cases": known,
        **verdict(rows),
        "provenance": "model_prediction (positions observed)",
        "rows": rows,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    if args.sheet:
        sheet(cat, rows, model, grid, args.out / "sheet.png")
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    return out


if __name__ == "__main__":
    main()

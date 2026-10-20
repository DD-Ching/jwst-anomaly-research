"""Abell 2744 system 700.1: fine source-redshift scan in Bergamini+23 or the UNCOVER v2.0 maps.

D-036 left open whether Bergamini+23 (``abell2744-bergamini23``, z of 700.1 free, best fit 1.217)
reproduces 700.1a/b/c as three separate images at any redshift; z = 1, 1.5, 2, 3, 5 did not.
Plans (frozen first): ``results/lens_consistency/a2744_s700_zscan/PLAN.md`` (Bergamini+23) and
``results/lens_consistency/a2744_s700_uncover/PLAN.md`` (``--model uncover2``: known cases are
UNCOVER's own list, at the radii ``UNCOVER_RADII``, gated by ``gated_match``).

At each z: one source (``a2744_uncover.wls_source`` of the three back-traced images),
``lensmodel.find_images``, then a one-to-one assignment of the catalogued images to predicted
images with the smallest total separation. A z *reproduces* the system when every catalogued
image gets its own predicted image within ``MATCH_ARCSEC`` (ASSUMPTION, the D-037 pair criterion).

Everything computed is ``model_prediction``; the image positions are ``observed``.

Usage: python scripts/a2744_s700.py [--model bergamini23|uncover2] [--out DIR] [--sheet]
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
# UNCOVER v2.0: the #412 radius and the D-041 match radius (1.8x its 0.84" own-list rms); the scan
# uses the smallest one that reproduces >= KNOWN_GATE of the own-list systems (ASSUMPTIONs)
UNCOVER_RADII = (1.0, 1.5)
KNOWN_GATE = 0.8
COUNTERPART_ARCSEC = (
    0.5  # an UNCOVER image this close to a catalogued one is the same image (ASSUMPTION)
)
GRID_STEP = 0.25
SHEET_BANDS = ("f150w", "f200w", "f444w")
DEFAULT_OUT = {
    "bergamini23": Path("results/lens_consistency/a2744_s700_zscan"),
    "uncover2": Path("results/lens_consistency/a2744_s700_uncover"),
}


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


def gated_match(known: dict[float, dict], gate: float = KNOWN_GATE) -> float | None:
    """Smallest match radius whose known-case fraction is >= ``gate`` (None: tool not validated)."""
    ok = [r for r, k in known.items() if k["fraction"] is not None and k["fraction"] >= gate]
    return min(ok) if ok else None


def nearest_index(pred: Table, cat: Table) -> tuple[list[int], np.ndarray]:
    """Index of, and separation (arcsec) to, each catalogued image's nearest predicted image."""
    if len(pred) == 0:
        return [-1] * len(cat), np.full(len(cat), np.nan)
    d = np.array(
        [au._sep_arcsec(pred["ra"], pred["dec"], float(c["ra"]), float(c["dec"])) for c in cat]
    )
    return [int(k) for k in d.argmin(axis=1)], d.min(axis=1)


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


def failure_kind(n_cat: int, nearest_idx: list[int], nearest_sep, match: float, n_pred: int) -> str:
    """Why a system is not reproduced: "missing" (fewer predicted images than catalogued),
    "shared" (every catalogued image has a predicted image within ``match``, but two share it),
    or "offset" (some catalogued image has no predicted image within ``match``)."""
    if n_pred < n_cat:
        return "missing"
    if np.all(np.asarray(nearest_sep, float) <= match) and len(set(nearest_idx)) < n_cat:
        return "shared"
    return "offset"


def validate(
    model, arcs: Table, grid, match: float = MATCH_ARCSEC, exclude: frozenset = frozenset()
) -> dict:
    """Known cases: the same assignment rule on every other system (not ``SYSTEM``, not in
    ``exclude``) with a fixed (catalogued) z and at least two images. Returns the reproduced
    fraction, the failing systems and, per failing system, the image counts, each catalogued
    image's nearest predicted image (index, separation) and the ``failure_kind``."""
    systems = np.asarray(arcs["system"]).astype(str)
    zcol = np.asarray(arcs["z"], float)
    done, failed, detail = 0, [], []
    for s in sorted(set(systems) - {SYSTEM} - set(exclude)):
        cat = arcs[systems == s]
        z = float(zcol[systems == s][0])
        if len(cat) < 2 or not z > 0:
            continue
        pred, _ = solve(model, cat, grid, z)
        idx, sep = assign(pred, cat)
        done += 1
        if not (idx and np.all(sep <= match)):
            failed.append(s)
            near_idx, near = nearest_index(pred, cat)
            detail.append(
                {
                    "system": s,
                    "images": [str(i) for i in cat["image_id"]],
                    "n_catalogued": len(cat),
                    "n_predicted": len(pred),
                    "max_assigned_arcsec": _r(float(np.max(sep))) if idx else None,
                    "nearest_pred": near_idx,
                    "nearest_arcsec": [_r(float(x)) for x in near],
                    "kind": failure_kind(len(cat), near_idx, near, match, len(pred)),
                }
            )
    return {
        "n_systems": done,
        "n_reproduced": done - len(failed),
        "fraction": _r((done - len(failed)) / done) if done else None,
        "failed": failed,
        "failed_detail": detail,
        "n_failed": {
            k: sum(d["kind"] == k for d in detail) for k in ("missing", "shared", "offset")
        },
    }


def scan(model, cat: Table, grid, match: float = MATCH_ARCSEC) -> list[dict]:
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
                "reproduced": bool(idx) and bool(np.all(sep <= match)),
            }
        )
    return rows


def _stretch(img: np.ndarray) -> np.ndarray:
    """arcsinh of the stamp in units of its robust noise (median-subtracted, MAD sigma)."""
    med = np.nanmedian(img)
    sig = 1.4826 * np.nanmedian(np.abs(img - med))
    return np.arcsinh((img - med) / (sig if sig > 0 else 1.0))


def sheet(
    cat: Table, rows: list[dict], model, grid, path: Path, label: str = "Bergamini+23"
) -> None:
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
    fig.suptitle(f'Abell 2744 700.1: catalogued (o), {label} at z {best["z"]} (x); 8" stamps')
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def load_uncover2(maps: list[Path]):
    """UNCOVER v2.0 maps as a ``MapLensModel`` and the image-table path; each map path is
    appended to ``maps`` as soon as it is cached, so the caller can delete it on failure."""
    from jwst_anomaly.photometry import fetch_catalog

    files = {}
    for key, (url, sha) in au.UNCOVER_V2.items():
        files[key] = fetch_catalog(url, sha)
        if key != "images":
            maps.append(files[key])
    return au.uncover2_model(files), files["images"]


def counterparts(cat: Table, images: Path) -> list[dict]:
    """Nearest UNCOVER own-list image of each catalogued image, with its redshift metadata
    (``observed`` catalogue values; post hoc, see the UNCOVER RESULT.md)."""
    raw = Table.read(images)
    out = []
    for c in cat:
        d = au._sep_arcsec(raw["RA"], raw["Dec"], float(c["ra"]), float(c["dec"]))
        k = int(np.argmin(d))
        row = {"image": str(c["image_id"]), "uncover_id": None, "sep_arcsec": _r(float(d[k]))}
        if d[k] <= COUNTERPART_ARCSEC:
            row["uncover_id"] = repr(float(raw["ID"][k]))
            row["z"] = au._col(raw, "z", k, np.nan)
            for col in ("Redshift_type", "Redshift_reference", "Redshift_instrument"):
                row[col.lower()] = au._col(raw, col, k, "")
        out.append(row)
    return out


def pair_geometry(images: Table) -> list[dict]:
    """Per system: separation of its first two images and their midpoint, as offsets (arcsec,
    West = +x) from the first system's first image (``observed`` positions, derived numbers)."""
    systems = np.asarray(images["system"]).astype(str)
    ra0, dec0 = float(images["ra"][0]), float(images["dec"][0])
    cos = np.cos(np.deg2rad(dec0))
    out = []
    for s in dict.fromkeys(systems):
        t = images[systems == s]
        if len(t) < 2:
            continue
        ra, dec = np.asarray(t["ra"][:2], float), np.asarray(t["dec"][:2], float)
        out.append(
            {
                "system": s,
                "images": [str(i) for i in t["image_id"][:2]],
                "sep_arcsec": _r(float(au._sep_arcsec(ra[1:], dec[1:], ra[0], dec[0])[0])),
                "mid_dx_arcsec": _r(-(ra.mean() - ra0) * cos * 3600.0),
                "mid_dy_arcsec": _r((dec.mean() - dec0) * 3600.0),
            }
        )
    return out


def counterpart_sources(rows: list[dict]) -> frozenset:
    """UNCOVER systems of the target's own source (every knot "S.K" of each matched "S")."""
    return frozenset(
        au.uncover_system(r["uncover_id"]).split(".")[0] for r in rows if r["uncover_id"]
    )


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", choices=sorted(DEFAULT_OUT), default="bergamini23")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--sheet", action="store_true", help="also write sheet.png (S3 range reads)")
    ap.add_argument("--keep-maps", action="store_true", help="keep the cached UNCOVER maps")
    args = ap.parse_args(argv)
    out_dir = args.out or DEFAULT_OUT[args.model]
    berg, files, par = lc.load_model(MODEL)
    arcs, _ = lc.image_list(MODEL, files, par)
    cat = arcs[np.asarray(arcs["system"]).astype(str) == SYSTEM]
    cat.sort("image_id")  # a, b, c: the sheet and the summary rely on this order
    maps: list[Path] = []
    try:
        exclude: frozenset = frozenset()
        if args.model == "bergamini23":
            model, known_list, label = berg, arcs, "Bergamini+23"
            radii: tuple[float, ...] = (MATCH_ARCSEC,)
        else:
            model, images = load_uncover2(maps)
            known_list = au.uncover_images(images)
            label, radii = "UNCOVER v2.0", UNCOVER_RADII
            found = counterparts(cat, images)
            # the target's own source is not a known case (post hoc, RESULT.md)
            src = counterpart_sources(found)
            sys_src = np.char.partition(np.asarray(known_list["system"]).astype(str), ".")[:, 0]
            exclude = frozenset(
                np.asarray(known_list["system"]).astype(str)[np.isin(sys_src, list(src))]
            )
        both = Table({"ra": np.r_[arcs["ra"], known_list["ra"]]})
        both["dec"] = np.r_[arcs["dec"], known_list["dec"]]
        grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, both), GRID_STEP)
        known = {r: validate(model, known_list, grid, r, exclude) for r in radii}
        print("known cases:", known, flush=True)
        # the gate is the UNCOVER plan's; the Bergamini+23 plan has none
        match = MATCH_ARCSEC if args.model == "bergamini23" else gated_match(known)
        rows = scan(model, cat, grid, match) if match is not None else []
        out = {
            "model": MODEL if args.model == "bergamini23" else "abell2744-uncover2",
            "system": SYSTEM,
            "images": [str(i) for i in cat["image_id"]],
            "match_arcsec": match,
            "grid_step_arcsec": GRID_STEP,
            "z_grid": [float(Z_GRID[0]), float(Z_GRID[-1]), 0.025],
        }
        if args.model == "bergamini23":
            out["known_cases"] = known[MATCH_ARCSEC]
        else:
            out["known_cases"] = {str(r): k for r, k in known.items()}
            out["known_gate"] = KNOWN_GATE
            out["uncover_counterparts"] = found
            out["known_excluded"] = sorted(exclude)
            # UNCOVER's own fit of the target's source (a constraint of that model, not a test)
            own = known_list[np.isin(np.asarray(known_list["system"]).astype(str), list(exclude))]
            out["target_systems_fit"] = {str(r): validate(model, own, grid, r) for r in radii}
            out["target_pairs"] = pair_geometry(own)
        if match is None:
            out["outcome"] = "tool not validated"
        else:
            out.update(verdict(rows))
        out["provenance"] = "model_prediction (positions observed)"
        out["rows"] = rows
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
        if args.sheet and rows:
            sheet(cat, rows, model, grid, out_dir / "sheet.png", label)
    finally:
        if not args.keep_maps:  # cloud-disk rule: ~250 MB of maps, deleted after use
            for p in maps:
                p.unlink(missing_ok=True)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    return out


if __name__ == "__main__":
    main()

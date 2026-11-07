"""Third epoch for SMACS ``f200w_384``: does 12774 o002 reproduce the 2736 -> 6882 shift?

Plan: ``results/smacs_jwst3/epoch3_384/PLAN.md`` (frozen in the claim commit). Reuses the centroid,
peer, injection and residual functions of ``smacs_epochs_img.py`` unchanged; epoch 2 is replaced
by 12774 (``e3``). ``decide`` encodes the plan's outcome rules.

    python scripts/smacs_epoch3_384.py <scratch dir> --cat-dir <dir with the CATs> \\
        --out results/smacs_jwst3/epoch3_384
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import catalog, cutouts

sys.path.insert(0, str(Path(__file__).resolve().parent))
from epoch_compare import offsets_mas  # noqa: E402
from smacs_epochs_1816_384 import _angle_deg as angle_deg  # noqa: E402
from smacs_epochs_img import (  # noqa: E402
    ANGLE_DEG,
    BANDS,
    INJECT_MAS,
    INJECT_TOL_MAS,
    MATCH_ARCSEC,
    MIN_PEERS,
    SIZE_ARCSEC,
    _read,
    band_stats,
    centroid_sky,
    select_peers,
    shifted,
)

NAME = "f200w_384"
CATS = {
    "e1": "jw02736-o001_t001_nircam_clear-{b}_cat.ecsv",
    "e2": "jw06882-o057_t057_nircam_clear-{b}_cat.ecsv",
    "e3": "jw12774-o002_t032_nircam_clear-{b}_cat.ecsv",
}
S3 = {"e1": "jw02736/L3/t/o001/", "e2": "jw06882/L3/t/o057/", "e3": "jw12774/L3/t/o002/"}
PREV = Path("results/smacs_jwst3/epochs_1816_384_img/summary.json")


def decide(per_band: dict[str, dict[str, Any] | None], prev_vec: dict[str, list[float]]) -> str:
    """Plan outcome rules (frozen): ``per_band[b]`` is ``band_stats`` (2736 -> 12774) or None;
    ``prev_vec[b]`` is the same band's 2736 -> 6882 residual vector (mas)."""
    res = [per_band.get(b) for b in BANDS]
    if any(r is None for r in res):
        return "undetermined"
    if all(
        r["target_residual_total_mas"] > r["peer_residual_p95_mas"]
        and angle_deg(np.array(r["target_residual_mas"]), np.array(prev_vec[b])) <= ANGLE_DEG
        for b, r in zip(BANDS, res, strict=True)
    ):
        return "shift reproduced"
    if all(r["target_residual_total_mas"] <= r["peer_residual_median_mas"] for r in res):
        return "shift not reproduced"
    return "undetermined"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("work", type=Path, help="scratch directory for cutout FITS files")
    ap.add_argument("--cat-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    prev = json.loads(PREV.read_text())["targets"][NAME]
    prev_vec = {b: prev["bands"][b]["target_residual_mas"] for b in BANDS}
    cat1 = catalog.load_pipeline_catalog(args.cat_dir / CATS["e1"].format(b="f150w"))
    cat3 = catalog.load_pipeline_catalog(args.cat_dir / CATS["e3"].format(b="f150w"))
    c1 = SkyCoord(cat1["ra"], cat1["dec"], unit="deg")
    t = SkyCoord(*prev["epoch1_f150w_position"], unit="deg")
    sep = t.separation(c1).arcsec
    k = int(np.argmin(sep))
    if sep[k] > MATCH_ARCSEC:
        raise SystemExit(f"{NAME}: no epoch-1 F150W row within {MATCH_ARCSEC} arcsec")
    peers = select_peers(cat1, c1[k], k, cat3)
    uids = [NAME] + [f"{NAME}_p{i}" for i in peers]
    rows = [k, *peers]
    targets = Table(
        {
            "source_uid": uids,
            "ra": [float(cat1["ra"][i]) for i in rows],
            "dec": [float(cat1["dec"][i]) for i in rows],
        }
    )

    def fetch(key: tuple[str, str]) -> dict[str, str]:
        b, ep = key
        i2d = CATS[ep].format(b=b).replace("_cat.ecsv", "_i2d.fits")
        tab = cutouts.make_cutouts(
            f"s3://stpubdata/jwst/public/{S3[ep]}{i2d}",
            targets,
            size_arcsec=SIZE_ARCSEC,
            out_dir=args.work / ep,
            storage_options={"anon": True},
        )
        return dict(zip(tab["source_uid"], tab["path"], strict=True))

    keys = [(b, ep) for b in BANDS for ep in ("e1", "e2", "e3")]
    with ThreadPoolExecutor(max_workers=len(keys)) as pool:  # six independent S3 reads
        files = dict(zip(keys, pool.map(fetch, keys), strict=True))
    cents = {
        key: {uid: (centroid_sky(*_read(p)) if p else None) for uid, p in fl.items()}
        for key, fl in files.items()
    }
    out: dict[str, Any] = {
        "plan": "results/smacs_jwst3/epoch3_384/PLAN.md",
        "epochs": {
            "e1": "2736 o001 (2022)",
            "e2": "6882 o057 (2026-06-05)",
            "e3": "12774 o002 (2026-07-06)",
        },
        "n_peers_selected": len(peers),
        "epoch1_f150w_row_sep_arcsec": float(sep[k]),
        "prev_2736_6882_residual_mas": prev_vec,
        "validation": {},
        "bands": {},
        "post_hoc_6882_12774": {},
    }
    valid = True
    per_band: dict[str, dict[str, Any] | None] = {}
    for b in BANDS:
        got = []
        for uid in uids[1:]:
            p, base = files[(b, "e3")].get(uid), cents[(b, "e3")].get(uid)
            if not p or base is None:
                continue
            data, wcs = _read(p)
            row = []
            for d in ((INJECT_MAS, 0.0), (0.0, INJECT_MAS)):
                c = centroid_sky(shifted(data, wcs, base, d), wcs)
                row.extend(offsets_mas(base, c)[0].tolist() if c is not None else [np.nan] * 2)
            got.append(row)
        g = np.array(got, dtype=float).reshape(-1, 4)
        ok = np.all(np.isfinite(g), axis=1)
        med = np.median(g[ok], axis=0) if ok.any() else np.full(4, np.nan)
        want = np.array([INJECT_MAS, 0.0, 0.0, INJECT_MAS])
        passed = bool(ok.sum() >= MIN_PEERS and np.all(np.abs(med - want) <= INJECT_TOL_MAS))
        valid &= passed
        out["validation"][b] = {
            "n": int(ok.sum()),
            "median_recovered_mas": med.tolist(),
            "max_abs_error_of_median_mas": float(np.max(np.abs(med - want))),
            "passed": passed,
            "label": "simulation",
        }

        def off(uid: str, a: str, z: str, b: str = b) -> np.ndarray | None:
            p, q = cents[(b, a)].get(uid), cents[(b, z)].get(uid)
            return None if p is None or q is None else offsets_mas(p, q)[0]

        def stats(a: str, z: str, members: list[str]) -> dict[str, Any] | None:
            po = [o for o in (off(m, a, z) for m in members) if o is not None]
            return band_stats(off(NAME, a, z), np.array(po).reshape(-1, 2))

        st = stats("e1", "e3", uids[1:])
        per_band[b] = st
        if st is not None:
            st["angle_to_2736_6882_residual_deg"] = angle_deg(
                np.array(st["target_residual_mas"]), np.array(prev_vec[b])
            )
        missing = {"not_measured": "target or < 10 peer centroids", "label": "derived"}
        out["bands"][b] = st or missing
        # Post hoc: 6882 -> 12774 with peers centroided in all three epochs.
        three = [m for m in uids[1:] if all(cents[(b, e)].get(m) is not None for e in CATS)]
        ph = stats("e2", "e3", three)
        out["post_hoc_6882_12774"][b] = ph or {"not_measured": f"{len(three)} peers in all epochs"}
    out["decision"] = decide(per_band, prev_vec) if valid else "not evaluated: validation failed"
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    sheet(args.out / "contact_sheet.png", files, cents)
    print(json.dumps({"decision": out["decision"], "validation": out["validation"]}))
    return 0


def sheet(path: Path, files: dict, cents: dict) -> None:
    """Target cutouts per band and epoch (0.6 arcsec), each epoch's own centroid marked."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy import units as u

    fig, axs = plt.subplots(len(BANDS), 3, figsize=(8, 2.8 * len(BANDS)), squeeze=False)
    for r, b in enumerate(BANDS):
        for c, ep in enumerate(("e1", "e2", "e3")):
            ax = axs[r, c]
            ax.set_xticks([])
            ax.set_yticks([])
            p = files[(b, ep)].get(NAME)
            if not p:
                ax.set_title(f"{NAME} {b} {ep}: no cutout", fontsize=7)
                continue
            data, wcs = _read(p)
            lo, hi = np.nanpercentile(data, [5, 99.5])
            ax.imshow(data, origin="lower", cmap="gray", vmin=lo, vmax=hi)
            sky = cents[(b, ep)].get(NAME)
            if sky is not None:
                x, y = wcs.world_to_pixel(sky)
                ax.plot(x, y, "+", color="r", ms=10)
            scale = np.mean([q.to_value(u.arcsec) for q in wcs.proj_plane_pixel_scales()])
            half = 0.3 / scale
            cy, cx = (np.array(data.shape) - 1) / 2
            ax.set_xlim(cx - half, cx + half)
            ax.set_ylim(cy - half, cy + half)
            ax.set_title(f"{NAME} {b} {ep} (stored orientation)", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

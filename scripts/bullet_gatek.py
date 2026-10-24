"""Bullet Cluster: gate K of the CANUCS DR1 Lenstool best-fit maps on the gold image list.

Frozen plan: ``results/lens_consistency/bullet_gatek/PLAN.md``. Gate K of ``m0417_images`` (same
code as MACS0417 / MACS1423 / MACS1149) on the 135 gold images at their system z; Lenstool's own
image-plane rms from ``chires.dat`` is reported next to ours. No extra-image or flag reading.

    python scripts/bullet_gatek.py [--grid-step 0.25] [--out DIR] [--posthoc-k1]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import numpy as np

from jwst_anomaly import lensmodel

MODEL = "bullet-canucs-maps"
OUT = Path("results/lens_consistency/bullet_gatek")
GOLD = Path("results/lens_consistency/bullet_inputs/bullet_gold_images.ecsv")
CHIRES = (
    "https://www.canfar.net/storage/vault/file/canucs/BULLET_DR1/lensing/best_fit_model/chires.dat"
)
GRID_STEP = 0.25  # m0417_images default (PLAN.md rule 1)


def parse_chires(text: str) -> dict:
    """Lenstool ``chires.dat``: total ``chimul`` row (chi2, source- and image-plane rms) and the
    per-image rows (Narcs = 1) with their image-plane offsets dx, dy (arcsec)."""
    images, total = {}, None
    for line in text.splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == "chimul":
            total = {"chi2": float(f[1]), "rms_source": float(f[5]), "rms_image": float(f[6])}
        elif len(f) == 13 and f[0].isdigit() and f[3] == "1":
            images[f[1]] = (float(f[10]), float(f[11]))
    if total is None:
        raise ValueError("chires.dat: no chimul row")
    d = np.hypot(*np.asarray(list(images.values()), float).T) if images else np.array([])
    total["n_images"] = len(images)
    total["rms_from_rows"] = float(np.sqrt(np.mean(d**2))) if len(d) else float("nan")
    return {"total": total, "images": images}


def per_system(ip) -> dict:
    """Reported, not a gate: rms (arcsec) and catalogued / predicted counts per system."""
    out = {}
    sysv = np.asarray(ip["system"], str)
    d = np.asarray(ip["dtheta_arcsec"], float)
    for s in dict.fromkeys(sysv):
        m = sysv == s
        out[s] = {
            "rms_arcsec": mi.rms_of(d[m]),
            "n_catalogued": int(m.sum()),
            "n_predicted": int(np.max(ip["n_predicted"][m])),
        }
    return out


def run(out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    images = mi.load_images(GOLD)
    solver = mi.Solver(MODEL, images, grid_step)
    ip = solver.image_plane(images)
    lc._write(ip, out / f"imageplane_bullet_{MODEL}.ecsv")
    systems = per_system(ip)
    with urllib.request.urlopen(CHIRES, timeout=60) as r:  # noqa: S310 (fixed https URL)
        raw = r.read()
    lt = parse_chires(raw.decode())
    summary = {
        "plan": str(out / "PLAN.md"),
        "model": MODEL,
        "model_sha256": solver.model.sha256,
        "grid_step_arcsec": grid_step,
        "gate_K": mi.gate_k(ip),
        "extra_predicted": {
            s: v for s, v in systems.items() if v["n_predicted"] > v["n_catalogued"]
        },
        "fewer_predicted": {
            s: v for s, v in systems.items() if v["n_predicted"] < v["n_catalogued"]
        },
        "per_system": systems,
        "lenstool_chires": {
            "url": CHIRES,
            "sha256": hashlib.sha256(raw).hexdigest(),
            **lt["total"],
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def posthoc_k1(out: Path = OUT, families=("K1b", "K1e", "K1c", "K1d")) -> dict:
    """Post hoc (not in PLAN.md): per family, the images each catalogued image's own back-traced
    source predicts, and the back-traced sources' largest offset from their mean (arcsec)."""
    images = mi.load_images(GOLD)
    s = mi.Solver(MODEL, images, GRID_STEP)
    res = {}
    for fam in families:
        t = images[np.asarray(images["system"], str) == fam]
        bt = lensmodel.backtrace_images(s.model, t, {})
        ox, oy = s.model.to_frame(np.asarray(t["ra"], float), np.asarray(t["dec"], float))
        rows = {}
        for r in bt:
            p = lensmodel.find_images(
                s.model, s.grid, float(r["beta_x"]), float(r["beta_y"]), float(r["z_used"])
            )
            d = [
                float(np.min(np.hypot(ox[j] - p["x"], oy[j] - p["y"]))) if len(p) else None
                for j in range(len(t))
            ]
            rows[str(r["image_id"])] = {
                "n_predicted": len(p),
                "dtheta_arcsec_per_catalogued_image": d,
            }
        b = np.c_[bt["beta_x"], bt["beta_y"]].astype(float)
        res[fam] = {
            "own_source_solutions": rows,
            "source_scatter_arcsec": float(np.max(np.hypot(*(b - b.mean(0)).T))),
        }
    (out / "posthoc_k1_sources.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def posthoc_compare(out: Path = OUT, skip=("K1b", "K1e")) -> dict:
    """Post hoc (not in PLAN.md): our per-image residuals against Lenstool's ``chires.dat`` rows,
    and gate K's numbers without the ``skip`` systems (reported, not a gate)."""
    from astropy.table import Table

    ip = Table.read(out / f"imageplane_bullet_{MODEL}.ecsv")
    with urllib.request.urlopen(CHIRES, timeout=60) as r:  # noqa: S310 (fixed https URL)
        lt = parse_chires(r.read().decode())["images"]
    ids = np.asarray(ip["image_id"], str)
    d = np.asarray(ip["dtheta_arcsec"], float)
    ltd = np.array([np.hypot(*lt[i]) if i in lt else np.nan for i in ids])
    keep = ~np.isin(np.asarray(ip["system"], str), skip)
    res = {
        "n_images_in_both": int(np.sum(np.isfinite(d) & np.isfinite(ltd))),
        "median_ours_arcsec": float(np.nanmedian(d)),
        "median_lenstool_arcsec": float(np.nanmedian(ltd)),
        "per_image_correlation": float(np.corrcoef(d[keep], ltd[keep])[0, 1]),
        "skipped_systems": list(skip),
        "without_skipped": {
            "n_images": int(keep.sum()),
            "rms_arcsec": mi.rms_of(d[keep]),
            "rms_lenstool_arcsec": mi.rms_of(ltd[keep]),
            "matched_fraction": float(np.mean(d[keep] <= mi.MATCH_ARCSEC)),
        },
        "worst_ours": {str(ids[k]): [float(d[k]), float(ltd[k])] for k in np.argsort(-d)[:6]},
    }
    (out / "posthoc_compare.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def run_wls(out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    """``PLAN_1b.md`` (deviation after gate K failed): gate K with the WLS family source of
    ``a2744_uncover.imageplane``; everything else as :func:`run`."""
    import a2744_uncover as au

    images = mi.load_images(GOLD)
    solver = mi.Solver(MODEL, images, grid_step)
    bt = lensmodel.backtrace_images(solver.model, images, {})
    ip = au.imageplane(solver.model, solver.grid, bt)
    lc._write(ip, out / f"imageplane_wls_bullet_{MODEL}.ecsv")
    systems = per_system(ip)
    res = {
        "plan": str(out / "PLAN_1b.md"),
        "label": "post-failure deviation (D-034 addendum k precedent)",
        "model_sha256": solver.model.sha256,
        "grid_step_arcsec": grid_step,
        "gate_K_wls": mi.gate_k(ip),
        "extra_predicted": {s: v for s, v in systems.items() if v["n_predicted"] > v["n_catalogued"]},
        "fewer_predicted": {s: v for s, v in systems.items() if v["n_predicted"] < v["n_catalogued"]},
    }
    (out / "summary_wls.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid-step", type=float, default=GRID_STEP)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--posthoc-k1", action="store_true", help="post hoc K1 source check only")
    ap.add_argument("--wls", action="store_true", help="PLAN_1b.md: gate K with the WLS source")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.posthoc_k1:
        posthoc_k1(args.out)
        posthoc_compare(args.out)
        return 0
    if args.wls:
        print(json.dumps(run_wls(args.out, args.grid_step)["gate_K_wls"], indent=1))
        return 0
    s = run(args.out, args.grid_step)
    print(json.dumps({k: s[k] for k in ("gate_K", "lenstool_chires")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

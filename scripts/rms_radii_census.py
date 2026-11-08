"""Census of stored image-plane tables for rms-scaled forced_check radii.

Plan: results/lens_consistency/rms_radii/PLAN.md. Offline: reads the stored
results/lens_consistency/*/imageplane*.ecsv tables."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("lc", ROOT / "scripts" / "lens_consistency.py")
lc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lc)

MATCH = 1.5  # predict_counter_images default match radius (arcsec)
SCALABLE_MAX_RMS = 2.0  # above this the scaled match radius exceeds the 3" cap (PLAN; ASSUMPTION)
MIN_N = 3


def census_row(path: Path) -> dict | None:
    d = np.asarray(Table.read(path, format="ascii.ecsv")["dtheta_arcsec"], dtype=float)
    d = d[np.isfinite(d)]
    if len(d) < MIN_N:
        return None
    rms = lc.image_plane_rms(d)
    match_s, search_s = lc.rms_scaled_radii(rms, match=MATCH)
    srt = np.sort(d)
    return {
        "table": str(path.relative_to(ROOT)),
        "n": int(len(d)),
        "rms_arcsec": round(rms, 4),
        "rms_leave_max_out_arcsec": round(lc.image_plane_rms(srt[:-1]), 4),
        "robust_scale_arcsec": round(1.4826 * float(np.median(np.abs(d))), 4),
        "max_dtheta_arcsec": round(float(srt[-1]), 4),
        "match_scaled_arcsec": round(match_s, 4),
        "search_scaled_arcsec": round(search_s, 4),
        "n_unmatched_fixed": int(np.sum(d > MATCH)),
        "n_unmatched_scaled": int(np.sum(d > match_s)),
        "outlier_driven": bool(lc.rms_outlier_driven(d)),
        "group": (
            "rms<=1"
            if rms <= lc.RMS_SCALE_THRESHOLD
            else ("scalable" if rms <= SCALABLE_MAX_RMS else "gate_k_failure")
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results/lens_consistency/rms_radii")
    args = ap.parse_args()
    paths = sorted((ROOT / "results/lens_consistency").glob("*/imageplane*.ecsv"))
    rows = [r for r in (census_row(p) for p in paths) if r]
    Table(rows=rows).write(args.out / "census.ecsv", format="ascii.ecsv", overwrite=True)
    scal = [r for r in rows if r["group"] == "scalable"]
    n_out = sum(r["outlier_driven"] for r in scal)
    reading = (
        "no scalable table"
        if not scal
        else ("rms scale outlier-driven" if 2 * n_out >= len(scal) else "rms scale broad")
    )
    summary = {
        "n_tables_read": len(paths),
        "n_tables_used": len(rows),
        "groups": {
            g: sum(r["group"] == g for r in rows) for g in ("rms<=1", "scalable", "gate_k_failure")
        },
        "scalable_outlier_driven": n_out,
        "reading": reading,
        "gate_k_failure_outlier_driven": sum(
            r["outlier_driven"] for r in rows if r["group"] == "gate_k_failure"
        ),
        "scalable_unmatched_fixed": sum(r["n_unmatched_fixed"] for r in scal),
        "scalable_unmatched_scaled": sum(r["n_unmatched_scaled"] for r in scal),
        "provenance": "derived",
        "assumptions": {
            "match_arcsec": MATCH,
            "rms_threshold_arcsec": lc.RMS_SCALE_THRESHOLD,
            "cap_arcsec": lc.MAX_SEARCH_ARCSEC,
            "scalable_max_rms_arcsec": SCALABLE_MAX_RMS,
            "min_n": MIN_N,
        },
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

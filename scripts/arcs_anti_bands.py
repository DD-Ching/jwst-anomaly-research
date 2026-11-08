"""Second-band (F150W, F444W / F410M) shapes of the 3 single-component pooled ``arcs`` anti rows.

Plan (frozen before the run): ``results/arcs/anti-bands/PLAN.md``. Matching, usability
and re-classing reuse ``scripts/arcs_bands.py`` unchanged; this script adds the band gate,
the per-row label and the frozen outcome. Catalogues are read from S3 into the cache, pinned
by sha256 in ``data/manifests/arcs_anti_bands.ecsv``, and deleted after reading (``--keep``
keeps them). Shapes ``observed``; classes and labels ``derived``.

    python scripts/arcs_anti_bands.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import paths
from jwst_anomaly.photometry import fetch_catalog

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ab = _load("arcs_bands")
lc = ab.lc

# row label, field directory, primary run, observation, SW and LW band of the same observation
ROWS = (
    (2194, "whl0137-relics", "whl0137-relics-lenstool", "jw02282-o010_t001", "f150w", "f444w"),
    (581, "rxc0600-relics", "rxc0600-relics-glafic", "jw06796-o010_t001", "f150w", "f444w"),
    (2157, "macs0308-relics", "macs0308-relics-glafic", "jw05293-o005_t001", "f150w", "f410m"),
)
BANDS = ("SW", "LW")

# frozen thresholds (ASSUMPTIONS, PLAN.md)
GATE_MIN_ROWS = 10
GATE_MIN_FRAC = 0.8
AGREE_DEG = 15.0


def band_gate(n_usable: int, n_stay: int) -> bool:
    """Rule 1: the band passes with >= 10 usable control rows and >= 0.8 of them still aligned."""
    return n_usable >= GATE_MIN_ROWS and n_stay / n_usable >= GATE_MIN_FRAC


def band_verdict(usable: bool, dpa: float) -> str:
    """Rule 2: ``agree`` if |dpa| <= 15 deg, ``differs`` otherwise, ``unusable`` without a shape."""
    if not usable or not np.isfinite(dpa):
        return "unusable"
    return "agree" if abs(dpa) <= AGREE_DEG else "differs"


def row_label(verdicts: dict[str, str], passed: dict[str, bool]) -> str:
    """Rule 3 (passing bands): ``differs`` -> blend_flag; else any ``agree`` -> same_axis."""
    used = [v for b, v in verdicts.items() if passed.get(b, False)]
    if "differs" in used:
        return "blend_flag"
    if "agree" in used:
        return "same_axis"
    return "undetermined"


def outcome(labels: list[str]) -> str:
    """Rule 4: the frozen outcome over the 3 rows."""
    n_blend, n_same = labels.count("blend_flag"), labels.count("same_axis")
    if n_blend >= 2:
        return (
            "open part of the pooled anti excess explained by colour blends or colour structure"
            " (exploratory)"
        )
    if n_same >= 2 and n_blend == 0:
        return "blend veto passed; the rows stay open (not a lensing test)"
    return "mixed / undetermined"


MANIFEST = "arcs_anti_bands"


def fetch(obs: str, band: str, cache: Path) -> tuple[Path, str]:
    """Pipeline catalogue of ``obs`` in ``band`` from S3, sha256-pinned by the manifest."""
    name = f"{obs}_nircam_clear-{band}_cat.ecsv"
    sha = ab.manifest_sha(MANIFEST, name)
    return fetch_catalog(ab.s3_url(name), sha, cache_dir=cache), sha


def row_index(t: Table, label: int) -> int:
    """Index of the single row ``label`` in ``t``; refuses a missing or duplicated label."""
    hits = np.flatnonzero(np.asarray(t["label"]) == label)
    if len(hits) != 1:
        raise SystemExit(f"error: label {label} found {len(hits)} times")
    return int(hits[0])


def band_args():
    """``arcs_bands`` defaults (match 0.3", e >= 0.3, S/N >= 10, class thresholds)."""
    p = argparse.Namespace(
        match_arcsec=0.3,
        min_e444=0.3,
        min_snr444=10.0,
        aligned_deg=30.0,
        anti_deg=60.0,
        z_margin=0.1,
    )
    return p


def run(out: Path, cache: Path, keep: bool) -> dict:
    args = band_args()
    if len({r[1] for r in ROWS}) != len(ROWS):
        raise SystemExit("error: one row per field (rule 1 pools each field's control once)")
    per_band = {b: [] for b in BANDS}  # control tables per band
    rows, inputs = {}, {}
    for label, field, model_name, obs, sw, lw in ROWS:
        stored = Table.read(
            paths.repo_root() / "results" / "arcs" / field / model_name / "arcs.ecsv"
        )
        model, _, _ = lc.load_model(model_name)
        lc.apply_frame_offset(model_name, model)
        rows[label] = {"field": field, "model": model_name, "bands": {}}
        for band_key, band in zip(BANDS, (sw, lw), strict=True):
            path, sha = fetch(obs, band, cache)
            inputs[path.name.split("_", 1)[1]] = sha
            try:
                shapes = lc.load_shapes(path)
            finally:
                if not keep:
                    path.unlink()
            t = ab.reclassify(model, ab.match_second_band(stored, shapes, args.match_arcsec), args)
            usable = ab.usable_mask(t, args)
            ctrl = (
                t["strong_shear"]
                & (t["population"] == "background")
                & (t["orientation_class"] == "aligned")
            )
            per_band[band_key].append((t[ctrl & usable], field, band))
            i = row_index(t, label)
            r = t[i]
            if str(r["orientation_class"]) != "anti":
                raise SystemExit(f"error: row {label} is not anti in {model_name}")
            dpa = float(r["dpa_bands"])
            rows[label]["bands"][band_key] = {
                "band": band,
                "status": ab.status_444(r, bool(usable[i])),
                "verdict": band_verdict(bool(usable[i]), dpa),
                "sep_nearest": round(float(r["sep_nearest"]), 2),
                "dpa": None if not np.isfinite(dpa) else round(dpa, 1),
                "e": None if not np.isfinite(r["e_444"]) else round(float(r["e_444"]), 3),
                "snr": None if not np.isfinite(r["snr_444"]) else round(float(r["snr_444"]), 1),
                "pa_200": round(float(r["pa_obs"]), 1),
                "g_min_200": round(float(r["g_min"]), 2),
                "e_200": round(float(r["ellipticity"]), 3),
            }
    gate = {}
    for b in BANDS:
        parts = per_band[b]
        n_use = sum(len(c) for c, _, _ in parts)
        n_stay = sum(int(np.sum(c["class_444"] == "aligned")) for c, _, _ in parts)
        dpa = np.concatenate([np.abs(np.asarray(c["dpa_bands"], float)) for c, _, _ in parts])
        gate[b] = {
            "per_field": {
                f"{f} {band}": [len(c), int(np.sum(c["class_444"] == "aligned"))]
                for c, f, band in parts
            },
            "usable": n_use,
            "stay_aligned": n_stay,
            "passes": band_gate(n_use, n_stay),
            "dpa_median": round(float(np.median(dpa)), 1) if len(dpa) else None,
            "dpa_p95": round(float(np.percentile(dpa, 95)), 1) if len(dpa) else None,
        }
    passed = {b: gate[b]["passes"] for b in BANDS}
    for r in rows.values():
        r["label"] = row_label({b: v["verdict"] for b, v in r["bands"].items()}, passed)
    labels = [r["label"] for r in rows.values()]
    summary = {
        "plan": "results/arcs/anti-bands/PLAN.md",
        "params": {
            "match_arcsec": args.match_arcsec,
            "min_e": args.min_e444,
            "min_snr": args.min_snr444,
            "aligned_deg": args.aligned_deg,
            "anti_deg": args.anti_deg,
            "gate_min_rows": GATE_MIN_ROWS,
            "gate_min_frac": GATE_MIN_FRAC,
            "agree_deg": AGREE_DEG,
        },
        "inputs_sha256": inputs,
        "gate": gate,
        "rows": {str(k): v for k, v in rows.items()},
        "outcome": outcome(labels),
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--out", type=Path, default=paths.repo_root() / "results" / "arcs" / "anti-bands"
    )
    p.add_argument("--cache", type=Path, default=paths.cache_dir() / "anti-bands")
    p.add_argument("--keep", action="store_true", help="keep the downloaded catalogues")
    a = p.parse_args(argv)
    print(json.dumps(run(a.out, a.cache, a.keep), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

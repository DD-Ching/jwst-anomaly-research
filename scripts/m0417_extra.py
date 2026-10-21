"""MACS0417 / MACS1423: forced photometry at every image the CANUCS v1 maps predict for the published
families (D-029 route: ``lens_consistency.predict_counter_images`` + ``forced_check``, unchanged).

The published image lists of ``m0417_images.py`` take the place of a Lenstool ``arcs.dat``; each
family's z is the ``z_used`` of that run's CANUCS image-plane table (z_spec, else its fitted z).
Gate C runs the catalogued images through the same forced photometry as a positive control.
Frozen plan: ``results/lens_consistency/macs0417_extra/PLAN.md``.

    python scripts/m0417_extra.py [--fields macs0417 macs1423]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m1149_sharon as ms
import numpy as np
from astropy.table import Table

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/macs0417_extra")
IMAGES = Path("results/lens_consistency/macs0417_images")
S3 = "s3://stpubdata/jwst/public/jw01208/L3/t/{0}/jw01208-{0}_{1}_nircam_clear-f200w_{2}"
FIELDS = {
    "macs0417": ("jauzac19_images.ecsv", "macs0417-canucs-maps", "o006", "t003"),
    "macs1423": ("kaleidoscope_m1423_images.ecsv", "macs1423-canucs-maps", "o010", "t005"),
}
# Frozen rules (PLAN.md; ASSUMPTIONs)
GATE_C_MIN = 0.8  # recovered (or confused) fraction of the testable catalogued images
TESTABLE = ("recovered", "absent", "confused", "ambiguous")
NOT_TARGETS = ("observed", "demagnified", "outside")
GRID_STEP = 0.1  # arcsec, as ``lens_consistency images``


def flag(image_class: str, forced_class: str) -> bool:
    """Frozen per-image rule: an unmatched, magnified, in-footprint predicted image that is absent."""
    return image_class not in NOT_TARGETS and forced_class == "absent"


def gate_c(forced_classes) -> dict:
    """Frozen per-field gate on the catalogued images run as forced targets."""
    fc = [str(c) for c in forced_classes]
    testable = [c for c in fc if c in TESTABLE]
    ok = sum(c in ("recovered", "confused") for c in testable)
    frac = ok / len(testable) if testable else float("nan")
    return {
        "n_targets": len(fc),
        "n_testable": len(testable),
        "n_recovered": ok,
        "fraction": frac,
        "passes": bool(testable) and frac >= GATE_C_MIN,
        "classes": {c: fc.count(c) for c in sorted(set(fc))},
    }


def load_images(field: str) -> Table:
    imgfile, model_name, _, _ = FIELDS[field]
    t = ms.gold_images(IMAGES / imgfile)
    ip = Table.read(IMAGES / f"imageplane_{field}_{model_name}.ecsv")
    z = dict(zip(np.asarray(ip["image_id"], str), np.asarray(ip["z_used"], float), strict=True))
    t["z"] = [z.get(str(i), np.nan) for i in t["image_id"]]
    return t


def run_field(field: str, cache: Path) -> dict:
    _, name, obs, tile = FIELDS[field]
    images = load_images(field)
    model, _, _ = lc.load_model(name)
    lc.apply_frame_offset(name, model)
    bt = lensmodel.backtrace_images(model, images, {})
    traced = images[np.isfinite(bt["beta_x"])]
    half = lc.grid_half_width(model, traced)
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, GRID_STEP), half, GRID_STEP
    )
    shapes = lc.load_shapes(cache / Path(S3.format(obs, tile, "cat.ecsv")).name)
    depth = lc.depth_mag(shapes)
    table, unpredicted = lc.predict_counter_images(model, grid, bt, shapes, depth=depth)
    stamp, close = lc.image_stamper(
        S3.format(obs, tile, "i2d.fits"), lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.1
    )
    try:
        lc.forced_check(table, bt, stamp, 1.0, unpredicted=unpredicted)
        control = table[np.asarray(table["image_class"], str) == "observed"]
        control["image_class"] = "control"  # catalogued images as forced targets (gate C)
        lc.forced_check(control, bt, stamp, 1.0, unpredicted=unpredicted)
    finally:
        close()
    table["reading"] = [
        "lensing_residual_input" if flag(str(r["image_class"]), str(r["forced_class"])) else "no_flag"
        for r in table
    ]
    gc = gate_c(control["forced_class"])
    lc._write(table, OUT / f"predicted_{field}.ecsv")
    lc._write(control, OUT / f"control_{field}.ecsv")
    cls = np.asarray(table["image_class"], str)
    return {
        "model": name,
        "model_sha256": model.sha256,
        "image": S3.format(obs, tile, "i2d.fits"),
        "depth_mag": depth,
        "unpredicted": unpredicted,
        "image_classes": {c: int((cls == c).sum()) for c in sorted(set(cls))},
        "forced_classes_targets": {
            c: int(((np.asarray(table["forced_class"], str) == c) & ~np.isin(cls, NOT_TARGETS)).sum())
            for c in sorted(set(np.asarray(table["forced_class"], str)))
        },
        "gate_C": gc,
        "flags": (
            [
                f"{r['system']}@{r['ra']:.6f},{r['dec']:.6f}"
                for r in table
                if r["reading"] == "lensing_residual_input"
            ]
            if gc["passes"]
            else []
        ),
        "flags_if_gate_ignored": int((table["reading"] == "lensing_residual_input").sum()),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fields", nargs="+", default=list(FIELDS))
    ap.add_argument("--cache", type=Path, default=Path("data/cache/m0417x"))
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {f: run_field(f, args.cache) for f in args.fields}
    summary["rules"] = {"GATE_C_MIN": GATE_C_MIN, "TESTABLE": TESTABLE, "NOT_TARGETS": NOT_TARGETS}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({f: (summary[f]["gate_C"], summary[f]["flags"]) for f in args.fields}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

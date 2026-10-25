"""Paterno-Mahler et al. 2018 images of SPT-CL J0615-5746 vs the RELICS Lenstool v1 and GLAFIC v2
maps.

Frozen plan: ``results/lens_consistency/spt0615_images/PLAN.md`` (D-034 addendum). Image list:
the paper's Table 2, transcribed in ``TABLE2``. Per model: gates P, K0 and F on families 1-3,
reading B on families 1 and 3 at z_spec, reading A on 3.3 / 3.4 / 3.6 (family 3b), extra /
fewer predicted images and the amplitude class (rule 7, families 1-3). Rule helpers are shared
with ``slice_relics_macs0358``.

    python scripts/spt0615_images.py --build-images
    python scripts/lens_consistency.py --model spt0615-relics-lenstool \
        --out results/lens_consistency/spt0615_images/validate validate   # and -glafic
    python scripts/spt0615_images.py [--grid-step 0.25]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import macs0647_clash_heldout as ho
import macs2129_clash_heldout as h29
import numpy as np
import slice_relics_glafic as sg
import slice_relics_heldout as srh
import slice_relics_rxj0152 as base
from astropy.coordinates import SkyCoord
from astropy.table import Table

OUT = Path("results/lens_consistency/spt0615_images")
IMAGES = OUT / "spt_images.ecsv"
CLUSTER = "spt0615"
LENSTOOL = "spt0615-relics-lenstool"
MODELS = (LENSTOOL, "spt0615-relics-glafic")
HELD = ("3b",)  # family 3b: all six system 3 images (exploratory reading A)
READ_A = ("3b.3", "3b.4", "3b.6")  # the images the paper's model 1 omits
AMP_FAMILIES = ("1", "2", "3")  # rule 7 median (3b repeats family 3)
TEX_SHA256 = "25b66b8d8a436fafe7990f7de4a76f80ef4d5e4e03d8d426cb189e388b2dc807"
EPRINT_SHA256 = "83c2bc6cab221cb2893d3719fd75f17dc924d4d471bc6dc83f642f529f6cdb26"

# Paterno-Mahler et al. 2018, Table 2 (ID, RA, Dec J2000); observed
TABLE2 = {
    "1.1": ("06 15 52.22", "-57 46 49.9"),
    "1.2": ("06 15 51.87", "-57 46 46.9"),
    "1.3": ("06 15 51.05", "-57 46 44.7"),
    "2.1": ("06 15 49.37", "-57 46 52.8"),
    "2.2": ("06 15 49.70", "-57 46 57.1"),
    "3.1": ("06 15 49.96", "-57 46 53.5"),
    "3.2": ("06 15 49.38", "-57 46 33.9"),
    "3.3": ("06 15 57.45", "-57 47 38.6"),
    "3.4": ("06 15 51.78", "-57 46 33.7"),
    "3.5": ("06 15 51.39", "-57 46 45.9"),
    "3.6": ("06 15 51.95", "-57 46 29.2"),
}
# family -> (images, table z, spectroscopic); family 2's z is the paper's model 1 value
FAMILIES = {
    "1": (("1.1", "1.2", "1.3"), 1.358, True),
    "2": (("2.1", "2.2"), 2.43, False),
    "3": (("3.1", "3.2", "3.5"), 4.013, True),
    "3b": (("3.1", "3.2", "3.3", "3.4", "3.5", "3.6"), 4.013, False),
}


def q1_answer(classes: dict[str, str | None]) -> str:
    """Rule 7 as in ``slice_relics_rxj0152.q1_answer``, with this cluster's Lenstool v1 model."""
    return base.q1_answer({(base.LENSTOOL if m == LENSTOOL else m): c for m, c in classes.items()})


def amplitude_median(families: dict[str, dict]) -> float | None:
    """Rule 7: median factor over ``AMP_FAMILIES`` that received one (edge families get none)."""
    f = [r["scale"] for k, r in families.items() if k in AMP_FAMILIES]
    return float(np.median(f)) if f else None


def build_images() -> Table:
    rows = []
    for fam, (ids, z, spec) in FAMILIES.items():
        for i in ids:
            c = SkyCoord(*TABLE2[i], unit=("hourangle", "deg"))
            iid = i if fam != "3b" else "3b." + i.split(".")[1]
            rows.append((iid, fam, i, c.ra.deg, c.dec.deg, z, spec, CLUSTER, fam in HELD, spec))
    names = ("image_id", "system", "paper_id", "ra", "dec", "z", "z_spec", "cluster")
    t = Table(rows=rows, names=(*names, "new_in_jwst", "use_spec"))
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        "Paterno-Mahler et al. 2018, ApJ 863, 154 (arXiv:1805.09834), Table 2, transcribed by hand "
        f"(spt0615_arxiv.tex sha256 {TEX_SHA256}; e-print sha256 {EPRINT_SHA256}); z_spec True: "
        "spectroscopic; family 2 z: the paper's model 1 value (model_prediction); family 3b "
        "repeats "
        "system 3 with the less-secure images (new_in_jwst here means held out, PLAN.md)"
    )
    return t


def load() -> Table:
    t = Table.read(IMAGES, format="ascii.ecsv")
    t["image_id"] = np.asarray(t["image_id"], str)
    t["system"] = np.asarray(t["system"], str)
    t["held_out"] = np.asarray(t["new_in_jwst"], bool)
    if set(t["system"][t["held_out"]]) != set(HELD):
        raise ValueError("stale image list (new_in_jwst vs HELD); rebuild it")
    t["z_spec"] = np.asarray(t["use_spec"], bool)
    t["a"] = np.ones(len(t))
    return t


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid-step", type=float, default=0.25)
    ap.add_argument("--build-images", action="store_true", help="write the image list")
    args = ap.parse_args(argv)
    if args.build_images:
        build_images().write(IMAGES, format="ascii.ecsv", overwrite=True)
        print("->", IMAGES)
        return 0
    t = load()
    images = Table.read(IMAGES, format="ascii.ecsv")
    ids = np.asarray(t["image_id"], str)
    read_a = np.isin(ids, READ_A)
    spec = np.asarray(t["z_spec"], bool)
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "rules": {
            "K0_MAX_RMS_ARCSEC": ho.K0_MAX_RMS_ARCSEC,
            "K0_MIN_MATCHED": ho.K0_MIN_MATCHED,
            "MATCH_ARCSEC": ho.MATCH_ARCSEC,
            "FRAME_MAX_ARCSEC": ho.FRAME_MAX_ARCSEC,
            "z_grid": [ho.Z_STEPS, ho.Z_MIN_OFFSET, ho.Z_MAX],
            "grid_step": args.grid_step,
            "near_one": sg.NEAR_ONE,
            "lenstool_like_max": sg.LENSTOOL_LIKE_MAX,
            "p_ratio": sg.P_RATIO,
            "amplitude_families": list(AMP_FAMILIES),
        },
        "models": {},
    }
    d_a, d_spec, usable, classes = {}, {}, [], {}
    for model in MODELS:
        validate = json.loads((OUT / "validate" / model / "validate.json").read_text())
        p_ok = sg.gate_p(validate)
        entry, ip, ds = h29.run_model(model, t, args.grid_step)
        ok = entry["usable"] and p_ok
        d = np.asarray(ip["dtheta_arcsec"], float)
        d_a[model], d_spec[model] = d[read_a], ds[spec]
        if ok:
            usable.append(model)
        amp = srh.posthoc_scale({CLUSTER: entry}, images, {CLUSTER: model})[CLUSTER]
        med = amplitude_median(amp["families"])
        classes[model] = sg.amplitude_class(med) if ok and med is not None else None
        entry.update(
            {
                "gate_P": {"passes": p_ok, "validate": str(OUT / "validate" / model)},
                "usable_with_P": ok,
                "misses_A": base.per_model_misses(ids[read_a], d[read_a]),
                "misses_B": base.per_model_misses(ids[spec], ds[spec]),
                "amplitude": {**amp, "median": med, "class": classes[model]},
            }
        )
        summary["models"][model] = entry
        ip.meta["model"] = model
        lc._write(ip, OUT / f"imageplane_{model}.ecsv")
    ra = ho.readings(d_a, usable)
    rb = h29.label_b(ho.readings(d_spec, usable))
    summary.update(
        {
            "images": len(t),
            "reading_A_images": int(read_a.sum()),
            "spec_z_images": int(spec.sum()),
            "usable": usable,
            "flags_A": [str(i) for i, v in zip(ids[read_a], ra, strict=True) if v != "no_flag"],
            "flags_B": [str(i) for i, v in zip(ids[spec], rb, strict=True) if v != "no_flag"],
            "amplitude_classes": classes,
            "q1_answer": q1_answer(classes),
        }
    )
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    for m, e in summary["models"].items():
        print(m, e["usable_with_P"], e["gate_P"]["passes"], e["gate_K0"], e["frame_F"], sep="\n  ")
        print("  A", e["misses_A"], "\n  B", e["spec_z_B"], e["misses_B"])
        print("  z_fit", e["z_fit"], "extra", e["extra_predicted"], "fewer", e["fewer_predicted"])
        print("  amplitude", e["amplitude"]["median"], e["amplitude"]["class"])
    print("joint", summary["flags_A"], summary["flags_B"], summary["q1_answer"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

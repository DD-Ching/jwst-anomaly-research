"""Post hoc (MACS1149 plan rule 6; not gates): MACS0417 and MACS1423 re-run with the
containing-source rule, each otherwise as its own frozen run (MACS0417: ``m0417_canucsref``,
no local-noise rule; MACS1423: ``m1423_canucsref`` with it). Both fields were seen while the rule
was designed. Also MACS1149 gate K at a 0.1" grid step. Writes
``results/lens_consistency/macs1149_canucsref/posthoc_contain.json``.

    python scripts/m1149_canucsref_posthoc.py --cat0417 <cat> --cat1423 <cat> --image1423 <img>
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import m0417_arcref as ar
import m0417_canucsref as cr
import m0417_images as mi
import m1149_canucsref as m11
import m1423_canucsref as m14
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from m0417_canucsref_posthoc import _clean

from jwst_anomaly.provenance import file_sha256

OUT = Path("results/lens_consistency/macs1149_canucsref/posthoc_contain")
KEYS = ("gate_C2", "control_classes", "target_classes", "flags")
# the files of the frozen MACS0417 / MACS1423 runs (their summary.json)
SHA = {
    "cat0417": "ca6b48aa6404444d2e0b8965a35376ff032917071ea64cbf6306d06ee57d90dd",
    "cat1423": "c257d8dff95319f4eea981cc838269c04aff0a79bac603c7ed64100deac622b0",
    "image1423": "b7c6739f6cbf25afaf80d10763eab1c6f33bf441ea69a63909fb386182736403",
}


GATE_K_STEP = 0.1  # post hoc: gate K at the stage-2 grid step instead of 0.25


def gate_k_at(step: float) -> dict:
    """Post hoc: MACS1149 gate K with another deflection-grid step (tables in a temporary dir)."""
    with tempfile.TemporaryDirectory() as tmp:
        m11.write_images(Path(tmp))
        m11.register()
        r = mi.run_field(m11.FIELD, Path(tmp), step)
    return {"grid_step_arcsec": step, "gate_K": r[m11.MODEL]["gate_K"]}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat0417", type=Path, required=True)
    ap.add_argument("--cat1423", type=Path, required=True)
    ap.add_argument("--image1423", type=Path, required=True)
    args = ap.parse_args(argv)
    for k, sha in SHA.items():
        if file_sha256(getattr(args, k)) != sha:
            raise ValueError(f"--{k} is not the file of the frozen run (sha256)")
    res = {"inputs_sha256": SHA}
    cat = cr.to_cat(Table.read(args.cat0417))
    depth = cr.depth_mag(cat["flux"], cat["snr"])
    s = ar.run(Path("."), cat=cat, depth=depth, cat_inputs={}, out=OUT / "macs0417", contain=True)
    res["macs0417"] = {k: s[k] for k in KEYS}

    t = Table.read(args.cat1423)
    cat = cr.to_cat(t)
    depth = cr.depth_mag(cat["flux"], cat["snr"])
    out = OUT / "macs1423"
    s = ar.run(
        Path("."),
        cat=cat,
        depth=depth,
        cat_inputs={},
        out=out,
        field=m14.FIELD,
        frame=m14.FRAME,
        offset=(0.0, 0.0),
        contain=True,
    )
    with fits.open(args.image1423) as hl:
        if hl[0].header.get("BUNIT") != "nJy":
            raise ValueError(f"unexpected BUNIT {hl[0].header.get('BUNIT')!r}")
        w, d = WCS(hl[0].header), np.asarray(hl[0].data, float)
    src_x, src_y = (
        np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(t["RA"], t["DEC"], unit="deg"))
    )
    rng = np.random.default_rng(m14.SEED)
    ctrl = Table.read(out / f"control_{m14.FIELD}.ecsv")
    pred = Table.read(out / f"predicted_{m14.FIELD}.ecsv")
    m14.apply_local(ctrl, d, w, src_x, src_y, rng)
    if len(pred):
        m14.apply_local(pred, d, w, src_x, src_y, rng)
    c2 = ar.gate_c2(list(ctrl["class_local"]), list(ctrl["arc_ref"]))
    res["macs1423"] = {
        "gate_C2": c2,
        "control_classes": m14._count(ctrl["class_local"]),
        "target_classes": m14._count(pred["class_local"]) if len(pred) else {},
        "flags": [
            f"{r['system']}@{r['ra']:.6f},{r['dec']:.6f}"
            for r in pred
            if c2["passes"] and str(r["class_local"]) == "absent"
        ],
        "controls": {
            str(r["images"]): str(r["class_local"]) for r in ctrl
        },  # per-control classes, for the 3.2 / 4.3 comparison
        "targets": [
            {
                "system": str(r["system"]),
                "image_class": str(r["image_class"]),
                "class": str(r["class_local"]),
            }
            for r in pred
        ],
        "without_local_noise_rule": {k: s[k] for k in KEYS},
    }
    res["macs1149_gate_K_grid_step"] = gate_k_at(GATE_K_STEP)
    path = OUT.parent / "posthoc_contain.json"
    path.write_text(json.dumps(_clean(res), indent=1, allow_nan=False) + "\n")
    print(json.dumps(_clean(res), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

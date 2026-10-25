"""El Gordo: the Caminha+2023 model vs Diego+2023 JWST-era images absent from its catalogue.

Frozen plan: ``results/lens_consistency/elgordo_heldout/PLAN.md``. Rules 1-5 of the
``mars_hff_heldout`` stage, unchanged, via ``m1149_hst_heldout.model_block``, with Caminha's
56 catalogued images as the "seen" list. The primary gate uses the held-out images of families
with a MUSE spectroscopic z; the gate over all held-out images (geometric z) is reported only.

    python scripts/elgordo_heldout.py --tex TableArcs_NEW3.tex [--out DIR]

``TableArcs_NEW3.tex`` is from the arXiv:2210.06514 source; the parsed list is stored as
``diego23_images.ecsv`` in ``--out``.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import m1149_hst_heldout as hh
import mars_hff_heldout as mo
import numpy as np
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

OUT = Path("results/lens_consistency/elgordo_heldout")
MODEL = "elgordo-caminha23"
SOURCE = "Diego et al. 2023, A&A 672, A3 (arXiv:2210.06514), Table tab_arcs (TableArcs_NEW3.tex)"
ROW = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?[a-z]+)\S*\s*&")


def _num(s: str) -> float | None:
    """A table number; ``\\no`` and -1 are missing, footnote marks (``\\ddag``) are dropped."""
    s = s.strip()
    if s in ("\\no", "-1", ""):
        return None
    return float(re.sub(r"\\[A-Za-z]+", "", s))


def parse_table(tex: str) -> Table:
    """PLAN.md Inputs: one row per image; family z = its z_spec (first row) when given, else
    z_geo (Diego's model; ASSUMPTION)."""
    rows = []
    for line in tex.splitlines():
        m = ROW.match(line)
        if not m:
            continue
        c = line.split("&")
        rows.append((m.group(1), float(c[1]), float(c[2]), _num(c[3]), _num(c[4])))
    fam_z: dict[str, tuple[float, bool]] = {}
    for iid, _, _, zs, zg in rows:
        f = re.match(r"[0-9]+", iid).group(0)
        if zs is not None:
            fam_z[f] = (zs, True)
        elif zg is not None and f not in fam_z:
            fam_z[f] = (zg, False)
    out = []
    for iid, ra, dec, _, _ in rows:
        z, spec = fam_z.get(re.match(r"[0-9]+", iid).group(0), (0.0, False))
        out.append((iid, ra, dec, z, spec))
    t = Table(rows=out, names=("image_id", "ra", "dec", "z", "z_spec"))
    t["ra"].unit = t["dec"].unit = "deg"
    t.meta = {
        "provenance": "observed (positions, z_spec); z_geo: model_prediction",
        "source": SOURCE,
    }
    return t


def primary_gate(r: Table, images: Table) -> dict:
    """PLAN.md rule 5: gate H on the held-out images of spectroscopic-z families."""
    spec = dict(
        zip(np.asarray(images["image_id"], str), np.asarray(images["z_spec"], bool), strict=True)
    )
    keep = np.array([spec[str(i)] for i in r["image_id"]], bool) if len(r) else np.array([], bool)
    g = mo.gate_h(r[keep])
    return {"gate_H": g, "outcome": mo.outcome(g), "n_rows": int(keep.sum())}


def run(tex: Path, out: Path = OUT) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    images = parse_table(tex.read_text())
    images.write(out / "diego23_images.ecsv", format="ascii.ecsv", overwrite=True)
    block = hh.model_block(MODEL, images, out)
    r = Table.read(out / f"heldout_{MODEL}.ecsv")
    prim = primary_gate(r, images)
    d = np.asarray(r["dtheta_arcsec"], float) if len(r) else np.array([])
    spec = set(np.asarray(images["image_id"][np.asarray(images["z_spec"], bool)], str))
    summary = {
        "plan": str(out / "PLAN.md"),
        "images": SOURCE,
        "n_images": len(images),
        "n_families_spec_z": len({re.match(r"[0-9]+", i).group(0) for i in spec}),
        "provenance": "positions and z_spec: observed; z_geo and residuals: model_prediction; "
        "thresholds: ASSUMPTION (mars_hff_heldout/PLAN.md)",
        "primary_spec_z": prim
        | {
            "over_3arcsec": [
                str(i) for i, x in zip(r["image_id"], d, strict=True) if x > 3.0 and str(i) in spec
            ]
        },
        "all_heldout (reported only)": block,
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tex", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    s = run(a.tex, a.out)
    b = s["all_heldout (reported only)"]
    keys = ("frame_offset", "seen", "held_out", "types", "outside_map", "gate_H", "outcome")
    print(
        json.dumps(
            _clean({"primary": s["primary_spec_z"], "all": {k: b.get(k) for k in keys}}), indent=1
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

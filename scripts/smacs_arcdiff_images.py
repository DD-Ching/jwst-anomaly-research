"""Distance of each SMACS 0723 ``arc_diff`` target to the nearest Mahler+22 ICLv2 multiple image.

Applies the rule frozen in ``results/transients/smacs_arcdiff/PLAN.md`` (ASSUMPTION): a hit
that passes ``arc_diff_morph.py`` and lies within ``ON_IMAGE_ARCSEC`` of a catalogued image is
"on an image" and goes to ``/vet-candidate``, unless it matches a known control (SN Atalanta).
The image list is the ICLv2 ``arcs.dat`` shifted by the model's JWST frame offset
(``lens_consistency.apply_frame_offset``), as in ``smacs_ref_fluxes``.

    python scripts/smacs_arcdiff_images.py results/transients/smacs_arcdiff/o001_o057/morph.json \
        --known results/transients/smacs_arcdiff/known_controls.csv \
        --out results/transients/smacs_arcdiff/o001_o057/image_distances.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402
import smacs_s17_vet as s17v  # noqa: E402

MODEL = s17v.MODEL
ON_IMAGE_ARCSEC = 1.0  # ASSUMPTION, frozen in PLAN.md before the run
KNOWN_MATCH_ARCSEC = 0.15  # ASSUMPTION: arc_diff's --match-arcsec default


def classify(verdict: str, dist_image: float, dist_known: float) -> str:
    """Frozen rule: ``known`` (within KNOWN_MATCH of a control, any verdict), else ``vet``
    (morph pass AND within ON_IMAGE of a catalogued image), else ``off_image`` (pass, farther)
    or ``vetoed`` (not a pass)."""
    if dist_known <= KNOWN_MATCH_ARCSEC:
        return "known"
    if verdict != "pass":
        return "vetoed"
    return "vet" if dist_image <= ON_IMAGE_ARCSEC else "off_image"


def nearest(targets: SkyCoord, refs: SkyCoord) -> tuple[np.ndarray, np.ndarray]:
    """Index of and distance (arcsec) to the nearest ``refs`` entry for each target."""
    idx, d2d, _ = targets.match_to_catalog_sky(refs)
    return np.asarray(idx), d2d.arcsec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("morph", type=Path)
    ap.add_argument("--known", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    morph = json.loads(a.morph.read_text())
    rows = morph["targets"]
    m0, files, par = lc.load_model(MODEL)
    images, _ = lc.image_list(MODEL, files, par)
    lc.apply_frame_offset(MODEL, m0, images)
    img = SkyCoord(np.asarray(images["ra"], float), np.asarray(images["dec"], float), unit="deg")
    with open(a.known, newline="") as f:
        known = list(csv.DictReader(f))
    kc = SkyCoord([float(k["ra"]) for k in known], [float(k["dec"]) for k in known], unit="deg")
    tc = SkyCoord([float(r["ra"]) for r in rows], [float(r["dec"]) for r in rows], unit="deg")
    ii, di = nearest(tc, img)
    ik, dk = nearest(tc, kc)
    out = []
    for r, i, d, j, e in zip(rows, ii, di, ik, dk, strict=True):
        out.append(
            {
                "id": r["id"],
                "verdict": r["verdict"],
                "nearest_image": str(images["image_id"][i]),
                "dist_arcsec": round(float(d), 2),
                "nearest_known": known[j]["name"],
                "known_dist_arcsec": round(float(e), 3),
                "class": classify(r["verdict"], float(d), float(e)),
            }
        )
    doc = {
        "provenance": "derived",
        "source": f"{MODEL} arcs.dat (Mahler+22 ICLv2) with its JWST frame offset; {a.known.name}",
        "rule": f'PLAN.md (frozen): known <= {KNOWN_MATCH_ARCSEC}"; '
        f'vet = pass and <= {ON_IMAGE_ARCSEC}"',
        "rows": out,
    }
    a.out.write_text(json.dumps(doc, indent=1) + "\n")
    for o in out:
        print(o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

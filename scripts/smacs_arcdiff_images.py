"""Distance of each SMACS 0723 ``arc_diff`` target to the nearest Mahler+22 ICLv2 multiple image.

Applies the rule frozen in ``results/transients/smacs_arcdiff/PLAN.md`` (ASSUMPTION): a hit
that passes ``arc_diff_morph.py`` and lies within ``ON_IMAGE_ARCSEC`` of a catalogued image is
"on an image" and goes to ``/vet-candidate``, unless it matches a known control (SN Atalanta).
The image list is the ICLv2 ``arcs.dat`` shifted by the model's JWST frame offset
(``lens_consistency.shift_images``). The 1.0″ radius fails the known case (Atalanta lies 1.6-1.7″
from its catalogued images), so rows also carry ``posthoc_2arcsec`` (pass and <= 2.0″, post hoc).

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

from jwst_anomaly.lensmodel import load_lenstool_images  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

MODEL = "smacs0723-iclv2"  # Mahler+22 ICLv2, as smacs_s17_vet.MODEL
ON_IMAGE_ARCSEC = 1.0  # ASSUMPTION, frozen in PLAN.md before the run
KNOWN_MATCH_ARCSEC = 0.15  # ASSUMPTION: arc_diff's --match-arcsec default
POST_HOC_ARCSEC = 2.0  # post hoc: smallest round radius holding both Atalanta images


def classify(verdict: str, dist_image: float, dist_known: float) -> str:
    """Frozen rule: ``known`` (within KNOWN_MATCH of a control, any verdict), else ``vet``
    (morph pass AND within ON_IMAGE of a catalogued image), else ``off_image`` (pass, farther)
    or ``vetoed`` (not a pass)."""
    if dist_known <= KNOWN_MATCH_ARCSEC:
        return "known" if verdict == "pass" else "known_not_pass"
    if verdict != "pass":
        return "vetoed"
    return "vet" if dist_image <= ON_IMAGE_ARCSEC else "off_image"


def nearest(targets: SkyCoord, refs: SkyCoord) -> tuple[np.ndarray, np.ndarray]:
    """Index of and distance (arcsec) to the nearest ``refs`` entry for each target."""
    idx, d2d, _ = targets.match_to_catalog_sky(refs)
    return np.asarray(idx), d2d.arcsec


def _control_dists(targets: SkyCoord, controls: SkyCoord) -> np.ndarray:
    """Distance (arcsec) from each control to its nearest target: > KNOWN_MATCH = control missed."""
    return controls.match_to_catalog_sky(targets)[1].arcsec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("morph", type=Path)
    ap.add_argument("--known", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    morph = json.loads(a.morph.read_text())
    rows = morph["targets"]
    with open(a.known, newline="") as f:
        known = list(csv.DictReader(f))
    if not rows or not known:
        raise SystemExit("no morph targets or no known controls: nothing to classify")
    url, sha = lc.MODELS[MODEL]["files"]["arcs.dat"]
    images = lc.shift_images(MODEL, load_lenstool_images(fetch_catalog(url, sha)))
    img = SkyCoord(np.asarray(images["ra"], float), np.asarray(images["dec"], float), unit="deg")
    kc = SkyCoord([float(k["ra"]) for k in known], [float(k["dec"]) for k in known], unit="deg")
    tc = SkyCoord([float(r["ra"]) for r in rows], [float(r["dec"]) for r in rows], unit="deg")
    ii, di = nearest(tc, img)
    ik, dk = nearest(tc, kc)
    out = []
    for r, i, d, j, e in zip(rows, ii, di, ik, dk, strict=True):
        d, e = round(float(d), 3), round(float(e), 3)  # classify the stored values
        out.append(
            {
                "id": r["id"],
                "verdict": r["verdict"],
                "nearest_image": str(images["image_id"][i]),
                "dist_arcsec": d,
                "nearest_known": known[j]["name"],
                "known_dist_arcsec": e,
                "class": classify(r["verdict"], d, e),
                "posthoc_2arcsec": r["verdict"] == "pass" and d <= POST_HOC_ARCSEC,
            }
        )
    doc = {
        "provenance": "derived",
        "source": f"{MODEL} arcs.dat (Mahler+22 ICLv2) with its JWST frame offset; {a.known.name}",
        "rule": f'PLAN.md (frozen): known <= {KNOWN_MATCH_ARCSEC}"; '
        f'vet = pass and <= {ON_IMAGE_ARCSEC}"',
        "controls_recovered": {
            # rounded like the rows, so both outputs judge the same stored value
            k["name"]: bool(round(float(c), 3) <= KNOWN_MATCH_ARCSEC)
            for k, c in zip(known, _control_dists(tc, kc), strict=True)
        },
        "rows": out,
    }
    a.out.write_text(json.dumps(doc, indent=1) + "\n")
    for o in out:
        print(o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

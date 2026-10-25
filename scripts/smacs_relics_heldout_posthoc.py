"""Post hoc (not in the frozen plan): sign readings of the RELICS SMACS 0723 maps on the seen
HST-era families, and their per-family residuals, to tell a map-convention problem apart from a
model miss (as ``smacs_chow23``'s rule U, after the fact).

    python scripts/smacs_relics_heldout_posthoc.py

Writes ``posthoc.json`` next to ``summary.json``; the maps are deleted after use.
"""

from __future__ import annotations

import json
import warnings

import lens_consistency as lc
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import numpy as np
import smacs_relics_heldout as sr
import smacs_s17_vet as sv
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

READINGS = {
    "A": (1.0, 1.0),
    "A_neg_x": (-1.0, 1.0),
    "A_neg_y": (1.0, -1.0),
    "A_neg_xy": (-1.0, -1.0),
    # lenstool v1 |alpha| is ~1.25x glafic v2 at the seen images (session check): a 1 / 1.25 scale
    "A_x0.8": (0.8, 0.8),
}


def _maps(name: str):
    paths = sr.fetch_maps(name)
    arrays = []
    try:
        for p in paths:
            with fits.open(p) as h, warnings.catch_warnings():
                warnings.simplefilter("ignore")
                arrays.append(np.asarray(h[0].data, np.float64))
                wcs = WCS(h[0].header)
    finally:
        for p in paths:
            p.unlink(missing_ok=True)
    return arrays, wcs


def per_family(model, grid, t: Table) -> dict[str, float]:
    """rms per seen family (mean source of its seen images), as ``sr.seen_rms``."""
    out = {}
    for fam in dict.fromkeys(np.asarray(t["system"], str)[np.asarray(t["seen"], bool)]):
        sub = t[(np.asarray(t["system"], str) == fam)]
        out[fam] = sr.seen_rms(model, grid, sub)["rms_arcsec"]
    return out


def run() -> dict:
    images, seen = sr.heldout_images(), Table.read(sr.SEEN)
    t = mo.split_seen(mo.shift(images, mo.frame_offset(images, seen)), seen)
    centre = (float(np.mean(seen["ra"])), float(np.mean(seen["dec"])))
    h0, om = sv.RELICS_COSMO
    out: dict = {"label": "post hoc; not in the frozen plan", "models": {}}
    for name in sr.MODELS:
        (ax, ay), wcs = _maps(name)
        row = {}
        for rd, (sx, sy) in READINGS.items():
            m = lensmodel.MapLensModel(
                ax * sx,
                ay * sy,
                wcs,
                sv.RELICS_Z_LENS,
                FlatLambdaCDM(H0=h0, Om0=om),
                name,
                "",
                centre,
            )
            g = lensmodel.DeflectionGrid.compute(m, lc.grid_half_width(m, t), mh.GRID_STEP)
            row[rd] = {"seen": sr.seen_rms(m, g, t)}
            if rd == "A":
                row[rd]["per_family"] = per_family(m, g, t)
        out["models"][name] = row
    return out


def main() -> int:
    s = _clean(run())
    (sr.OUT / "posthoc.json").write_text(json.dumps(s, indent=1, allow_nan=False) + "\n")
    print(json.dumps(s, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

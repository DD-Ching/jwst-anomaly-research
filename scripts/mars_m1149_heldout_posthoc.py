"""Post hoc (not in the frozen plan): the Schuldt+2024 families under the CANUCS v1 MACS1149 maps.

The CANUCS model was fitted to a gold list that includes Schuldt et al. 2024's systems (readme),
so this is a training-set check, not a held-out test: it asks whether the listed families are
self-consistent under a model that saw them, to separate a MARS-specific miss from an
identification or redshift problem. Same rule 3 residuals (all images held out, type N) and
grid step as the frozen stage.

    python scripts/mars_m1149_heldout_posthoc.py
"""

from __future__ import annotations

import json
from pathlib import Path

import lens_consistency as lc
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import mars_m1149_heldout as m1
import numpy as np
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

MODEL = "macs1149-canucs-maps"


def main() -> int:
    model, files, _ = lc.load_model(MODEL)
    t = Table.read(m1.IMAGES)
    t = mo.split_seen(t, Table({"ra": [0.0], "dec": [0.0]}))  # every image held out (type N)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, t), mh.GRID_STEP)
    r = mo.heldout_residuals(model, grid, t)
    r["dtheta_arcsec"].info.format = ".3f"
    r.write(m1.OUT / "posthoc_canucs.ecsv", format="ascii.ecsv", overwrite=True)
    mars = Table.read(m1.OUT / "heldout_MACSJ1149.ecsv")
    both = np.isin(r["image_id"], mars["image_id"][np.asarray(mars["inside"], bool)])
    out = {
        "label": "post hoc; CANUCS v1 saw these systems (training-set check, not held-out)",
        "model": MODEL,
        "files": {k: Path(v).name for k, v in files.items()},
        "all_inside": mo.gate_h(r),
        "same_images_as_mars_test": mo.gate_h(r[both]),
        "over_3arcsec": [str(i) for i in r["image_id"][np.asarray(r["dtheta_arcsec"]) > 3.0]],
    }
    (m1.OUT / "posthoc.json").write_text(json.dumps(_clean(out), indent=1) + "\n")
    print(json.dumps(_clean(out), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

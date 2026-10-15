"""SMACS 0723 two-epoch grid pair (2022 ERO vs 2026 VENUS o057): one lensed source?

``transient_grid.py`` found two point sources, 3.2" apart on either side of the system 5 arc,
that are present only in 2026. This traces both positions through the Mahler+22 ICLv2 model over a
redshift grid. Per z it reports the source-plane separation of the two back-traced positions, their
model magnifications, and every image of their mean source with its distance to the observed pair.
A pair from one source point gives a small separation, and a third image is predicted elsewhere.
Positions are observed; everything traced is model_prediction.

    python scripts/smacs_transient_pair.py --out results/transients/smacs_grid/pair_trace.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

MODEL = "smacs0723-iclv2"
PAIR = {"A": (110.825136, -73.452178), "B": (110.823747, -73.452975)}  # grid g03415, g03423
Z_GRID = (1.0, 1.2, 1.3, 1.4, 1.425, 1.45, 1.5, 1.6, 2.0, 3.0)  # 1.425: system 5 (arcs.dat)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    model, files, par = lc.load_model(MODEL)
    images, _ = lc.image_list(MODEL, files, par)
    lc.apply_frame_offset(MODEL, model, images)
    grid = lensmodel.DeflectionGrid.compute(model, 45.0, 0.25)
    pair = Table(
        {
            "image_id": list(PAIR),
            "system": ["T", "T"],
            "ra": [v[0] for v in PAIR.values()],
            "dec": [v[1] for v in PAIR.values()],
            "z": [0.0, 0.0],
        }
    )
    obs = SkyCoord(pair["ra"] * u.deg, pair["dec"] * u.deg)
    rows = []
    for z in Z_GRID:
        bt = lensmodel.backtrace_images(model, pair, {"T": z})
        sep = float(np.hypot(*np.diff(bt["beta_x"]), *np.diff(bt["beta_y"])))
        im = lensmodel.find_images(
            model, grid, float(np.mean(bt["beta_x"])), float(np.mean(bt["beta_y"])), z
        )
        pred = SkyCoord(im["ra"] * u.deg, im["dec"] * u.deg)
        d = np.array([pred.separation(o).arcsec for o in obs])  # (2, n_images)
        rows.append(
            {
                "z": z,
                "source_sep_arcsec": round(sep, 3),
                "mu_observed_positions": [round(float(m), 1) for m in bt["magnification"]],
                "closest_predicted_arcsec": [round(float(x), 2) for x in d.min(axis=1)],
                "other_images": [
                    {
                        "ra": round(float(r["ra"]), 6),
                        "dec": round(float(r["dec"]), 6),
                        "mu": round(float(r["magnification"]), 1),
                    }
                    for r, dd in zip(im, d.min(axis=0), strict=True)
                    if dd > 0.5
                ],
            }
        )
    out = {
        "provenance": "model_prediction",
        "source": f"{MODEL} (lens_consistency.load_model, frame offset applied)",
        "pair_observed": PAIR,
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(rows, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Magnitude dependence of the ``arcs`` classes on a stored table (post hoc sensitivity check).

The ``arcs`` shape-noise null (``anti_expected``) models intrinsic shapes only (σ_ε per
component). For faint, small sources the measured axis also carries photometric noise, so an
``anti`` excess concentrated at the faint end points to measurement noise rather than to the
lens model. For the strong-shear background rows of ``results/arcs/<model>/arcs.ecsv`` this
script reports:

* ``rank_test``: the chance that all ``anti`` rows fall among the faintest ``k`` rows, where
  ``k`` is the rank of the brightest ``anti`` row counted from the faint end (hypergeometric,
  random class labels; ``derived``);
* ``class_stats`` and ``anti_expected`` above and below ``--mag-cut`` (an ASSUMPTION).

Output: ``results/arcs/<model>/faint.json``.

    python scripts/arcs_faint.py --model abells1063-cats --mag-cut 27.0
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from math import comb
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import paths

_spec = importlib.util.spec_from_file_location(
    "lens_consistency", Path(__file__).resolve().parent / "lens_consistency.py"
)
lc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lc)


def faint_rank_p(mag: np.ndarray, is_anti: np.ndarray) -> dict:
    """P that ``n_anti`` random rows all fall among the ``k`` faintest, ``k`` as observed."""
    order = np.argsort(-np.asarray(mag, float))  # faintest first
    anti_sorted = np.asarray(is_anti, bool)[order]
    n, a = len(anti_sorted), int(anti_sorted.sum())
    if a == 0:
        return {"n": n, "anti": 0, "k_faintest": 0, "p": 1.0}
    k = int(np.flatnonzero(anti_sorted).max()) + 1
    return {"n": n, "anti": a, "k_faintest": k, "p": comb(k, a) / comb(n, a)}


def run(args) -> dict:
    d = args.arcs_dir / args.model
    t = Table.read(d / "arcs.ecsv")
    s = json.loads((d / "arcs.json").read_text())["assumptions"]
    t = t[np.asarray(t["strong_shear"], bool) & (t["population"] == "background")]
    anti = np.asarray(t["orientation_class"] == "anti")
    out = {
        "model": args.model,
        "label": "post hoc sensitivity check (derived)",
        "rows": "strong-shear background rows of arcs.ecsv",
        "mag_cut": args.mag_cut,
        "rank_test": faint_rank_p(t["mag"], anti),
        "anti_mags": sorted(round(float(m), 2) for m in t["mag"][anti]),
    }
    for name, m in (("bright", t["mag"] < args.mag_cut), ("faint", t["mag"] >= args.mag_cut)):
        sub = t[np.asarray(m, bool)]
        out[name] = {
            "classes": lc.class_stats(sub, s["aligned_deg"]),
            "anti_expected": lc.anti_expected(
                sub, s["sigma_eps"], s["min_ellipticity"], s["anti_deg"]
            ),
        }
    (d / "faint.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", required=True)
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--mag-cut", type=float, default=27.0, help="ASSUMPTION: F200W AB")
    args = p.parse_args(argv)
    print(json.dumps(run(args), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

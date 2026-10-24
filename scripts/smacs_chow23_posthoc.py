"""Post hoc (not in the frozen plan) readings of the SMACS 0723 Chow+2023 maps after gate K failed.

1. Per system, the smallest backtrace rms over a source-z grid and its z (does a free z fix
   system 7?).
2. Gate K numbers without system 7 (mean source and WLS), as the S1063 post hoc did without
   system 10.

    python scripts/smacs_chow23_posthoc.py --maps DIR

Writes ``results/lens_consistency/smacs_chow23/posthoc.json``; the maps are deleted unless
``--keep``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import requests
import s1063_wslap as sw
import smacs_chow23 as sc
import spt0546_inputs as sp

Z_GRID = np.round(np.arange(0.5, 10.0001, 0.01), 4)
DROP = "7"


def run(maps: Path, out: Path = sc.OUT) -> dict:
    maps.mkdir(parents=True, exist_ok=True)
    shas = {m: sc.fetch(requests.Session(), m, maps / m) for m in sc.MAPS}
    frozen = json.loads((out / "summary.json").read_text())["files"]
    if any(shas[m] != frozen[m]["sha256"] for m in sc.MAPS):
        raise ValueError("maps differ from the frozen run")
    (ax, hx), (ay, _) = (sw._read(maps / m) for m in sc.MAPS)
    sha = hashlib.sha256((shas[sc.MAPS[0]] + shas[sc.MAPS[1]]).encode()).hexdigest()
    scale, sign = sc.READINGS["A"]
    model = sc.build_model(ax, ay, hx, scale, sign, sha, "smacs0723-chow23 (A)")
    images = sc.fixed_z_images()
    kept = images[np.asarray(images["system"], str) != DROP]
    k = sw.k_block(model, kept, sc.GRID_STEP)
    res = {
        "label": "post hoc (not in PLAN.md)",
        "z_grid": [float(Z_GRID[0]), float(Z_GRID[-1]), 0.01],
        "free_z_minima": sp.unit_minima(model, images, Z_GRID),
        "catalogue_z": {
            s: float(images["z"][images["system"] == s][0]) for s in sorted(set(images["system"]))
        },
        f"without_system_{DROP}": {"K": k["gate_K"], "K_wls": k["K_wls"]},
    }
    (out / "posthoc.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args(argv)
    try:
        print(json.dumps(run(a.maps), indent=1))
    finally:
        if not a.keep:
            for f in sc.MAPS:
                (a.maps / f).unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

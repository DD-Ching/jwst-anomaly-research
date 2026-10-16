"""Real morphology passes of 4 Abell 2744 tiles: stored peer cut vs ``arc_diff_extended.py``.

Reads ``passes.ecsv`` (the #261 tiling), the real primary-pair maps
``<scratch>/maps/real_<tile>.npz`` and the injection runs' confirming-pair maps
``<scratch>/maps/<tile>_o001.npz`` (passes within 0.5″ of an injected star are marked ``inj_near``
and left out of the counts), and writes ``<scratch>/real_eval.ecsv``. Every number is ``derived``;
thresholds are ASSUMPTIONs.

Usage: ``arc_diff_extended_eval.py <scratch>`` (CHANGELOG 2026-10-16).
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

R = Path(__file__).resolve().parents[1]
TILES = ("t-1+0", "t+0+0", "t+0-1", "t+1-1")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, R / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main(scratch: Path) -> None:
    ade, ch = _load("arc_diff_extended"), _load("arc_diff_inject_chain")
    passes = Table.read(R / "results/transients/a2744_grid_alt/arcdiff_tiles/passes.ecsv")
    out = []
    for t in TILES:
        p = passes[passes["tile"] == t]
        e = ade.run(scratch / f"maps/real_{t}.npz", p, "F356W", "hit")
        with np.load(scratch / f"maps/{t}_o001.npz") as m:
            s1 = m["F356W_s"]
            w1 = WCS(fits.Header.fromstring(str(m["F356W_wcs"])))
        with open(scratch / f"inj/{t}_o001/injected.csv", newline="") as f:
            inj = np.array([(float(r["ra"]), float(r["dec"])) for r in csv.DictReader(f)])
        for r, x in zip(p, e, strict=True):
            row = dict(zip(e.colnames, x, strict=True))
            old = bool(abs(r["z_peers"]) >= 3 or r["n_peers"] < 1)
            new = old if row["compact"] else ade.ext_pass(row)
            veto = ade.ext_pass_vetoed(row, old)
            c = ch.confirm_value(s1, w1, r["ra"], r["dec"])
            dra = (inj[:, 0] - r["ra"]) * np.cos(np.radians(r["dec"]))
            dmin = float(np.min(np.hypot(dra, inj[:, 1] - r["dec"]))) * 3600
            out.append(
                {
                    "tile": t,
                    "hit": r["hit"],
                    "ra": r["ra"],
                    "dec": r["dec"],
                    "s_hit": r["s_hit"],
                    "z_peers": r["z_peers"],
                    "compact": row["compact"],
                    "contrast": row["contrast"],
                    "z_ann": row["z_ann"],
                    "core_off_contrast": row["core_off_contrast"],
                    "line_frac": row["line_frac"],
                    "line_med": row["line_med"],
                    "line_theta": row["line_theta"],
                    "old": old,
                    "new": new,
                    "veto": veto,
                    "conf_s": c,
                    "conf": bool(np.sign(c) == np.sign(r["s_hit"]) and abs(c) >= 3),
                    "inj_near": dmin < 0.5,
                }
            )
    tab = Table(rows=out)
    tab.meta = {"provenance": "derived", "thresholds": "ASSUMPTION: see arc_diff_extended.py"}
    tab.write(scratch / "real_eval.ecsv", overwrite=True)
    use = ~tab["inj_near"]
    for k in ("old", "new", "veto"):
        sel = tab[k] & use
        print(k, "peer-stage", int(sel.sum()), "confirmed", int((sel & tab["conf"]).sum()))
    added = tab["new"] & ~tab["old"] & use
    print("n", len(tab), "compact", int(tab["compact"].sum()), "newly passing", int(added.sum()))
    added_v = tab["veto"] & ~tab["old"] & use
    print("newly passing with vetoes", int(added_v.sum()), list(tab["tile", "hit"][added_v]))


if __name__ == "__main__":
    main(Path(sys.argv[1]))

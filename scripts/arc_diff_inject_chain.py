"""Whole-chain efficiency of ``arc_diff.py --inject-grid`` stars for a two-pair tiling (2026-10-16).

The Abell 2744 UNCOVER vs ALT F356W tiling (CHANGELOG 2026-10-16) keeps a hit when it is a hit,
not ``near_bright_core``, passes ``arc_diff_morph.py``, passes the peer-population cut
(``arc_diff_peers.py``, z against the stored peer statistics) and has S ≥ ``--confirm-sigma`` with
the same sign in the confirming pair (UNCOVER→ALT o001). For each tile this expects
``<inject>/<tile>_<pair>/`` outputs of ``arc_diff.py --inject-grid`` run on both pairs with the same
grid (same epoch 1 and seed, so the same stars go into both second epochs) and their ``--maps``
files ``<maps>/<tile>_<pair>.npz``; it runs morphology and peers on the primary pair and counts the
injected stars surviving each stage, per magnitude. Injected stars brighten epoch 2, so "same sign"
means S > 0. Every number is ``derived`` from simulated stars in real maps; the thresholds are
ASSUMPTIONs.

Usage: ``arc_diff_inject_chain.py --inject DIR --maps DIR --ref-passes passes.ecsv
--tiles t-1+0 ...``
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parent
STAGES = ("n", "hit", "unflagged", "morph", "peers", "confirm")


def confirm_value(s: np.ndarray, w: WCS, ra: float, dec: float, half_px: int = 2) -> float:
    """Signed extreme of S in a (2 half_px + 1)² box at the position (nan outside the map)."""
    x, y = (int(round(float(v))) for v in w.all_world2pix(ra, dec, 0))
    win = s[max(y - half_px, 0) : y + half_px + 1, max(x - half_px, 0) : x + half_px + 1]
    if not win.size or not np.isfinite(win).any():
        return float("nan")
    return float(win.flat[np.nanargmax(np.abs(np.nan_to_num(win)))])


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--inject", type=Path, required=True)
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--ref-passes", type=Path, required=True)
    ap.add_argument("--tiles", nargs="+", required=True)
    ap.add_argument("--pair", default="o002", help="primary pair suffix")
    ap.add_argument("--confirm", default="o001", help="confirming pair suffix")
    ap.add_argument("--confirm-sigma", type=float, default=3.0)
    ap.add_argument("--band", default="F356W")
    a = ap.parse_args(argv)
    py = sys.executable
    out = {
        "provenance": "derived (simulated stars in real maps)",
        "thresholds_assumption": {"confirm_sigma": a.confirm_sigma, "peer_min_z": 3.0},
        "tiles": {},
    }
    total: dict[str, dict[str, int]] = {}
    for t in a.tiles:
        d = a.inject / f"{t}_{a.pair}"
        m2 = a.maps / f"{t}_{a.pair}.npz"
        morph = d / "morph.json"
        subprocess.run(
            [py, str(_DIR / "arc_diff_morph.py"), "--hits", str(d / "hits.ecsv"), "--maps", str(m2)]
            + ["--min-rosette-bands", "1", "--min-dipole-bands", "1", "--min-ring-bands", "1"]
            + ["--out", str(morph)],
            check=True,
        )
        peers = d / "peers.ecsv"
        subprocess.run(
            [py, str(_DIR / "arc_diff_peers.py"), "--maps", str(m2), "--hits", str(d / "hits.ecsv")]
            + ["--tile", t, "--ref-passes", str(a.ref_passes), "--out", str(peers)],
            check=True,
        )
        verdict = {x["id"]: x["verdict"] for x in json.loads(morph.read_text())["targets"]}
        pt = {r["hit"]: r for r in Table.read(peers)}
        with np.load(a.maps / f"{t}_{a.confirm}.npz") as m:
            s1 = m[f"{a.band}_s"]
            w1 = WCS(fits.Header.fromstring(str(m[f"{a.band}_wcs"])))
        per: dict[str, dict] = {}
        with open(d / "injected.csv", newline="") as f:
            for r in csv.DictReader(f):
                c = per.setdefault(r["mag"], {k: 0 for k in STAGES} | {"z_stored": []})
                c["n"] += 1
                h = r["hit"]
                if not h:
                    continue
                c["hit"] += 1
                if r["near_bright_core"] in ("True", "true", "1"):
                    continue
                c["unflagged"] += 1
                if verdict.get(h) != "pass":
                    continue
                c["morph"] += 1
                c["z_stored"].append(round(float(pt[h]["z_stored"]), 1))
                if not bool(pt[h]["stored_pass"]):
                    continue
                c["peers"] += 1
                if confirm_value(s1, w1, float(r["ra"]), float(r["dec"])) >= a.confirm_sigma:
                    c["confirm"] += 1
        out["tiles"][t] = per
        for mag, c in per.items():
            tot = total.setdefault(mag, {k: 0 for k in STAGES})
            for k in STAGES:
                tot[k] += c[k]
    for c in total.values():
        c["eff_chain"] = round(c["confirm"] / c["n"], 3) if c["n"] else None
    out["total"] = total
    text = json.dumps(out, indent=1) + "\n"
    (a.inject / "summary.json").write_text(text)
    print(json.dumps(total, indent=1))


if __name__ == "__main__":
    main()

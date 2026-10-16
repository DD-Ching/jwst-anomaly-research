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

``--extended`` adds the alternative stage of ``arc_diff_extended.py``: the peer cut only for hits
at a compact epoch-1 peak, the local-annulus cut |z_ann| ≥ ``--ext-min-z`` otherwise (stages
``peers_ext`` and ``confirm_ext``). ``--vetoes`` (with ``--extended``) adds its second version with
the line-coherence and offset-core vetoes (stages ``peers_veto`` and ``confirm_veto``).

Usage: ``arc_diff_inject_chain.py --inject DIR --maps DIR --ref-passes passes.ecsv
--tiles t-1+0 ... [--extended]``
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("arc_diff_peers", _DIR / "arc_diff_peers.py")
adp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(adp)
_spec = importlib.util.spec_from_file_location("arc_diff_extended", _DIR / "arc_diff_extended.py")
ade = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ade)
STAGES = (
    "n",
    "hit",
    "unflagged",
    "morph",
    "peers",
    "confirm",
    "peers_ext",
    "confirm_ext",
    "peers_veto",
    "confirm_veto",
)


def confirm_value(s: np.ndarray, w: WCS, ra: float, dec: float, half_px: int = 2) -> float:
    """Signed extreme of S in a (2 half_px + 1)² box at the position (nan outside the map)."""
    x, y = (int(round(float(v))) for v in w.all_world2pix(ra, dec, 0))
    return adp.signed_extreme(adp.window(s, y, x, half_px))


def run_stages(py: str, t: str, d: Path, m2: Path, ref_passes: Path) -> None:
    """Morphology and peers for one tile (subprocesses; tiles run concurrently)."""
    subprocess.run(
        [py, str(_DIR / "arc_diff_morph.py"), "--hits", str(d / "hits.ecsv"), "--maps", str(m2)]
        + ["--min-rosette-bands", "1", "--min-dipole-bands", "1", "--min-ring-bands", "1"]
        + ["--out", str(d / "morph.json")],
        check=True,
    )
    # --all-hits: the chain applies near_bright_core itself, from injected.csv.
    subprocess.run(
        [py, str(_DIR / "arc_diff_peers.py"), "--maps", str(m2), "--hits", str(d / "hits.ecsv")]
        + ["--all-hits", "--tile", t, "--ref-passes", str(ref_passes)]
        + ["--out", str(d / "peers.ecsv")],
        check=True,
    )


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
    ap.add_argument("--procs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--extended", action="store_true", help="add the arc_diff_extended.py stage")
    ap.add_argument("--ext-min-z", type=float, default=ade.DEFAULTS["min_z"])
    ap.add_argument("--vetoes", action="store_true", help="with --extended: add the vetoed version")
    ap.add_argument("--skip-stages", action="store_true", help="reuse existing morph/peers outputs")
    a = ap.parse_args(argv)
    py = sys.executable
    out = {
        "provenance": "derived (simulated stars in real maps)",
        "thresholds_assumption": {"confirm_sigma": a.confirm_sigma, "peer_min_z": 3.0}
        | ({"extended": ade.DEFAULTS | {"min_z": a.ext_min_z}} if a.extended else {})
        | ({"vetoes": ade.VETO} if a.extended and a.vetoes else {}),
        "tiles": {},
    }
    total: dict[str, dict[str, int]] = {}
    with ThreadPoolExecutor(max_workers=a.procs) as ex:
        futs = (
            []
            if a.skip_stages
            else [
                ex.submit(
                    run_stages,
                    py,
                    t,
                    a.inject / f"{t}_{a.pair}",
                    a.maps / f"{t}_{a.pair}.npz",
                    a.ref_passes,
                )
                for t in a.tiles
            ]
        )
        for f in futs:
            f.result()
    for t in a.tiles:
        d = a.inject / f"{t}_{a.pair}"
        morph = d / "morph.json"
        peers = d / "peers.ecsv"
        verdict = {x["id"]: x["verdict"] for x in json.loads(morph.read_text())["targets"]}
        pt = {r["hit"]: r for r in Table.read(peers)}
        ext = {}
        if a.extended:
            et = ade.run(
                a.maps / f"{t}_{a.pair}.npz", Table.read(d / "hits.ecsv"), a.band, "hit_id"
            )
            ext = {r["hit_id"]: dict(zip(et.colnames, r, strict=True)) for r in et}
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
                c["z_stored"].append(float(pt[h]["z_stored"]))
                confirmed = (
                    confirm_value(s1, w1, float(r["ra"]), float(r["dec"])) >= a.confirm_sigma
                )
                if bool(pt[h]["stored_pass"]):
                    c["peers"] += 1
                    c["confirm"] += confirmed
                if a.extended:
                    e = ext[h]
                    stored = bool(pt[h]["stored_pass"])
                    ok = stored if e["compact"] else ade.ext_pass(e, a.ext_min_z)
                    c["peers_ext"] += ok
                    c["confirm_ext"] += ok and confirmed
                    if a.vetoes:
                        okv = ade.ext_pass_vetoed(e, stored, a.ext_min_z)
                        c["peers_veto"] += okv
                        c["confirm_veto"] += okv and confirmed
        out["tiles"][t] = per
        for mag, c in per.items():
            tot = total.setdefault(mag, {k: 0 for k in STAGES})
            for k in STAGES:
                tot[k] += c[k]
    for c in total.values():
        c["eff_chain"] = round(c["confirm"] / c["n"], 3) if c["n"] else None
        if a.extended:
            c["eff_chain_ext"] = round(c["confirm_ext"] / c["n"], 3) if c["n"] else None
        if a.extended and a.vetoes:
            c["eff_chain_veto"] = round(c["confirm_veto"] / c["n"], 3) if c["n"] else None
    out["total"] = total
    text = json.dumps(out, indent=1) + "\n"
    (a.inject / "summary.json").write_text(text)
    print(json.dumps(total, indent=1))


if __name__ == "__main__":
    main()

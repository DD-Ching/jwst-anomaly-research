"""Choose boxes for a bounded two-epoch batch search from pair metadata only (no differences).

``results/transients/batch_venus/PLAN.md`` fixes the rule. For every VENUS (6882) target with
another public NIRCam programme (``venus_pairs.json``), MAST CAOM level-3 F150W and F444W
observations of VENUS and of each other programme are listed. A pair = (VENUS obs, other obs)
whose F150W and F444W footprints both contain the target centre. Candidate boxes are a 3 x 3
grid of 24″ boxes centred on the target (the cluster core, where lensed images concentrate); a
box is kept if its corners and centre lie in all four footprints and it does not overlap an
already searched box. Pairs are ranked by expected sensitivity (the shallower F150W exposure of
the two epochs, longest first; one pair per target; equal depth within a target: the lowest
(other, VENUS) obs_id pair; across targets: input order), boxes by distance from the centre
(ties: lowest RA, then lowest Dec); up to ``--per-pair`` boxes per pair are taken in rank order
until ``--max-boxes``. Each kept box is then quality-checked (``read_box`` finite fraction >=
``--min-finite`` in all four images; header DATE-OBS, CAL_VER, CRDS context) and replaced by the
next candidate of the same rule if it fails.

    python scripts/batch_pairs.py --pairs venus_pairs.json --exclude RA DEC [...] --out boxes.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MPath

sys.path.insert(0, str(Path(__file__).resolve().parent))

BOX = 24.0  # arcsec
TOUCH_TOL = (
    0.01  # arcsec: boxes that only touch (centres BOX apart, rounded to 1e-6 deg) do not overlap
)


def overlaps(ra: float, dec: float, others: list[tuple[float, float]]) -> bool:
    """True if the BOX-wide box at (ra, dec) shares area with one centred on any of ``others``."""
    c = math.cos(math.radians(dec))
    return any(
        abs((ra - x) * c * 3600) < BOX - TOUCH_TOL and abs((dec - y) * 3600) < BOX - TOUCH_TOL
        for x, y in others
    )


def polygons(s_region: str) -> list[MPath]:
    out = []
    for chunk in str(s_region).split("POLYGON")[1:]:
        v = [float(x) for x in chunk.replace("ICRS", "").split()]
        out.append(MPath(np.array(v).reshape(-1, 2)))
    return out


def inside(polys: list[MPath], ra: float, dec: float) -> bool:
    return any(p.contains_point((ra, dec)) for p in polys)


def box_points(ra: float, dec: float, half: float = BOX / 2) -> list[tuple[float, float]]:
    c = math.cos(math.radians(dec))
    pts = [(ra, dec)]
    for sx in (-1, 1):
        for sy in (-1, 1):
            pts.append((ra + sx * half / 3600 / c, dec + sy * half / 3600))
    return pts


def uri(obs_id: str, pid: str) -> str:
    o = obs_id.split("-")[1].split("_")[0]
    return f"s3://stpubdata/jwst/public/jw{int(pid):05d}/L3/t/{o}/{obs_id}_i2d.fits"


def list_obs(pid: str, ra: float, dec: float) -> dict:
    from astroquery.mast import Observations

    t = Observations.query_criteria(
        obs_collection="JWST",
        instrument_name="NIRCAM/IMAGE",
        calib_level=3,
        proposal_id=pid,
        filters=["F150W", "F444W"],
        s_ra=[ra - 0.05, ra + 0.05],
        s_dec=[dec - 0.05, dec + 0.05],
        dataRights="PUBLIC",
    )
    groups: dict = {}
    for r in t:
        oid = str(r["obs_id"])
        if not oid.endswith(("-f150w", "-f444w")) or "clear" not in oid:
            continue
        key = oid.rsplit("_", 1)[0]  # jwPPPPP-oNNN_tNNN_nircam
        band = oid.rsplit("-", 1)[1].upper()
        groups.setdefault(key, {})[band] = {
            "obs_id": oid,
            "uri": uri(oid, pid),
            "mjd": float(r["t_min"]),
            "exptime": float(r["t_exptime"]),
            "polys": polygons(r["s_region"]),
        }
    return {k: v for k, v in groups.items() if {"F150W", "F444W"} <= set(v)}


def extra_epochs(p: dict, ra: float, dec: float) -> int:
    """Public F150W visits (other than the pair's two) whose footprint contains the box centre."""
    used = {p["e1"][0], p["e2"][0]}
    return sum(
        1
        for k, v in p.get("groups", {}).items()
        if k not in used and inside(v["F150W"]["polys"], ra, dec)
    )


def rank_key(depth: float, extra: int, rank: int, dist: float) -> tuple:
    """``--prefer-extra-epochs`` order: depth tier (factor-2 bins of the shallower F150W exposure,
    deeper first), then more extra epochs, then pair rank and distance from the centre."""
    return (-math.floor(math.log2(depth)), -extra, rank, dist)


def header_info(u: str) -> dict:
    import fsspec
    from astropy.io import fits

    fs, path = fsspec.core.url_to_fs(u, anon=True)
    with fs.open(path, "rb", block_size=2**20) as f, fits.open(f, lazy_load_hdus=True) as h:
        p, sci = h[0].header, h["SCI"].header
        out = {k: p.get(k) for k in ("DATE-OBS", "CAL_VER", "CRDS_CTX")}
        out["PA_V3"] = sci.get("PA_V3", p.get("PA_V3"))  # level-3: in the SCI extension
        return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pairs", type=Path, required=True)
    ap.add_argument("--exclude", type=float, nargs=2, action="append", default=[])
    ap.add_argument("--max-boxes", type=int, default=20)
    ap.add_argument("--per-pair", type=int, default=2)
    ap.add_argument("--min-finite", type=float, default=0.95)
    ap.add_argument(
        "--prefer-extra-epochs",
        action="store_true",
        help="among boxes of similar depth, prefer those with more extra public F150W epochs",
    )
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    from arc_diff import read_box

    targets = json.loads(a.pairs.read_text())["targets"]
    pairs, log = [], []
    for tg in targets:
        if not tg["other_programs_F150W_L3"]:
            continue
        ra, dec = tg["ra"], tg["dec"]
        venus_all = list_obs("6882", ra, dec)
        venus = {
            k: v for k, v in venus_all.items() if all(inside(v[b]["polys"], ra, dec) for b in v)
        }
        groups = dict(venus_all)  # every F150W+F444W visit of the listed programmes (extra epochs)
        best = None
        for pid in tg["other_programs_F150W_L3"]:
            obs = list_obs(pid, ra, dec)
            groups.update(obs)
            for k, v in obs.items():
                if not all(inside(v[b]["polys"], ra, dec) for b in v):
                    continue
                for vk, vv in venus.items():
                    depth = min(v["F150W"]["exptime"], vv["F150W"]["exptime"])
                    cand = {
                        "target": tg["target"],
                        "ra": ra,
                        "dec": dec,
                        "e1": (k, v),
                        "e2": (vk, vv),
                        "depth": depth,
                    }
                    # equal depth: lowest (other, VENUS) obs_id pair, not MAST row order
                    tie = (k, vk) < (best["e1"][0], best["e2"][0]) if best else True
                    if best is None or depth > best["depth"] or (depth == best["depth"] and tie):
                        best = cand
        log.append(
            {"target": tg["target"], "pair": [best["e1"][0], best["e2"][0]] if best else None}
        )
        if best:
            best["groups"] = groups
            pairs.append(best)
    pairs.sort(key=lambda p: -p["depth"])
    seq = []  # (rank, pair, distance, ra, dec, extra epochs) in visiting order
    for rank, p in enumerate(pairs):
        c = math.cos(math.radians(p["dec"]))
        grid = []
        for ix in (-1, 0, 1):
            for iy in (-1, 0, 1):
                bra, bdec = p["ra"] + ix * BOX / 3600 / c, p["dec"] + iy * BOX / 3600
                grid.append((math.hypot(ix, iy), bra, bdec))
        grid.sort()
        for dist, bra, bdec in grid:
            seq.append((rank, p, dist, bra, bdec, extra_epochs(p, bra, bdec)))
    if a.prefer_extra_epochs:
        seq.sort(key=lambda t: rank_key(t[1]["depth"], t[5], t[0], t[2]))
    chosen, rejected = [], []
    n_pair: dict = {}
    for rank, p, dist, bra, bdec, extra in seq:
        if len(chosen) >= a.max_boxes:
            break
        if n_pair.get(rank, 0) >= a.per_pair:
            continue
        e1, e2 = p["e1"][1], p["e2"][1]
        if overlaps(bra, bdec, a.exclude):
            rejected.append(
                {"target": p["target"], "ra": bra, "dec": bdec, "why": "overlaps searched box"}
            )
            continue
        pts = box_points(bra, bdec)
        if not all(
            inside(e[b]["polys"], x, y)
            for e in (e1, e2)
            for b in ("F150W", "F444W")
            for x, y in pts
        ):
            rejected.append(
                {
                    "target": p["target"],
                    "ra": bra,
                    "dec": bdec,
                    "why": "not inside all footprints",
                }
            )
            continue
        fin = {}
        for e, tag in ((e1, "e1"), (e2, "e2")):
            for b in ("F150W", "F444W"):
                try:
                    bx = read_box(e[b]["uri"], bra, bdec, BOX / 2)
                    full = (2 * math.ceil(BOX / 2 / bx["scale"]) + 1) ** 2
                    fin[f"{tag}_{b}"] = round(float(np.isfinite(bx["sci"]).sum()) / full, 4)
                except (ValueError, OSError) as err:
                    fin[f"{tag}_{b}"] = 0.0
                    fin[f"{tag}_{b}_error"] = str(err)[:120]
        if min(v for k, v in fin.items() if not k.endswith("error")) < a.min_finite:
            rejected.append(
                {
                    "target": p["target"],
                    "ra": bra,
                    "dec": bdec,
                    "why": "finite fraction",
                    "finite": fin,
                }
            )
            continue
        k_in = sum(x["target"] == p["target"] for x in chosen) + 1
        chosen.append(
            {
                "box_id": f"{p['target']}_{k_in}",
                "target": p["target"],
                "ra": round(bra, 6),
                "dec": round(bdec, 6),
                "pair_rank": rank,
                "dist_from_centre_boxes": round(dist, 3),  # in box widths
                "extra_f150w_epochs": extra,
                "finite": fin,
                "e1": {
                    b: {k: e1[b][k] for k in ("obs_id", "uri", "mjd", "exptime")}
                    for b in ("F150W", "F444W")
                },
                "e2": {
                    b: {k: e2[b][k] for k in ("obs_id", "uri", "mjd", "exptime")}
                    for b in ("F150W", "F444W")
                },
            }
        )
        n_pair[rank] = n_pair.get(rank, 0) + 1
    for bx in chosen:
        for tag in ("e1", "e2"):
            for b in ("F150W", "F444W"):
                try:
                    bx[tag][b]["header"] = header_info(bx[tag][b]["uri"])
                except Exception as err:  # noqa: BLE001 (metadata only; recorded)
                    bx[tag][b]["header"] = {"error": str(err)[:120]}
    out = {
        "label": "observed metadata + derived box choice (no difference images used)",
        "rule": vars(a) | {"pairs": str(a.pairs), "out": str(a.out)},
        "pairs_considered": log,
        "chosen": chosen,
        "rejected": rejected,
        "effective_area_arcmin2": round(len(chosen) * (BOX / 60) ** 2, 3),
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(len(chosen), "boxes;", sum(1 for x in log if x["pair"]), "pairs")
    return 0


if __name__ == "__main__":
    sys.exit(main())

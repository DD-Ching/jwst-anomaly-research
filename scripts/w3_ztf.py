"""W3 (negative-mass point lens, n = 1, eps < 0) feasibility trial on ZTF DR24 light curves.

Plan: ``results/w3_ztf/PLAN.md`` (fixed before any light curve was read). Data: IRSA's ZTF DR24 HATS
collection on S3 (``s3://ipac-irsa-ztf/ztf/enhanced/dr24/``, anonymous), read column-wise with
pyarrow; nothing is stored except the small tracked tables under ``results/w3_ztf/`` and a
session cache.

Stages (``python scripts/w3_ztf.py <stage>``):

- ``sample``: metadata-only star selection, groups and the development/holdout split;
- ``fetch``: the light-curve partition's rows for the sample and for the comparison stars (cache);
- ``coverage``: good epochs per (field, rcid, filter) and per-epoch limits from comparison stars;
- ``inject``: W3 and ordinary-control injections on development carriers, through pre-screen and
  vetting;
- ``search``: pre-screen + vetting on one half (``--half dev|holdout``).

Every number is ``derived`` (``simulated`` for injections); thresholds are ASSUMPTIONs in the plan.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import zlib
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

BUCKET = "ipac-irsa-ztf/ztf/enhanced/dr24"
LC_PART = f"{BUCKET}/lc/hats/ztf_dr24_lc-hats/dataset/Norder=5/Dir=0/Npix=3839/part0.snappy.parquet"
OBJ_PART = (
    f"{BUCKET}/objects/hats/ztf_dr24_objects-hats/dataset/Norder=4/Dir=0/Npix=959"
    "/part0.snappy.parquet"
)
LC_ORDER, LC_PIX = 5, 3839
OUT = Path(__file__).resolve().parents[1] / "results" / "w3_ztf"
CACHE = Path(os.environ.get("JWST_ANOMALY_DATA", "data")) / "cache" / "w3_ztf"

# plan section 4
MAG_RANGE = (14.5, 17.5)
SHARP_MAX = 0.3
MIN_GOODOBS = 300
ISO_RADIUS = 5.0  # arcsec
ISO_DMAG = 2.5
GROUP_RADIUS = 1.5  # arcsec
N_MAX = 1000
SEED = 7
COMP_MAG = (14.0, 18.0)
# plan sections 3, 5, 6
GOOD_FRAC = 0.5
EPOCH_TOL = 0.0005  # d
DIM_MAG = 1.0
MIN_RUN = 3
MIN_NIGHTS = 3
MIN_BRACKET = 20
COMMON_FRAC = 0.2
INTER = (0.3, 1.0)
MAX_INTER = 1


def s3():
    import s3fs

    return s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-west-2"})


def holdout(key: int) -> bool:
    return zlib.crc32(str(int(key)).encode()) % 2 == 0


def in_partition(ra, dec) -> np.ndarray:
    import astropy.units as u
    from astropy_healpix import HEALPix

    hp = HEALPix(nside=2**LC_ORDER, order="nested")
    return hp.lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg) == LC_PIX


# ----------------------------------------------------------------------------------------- sample


def select(obj) -> tuple[list[dict], dict]:
    """Groups of the sample from an objects table (pandas), plan section 4."""
    from scipy.spatial import cKDTree

    obj = obj[in_partition(obj["ra"].to_numpy(), obj["dec"].to_numpy())].reset_index(drop=True)
    dec0 = float(np.median(obj["dec"]))
    xy = np.column_stack(
        [obj["ra"].to_numpy() * np.cos(np.radians(dec0)) * 3600, obj["dec"].to_numpy() * 3600]
    )
    tree = cKDTree(xy)
    ok = (
        (obj["fid"] == 2)
        & obj["refmag"].between(*MAG_RANGE)
        & (obj["refsharp"].abs() <= SHARP_MAX)
        & (obj["transient"] == 0)
        & (obj["ngoodobs"] >= MIN_GOODOBS)
    ).to_numpy()
    refmag = obj["refmag"].to_numpy()
    cand = np.flatnonzero(ok)
    iso = []
    for i in (
        cand
    ):  # isolated: no other star (beyond GROUP_RADIUS) within ISO_RADIUS brighter than m + ISO_DMAG
        nb = tree.query_ball_point(xy[i], ISO_RADIUS)
        if not any(
            np.hypot(*(xy[j] - xy[i])) > GROUP_RADIUS
            and np.isfinite(refmag[j])
            and refmag[j] < refmag[i] + ISO_DMAG
            for j in nb
        ):
            iso.append(i)
    iso = sorted(iso, key=lambda i: -int(obj["ngoodobs"][i]))
    used = set()
    groups = []
    for i in iso:
        if i in used:
            continue
        mem = [j for j in tree.query_ball_point(xy[i], GROUP_RADIUS) if j not in used]
        used.update(mem)
        groups.append(
            {
                "key": int(obj["oid"][i]),
                "ra": float(obj["ra"][i]),
                "dec": float(obj["dec"][i]),
                "refmag_r": float(refmag[i]),
                "members": [int(obj["oid"][j]) for j in mem],
            }
        )
    n_all = len(groups)
    rng = np.random.default_rng(SEED)
    pick = rng.permutation(n_all)[:N_MAX]
    groups = [groups[k] for k in sorted(pick)]
    for g in groups:
        g["holdout"] = holdout(g["key"])
    info = {
        "objects_in_area": int(len(obj)),
        "r_band_passing_cuts": int(ok.sum()),
        "isolated": len(iso),
        "groups_before_draw": int(n_all),
        "groups": len(groups),
        "holdout": int(sum(g["holdout"] for g in groups)),
    }
    return groups, info


def stage_sample(_a) -> None:
    import pyarrow.parquet as pq

    cols = ["oid", "ra", "dec", "field", "ccdid", "qid", "fid", "refmag", "refsharp"]
    cols += ["ngoodobs", "transient"]
    t0 = time.time()
    obj = pq.ParquetFile(s3().open(OBJ_PART)).read(columns=cols).to_pandas()
    groups, info = select(obj)
    info["read_s"] = round(time.time() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "sample.json").write_text(
        json.dumps({"label": "derived; plan section 4", "info": info, "groups": groups}) + "\n"
    )
    # comparison stars for coverage (all objects in the area in the comparison magnitude range)
    comp = obj[in_partition(obj["ra"].to_numpy(), obj["dec"].to_numpy())]
    comp = comp[comp["refmag"].between(*COMP_MAG) & (comp["transient"] == 0)]
    CACHE.mkdir(parents=True, exist_ok=True)
    np.save(CACHE / "comparison_oids.npy", comp["oid"].to_numpy())
    print(json.dumps(info), len(comp))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=["sample"])
    a = ap.parse_args(argv)
    {"sample": stage_sample}[a.stage](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

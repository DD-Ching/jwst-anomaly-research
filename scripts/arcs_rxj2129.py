"""`arcs` on RX J2129.7+0005: a further independent test of the Abell S1063 anti excess.

Plan (frozen before the run): ``results/arcs/rxj2129/PLAN.md``. This script stages the pinned map
members of the Zenodo 21460864 zip (HTTP range reads, no full download) into the cache, runs
``lens_consistency.py validate`` and ``arcs`` (D-024 defaults) with the Lenstool (primary) and
GLAFIC (second) maps on the program 2767 F200W catalogue, applies the frozen outcome rule
(``arcs_macs0647.outcome``) and lists the anti rows of the deciding run with the other model's
class. Shapes are ``observed``, predicted axes ``model_prediction``, classes, counts and P values
``derived``; thresholds are ASSUMPTIONs.

    python scripts/arcs_rxj2129.py --catalog <jw02767-o002_t001 F200W _cat.ecsv> \
        --photoz <rxj2129-grizli-v7.0-fix.eazypy.zout.fits>
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import lensmodel, paths
from jwst_anomaly.photometry import sha256_file

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402
from arcs_macs0647 import (  # noqa: E402
    _arcs,
    _background_by_gmin,
    _relative_paths,
    match_class,
    outcome,
    validate_passes,
)

PRIMARY = "rxj2129-lenstool-maps"
SECOND = "rxj2129-glafic-maps"
RETRIES = 5  # range-read attempts per block, with 1-16 s back-off
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "rxj2129"


class _RangeFile(io.RawIOBase):
    """A seekable read-only view of a remote file through HTTP range requests."""

    def __init__(self, url: str):
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:  # noqa: S310
            self.url, self.size = r.url, int(r.headers["Content-Length"])
        self.pos = 0

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def seek(self, offset: int, whence: int = 0) -> int:
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def tell(self) -> int:
        return self.pos

    def readinto(self, buf) -> int:
        if self.pos >= self.size:
            return 0
        end = min(self.pos + len(buf), self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        for attempt in range(RETRIES):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
                    if r.status != 206:  # a 200 would stream the whole archive
                        raise ValueError(f"{self.url}: range request answered {r.status}")
                    data = r.read(len(buf) + 1)
                break
            except OSError:  # URLError and dropped TLS connections; retried, never skipped
                if attempt == RETRIES - 1:
                    raise
                time.sleep(2**attempt)
        if not 0 < len(data) <= len(buf):
            raise ValueError(f"{self.url}: range read returned {len(data)} bytes for {len(buf)}")
        buf[: len(data)] = data
        self.pos += len(data)
        return len(data)


def stage_members(files: dict[str, tuple[str, str]]) -> None:
    """Put each pinned zip member (``<zip url>#<member>``) where ``fetch_catalog`` looks for it,
    verified against its sha256; members already cached are left alone."""
    cache = paths.cache_dir() / "external"
    cache.mkdir(parents=True, exist_ok=True)
    archive = None
    for url, sha in files.values():
        zip_url, member = url.split("#", 1)
        target = cache / f"{sha[:12]}_{member.rsplit('/', 1)[-1]}"
        if target.exists() and sha256_file(target) == sha:
            continue
        if archive is None:
            archive = zipfile.ZipFile(io.BufferedReader(_RangeFile(zip_url), buffer_size=4 << 20))
        tmp = target.with_suffix(".part")
        with archive.open(member) as src, open(tmp, "wb") as out:
            while chunk := src.read(1 << 20):
                out.write(chunk)
        if sha256_file(tmp) != sha:
            tmp.unlink()
            raise ValueError(f"{member}: sha256 mismatch")
        tmp.replace(target)


def _run(model: str, out: Path, catalog: Path, photoz: Path) -> dict:
    """``validate`` then ``arcs``; a model whose validate fails gets no ``arcs`` run and an empty
    control (n = 0 fails it), as PLAN.md says."""
    lc.main(["--model", model, "--out", str(out), "validate"])
    _relative_paths(out / model / "validate.json")
    if not validate_passes(json.loads((out / model / "validate.json").read_text())):
        return {"model": model, "validate_passes": False, "strong_shear_background": {"n": 0}}
    run = _arcs(model, out, catalog, photoz)
    return {"model": model, "validate_passes": True, **run}


def decide(primary: str, second: str, out: Path, catalog: Path, photoz: Path, plan: str) -> dict:
    """Run both models, apply the frozen outcome rule and list the deciding run's anti rows with
    the other model's class; writes ``out/summary.json`` and returns it."""
    # dict.fromkeys: a one-model field passes the same name twice and is run once
    runs = {m: _run(m, out, catalog, photoz) for m in dict.fromkeys((primary, second))}
    res = outcome(runs[primary], runs[second])
    summary = {
        "plan": plan,
        "inputs": {
            "catalog": catalog.name,
            "catalog_sha256": sha256_file(catalog),
            "photoz": photoz.name,
            "photoz_sha256": sha256_file(photoz),
        },
        "runs": {
            m: {
                "validate_passes": r["validate_passes"],
                "strong_shear_background": r["strong_shear_background"],
                "anti_expected_shape_noise_background": r.get(
                    "anti_expected_shape_noise_background"
                ),
            }
            for m, r in runs.items()
        },
        **res,
    }
    if res["deciding_model"] is not None:
        dec = res["deciding_model"]
        other = second if dec == primary else primary
        t = Table.read(out / dec / "arcs.ecsv")
        summary["background_by_gmin"] = _background_by_gmin(t)
        anti = t[
            t["strong_shear"]
            & (t["population"] == "background")
            & (t["orientation_class"] == "anti")
        ]
        if other == dec:  # one-model field: there is no second run to compare with
            other_cls = ["no second model"] * len(anti)
        elif runs[other]["validate_passes"]:
            other_cls = match_class(anti, Table.read(out / other / "arcs.ecsv"), None)
        else:
            other_cls = ["no run"] * len(anti)
        summary["anti_background"] = [
            {
                "label": int(r["label"]),
                "ra": round(float(r["ra"]), 6),
                "dec": round(float(r["dec"]), 6),
                "mag": round(float(r["mag"]), 2),
                "g_min": round(float(r["g_min"]), 3),
                "z_phot": round(float(r["z_phot"]), 2),
                "ellipticity": round(float(r["ellipticity"]), 3),
                f"class_{other}": c,
            }
            for r, c in zip(anti, other_cls, strict=True)
        ]
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("deciding_model", "outcome", "p_ge_observed")}))
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--photoz", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    stage_members(lensmodel.RXJ2129_LENSTOOL_MAPS)
    stage_members(lensmodel.RXJ2129_GLAFIC_MAPS)
    decide(PRIMARY, SECOND, args.out, args.catalog, args.photoz, "results/arcs/rxj2129/PLAN.md")
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())

"""W3 / Kepler effective search volume of the #542 sample (injection-recovery, searched stars).

Plan: ``results/w3_kepler/effsearch/PLAN.md`` (frozen before any data were restored or
injected). The #529 / #542 screen is used unchanged (``w3_kepler.screen_star(local=True)``,
survivor ``pass_local``, common mode from ``w3_kepler_expand``). Nothing here changes a threshold
or adds a test.

Stages (``python scripts/w3_kepler_effsearch.py <stage>``):

- ``fetch``: restore the 500 stars' Q1-Q17 files by the names and sha256 in
  ``expand/files_expand.json.gz`` through the download ledger (Q0 not fetched);
- ``check``: re-run the uninjected #542 screen and compare it with ``expand/screen_expand.json.gz``
  and ``expand/runs_expand.json.gz``; per-star sampling windows;
- ``inject``: 2 injections per star and W3 cell, each from the original SAP series;
- ``summarise``: effective search volume per cell, uncertainties, loss stages, conditional limits
  and the expansion estimate.

Labels: injections simulated; light curves observed; recoveries and exposures derived.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w3_kepler as wk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w3_kepler" / "effsearch"
EXP = ROOT / "results" / "w3_kepler" / "expand"
CACHE = Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data")) / "cache" / "w3_kepler_effsearch"
CAP = 6_000_000_000  # owner limit for this round, bytes


def configure() -> None:
    wk.OUT, wk.CACHE, wk.LEDGER_CAP = OUT, CACHE, CAP


def usage(stage: str) -> None:
    with open(OUT / "source_usage.jsonl", "a") as fh:
        fh.write(json.dumps({"t": round(time.time(), 1), "stage": stage}) + "\n")


def stars() -> list[int]:
    return json.loads((EXP / "sample_expand.json").read_text())["stars"]


def file_list() -> dict[int, list[dict]]:
    d = json.loads(gzip.open(EXP / "files_expand.json.gz", "rt").read())
    return {s["kepid"]: s["files"] for s in d["stars"]}


# ------------------------------------------------------------------------------------------- fetch


def columns(b: bytes) -> dict | None:
    """The columns ``w3_kepler._fetch_star`` caches from one ``_llc.fits`` file (None for Q0)."""
    from astropy.io import fits

    with fits.open(io.BytesIO(b)) as h:
        h0, d = h[0].header, h[1].data
        q = int(h0["QUARTER"])
        if not 1 <= q <= 17:
            return None
        arr = {c: np.asarray(d[c], float) for c in wk.COLS}
        arr["QUARTER"] = np.full(len(d), q, float)
        arr["CHANNEL"] = np.full(len(d), int(h0["CHANNEL"]), float)
        return arr


def restore_star(kepid: int, files: list[dict]) -> dict:
    """Fetch the recorded Q1-Q17 files of one star, check each sha256, cache as #542 did."""
    configure()
    kid = f"{kepid:09d}"
    if (CACHE / f"{kid}.npz").exists():
        return {"kepid": kepid, "cached": True}
    lg = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP)
    base = f"{wk.ARCHIVE}/{kid[:4]}/{kid}/"
    parts, got = [], []
    try:
        for f in sorted(files, key=lambda f: f["name"]):
            if not 1 <= f["quarter"] <= 17:
                continue  # Q0 is not fetched (PLAN section 1)
            b = wk.nl.fetch(base + f["name"], lg, est_bytes=f["bytes"], retries=2)
            ok = hashlib.sha256(b).hexdigest() == f["sha256"]
            got.append({"name": f["name"], "bytes": len(b), "sha256_match": ok})
            if not ok:
                return {"kepid": kepid, "error": f"sha256 mismatch {f['name']}", "files": got}
            parts.append(columns(b))
    except Exception as e:  # recorded; the star is not replaced
        return {"kepid": kepid, "error": f"{type(e).__name__}: {e}", "files": got}
    cat = {c: np.concatenate([p[c] for p in parts]) for c in parts[0]}
    o = np.argsort(cat["TIME"], kind="stable")
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / f"{kid}.npz", **{c: v[o] for c, v in cat.items()})
    return {"kepid": kepid, "files": got}


def stage_fetch(_a) -> None:
    usage("fetch")
    fl = file_list()
    order = stars()
    t0 = time.time()
    with ThreadPoolExecutor(8) as ex:
        info = list(ex.map(lambda k: restore_star(k, fl[k]), order))
    out = {"label": "observed", "elapsed_s": round(time.time() - t0, 1), "stars": info}
    with gzip.open(OUT / "restore_files.json.gz", "wt") as fh:
        fh.write(json.dumps(out) + "\n")
    used = wk.nl.Ledger(OUT / "net_ledger.jsonl", CAP).used()
    print(len(info), "stars;", sum("error" in x for x in info), "errors; ledger", used)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("fetch")
    p.set_defaults(fn=stage_fetch)
    a = ap.parse_args(argv)
    configure()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

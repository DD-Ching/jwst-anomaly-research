"""Compare the top k of two pipeline runs of one sample (MIRI re-rank, ``results/miri_confirm/PLAN_rerank.md``).

Reads ``quality.ecsv``, ``scores.ecsv`` and ``sources.ecsv`` of ``<run>/<sample>`` for a run with
the band-restricted NIRCam confirmation (``new``) and one without (``base``), and writes the gate
counts, the top-k Jaccard overlap and the new / dropped top-k entries with their detection-band
count and top features. Provenance: scores are ``model_prediction``; the comparison is ``derived``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

TOP_K = 20  # ASSUMPTION (plan): the cutout stage's top_k


def top_uids(scores: Table, k: int = TOP_K) -> list[str]:
    """Source uids of the ``k`` best ranks (rank 1 first)."""
    order = np.argsort(np.asarray(scores["rank"]), kind="stable")
    return [str(scores["source_uid"][i]) for i in order[:k]]


def jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if sa | sb else float("nan")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("new_run", type=Path)
    p.add_argument("base_run", type=Path)
    p.add_argument("--sample", default="smacs0723_miri")
    p.add_argument("--out", type=Path, default=Path("results/miri_confirm/summary_rerank.json"))
    a = p.parse_args(argv)
    new_dir, base_dir = a.new_run / a.sample, a.base_run / a.sample
    out: dict = {"new_run": a.new_run.name, "base_run": a.base_run.name, "top_k": TOP_K}
    for name, d in (("new", new_dir), ("base", base_dir)):
        q = Table.read(d / "quality.ecsv")
        out[f"{name}_gate_passed"] = int(np.sum(np.asarray(q["quality_ok"], bool)))
    new_scores = Table.read(new_dir / "scores.ecsv")
    tn, tb = top_uids(new_scores), top_uids(Table.read(base_dir / "scores.ecsv"))
    src = Table.read(new_dir / "sources.ecsv")
    n_bands = dict(zip(np.asarray(src["source_uid"]).astype(str), np.asarray(src["n_bands"])))
    feats = dict(
        zip(
            np.asarray(new_scores["source_uid"]).astype(str),
            np.asarray(new_scores["top_features"]).astype(str),
        )
    )
    out["jaccard"] = round(jaccard(tn, tb), 3)
    out["new_top"] = [
        {"rank": tn.index(u) + 1, "source_uid": u, "n_bands": int(n_bands[u]), "top": feats[u]}
        for u in tn
        if u not in tb
    ]
    out["dropped_top"] = sorted(set(tb) - set(tn))
    out["new_top_single_band"] = sum(e["n_bands"] <= 1 for e in out["new_top"])
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k not in ("new_top",)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

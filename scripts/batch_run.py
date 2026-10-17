"""Run the fixed ``arc_diff`` chain over independent boxes (``batch_pairs.py`` output) in parallel.

``results/transients/batch_venus/PLAN.md`` §3–§5. For every box: the injection sample
(``fade_calib.py make`` on the epoch-1 F150W mosaic, seed ``--seed0`` + box index), one screen
and one pre-matching injection run of ``arc_diff.py`` with the fixed settings,
``arc_diff_morph.py --all-hits`` on both, ``fade_calib.py eval`` on the injection run, and the
per-box pass rule (vanish 26.0 + 27.0 AB on the review list with the correct sign >= ``--pass-
vanish``, unchanged-control false hits <= ``--max-false``). Tasks are subprocesses in a
``--workers`` pool (one process each, ``OMP_NUM_THREADS=1``, shared ``--cache-dir``); each
task's wall time and peak memory are recorded. ``summarize`` then writes ``batch.json`` (per-box
calibration, screen counts, status ``passed`` / ``undeterminable``) and ``review_queue.json``
(review-list entries ordered: passing boxes first, then max|S| descending).

    python scripts/batch_run.py run --boxes boxes.json --out DIR --scratch SCRATCH [--workers 4]
    python scripts/batch_run.py summarize --boxes boxes.json --out DIR
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable
SETTINGS = [
    "--half-arcsec", "12", "--kernel-px", "4", "--kernel-stamps", "--min-bands", "2",
    "--reg-tile-arcsec", "6", "--procs", "1",
]  # fmt: skip


def _run(cmd: list[str], log: Path) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    t0 = time.time()
    with open(log, "w") as f:
        p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
        _, status, ru = os.wait4(p.pid, 0)
    return {
        "rc": os.waitstatus_to_exitcode(status),
        "wall_s": round(time.time() - t0, 1),
        "peak_rss_mb": round(ru.ru_maxrss / 1024, 1),
    }


def _bands(box: dict) -> list[str]:
    args = []
    for b in ("F150W", "F444W"):
        args += ["--band", b.lower(), box["e1"][b]["uri"], box["e2"][b]["uri"]]
    return args


def task(kind: str, i: int, box: dict, out: Path, scratch: Path, cache: Path, seed0: int) -> dict:
    d = out / box["box_id"]
    d.mkdir(parents=True, exist_ok=True)
    sample = d / "sample.csv"
    rec = {"box_id": box["box_id"], "kind": kind}
    if kind == "inj" and not sample.exists():
        rec["sample"] = _run(
            [PY, str(HERE / "fade_calib.py"), "make", "--uri", box["e1"]["F150W"]["uri"],
             "--ra", str(box["ra"]), "--dec", str(box["dec"]), "--seed", str(seed0 + i),
             "--out", str(sample)],
            d / "sample.log",
        )  # fmt: skip
    run = d / kind
    maps = scratch / box["box_id"] / kind / "maps.npz"
    maps.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        PY,
        str(HERE / "arc_diff.py"),
        *_bands(box),
        "--ra",
        str(box["ra"]),
        "--dec",
        str(box["dec"]),
    ]
    cmd += [*SETTINGS, "--cache-dir", str(cache), "--out", str(run), "--maps", str(maps)]
    if kind == "inj":
        cmd += ["--inject-pre", str(sample), "--inject-pre-psf", "empirical"]
    rec["arc_diff"] = _run(cmd, d / f"{kind}.log")
    if rec["arc_diff"]["rc"] == 0:
        rec["morph"] = _run(
            [PY, str(HERE / "arc_diff_morph.py"), "--hits", str(run / "hits.ecsv"),
             "--maps", str(maps), "--all-hits", "--out", str(run / "morph_all.json")],
            d / f"{kind}_morph.log",
        )  # fmt: skip
    (d / f"{kind}_task.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def summarize(boxes: list[dict], out: Path, pass_vanish: float, max_false: int) -> None:
    from astropy.table import Table

    sys.path.insert(0, str(HERE))
    import fade_calib as fc

    rows, queue = [], []
    for box in boxes:
        d = out / box["box_id"]
        r = {"box_id": box["box_id"], "target": box["target"], "ra": box["ra"], "dec": box["dec"]}
        for kind in ("screen", "inj"):
            tj = d / f"{kind}_task.json"
            r[f"{kind}_task"] = json.loads(tj.read_text()) if tj.exists() else None
        inj = d / "inj"
        if (inj / "morph_all.json").exists():
            st = fc.summarize(fc.stage_rows(inj, inj / "morph_all.json"))
            (d / "inj_eval.json").write_text(json.dumps(st, indent=1) + "\n")
            v = [st["vanish"][m] for m in ("26.00", "27.00") if m in st.get("vanish", {})]
            vk, vn = sum(c["review_sign"] for c in v), sum(c["n"] for c in v)
            fh = st["unchanged_false_hits"]
            r["calibration"] = {
                "vanish_26_27_review_sign": [vk, vn],
                "unchanged_false_hits": [fh["false_hits"], fh["n"]],
                "per_cell": {
                    c: {
                        m: [x["review_sign"], x["n"], x["review_sign_68"]] for m, x in st[c].items()
                    }
                    for c in ("vanish", "partial", "appear")
                    if c in st
                },
            }
            ok = vn > 0 and vk / vn >= pass_vanish and fh["false_hits"] <= max_false
            r["status"] = "passed" if ok else "undeterminable"
        else:
            r["status"] = "undeterminable"
            r["why"] = "injection run missing"
        scr = d / "screen"
        if (scr / "morph_all.json").exists():
            h = Table.read(scr / "hits.ecsv")
            verdict = {
                t["id"]: t["verdict"]
                for t in json.loads((scr / "morph_all.json").read_text())["targets"]
            }
            rev = [x for x in h if not x["near_bright_core"] and verdict.get(x["hit_id"]) == "pass"]
            r["screen"] = {
                "hits": len(h),
                "unflagged": int(sum(not x["near_bright_core"] for x in h)),
                "review": len(rev),
                "summary_shift_px": {
                    k: v["shift_px"]
                    for k, v in json.loads((scr / "summary.json").read_text())["bands"].items()
                },
            }
            for x in rev:
                queue.append(
                    {
                        "box_id": box["box_id"],
                        "box_status": r["status"],
                        "hit_id": str(x["hit_id"]),
                        "ra": float(x["ra"]),
                        "dec": float(x["dec"]),
                        "sign": str(x["sign"]),
                        "bands": str(x["bands"]),
                        "max_abs_sigma": float(x["max_abs_sigma"]),
                    }
                )
        else:
            r["screen"] = None
            r["status"] = "undeterminable"
            r["why"] = "screen missing"
        rows.append(r)
    queue.sort(key=lambda q: (q["box_status"] != "passed", -q["max_abs_sigma"]))
    for k, q in enumerate(queue):
        q["rank"] = k + 1
    (out / "batch.json").write_text(
        json.dumps({"label": "derived; PLAN §3-§5", "boxes": rows}, indent=1) + "\n"
    )
    (out / "review_queue.json").write_text(
        json.dumps({"label": "derived; PLAN §5 order", "entries": queue}, indent=1) + "\n"
    )
    print(
        json.dumps(
            {
                "passed": sum(r["status"] == "passed" for r in rows),
                "boxes": len(rows),
                "queue": len(queue),
            }
        )
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=("run", "summarize"))
    ap.add_argument("--boxes", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--scratch", type=Path)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed0", type=int, default=100)
    ap.add_argument("--pass-vanish", type=float, default=0.70)
    ap.add_argument("--max-false", type=int, default=2)
    a = ap.parse_args(argv)
    boxes = json.loads(a.boxes.read_text())["chosen"]
    if a.cmd == "run":
        cache = a.scratch / "cache"
        jobs = [(k, i, b) for i, b in enumerate(boxes) for k in ("screen", "inj")]
        with ThreadPoolExecutor(a.workers) as ex:
            futs = [ex.submit(task, k, i, b, a.out, a.scratch, cache, a.seed0) for k, i, b in jobs]
            for f in futs:
                r = f.result()
                print(r["box_id"], r["kind"], r.get("arc_diff"), flush=True)
    summarize(boxes, a.out, a.pass_vanish, a.max_false)
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""D115b Lenstool re-fit, stage V: Lenstool's chi2 of the pinned CANUCS MACS0416 model.

Nothing is freed.

Plan: ``results/lens_consistency/d115b_lenstool_refit/PLAN.md`` (gate V, frozen before the run).

    python scripts/d115b_lenstool_refit.py --lenstool <prefix>/bin/lenstool --workdir <scratch dir>
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import subprocess
import time
from pathlib import Path

import lens_consistency as lc

MODEL = "macs0416-canucs"
REF_CHI2POS = 344.298308  # best.par header "#Chi2pos" (CANUCS optimisation, 2024-09-05)
TOL_FRAC = 0.01  # ASSUMPTION (plan)
OUT = Path(__file__).resolve().parents[1] / "results/lens_consistency/d115b_lenstool_refit"


def _sections(text: str) -> list[tuple[str, list[str]]]:
    """Split a Lenstool .par into (header line, body lines) blocks.

    Lines before the first block go to a leading ('', ...) block."""
    blocks: list[tuple[str, list[str]]] = [("", [])]
    in_block = False
    for line in text.splitlines():
        s = line.strip()
        if not in_block and s and not s.startswith("#") and not line[:1].isspace():
            blocks.append((line, []))
            in_block = s.split()[0] not in ("fini", "finish")
            continue
        blocks[-1][1].append(line)
        if s.split()[:1] == ["end"]:
            in_block = False
    return blocks


def _block(text: str, name: str) -> tuple[str, list[str]]:
    for head, body in _sections(text):
        if head.split()[:1] == [name]:
            return head, body
    raise ValueError(f"no '{name}' section")


def build_v0(best_par: str, input_par: str, multfile: str = "arcs.dat") -> str:
    """``best.par`` with runmode ``inverse 0`` + ``chi2 1`` and ``input.par``'s ``image`` section.

    The image section's ``multfile`` is pointed at ``multfile`` (the published image list): Lenstool
    8.15.5 silently reads no image, and reports chi2 0, when the named file is absent. ``chi2 1`` is
    what makes Lenstool 8.15.5 print the chi2 and write ``chires.dat`` (deviation, see RESULT)."""
    _, rm = _block(best_par, "runmode")
    ref = [ln for ln in rm if ln.split()[:1] == ["reference"]]
    if len(ref) != 1:
        raise ValueError("best.par runmode needs exactly one reference line")
    _, img = _block(input_par, "image")
    img = [ln for ln in img if ln.strip()]
    mf = [i for i, ln in enumerate(img) if ln.split()[:1] == ["multfile"]]
    if len(mf) != 1:
        raise ValueError("input.par image section needs exactly one multfile line")
    img[mf[0]] = f"    multfile    {img[mf[0]].split()[1]} {multfile}"
    out = []
    for head, body in _sections(best_par):
        if head.split()[:1] == ["runmode"]:
            out += ["runmode", ref[0], "\tinverse 0", "\tchi2 1", "\tend", "image", *img]
            continue
        if head:
            out.append(head)
        out += body
    return "\n".join(out) + "\n"


_CHI2 = re.compile(r"INFO: chi2 ([-+0-9.eE]+)")


def read_chi2pos(log: str) -> float:
    """Last ``INFO: chi2 <value>`` Lenstool printed (``forme -1``: image-plane position chi2)."""
    vals = _CHI2.findall(log)
    if not vals:
        raise ValueError("no 'INFO: chi2' in the Lenstool output")
    return float(vals[-1])


def gate_v(chi2: float) -> bool:
    """Frozen gate V: Lenstool's position chi2 within TOL_FRAC of the best.par header value."""
    return abs(chi2 - REF_CHI2POS) <= TOL_FRAC * REF_CHI2POS


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lenstool", required=True)
    ap.add_argument("--workdir", required=True, type=Path)
    a = ap.parse_args()
    files = lc.model_files(MODEL)
    wd = a.workdir
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(files["arcs.dat"], wd / "arcs.dat")
    (wd / "v0.par").write_text(
        build_v0(files["best.par"].read_text(), files["input.par"].read_text())
    )
    t0 = time.time()
    p = subprocess.run(
        [a.lenstool, "v0.par", "-n"], cwd=wd, capture_output=True, text=True, timeout=3600
    )
    dt = time.time() - t0
    log = p.stdout + p.stderr
    OUT.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / "lenstool_v0.log.gz", "wt") as f:
        f.write(log)
    chi2 = read_chi2pos(log)
    out = {
        "plan": "results/lens_consistency/d115b_lenstool_refit/PLAN.md",
        "provenance": "model_prediction",
        "lenstool": "conda-forge lenstool 8.15.5 (linux-64), runmode inverse 0, pinned best.par",
        "returncode": p.returncode,
        "runtime_s": round(dt, 1),
        "chi2pos_lenstool": chi2,
        "chi2pos_best_par_header": REF_CHI2POS,
        "rel_diff": (chi2 - REF_CHI2POS) / REF_CHI2POS,
        "tol_frac_assumption": TOL_FRAC,
        "gate_V": gate_v(chi2),
    }
    (OUT / "result.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

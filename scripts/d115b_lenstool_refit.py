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
REF_CHI2POS = (
    344.298308  # best.par header "#Chi2pos" (CANUCS optimisation); checked against the file
)
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


def read_chires(text: str) -> tuple[float, dict[str, float], int]:
    """``chires.dat``: the ``chimul`` (position) chi2, per-family chi2 and the summed warnings.

    A family's summary row is the one whose ``dx`` is ``N/A``; it is keyed by the family name (the
    row's image id without its ``.n`` suffix)."""
    chimul, fam, nwarn = None, {}, 0
    for line in text.splitlines():
        p = line.split()
        if p[:1] == ["chimul"]:
            chimul = float(p[1])
        elif len(p) >= 13 and p[0].isdigit():
            nwarn += int(p[12])
            if p[10] == "N/A":
                fam[p[1].rsplit(".", 1)[0]] = float(p[4])
    if chimul is None:
        raise ValueError("no chimul row in chires.dat")
    return chimul, fam, nwarn


def gate_v(chi2: float) -> bool:
    """Frozen gate V: Lenstool's position chi2 within TOL_FRAC of the best.par header value."""
    return abs(chi2 - REF_CHI2POS) <= TOL_FRAC * REF_CHI2POS


def run(lenstool: Path, wd: Path, par: str, arcs: Path) -> tuple[subprocess.CompletedProcess, str]:
    """Run Lenstool on ``par`` in ``wd`` (a nonzero exit raises); return the run and chires.dat."""
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(arcs, wd / "arcs.dat")
    (wd / "v.par").write_text(par)
    (wd / "chires.dat").unlink(missing_ok=True)
    p = subprocess.run([str(lenstool), "v.par", "-n"], cwd=wd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"lenstool exited {p.returncode} in {wd}")
    return p, (wd / "chires.dat").read_text()


def version(lenstool: Path) -> str:
    """Package version and build from the conda-meta record next to the binary."""
    meta = sorted((lenstool.parent.parent / "conda-meta").glob("lenstool-*.json"))
    if len(meta) != 1:
        return "unknown"
    d = json.loads(meta[0].read_text())
    return f"{d['name']} {d['version']} ({d['build']}, {d.get('channel', '?')})"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lenstool", required=True, type=Path)
    ap.add_argument("--workdir", required=True, type=Path)
    a = ap.parse_args()
    lt, wd = a.lenstool.resolve(), a.workdir.resolve()
    files = lc.model_files(MODEL)
    if lc.chi2pos_from_par(files["best.par"]) != REF_CHI2POS:
        raise ValueError("best.par #Chi2pos differs from REF_CHI2POS")
    best, inp = files["best.par"].read_text(), files["input.par"].read_text()
    v0 = build_v0(best, inp)
    t0 = time.time()
    p, chires = run(lt, wd / "v0", v0, files["arcs.dat"])
    dt = time.time() - t0
    chimul, fam, nwarn = read_chires(chires)
    OUT.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / "lenstool_v0.log.gz", "wt") as f:
        f.write(p.stdout + p.stderr)
    with gzip.open(OUT / "chires_v0.dat.gz", "wt") as f:
        f.write(chires)
    out = {
        "plan": "results/lens_consistency/d115b_lenstool_refit/PLAN.md",
        "provenance": "model_prediction",
        "lenstool": version(lt),
        "runmode": "inverse 0, chi2 1 (image plane, forme -1, sigposArcsec 0.49)",
        "runtime_s": round(dt, 1),
        "chi2_info_line": read_chi2pos(p.stdout + p.stderr),
        "chi2pos_lenstool": chimul,
        "chi2pos_best_par_header": REF_CHI2POS,
        "rel_diff": (chimul - REF_CHI2POS) / REF_CHI2POS,
        "tol_frac_assumption": TOL_FRAC,
        "n_families": len(fam),
        "n_warnings": nwarn,
        "gate_V": gate_v(chimul),
    }
    (OUT / "result.json").write_text(json.dumps(out, indent=1) + "\n")
    # post hoc (after gate V failed): solver grid size and the full-precision input.par reference
    post = {"label": "post hoc", "chi2pos_by_grille_nombre": {}}
    for n in (128, 256, 512):
        par = v0.replace("\tnombre      128\n", f"\tnombre      {n}\n")
        if n != 128 and par == v0:
            raise ValueError("grille nombre line not found")
        post["chi2pos_by_grille_nombre"][n] = read_chires(
            run(lt, wd / f"grid{n}", par, files["arcs.dat"])[1]
        )[0]
    _, rm_in = _block(inp, "runmode")
    ref_in = [ln for ln in rm_in if ln.split()[:1] == ["reference"]][0]
    _, rm_best = _block(best, "runmode")
    ref_best = [ln for ln in rm_best if ln.split()[:1] == ["reference"]][0]
    post["reference_input_par"] = ref_in.split()[1:4]
    post["chi2pos_input_par_reference"] = read_chires(
        run(lt, wd / "refin", v0.replace(ref_best, ref_in), files["arcs.dat"])[1]
    )[0]
    top = sorted(fam.items(), key=lambda kv: -kv[1])
    post["top8_families_chi2"] = top[:8]
    post["family_chi2_D115a_D115b"] = [fam.get("D115a"), fam.get("D115b")]
    post["sum_family_chi2"] = round(sum(fam.values()), 2)
    (OUT / "posthoc.json").write_text(json.dumps(post, indent=1) + "\n")
    print(json.dumps(out, indent=1), json.dumps(post, indent=1))


if __name__ == "__main__":
    main()

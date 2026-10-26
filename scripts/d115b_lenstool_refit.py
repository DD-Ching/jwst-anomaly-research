"""D115b Lenstool re-fit: Lenstool's chi2 of the pinned CANUCS MACS0416 model.

Stage V (default): nothing is freed; plan ``results/lens_consistency/d115b_lenstool_refit/PLAN.md``.
Stage V' and F (``--stage vprime-f``): the full-precision reference, then a chi2 profile in member
8515's sigma; plan ``PLAN_VPRIME_F.md`` (gates V' and rule F, frozen before the run).

    python scripts/d115b_lenstool_refit.py --lenstool <prefix>/bin/lenstool --workdir <scratch dir>
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import lens_consistency as lc
import numpy as np
from lenstool_cloud_check import FLAG, FLAG_IMAGES, near_flag, read_image_all

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


SIGMA_8515 = 120.600977  # best.par "potentiel 8515" v_disp (km/s)
F_GRID = (1.00, 0.90, 0.85, 0.80, 0.75, 0.70)  # plan PLAN_VPRIME_F.md
DCHI2_MAX = 6.0  # ASSUMPTION (plan): system-51 precedent, as in the port scans
F2_WALL_S = 600  # plan: wall cap of the stage F2 cost probe


def _ref_line(text: str) -> str:
    _, rm = _block(text, "runmode")
    return [ln for ln in rm if ln.split()[:1] == ["reference"]][0]


def build_vprime(v0: str, best_par: str, input_par: str) -> str:
    """Stage V': the stage V run file with ``input.par``'s full-precision reference line."""
    old, new = _ref_line(best_par), _ref_line(input_par)
    if v0.count(old + "\n") != 1:
        raise ValueError("best.par reference line not found once in the run file")
    return v0.replace(old + "\n", new + "\n")


def _pot_body(text: str, name: str) -> tuple[int, int]:
    """Line span (start, end) of ``potentiel <name>``'s body; ``end`` is its ``end`` line."""
    lines = text.split("\n")
    heads = [i for i, ln in enumerate(lines) if ln.split()[:2] == ["potentiel", name]]
    if len(heads) != 1:
        raise ValueError(f"potentiel {name} not found once")
    j = heads[0] + 1
    while lines[j].split()[:1] != ["end"]:
        j += 1
    return heads[0], j


def set_vdisp(par: str, name: str, vdisp: float) -> str:
    """Set ``v_disp`` of ``potentiel <name>`` (exactly one such line in its body)."""
    lines = par.split("\n")
    a, b = _pot_body(par, name)
    hits = [i for i in range(a + 1, b) if lines[i].split()[:1] == ["v_disp"]]
    if len(hits) != 1:
        raise ValueError(f"potentiel {name}: need one v_disp line")
    lines[hits[0]] = f"\tv_disp     {vdisp:.6f}"
    return "\n".join(lines)


def add_image_mode(par: str) -> str:
    """Add ``image 1 arcs.dat`` to the runmode (Lenstool then writes ``image.all``)."""
    if par.count("\tinverse 0\n") != 1:
        raise ValueError("runmode inverse 0 line not found once")
    return par.replace("\tinverse 0\n", "\tinverse 0\n\timage 1 arcs.dat\n")


def rule_f(dchi2: dict[float, float], flag: dict[float, bool]) -> tuple[str, float | None]:
    """Frozen rule F on the sigma(8515) profile: (outcome, f_L*)."""
    fs = sorted(flag, reverse=True)
    if not flag[1.0] or all(flag[f] for f in fs):
        return "undetermined", None
    f_star = next(f for f in fs if not flag[f])
    return ("cheap" if dchi2[f_star] <= DCHI2_MAX else "costly"), f_star


def flag_images(image_all: Path) -> dict[str, bool]:
    """Per FLAG_IMAGES id: any ``image.all`` image within FLAG_RADIUS of the #458 flag."""
    (ra0, dec0), rows = read_image_all(image_all)
    # Lenstool's frame, as lensmodel.LensModel.to_frame: x = -dRA cos(dec0), y = dDec (arcsec)
    fx = -((FLAG[0] - ra0 + 180.0) % 360.0 - 180.0) * np.cos(np.radians(dec0)) * 3600.0
    fy = (FLAG[1] - dec0) * 3600.0
    return {i: i in rows and near_flag(rows[i], fx, fy) for i in FLAG_IMAGES}


def build_f2(vprime: str, best_par: str, input_par: str) -> str:
    """Stage F2 cost probe: ``inverse 3`` (input.par's rate and sample count), O1-O4 freed with
    ``input.par``'s ``limit 1``-``limit 4`` and 8515's v_disp free in [60, 200] km/s."""
    names = [ln.split()[1] for ln in best_par.splitlines() if ln.split()[:1] == ["potentiel"]]
    if names[:4] != ["O1", "O2", "O3", "O4"]:
        raise ValueError("best.par does not start with potentials O1-O4")
    lim = []
    for i in (1, 2, 3, 4):
        _, body = _limit(input_par, i)
        lim += [f"limit {i}", *[ln for ln in body if ln.strip()]]
    lim += [f"limit {names.index('8515') + 1}", "    v_disp      1 60. 200. 0.1", "    end"]
    _, rm = _block(input_par, "runmode")
    inv = [ln for ln in rm if ln.split()[:1] == ["inverse"]]
    if len(inv) != 1 or vprime.count("\tinverse 0\n") != 1:
        raise ValueError("inverse line not found once")
    par = vprime.replace("\tinverse 0\n", "\t" + inv[0].split("#")[0].strip() + "\n")
    lines = par.rstrip("\n").split("\n")
    if lines[-1].strip() not in ("fini", "finish"):
        raise ValueError("run file does not end with fini")
    return "\n".join(lines[:-1] + lim + lines[-1:]) + "\n"


def _limit(text: str, i: int) -> tuple[str, list[str]]:
    for head, body in _sections(text):
        if head.split()[:2] == ["limit", str(i)]:
            return head, body
    raise ValueError(f"no 'limit {i}' section")


def run(
    lenstool: Path, wd: Path, par: str, arcs: Path, timeout: float | None = None
) -> tuple[subprocess.CompletedProcess, str]:
    """Run Lenstool on ``par`` in ``wd`` (a nonzero exit raises); return the run and chires.dat."""
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(arcs, wd / "arcs.dat")
    (wd / "v.par").write_text(par)
    (wd / "chires.dat").unlink(missing_ok=True)
    env = {**os.environ, "OMP_NUM_THREADS": str(os.cpu_count() or 1)}  # undo the import's 1
    p = subprocess.run(
        [str(lenstool), "v.par", "-n"],
        cwd=wd,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
    )
    if p.returncode != 0:
        raise RuntimeError(f"lenstool exited {p.returncode} in {wd}")
    chires = wd / "chires.dat"
    return p, chires.read_text() if chires.exists() else ""


def version(lenstool: Path) -> str:
    """Package version and build from the conda-meta record next to the binary."""
    meta = sorted((lenstool.parent.parent / "conda-meta").glob("lenstool-*.json"))
    if len(meta) != 1:
        return "unknown"
    d = json.loads(meta[0].read_text())
    return f"{d['name']} {d['version']} ({d['build']}, {d.get('channel', '?')})"


def main_vprime_f(lt: Path, wd: Path, files: dict[str, Path]) -> None:
    """Stages V' and F, and the F2 cost probe (plan PLAN_VPRIME_F.md)."""
    best, inp = files["best.par"].read_text(), files["input.par"].read_text()
    vp = build_vprime(build_v0(best, inp), best, inp)
    p, chires = run(lt, wd / "vprime", vp, files["arcs.dat"])
    chimul, fam, nwarn = read_chires(chires)
    out = {
        "plan": "results/lens_consistency/d115b_lenstool_refit/PLAN_VPRIME_F.md",
        "provenance": "model_prediction",
        "lenstool": version(lt),
        "reference_input_par": _ref_line(inp).split()[1:4],
        "chi2pos_vprime": chimul,
        "chi2_info_line_vprime": read_chi2pos(p.stdout + p.stderr),
        "rel_diff_vprime": (chimul - REF_CHI2POS) / REF_CHI2POS,
        "n_families": len(fam),
        "n_warnings": nwarn,
        "gate_Vprime": gate_v(chimul),
    }
    if out["gate_Vprime"]:
        prof = {}
        for f in F_GRID:
            par = add_image_mode(set_vdisp(vp, "8515", f * SIGMA_8515))
            d = wd / f"f{f:.2f}"
            prof[f] = (
                read_chires(run(lt, d, par, files["arcs.dat"])[1])[0],
                flag_images(d / "image.all"),
            )
        dchi2 = {f: round(c - prof[1.0][0], 4) for f, (c, _) in prof.items()}
        flag = {f: any(fl.values()) for f, (_, fl) in prof.items()}
        outcome, f_star = rule_f(dchi2, flag)
        out["stage_F"] = {
            "sigma_8515_fiducial_kms": SIGMA_8515,
            "dchi2_max_assumption": DCHI2_MAX,
            "profile": [
                {
                    "f": f,
                    "v_disp_kms": round(f * SIGMA_8515, 4),
                    "chi2pos": c,
                    "dchi2": dchi2[f],
                    "flag_by_image": fl,
                    "flag": flag[f],
                }
                for f, (c, fl) in prof.items()
            ],
            "f_L_star": f_star,
            "rule_F": outcome,
        }
        t0 = time.time()
        try:
            p2, _ = run(lt, wd / "f2", build_f2(vp, best, inp), files["arcs.dat"], F2_WALL_S)
            done = True
        except subprocess.TimeoutExpired:
            done = False
        out["stage_F2_probe"] = {
            "label": "reported, not gated",
            "wall_cap_s": F2_WALL_S,
            "finished": done,
            "runtime_s": round(time.time() - t0, 1),
        }
    (OUT / "result_vprime_f.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lenstool", required=True, type=Path)
    ap.add_argument("--workdir", required=True, type=Path)
    ap.add_argument("--stage", choices=("v", "vprime-f"), default="v")
    a = ap.parse_args()
    lt, wd = a.lenstool.resolve(), a.workdir.resolve()
    files = lc.model_files(MODEL)
    if lc.chi2pos_from_par(files["best.par"]) != REF_CHI2POS:
        raise ValueError("best.par #Chi2pos differs from REF_CHI2POS")
    if a.stage == "vprime-f":
        main_vprime_f(lt, wd, files)
        return
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

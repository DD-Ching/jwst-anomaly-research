"""D115b Lenstool re-fit: Lenstool's chi2 of the pinned CANUCS MACS0416 model.

Stage V (default): nothing is freed; plan ``results/lens_consistency/d115b_lenstool_refit/PLAN.md``.
Stage V' and F (``--stage vprime-f``): the full-precision reference, then a chi2 profile in member
8515's sigma; plan ``PLAN_VPRIME_F.md`` (gates V' and rule F, frozen before the run).
Stage D14a (``--stage d14a``): the stage F runs again, read for Lenstool's image-proximity penalty;
plan ``PLAN_D14A.md`` (rule P, frozen before the run).
Stage F2-B (``--stage f2b``): a short ``inverse 3`` run with O1-O6 free and 8515 fixed, a budget
probe of Lenstool's optimiser; plan ``PLAN_F2.md`` (gate V2, frozen before the run).
Stage F2-L V3 (``--stage f2l-v3``): bounded Nelder-Mead over O1-O6 driving Lenstool
``inverse 0``, validated from a perturbed best.par at f 1.00; plan ``PLAN_F2L.md`` (gate V3).
Stage V3' (``--stage v3p``): the same with a small-step start ladder; plan ``PLAN_V3P.md``.
Stage VS (``--stage vs``): continuity of the source-plane chi2 (``forme 0``) at best.par under
single-parameter offsets; plan ``PLAN_VS.md`` (gate VS, frozen before the run).

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

_OMP = os.environ.get("OMP_NUM_THREADS")  # the caller's setting, before lenstool_cloud_check sets 1

import lens_consistency as lc  # noqa: E402
import numpy as np  # noqa: E402
from lenstool_cloud_check import FLAG, FLAG_IMAGES, near_flag, read_image_all  # noqa: E402

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
    """Last ``INFO: chi2 <value>`` Lenstool printed (image-plane position chi2 with ``forme -1``,
    source-plane chi2 with ``forme 0``: stage VS, where it matched chires.dat's ``chitot``)."""
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
# Lenstool 8.15.5 src/o_chires.c chi2_img: (chip_penalty_dist 30")^2 / sigma^2,
# sigma = sigposArcsec 0.49
PENALTY = (30.0 / 0.49) ** 2
_PEN_WARN = re.compile(r"WARNING: chi2 penalty for image (\S+)")
F2_WALL_S = 600  # plan: wall cap of the stage F2 cost probe
F2B_INVERSE = "inverse 3 0.5 20"  # ASSUMPTION (PLAN_F2.md): smallest budget that writes bayes.dat
F2B_NLENS_OPT = 6  # input.par's limit 1-6 = best.par's first six potentiels O1-O6
F2B_WALL_S = 900  # ASSUMPTION (PLAN_F2.md): wall cap of the probe
V2_REF = 344.42  # gate V' value (result_vprime_f.json chi2pos_vprime)


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


def build_f2b(vprime: str, best_par: str, input_par: str) -> str:
    """Stage F2-B: ``F2B_INVERSE``, ``nlens_opt`` in ``grille``, ``input.par``'s ``limit 1``-``6``
    appended unchanged; 8515 and the potfile stay fixed."""
    names = [ln.split()[1] for ln in best_par.splitlines() if ln.split()[:1] == ["potentiel"]]
    if names[:F2B_NLENS_OPT] != [f"O{i}" for i in range(1, F2B_NLENS_OPT + 1)]:
        raise ValueError("best.par does not start with potentials O1-O6")
    lim = []
    for i in range(1, F2B_NLENS_OPT + 1):
        _, body = _limit(input_par, i)
        lim += [f"limit {i}", *[ln for ln in body if ln.strip()]]
    if vprime.count("\tinverse 0\n") != 1:
        raise ValueError("inverse line not found once")
    par = vprime.replace("\tinverse 0\n", "\t" + F2B_INVERSE + "\n")
    lines = par.rstrip("\n").split("\n")
    g = [i for i, ln in enumerate(lines) if ln.strip() == "grille"]
    if len(g) != 1:
        raise ValueError("grille section not found once")
    lines.insert(g[0] + 1, f"\tnlens_opt   {F2B_NLENS_OPT}")
    if lines[-1].strip() not in ("fini", "finish"):
        raise ValueError("run file does not end with fini")
    return "\n".join(lines[:-1] + lim + lines[-1:]) + "\n"


def read_bayes(text: str) -> tuple[int, float | None]:
    """(number of samples, lowest chi2) from Lenstool's ``bayes.dat``: ``#`` header lines name the
    columns in order; the chi2 column is the one named ``Chi2`` (case-insensitive)."""
    cols = [ln[1:].strip() for ln in text.splitlines() if ln.startswith("#")]
    rows = [ln.split() for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    idx = [i for i, c in enumerate(cols) if c.lower() == "chi2"]
    if len(idx) != 1:
        raise ValueError("bayes.dat: need one Chi2 column")
    vals = [float(r[idx[0]]) for r in rows if len(r) == len(cols)]
    return len(vals), (min(vals) if vals else None)


def gate_v2(status: str, chi2_min: float | None) -> str:
    """Frozen gate V2 (PLAN_F2.md): ``pass`` iff the run finished and its lowest chi2 is within 1 %
    above the gate V' value; a run that did not finish (or wrote no sample) is ``undetermined``."""
    if status != "finished" or chi2_min is None:
        return "undetermined"
    return "pass" if chi2_min <= V2_REF * (1.0 + TOL_FRAC) else "fail"


def main_f2b(lt: Path, wd: Path, files: dict[str, Path]) -> None:
    """Stage F2-B budget probe (plan PLAN_F2.md)."""
    best, inp = files["best.par"].read_text(), files["input.par"].read_text()
    par = build_f2b(build_vprime(build_v0(best, inp), best, inp), best, inp)
    d = wd / "f2b"
    for stale in ("bayes.dat", "burnin.dat"):
        (d / stale).unlink(missing_ok=True)
    t0 = time.time()
    log, chires = "", ""
    try:
        p, chires = run(lt, d, par, files["arcs.dat"], F2B_WALL_S, False)
        status, log = "finished", p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        status = "wall cap reached"
    except LenstoolError as e:
        status, log = str(e), e.log
    dt = time.time() - t0
    bayes = d / "bayes.dat"
    n, cmin = read_bayes(bayes.read_text()) if bayes.exists() else (0, None)
    out = {
        "plan": "results/lens_consistency/d115b_lenstool_refit/PLAN_F2.md",
        "provenance": "model_prediction",
        "lenstool": version(lt),
        "inverse": F2B_INVERSE,
        "nlens_opt": F2B_NLENS_OPT,
        "wall_cap_s": F2B_WALL_S,
        "omp_threads": _OMP or str(os.cpu_count() or 1),
        "status": status,
        "runtime_s": round(dt, 1),
        "lenstool_error_lines": [ln for ln in log.splitlines() if "ERROR" in ln],
        "n_samples_bayes": n,
        "chi2_min_bayes": cmin,
        "chi2_ref_vprime": V2_REF,
        "tol_frac_assumption": TOL_FRAC,
        "gate_V2": gate_v2(status, cmin),
    }
    if chires:
        out["chimul_final"] = read_chires(chires)[0]
        out["chipenalty_final"] = read_penalty(chires)
    if log:
        with gzip.open(OUT / "lenstool_f2b.log.gz", "wt") as g:
            g.write(log)
    if bayes.exists():
        with gzip.open(OUT / "bayes_f2b.dat.gz", "wt") as g:
            g.write(bayes.read_text())
    _write(out, "result_f2b.json")
    print(json.dumps(out, indent=1))


def posthoc_f2b(wd: Path) -> dict:
    """Post hoc (after gate V2 was undetermined): BayeSys burn-in lines and their best ln(Lhood)."""
    text = (wd / "f2b" / "burnin.dat").read_text()
    with gzip.open(OUT / "burnin_f2b.dat.gz", "wt") as g:
        g.write(text)
    rows = [ln.split() for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    lnl = [float(r[1]) for r in rows]
    runtime = json.loads((OUT / "result_f2b.json").read_text())["runtime_s"]
    return {
        "label": "post hoc",
        "n_burnin_lines": len(rows),
        "n_columns": sorted({len(r) for r in rows}),
        "lnLhood_max_burnin": max(lnl) if lnl else None,
        # for scale only: stage V's INFO line (chi2 356.67, rounded reference)
        "lnLhood_stage_v_info_line": float(
            re.findall(r"Lhood (\S+)", gzip.open(OUT / "lenstool_v0.log.gz", "rt").read())[-1]
        ),
        "seconds_per_burnin_line": round(runtime / len(rows), 1) if rows else None,
    }


def _limit(text: str, i: int) -> tuple[str, list[str]]:
    for head, body in _sections(text):
        if head.split()[:2] == ["limit", str(i)]:
            return head, body
    raise ValueError(f"no 'limit {i}' section")


class LenstoolError(RuntimeError):
    """A nonzero Lenstool exit; ``log`` holds its stdout + stderr."""

    def __init__(self, msg: str, log: str):
        super().__init__(msg)
        self.log = log


def run(
    lenstool: Path,
    wd: Path,
    par: str,
    arcs: Path,
    timeout: float | None = None,
    need_chires: bool = True,
) -> tuple[subprocess.CompletedProcess, str]:
    """Run Lenstool on ``par`` in ``wd`` (a nonzero exit raises LenstoolError); return the run and
    chires.dat ("" only when ``need_chires`` is false and Lenstool wrote none)."""
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(arcs, wd / "arcs.dat")
    (wd / "v.par").write_text(par)
    for stale in ("chires.dat", "image.all"):
        (wd / stale).unlink(missing_ok=True)
    env = {**os.environ, "OMP_NUM_THREADS": _OMP or str(os.cpu_count() or 1)}
    p = subprocess.run(
        [str(lenstool), "v.par", "-n"],
        cwd=wd,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
    )
    if p.returncode != 0:
        how = f"signal {-p.returncode}" if p.returncode < 0 else f"exit {p.returncode}"
        raise LenstoolError(f"lenstool stopped ({how})", p.stdout + p.stderr)
    chires = wd / "chires.dat"
    if need_chires or chires.exists():
        return p, chires.read_text()
    return p, ""


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
        _write(out)  # stage V' and F are stored before the probe runs
        for f in (1.00, 0.90):  # the fiducial and the first factor of the post hoc reading
            with gzip.open(OUT / f"chires_f{f:.2f}.dat.gz", "wt") as g:
                g.write((wd / f"f{f:.2f}" / "chires.dat").read_text())
        t0 = time.time()
        try:
            run(lt, wd / "f2", build_f2(vp, best, inp), files["arcs.dat"], F2_WALL_S, False)
            status, errors = "finished", []
        except subprocess.TimeoutExpired:
            status, errors = "wall cap reached", []
        except LenstoolError as e:
            status, errors = str(e), [ln for ln in e.log.splitlines() if "ERROR" in ln]
        out["stage_F2_probe"] = {
            "label": "reported, not gated",
            "wall_cap_s": F2_WALL_S,
            "status": status,
            "lenstool_error_lines": errors,
            "runtime_s": round(time.time() - t0, 1),
        }
        _write(posthoc_vprime_f(wd), "posthoc_vprime_f.json")
    _write(out)
    print(json.dumps(out, indent=1))


def _write(out: dict, name: str = "result_vprime_f.json") -> None:
    (OUT / name).write_text(json.dumps(out, indent=1) + "\n")


def _family_rows(text: str, fam: str) -> tuple[float, float]:
    """``chires.dat``: (family summary-row chi2, sum of the family's image-row chi2)."""
    rows = [ln.split() for ln in text.splitlines()]
    rows = [r for r in rows if len(r) >= 13 and r[0].isdigit() and r[1].rsplit(".", 1)[0] == fam]
    summ = [float(r[4]) for r in rows if r[10] == "N/A"]
    if len(summ) != 1:
        raise ValueError(f"{fam}: need one summary row")
    return summ[0], round(sum(float(r[4]) for r in rows if r[10] != "N/A"), 2)


def posthoc_vprime_f(wd: Path) -> dict:
    """Post hoc reading of the stage F runs (after rule F was applied): per-family chi2 changes,
    families with warnings, and D14a's summary row against its image rows."""
    texts = {f: (wd / f"f{f:.2f}" / "chires.dat").read_text() for f in F_GRID}
    c0, base, _ = read_chires(texts[1.0])
    per_f = {}
    for f, t in texts.items():
        c, fam, nwarn = read_chires(t)
        d = {k: round(fam[k] - base.get(k, 0.0), 2) for k in fam}
        warn = {}
        for r in (ln.split() for ln in t.splitlines()):
            if len(r) >= 13 and r[0].isdigit() and r[10] == "N/A" and int(r[12]) > 0:
                warn[r[1].rsplit(".", 1)[0]] = int(r[12])
        row, img = _family_rows(t, "D14a")
        per_f[f"{f:.2f}"] = {
            "chi2pos": c,
            "n_warnings": nwarn,
            "families_with_warnings": warn,
            "largest_family_dchi2": sorted(d.items(), key=lambda kv: -kv[1])[:4],
            "smallest_family_dchi2": sorted(d.items(), key=lambda kv: kv[1])[:3],
            "D14a_family_row_chi2": row,
            "D14a_image_rows_chi2_sum": img,
            "dchi2_without_D14a": round(c - c0 - d["D14a"], 2),
        }
    return {"label": "post hoc (after rule F was applied)", "per_f": per_f}


def read_penalty(chires: str) -> float:
    """Total of Lenstool's ``chipenalty`` line (it prints one only when a penalty was added)."""
    vals = [float(ln.split()[1]) for ln in chires.splitlines() if ln.split()[:1] == ["chipenalty"]]
    if len(vals) > 1:
        raise ValueError("more than one chipenalty line")
    return vals[0] if vals else 0.0


def penalty_images(log: str) -> list[str]:
    """Image ids named in Lenstool's ``WARNING: chi2 penalty for image <id>`` lines."""
    return _PEN_WARN.findall(log)


def rule_p(rows: dict[float, dict]) -> tuple[str, list[str]]:
    """Frozen rule P (PLAN_D14A.md): (outcome, reasons it failed). Each row holds ``penalty``,
    ``penalty_images``, ``D14a_row`` and ``D14a_images``."""
    bad = []
    for f, r in sorted(rows.items(), reverse=True):
        k = round(r["penalty"] / PENALTY)
        if abs(r["penalty"] - k * PENALTY) > 0.02:
            bad.append(f"f {f:.2f}: penalty {r['penalty']} is not a multiple of {PENALTY:.2f}")
        n = sum(i.rsplit(".", 1)[0] == "D14a" for i in r["penalty_images"])
        if abs(r["D14a_row"] - r["D14a_images"] - n * PENALTY) > 0.05:
            ex = r["D14a_row"] - r["D14a_images"]
            bad.append(f"f {f:.2f}: D14a row - images {ex:.2f} vs {n} penalties")
    return ("rejected" if bad else "confirmed"), bad


def main_d14a(lt: Path, wd: Path, files: dict[str, Path]) -> None:
    """Stage D14a: the stage F runs, read for the proximity penalty (plan PLAN_D14A.md)."""
    best, inp = files["best.par"].read_text(), files["input.par"].read_text()
    vp = build_vprime(build_v0(best, inp), best, inp)
    rows = {}
    for f in F_GRID:
        p, chires = run(
            lt, wd / f"p{f:.2f}", set_vdisp(vp, "8515", f * SIGMA_8515), files["arcs.dat"]
        )
        row, img = _family_rows(chires, "D14a")
        rows[f] = {
            "chi2pos": read_chires(chires)[0],
            "penalty": read_penalty(chires),
            "penalty_images": penalty_images(p.stdout + p.stderr),
            "D14a_row": row,
            "D14a_images": img,
        }
    outcome, bad = rule_p(rows)
    c0 = rows[1.0]["chi2pos"] - rows[1.0]["penalty"]
    out = {
        "plan": "results/lens_consistency/d115b_lenstool_refit/PLAN_D14A.md",
        "provenance": "model_prediction",
        "lenstool": version(lt),
        "penalty_unit": round(PENALTY, 4),
        "profile": [
            {"f": f, **r, "dchi2_no_penalty": round(r["chi2pos"] - r["penalty"] - c0, 2)}
            for f, r in rows.items()
        ],
        "rule_P": outcome,
        "rule_P_failures": bad,
    }
    _write(out, "result_d14a.json")
    _write(posthoc_d14a(wd), "posthoc_d14a.json")
    print(json.dumps(out, indent=1))


def posthoc_d14a(wd: Path) -> dict:
    """Post hoc (after rule P): Lenstool-predicted positions of D14a.1 and D14a.2 per f, from the
    observed positions in ``arcs.dat`` plus the chires offsets (chires prints -dx and dy)."""
    out = {}
    for f in F_GRID:
        d = wd / f"p{f:.2f}"
        obs = {}
        for ln in (d / "arcs.dat").read_text().splitlines():
            p = ln.split()
            if p[:1] and p[0] in ("D14a.1", "D14a.2"):
                obs[p[0]] = (float(p[1]), float(p[2]))
        ra0, dec0 = obs["D14a.1"]
        pred = {}
        for r in (ln.split() for ln in (d / "chires.dat").read_text().splitlines()):
            if len(r) >= 13 and r[1] in obs and r[10] != "N/A":
                ra, dec = obs[r[1]]
                x = -(ra - ra0) * np.cos(np.radians(dec0)) * 3600.0
                y = (dec - dec0) * 3600.0
                pred[r[1]] = (x - float(r[10]), y + float(r[11]))
        a, b = pred.get("D14a.1"), pred.get("D14a.2")
        out[f"{f:.2f}"] = {
            "pred_D14a.1_xy_arcsec": a and [round(v, 3) for v in a],
            "pred_D14a.2_xy_arcsec": b and [round(v, 3) for v in b],
            "pred_separation_arcsec": (
                round(float(np.hypot(a[0] - b[0], a[1] - b[1])), 3) if a and b else None
            ),
        }
    return {
        "label": "post hoc (after rule P); frame centred on observed D14a.1,"
        " x = -dRA cos(dec), y = dDec; offsets are chires' 2-decimal values",
        "per_f": out,
    }


F2L_NLENS = 6  # PLAN_F2L.md: O1-O6 = input.par's limit 1-6
F2L_PERTURB = 0.03  # ASSUMPTION (plan): start offset, fraction of each parameter's range
F2L_SIMPLEX = 0.05  # ASSUMPTION (plan): initial simplex step in u
F2L_MAXFEV = 1500  # ASSUMPTION (plan)
F2L_WALL_S = 1500  # ASSUMPTION (plan)
F2L_INFEASIBLE = 1.0e6  # ASSUMPTION (plan): objective of an infeasible point
_KEY_BEST = {"ellipticity": "ellipticite"}  # input.par limit keyword -> best.par keyword


def free_params(input_par: str, n: int = F2L_NLENS) -> list[tuple[str, str, float, float]]:
    """(potential name O<i>, keyword, lo, hi) for every line of ``limit 1``-``n``; type 1 only."""
    out = []
    for i in range(1, n + 1):
        _, body = _limit(input_par, i)
        for ln in body:
            p = ln.split()
            if not p or p[0] == "end" or p[0].startswith("#"):
                continue
            if len(p) < 4:
                raise ValueError(f"limit {i}: short line {ln.strip()!r}")
            if p[1] != "1":
                raise ValueError(f"limit {i} {p[0]}: only type 1 (uniform) is handled")
            out.append((f"O{i}", p[0], float(p[2]), float(p[3])))
    return out


def get_param(par: str, name: str, key: str) -> float:
    """Value of ``key`` (input.par keyword) in ``potentiel <name>``."""
    lines = par.split("\n")
    a, b = _pot_body(par, name)
    k = _KEY_BEST.get(key, key)
    hits = [lines[i].split()[1] for i in range(a + 1, b) if lines[i].split()[:1] == [k]]
    if len(hits) != 1:
        raise ValueError(f"potentiel {name}: need one {k} line")
    return float(hits[0])


def set_param(par: str, name: str, key: str, value: float) -> str:
    """Set ``key`` (input.par keyword) in ``potentiel <name>``; a ``<key>_kpc`` companion line is
    removed so the arcsec value is the only one Lenstool reads."""
    lines = par.split("\n")
    a, b = _pot_body(par, name)
    k = _KEY_BEST.get(key, key)
    hits = [i for i in range(a + 1, b) if lines[i].split()[:1] == [k]]
    if len(hits) != 1:
        raise ValueError(f"potentiel {name}: need one {k} line")
    lines[hits[0]] = f"\t{k}     {value:.6f}"
    drop = {i for i in range(a + 1, b) if lines[i].split()[:1] == [k + "_kpc"]}
    return "\n".join(ln for i, ln in enumerate(lines) if i not in drop)


def objective(chires: str) -> tuple[float, bool]:
    """(objective, feasible): feasible iff no warning (lost image) and no proximity penalty."""
    chimul, _, nwarn = read_chires(chires)
    ok = nwarn == 0 and read_penalty(chires) == 0.0
    return (chimul if ok else F2L_INFEASIBLE), ok


def perturbed_start(u_best: np.ndarray, frac: float = F2L_PERTURB) -> np.ndarray:
    """u_best + frac * (+1, -1, +1, ...), clipped to [0.001, 0.999] (PLAN_F2L.md)."""
    s = np.where(np.arange(len(u_best)) % 2 == 0, 1.0, -1.0)
    return np.clip(u_best + frac * s, 0.001, 0.999)


def gate_v3(
    start: tuple[float, bool], best: float | None, trace: list[float], stopped: bool
) -> str:
    """Frozen gate V3 (PLAN_F2L.md). ``start`` = (chi2, feasible) of the perturbed start, ``best``
    the lowest feasible chi2 found, ``trace`` the running best objective per evaluation, ``stopped``
    true when the budget or the wall cap ended the run."""
    lim = V2_REF * (1.0 + TOL_FRAC)
    if not start[1] or start[0] <= lim:
        return "uninformative"
    if best is not None and best <= lim:
        return "pass"
    if stopped and len(trace) > 300 and (trace[-301] - trace[-1]) > 3.0:
        return "undetermined"
    return "fail"


F2L_RUN_TIMEOUT_S = (
    60  # deviation (review fix): one Lenstool evaluation slower than this is infeasible
)


def _f2l_setup(files: dict[str, Path]) -> tuple[str, list, np.ndarray, np.ndarray, np.ndarray]:
    """Stage V' run file, the free parameters, their bounds and best.par's u (PLAN_F2L.md)."""
    best, inp = files["best.par"].read_text(), files["input.par"].read_text()
    vp = build_vprime(build_v0(best, inp), best, inp)
    free = free_params(inp)
    lo = np.array([f[2] for f in free])
    hi = np.array([f[3] for f in free])
    p_best = np.array([get_param(vp, n, k) for n, k, _, _ in free])
    if np.any(p_best < lo) or np.any(p_best > hi):
        raise ValueError("best.par value outside its limit range")
    return vp, free, lo, hi, (p_best - lo) / (hi - lo)


def _par_at(vp: str, free: list, lo: np.ndarray, hi: np.ndarray, u: np.ndarray) -> str:
    """The run file with every free parameter at u (clipped to [0, 1])."""
    par = vp
    for (n, k, _, _), v in zip(free, lo + np.clip(u, 0, 1) * (hi - lo), strict=True):
        par = set_param(par, n, k, float(v))
    return par


def _f2l_eval(
    lt: Path, d: Path, par: str, arcs: Path
) -> tuple[subprocess.CompletedProcess | None, str]:
    """One Lenstool run; (None, "") when it fails, stalls or writes no chires.dat."""
    try:
        return run(lt, d, par, arcs, F2L_RUN_TIMEOUT_S)
    except (LenstoolError, FileNotFoundError, subprocess.TimeoutExpired):
        return None, ""


def _nelder_mead(
    lt: Path, wd: Path, files: dict[str, Path], u0: np.ndarray, step: float
) -> tuple[float, bool, list, bool, str, float]:
    """Bounded Nelder-Mead from u0 with initial simplex ``step`` (PLAN_F2L.md settings).
    Returns (start objective, start feasible, history, stopped, message, t0)."""
    from scipy.optimize import minimize

    vp, free, lo, hi, _ = _f2l_setup(files)
    hist: list[tuple[float, bool, np.ndarray]] = []
    cache: dict[bytes, float] = {}
    t0 = time.time()

    class _Wall(Exception):
        pass

    def f(u: np.ndarray) -> float:
        key = np.clip(u, 0, 1).tobytes()
        if key in cache:  # Nelder-Mead's first vertex repeats the start
            return cache[key]
        if time.time() - t0 > F2L_WALL_S:
            raise _Wall
        _, chires = _f2l_eval(lt, wd / "f2l", _par_at(vp, free, lo, hi, u), files["arcs.dat"])
        try:
            val, ok = objective(chires)
        except ValueError:
            val, ok = F2L_INFEASIBLE, False
        hist.append((val, ok, np.clip(u, 0, 1)))
        cache[key] = val
        return val

    start = f(u0)
    start_ok = hist[0][1]
    sim = np.vstack(
        [u0] + [u0 + step * np.eye(len(u0))[i] * (1 if u0[i] < 0.9 else -1) for i in range(len(u0))]
    )
    stopped = False
    try:
        r = minimize(
            f,
            u0,
            method="Nelder-Mead",
            bounds=[(0.0, 1.0)] * len(u0),
            options={"maxfev": F2L_MAXFEV, "xatol": 1e-4, "fatol": 0.01, "initial_simplex": sim},
        )
        stopped = not r.success
        message = str(r.message)
    except _Wall:
        stopped, message = True, "wall cap reached"
    return start, start_ok, hist, stopped, message, t0


def _nm_summary(lt: Path, files: dict[str, Path], hist: list, u_best: np.ndarray) -> dict:
    """Fields shared by the V3 and V3' result files."""
    free = free_params(files["input.par"].read_text())
    vals = [v for v, _, _ in hist]
    trace = list(np.minimum.accumulate(vals))
    feas = [(v, u) for v, ok, u in hist if ok]
    bestv, xbest = min(feas, key=lambda vu: vu[0]) if feas else (None, None)
    out = {
        "provenance": "model_prediction",
        "lenstool": version(lt),
        "omp_threads": _OMP or str(os.cpu_count() or 1),
        "n_free": len(free),
        "free_params": [f"{n}:{k}" for n, k, _, _ in free],
        "maxfev_assumption": F2L_MAXFEV,
        "wall_cap_s": F2L_WALL_S,
        "run_timeout_s": F2L_RUN_TIMEOUT_S,
        "chi2_best_feasible": bestv,
        "chi2_ref_vprime": V2_REF,
        "tol_frac_assumption": TOL_FRAC,
        "n_evaluations": len(hist),
        "n_infeasible": sum(not ok for _, ok, _ in hist),
        "running_best_every_100": [round(float(trace[i]), 3) for i in range(0, len(trace), 100)],
        "_trace": trace,
    }
    if xbest is not None:  # only a feasible best point has a deviation worth reporting
        dev = np.abs(xbest - u_best)
        j = int(np.argmax(dev))
        out["max_abs_du_from_best_par"] = round(float(dev[j]), 4)
        out["max_abs_du_param"] = out["free_params"][j]
    return out


def main_f2l_v3(lt: Path, wd: Path, files: dict[str, Path]) -> None:
    """Stage F2-L gate V3 (plan PLAN_F2L.md): perturbed start at f 1.00."""
    *_, u_best = _f2l_setup(files)
    start, start_ok, hist, stopped, message, t0 = _nelder_mead(
        lt, wd, files, perturbed_start(u_best), F2L_SIMPLEX
    )
    out = {"plan": "results/lens_consistency/d115b_lenstool_refit/PLAN_F2L.md"}
    out |= _nm_summary(lt, files, hist, u_best)
    trace = out.pop("_trace")
    out |= {
        "perturb_frac_assumption": F2L_PERTURB,
        "chi2_start": start if start_ok else None,
        "start_feasible": start_ok,
        "optimizer_message": message,
        "stopped_by_budget": stopped,
        "runtime_s": round(time.time() - t0, 1),
        "gate_V3": gate_v3((start, start_ok), out["chi2_best_feasible"], trace, stopped),
    }
    _write(out, "result_f2l_v3.json")
    print(json.dumps(out, indent=1))


V3P_LADDER = (0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02)  # ASSUMPTION (PLAN_V3P.md)


def choose_start(ladder_rows: list[tuple[float, float, bool]]) -> float | None:
    """PLAN_V3P.md: the first ladder frac whose start is feasible with chi2 above the tolerance.
    ``ladder_rows`` = (frac, objective, feasible) in ladder order."""
    lim = V2_REF * (1.0 + TOL_FRAC)
    for frac, val, ok in ladder_rows:
        if ok and val > lim:
            return frac
    return None


def gate_v3p(chosen: float | None, best: float | None, trace: list[float], stopped: bool) -> str:
    """Frozen gate V3' (PLAN_V3P.md): ``uninformative`` without a qualifying start, else the
    PLAN_F2L.md rule on the optimiser run."""
    if chosen is None:
        return "uninformative"
    # the chosen start is feasible and above the tolerance, so gate_v3's start test passes
    return gate_v3((float("inf"), True), best, trace, stopped)


def main_v3p(lt: Path, wd: Path, files: dict[str, Path]) -> None:
    """Stage V3' (plan PLAN_V3P.md): start ladder, then Nelder-Mead with step = chosen frac."""
    vp, free, lo, hi, u_best = _f2l_setup(files)
    rows = []
    for frac in V3P_LADDER:
        _, chires = _f2l_eval(
            lt,
            wd / "v3p_ladder",
            _par_at(vp, free, lo, hi, perturbed_start(u_best, frac)),
            files["arcs.dat"],
        )
        try:
            val, ok = objective(chires)
        except ValueError:
            val, ok = F2L_INFEASIBLE, False
        rows.append((frac, val, ok))
        if choose_start(rows) is not None:
            break
    chosen = choose_start(rows)
    out = {
        "plan": "results/lens_consistency/d115b_lenstool_refit/PLAN_V3P.md",
        "ladder": [{"frac": fr, "chi2": v if ok else None, "feasible": ok} for fr, v, ok in rows],
        "chosen_frac": chosen,
    }
    if chosen is None:
        out |= {"provenance": "model_prediction", "lenstool": version(lt)}
        out["gate_V3p"] = "uninformative"
    else:
        start, _, hist, stopped, message, t0 = _nelder_mead(
            lt, wd, files, perturbed_start(u_best, chosen), chosen
        )
        out |= _nm_summary(lt, files, hist, u_best)
        trace = out.pop("_trace")
        out |= {
            "chi2_start": start,
            "optimizer_message": message,
            "stopped_by_budget": stopped,
            "runtime_s": round(time.time() - t0, 1),
            "gate_V3p": gate_v3p(chosen, out["chi2_best_feasible"], trace, stopped),
        }
    _write(out, "result_v3p.json")
    print(json.dumps(out, indent=1))


def _warn_families(chires: str) -> dict[str, int]:
    """Family summary rows with nwarn > 0 (lost images)."""
    out = {}
    for r in (ln.split() for ln in chires.splitlines()):
        if len(r) >= 13 and r[0].isdigit() and r[10] == "N/A" and int(r[12]) > 0:
            out[r[1].rsplit(".", 1)[0]] = int(r[12])
    return out


def posthoc_f2l_v3(lt: Path, wd: Path, files: dict[str, Path]) -> dict:
    """Post hoc (after gate V3 was ``uninformative``): why the perturbed start is infeasible, and
    which single-parameter offsets of +-F2L_PERTURB from best.par stay feasible."""
    vp, free, lo, hi, u_best = _f2l_setup(files)

    def ev(u: np.ndarray) -> dict:
        p, chires = _f2l_eval(lt, wd / "f2l_post", _par_at(vp, free, lo, hi, u), files["arcs.dat"])
        if p is None or "chimul" not in chires:
            return {"feasible": False, "error": "lenstool failed or wrote no chimul"}
        _, ok = objective(chires)
        return {
            "feasible": ok,
            "chimul": read_chires(chires)[0],
            "penalty": read_penalty(chires),
            "penalty_images": penalty_images(p.stdout + p.stderr),
            "warn_families": _warn_families(chires),
        }

    single = {}
    for j, (n, k, _, _) in enumerate(free):
        for s in (1.0, -1.0):
            u = u_best.copy()
            u[j] = np.clip(u[j] + s * F2L_PERTURB, 0.001, 0.999)
            r = ev(u)
            r.pop("penalty_images", None)
            single[f"{n}:{k}:{'+' if s > 0 else '-'}"] = r
    return {
        "label": "post hoc (after gate V3 was uninformative)",
        "start": ev(perturbed_start(u_best)),
        "single_offsets": single,
        "n_single_feasible": sum(v["feasible"] for v in single.values()),
        "n_single": len(single),
    }


V3P_POST_STEPS = (0.001, 0.0003, 0.0001, 0.00003)  # post hoc: offset sizes, fraction of range
V3P_POST_PARAMS = ("O1:ellipticity", "O2:angle_pos", "O3:v_disp", "O5:cut_radius", "O6:v_disp")


def posthoc_v3p(lt: Path, wd: Path, files: dict[str, Path]) -> dict:
    """Post hoc (after gate V3' was ``uninformative``): the all-parameter start at 1e-4, single
    +-0.001 offsets, and the family chi2 changes of + offsets shrinking to 3e-5 of the range."""
    vp, free, lo, hi, u_best = _f2l_setup(files)
    names = [f"{n}:{k}" for n, k, _, _ in free]

    def ev(u: np.ndarray) -> dict:
        p, chires = _f2l_eval(lt, wd / "v3p_post", _par_at(vp, free, lo, hi, u), files["arcs.dat"])
        if p is None or "chimul" not in chires:
            return {"feasible": False, "error": "lenstool failed or wrote no chimul"}
        chimul, fam, nwarn = read_chires(chires)
        return {
            "feasible": objective(chires)[1],
            "chimul": chimul,
            "nwarn": nwarn,
            "penalty": read_penalty(chires),
            "fam": fam,
        }

    ref = ev(u_best)
    fam0 = ref.pop("fam")
    single = {}
    for j, nm in enumerate(names):
        for s in (1.0, -1.0):
            u = u_best.copy()
            u[j] = np.clip(u[j] + s * 0.001, 0.001, 0.999)
            r = ev(u)
            r.pop("fam", None)
            single[f"{nm}:{'+' if s > 0 else '-'}"] = r
    shrink = {}
    for nm in V3P_POST_PARAMS:
        for d in V3P_POST_STEPS:
            u = u_best.copy()
            j = names.index(nm)
            u[j] = np.clip(u[j] + d, 0.001, 0.999)
            r = ev(u)
            fam = r.pop("fam", {})
            r["family_dchi2_gt1"] = {
                k: round(v - fam0.get(k, 0.0), 2)
                for k, v in fam.items()
                if abs(v - fam0.get(k, 0.0)) > 1
            }
            shrink[f"{nm}:+{d}"] = r
    start = ev(perturbed_start(u_best, 0.0001))
    start.pop("fam", None)
    return {
        "label": "post hoc (after gate V3' was uninformative)",
        "best_par": ref,
        "start_frac_1e-4": start,
        "single_offsets_0.001": single,
        "n_single_feasible": sum(v["feasible"] for v in single.values()),
        "n_single": len(single),
        "shrinking_plus_offsets": shrink,
    }


VS_STEPS = (
    3e-5,
    1e-4,
    3e-4,
    1e-3,
    3e-3,
    1e-2,
)  # ASSUMPTION (PLAN_VS.md): offsets, fraction of range
VS_JUMP_ABS = 2.0  # ASSUMPTION (PLAN_VS.md): |delta chi2| floor of a jump
VS_JUMP_RATIO = 0.6  # ASSUMPTION (PLAN_VS.md): |D_k| >= ratio * |D_k+1| marks a step


def build_vs(vp: str) -> str:
    """PLAN_VS.md: the stage V' run file with its one ``forme -1`` line set to ``forme 0``."""
    lines = vp.split("\n")
    hits = [i for i, ln in enumerate(lines) if ln.split()[:2] == ["forme", "-1"]]
    if len(hits) != 1:
        raise ValueError("need exactly one 'forme -1' line")
    head = [ln for ln in lines[: hits[0]] if ln.strip() and not ln[0].isspace()][-1:]
    if head != ["image"]:
        raise ValueError("the 'forme -1' line is not in the image section")
    lines[hits[0]] = "    forme       0"
    return "\n".join(lines)


def series_vs(deltas: list[float | None]) -> str:
    """One (parameter, sign) series of chi2 changes in increasing step order (PLAN_VS.md):
    ``jump``, ``incomplete`` (a value missing, no jump) or ``smooth``."""
    for k in range(min(4, len(deltas) - 1)):  # k <= 4, i.e. steps up to 1e-3
        a, b = deltas[k], deltas[k + 1]
        if a is None or b is None:
            continue
        if abs(a) > VS_JUMP_ABS and abs(a) >= VS_JUMP_RATIO * abs(b):
            return "jump"
    if len(deltas) != len(VS_STEPS) or any(d is None for d in deltas):
        return "incomplete"
    return "smooth"


def gate_vs(
    vs0: float | None, series: dict[tuple[str, int], list[float | None]], params: list[str]
) -> str:
    """Frozen gate VS (PLAN_VS.md). ``series`` maps (parameter, sign) to its six chi2 changes
    (skipped signs are absent); ``params`` lists all free parameters."""
    if vs0 is None:
        return "undetermined"
    kinds = [series_vs(d) for d in series.values()]
    if "jump" in kinds:
        return "fail"
    covered = {p for p, _ in series}
    if "incomplete" in kinds or any(p not in covered for p in params):
        return "undetermined"
    return "pass"


def main_vs(lt: Path, wd: Path, files: dict[str, Path]) -> None:
    """Stage VS (plan PLAN_VS.md): VS0 at best.par, then single offsets of every O1-O6 parameter."""
    vp, free, lo, hi, u_best = _f2l_setup(files)
    vs = build_vs(vp)
    names = [f"{n}:{k}" for n, k, _, _ in free]
    t0 = time.time()

    def ev(u: np.ndarray, tag: str) -> float | None:
        par = _par_at(vs, free, lo, hi, u)
        try:  # chires.dat is optional here: the objective is the printed chi2
            p, _ = run(lt, wd / "vs" / tag, par, files["arcs.dat"], F2L_RUN_TIMEOUT_S, False)
            return read_chi2pos(p.stdout + p.stderr)
        except (LenstoolError, subprocess.TimeoutExpired, ValueError):
            return None

    vs0 = ev(u_best, "vs0")
    out = {
        "plan": "results/lens_consistency/d115b_lenstool_refit/PLAN_VS.md",
        "provenance": "model_prediction",
        "lenstool": version(lt),
        "chi2_vs0": vs0,
        "steps": list(VS_STEPS),
    }
    series: dict[tuple[str, int], list[float | None]] = {}
    if vs0 is not None:
        for j, nm in enumerate(names):
            for s in (1, -1):
                if not 0.0 <= u_best[j] + s * VS_STEPS[-1] <= 1.0:
                    continue
                vals = []
                for d in VS_STEPS:
                    u = u_best.copy()
                    u[j] += s * d
                    v = ev(u, "off")
                    vals.append(None if v is None else round(v - vs0, 4))
                series[(nm, s)] = vals
    out["series"] = {f"{p}:{'+' if s > 0 else '-'}": v for (p, s), v in series.items()}
    out["kind"] = {k: series_vs(v) for k, v in out["series"].items()}
    out["runtime_s"] = round(time.time() - t0, 1)
    out["gate_VS"] = gate_vs(vs0, series, names)
    _write(out, "result_vs.json")
    print(json.dumps(out, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lenstool", required=True, type=Path)
    ap.add_argument("--workdir", required=True, type=Path)
    ap.add_argument(
        "--stage",
        choices=(
            "v",
            "vprime-f",
            "d14a",
            "f2b",
            "f2b-posthoc",
            "f2l-v3",
            "f2l-v3-posthoc",
            "v3p",
            "v3p-posthoc",
            "vs",
        ),
        default="v",
    )
    a = ap.parse_args()
    lt, wd = a.lenstool.resolve(), a.workdir.resolve()
    files = lc.model_files(MODEL)
    if lc.chi2pos_from_par(files["best.par"]) != REF_CHI2POS:
        raise ValueError("best.par #Chi2pos differs from REF_CHI2POS")
    if a.stage == "vprime-f":
        main_vprime_f(lt, wd, files)
        return
    if a.stage == "d14a":
        main_d14a(lt, wd, files)
        return
    if a.stage == "f2b":
        main_f2b(lt, wd, files)
        return
    if a.stage == "f2l-v3-posthoc":
        _write(posthoc_f2l_v3(lt, wd, files), "posthoc_f2l_v3.json")
        return
    if a.stage == "f2l-v3":
        main_f2l_v3(lt, wd, files)
        return
    if a.stage == "v3p-posthoc":
        _write(posthoc_v3p(lt, wd, files), "posthoc_v3p.json")
        return
    if a.stage == "vs":
        main_vs(lt, wd, files)
        return
    if a.stage == "v3p":
        main_v3p(lt, wd, files)
        return
    if a.stage == "f2b-posthoc":
        _write(posthoc_f2b(wd), "posthoc_f2b.json")
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

"""SLICE JWST-era image lists vs HST-era RELICS Lenstool v1 maps on three clusters.

Frozen plan: ``results/lens_consistency/slice_relics_heldout/PLAN.md``. Images: SLICE Table 8
(clump 1 of each source; ``slice_images.ecsv``, observed). Per cluster (one model each): every
family at its fitted z (rule 1), gate K0 and frame check F on the families not held out (rules
2-3), reading A on the JWST-new families (rule 4) and reading B on the spectroscopic families at
their spectroscopic z (rule 5). Rule functions are those of ``macs0647_clash_heldout`` and
``macs2129_clash_heldout`` (same frozen rules). Writes ``summary.json`` and per-model tables.

    python scripts/slice_relics_heldout.py --build-images main.tex   # once, from the e-print
    python scripts/slice_relics_heldout.py [--grid-step 0.25]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import lens_consistency as lc
import macs0647_clash_heldout as ho
import macs2129_clash_heldout as h29
import numpy as np
from astropy.table import Table, vstack

OUT = Path("results/lens_consistency/slice_relics_heldout")
IMAGES = OUT / "slice_images.ecsv"
EPRINT_SHA256 = "f7120d94a5f1c302596838baa5353128ded977119923f53eddab5fa089dea9a7"
# cluster -> (Table 8 side-head, model, held-out families, spectroscopic families for reading B)
CLUSTERS: dict[str, tuple[str, str, tuple[str, ...], tuple[str, ...]]] = {
    "rxc0232": ("RXC\\,J0232.2$-$4420", "rxc0232-relics-lenstool", (), ("5", "7")),
    # family 1: Table 8 z 1.051, text 1.501 -> not used in reading B (PLAN.md)
    "rxc2211": (
        "RXC\\,J2211.7$-$0349",
        "rxc2211-relics-lenstool",
        tuple(str(i) for i in range(4, 11)),
        (),
    ),
    "macs0553": (
        "MACS\\,J0553.4$-$3342",
        "macs0553-relics-lenstool",
        tuple(str(i) for i in range(11, 22)),
        ("1", "5", "10", "13"),
    ),
}
_ROW = re.compile(r"^\s*(\d+)\.(\d+)\.(\d+)\s*&\s*([-\d.]+)\s*&\s*([-\d.]+)\s*&(.*?)&")
_ZSPEC = re.compile(r"^\s*(\d+\.\d+)\s*$")
_ZMODEL = re.compile(r"\$\s*(\d+\.\d+)_")


def parse_block(tex: str, sidehead: str) -> Table:
    """Clump-1 rows (``S.1.k``) of one Table 8 block; z from the clump's first row; ``z_spec``
    True for a plain number (spectroscopic), False for a model value with errors."""
    start = tex.index(f"\\sidehead {{{{{sidehead}}}}}")
    end = tex.index("\\hline\\hline", start)
    rows, z_of = [], {}
    for line in tex[start:end].splitlines():
        m = _ROW.match(line)
        if not m:
            continue
        src, clump, img, ra, dec, zcell = m.groups()
        key = (src, clump)
        if zcell.strip():
            zs, zm = _ZSPEC.match(zcell), _ZMODEL.search(zcell)
            if zs:
                z_of[key] = (float(zs.group(1)), True)
            elif zm:
                z_of[key] = (float(zm.group(1)), False)
            else:
                raise ValueError(f"unparsed z cell: {line!r}")
        if clump != "1":
            continue
        z, spec = z_of[key]
        rows.append((f"{src}.{img}", src, f"{src}.{clump}.{img}", float(ra), float(dec), z, spec))
    names = ("image_id", "system", "paper_id", "ra", "dec", "z", "z_spec")
    return Table(rows=rows, names=names)


def build_images(tex_path: Path) -> Table:
    raw = tex_path.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    tex = raw.decode()
    parts = []
    for cluster, (head, _, held, spec_fams) in CLUSTERS.items():
        t = parse_block(tex, head)
        t["cluster"] = cluster
        t["new_in_jwst"] = np.isin(np.asarray(t["system"], str), held)
        t["use_spec"] = np.asarray(t["z_spec"], bool) & np.isin(
            np.asarray(t["system"], str), spec_fams
        )
        parts.append(t)
    t = vstack(parts)
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        "Cerny et al. 2026, ApJ 1001, 60 (SLICE; arXiv:2503.17498v2), Table 8 blocks of "
        "RXC J0232.2-4420, "
        f"RXC J2211.7-0349 and MACS J0553.4-3342 (main.tex sha256 {sha}; e-print tarball sha256 "
        f"{EPRINT_SHA256}); clump 1 of each source only; z_spec True: spectroscopic, False: SLICE "
        "model-optimized; new_in_jwst and use_spec per the paper's text (PLAN.md)"
    )
    return t


def load_cluster(images: Table, cluster: str) -> Table:
    t = images[np.asarray(images["cluster"], str) == cluster]
    t["image_id"] = np.asarray(t["image_id"], str)
    t["system"] = np.asarray(t["system"], str)
    t["held_out"] = np.asarray(t["new_in_jwst"], bool)
    held = CLUSTERS[cluster][2]
    if set(t["system"][t["held_out"]]) != set(held):
        raise ValueError(f"{cluster}: new_in_jwst disagrees with the frozen held-out families")
    t["z_spec"] = np.asarray(t["use_spec"], bool)  # rule 5 rows (h29.run_model reads z_spec)
    t["a"] = np.ones(len(t))
    return t


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid-step", type=float, default=0.25)
    ap.add_argument(
        "--build-images", type=Path, help="SLICE e-print main.tex: write the image list"
    )
    args = ap.parse_args(argv)
    if args.build_images:
        t = build_images(args.build_images)
        t.write(IMAGES, format="ascii.ecsv", overwrite=True)
        print(len(t), "images ->", IMAGES)
        return 0
    images = Table.read(IMAGES, format="ascii.ecsv")
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "rules": {
            "K0_MAX_RMS_ARCSEC": ho.K0_MAX_RMS_ARCSEC,
            "K0_MIN_MATCHED": ho.K0_MIN_MATCHED,
            "MATCH_ARCSEC": ho.MATCH_ARCSEC,
            "FRAME_MAX_ARCSEC": ho.FRAME_MAX_ARCSEC,
            "z_grid": [ho.Z_STEPS, ho.Z_MIN_OFFSET, ho.Z_MAX],
            "grid_step": args.grid_step,
        },
    }
    for cluster, (_, model, _, _) in CLUSTERS.items():
        t = load_cluster(images, cluster)
        held = np.asarray(t["held_out"], bool)
        spec = np.asarray(t["z_spec"], bool)
        entry, ip, ds = h29.run_model(model, t, args.grid_step)
        usable = [model] if entry["usable"] else []
        ra = ho.readings({model: np.asarray(ip["dtheta_arcsec"], float)[held]}, usable)
        rb = h29.reading_b({model: ds}, spec, usable) if spec.any() else np.array([], str)
        entry.update(
            {
                "model": model,
                "images": len(t),
                "held_out_images": int(held.sum()),
                "spec_z_images": int(spec.sum()),
                "flags_A": [
                    str(i) for i, v in zip(t["image_id"][held], ra, strict=True) if v != "no_flag"
                ],
                "flags_B": [
                    str(i) for i, v in zip(t["image_id"][spec], rb, strict=True) if v != "no_flag"
                ],
            }
        )
        if not held.any():
            entry["held_out_A"] = None
        if not spec.any():
            entry["spec_z_B"] = entry["spec_fitted_z"] = None
        summary[cluster] = entry
        ip.meta["model"] = model
        lc._write(ip, OUT / f"imageplane_{cluster}.ecsv")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    for c in CLUSTERS:
        e = summary[c]
        print(
            c, e["usable"], e["gate_K0"], e["frame_F"], e["held_out_A"], e["spec_z_B"], sep="\n  "
        )
        print("  flags", e["flags_A"], e["flags_B"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""SMACS 0723: the reference-catalogue extra-image test (stages 1 and 2).

Frozen plan: ``results/lens_consistency/smacs_djaref/PLAN_1.md``. Rule 0 takes the systems whose
redshift ``input.par`` keeps fixed (no ``z_m_limit``) as spectroscopic; rule 1 runs gate K of
``m0417_images`` on them (``smacs0723-iclv2``, no frame offset); rule 2 reads only the FITS headers
of the DJA v7.4 photometry (HTTP range) and checks the columns El Gordo stage 2b used. Writes
``summary.json``. ``PLAN_2.md`` (``--stage 2``): El Gordo stage 2b's rules on GO 2736 F200W and
the DJA v7.4 reference, with families at the ``fitted_z`` grid edge not flagged (deviation D2).

    python scripts/smacs_djaref.py [--stage 2] [--cache data/cache/smacsx]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import a2744_uncoverref as uref
import elgordo_djaref as eg
import lens_consistency as lc
import m0417_canucsref as cr
import m0417_extra as mx
import m0417_images as mi
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

FIELD = "smacs0723"
MODEL = "smacs0723-iclv2"
OUT = Path("results/lens_consistency/smacs_djaref")
IMAGES = "mahler22_images.ecsv"
DJA_PHOT = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/smacs0723-grizli-v7.4-fix_phot.fits"
GRID_STEP = 0.25  # m0417_images default (PLAN_1.md rule 1)
# PLAN_1.md rule 2: the columns El Gordo stage 2b read (D-034 addenda p-q)
STAGE2_COLS = (
    *eg.SHAPE_COLS,
    "flux_auto",
    "mag_auto",
    "tot_corr",
    "flux_aper_1",
    "f200w_flux_aper_1",
    "f200w_fluxerr_aper_1",
    "f200w_mask_aper_1",
    "f200w_tot_corr",
)
STAGE2_CARDS = ("APER_1",)
# stage 2 (PLAN_2.md; ASSUMPTIONs)
OBS, TILE = (
    "o001",
    "t001",
)  # GO 2736 F200W tile containing all 60 images (PLAN_2.md footprint check)
S3_2736 = "s3://stpubdata/jwst/public/jw02736/L3/t/{0}/jw02736-{0}_{1}_nircam_clear-f200w_{2}"
FRAME = OUT / "frame"
Z_EDGE = 10.0  # m0417_images.fitted_z grid's upper edge (deviation D2)


def fixed_z_systems(input_par: str, systems) -> list[str]:
    """PLAN_1.md rule 0: systems without a ``z_m_limit`` line in ``input.par``."""
    free = set()
    for line in input_par.splitlines():
        f = line.split()
        if f and f[0] == "z_m_limit" and len(f) > 2:
            free.add(str(int(float(f[2]))))
    return sorted({str(s) for s in systems} - free, key=lambda s: (len(s), s))


def stage2_inputs(blob: bytes, total_bytes: int) -> dict:
    """PLAN_1.md rule 2 on the header bytes of the DJA photometry table."""
    heads = eg.header_cards(blob)
    tab = next((h for h in heads if h.get("XTENSION") == "BINTABLE"), None)
    if tab is None:
        return {"reachable": True, "bytes": total_bytes, "table_header": False, "passes": False}
    cols = {str(tab[f"TTYPE{i}"]).lower() for i in range(1, tab["TFIELDS"] + 1)}
    missing = [c for c in STAGE2_COLS if c not in cols]
    missing += [k for k in STAGE2_CARDS if k not in tab]
    return {
        "reachable": True,
        "bytes": total_bytes,
        "rows": int(tab["NAXIS2"]),
        "n_columns": len(cols),
        "missing": missing,
        "aper_1_pix": float(tab["APER_1"]) if "APER_1" in tab else None,
        "passes": not missing,
    }


def write_images(out: Path = OUT) -> tuple[Path, list[str]]:
    """The model's ``arcs.dat`` with ``z_spec`` from rule 0."""
    _, files, _ = lc.load_model(MODEL)
    t = lensmodel.load_lenstool_images(files["arcs.dat"])
    fixed = fixed_z_systems(Path(files["input.par"]).read_text(), t["system"])
    t["z_spec"] = np.isin(np.asarray(t["system"], str), fixed)
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{MODEL} arcs.dat (sha256 in lensmodel.SMACS0723_MAHLER22_ICLV2); z_spec = systems "
        "without z_m_limit in input.par (PLAN_1.md rule 0, ASSUMPTION)"
    )
    out.mkdir(parents=True, exist_ok=True)
    path = out / IMAGES
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path, fixed


def split_z_edge(flags: list[str], z_used: dict, fixed: list[str]) -> dict:
    """PLAN_2.md D2 / D3: flags (``system@ra,dec``) of families with ``z_used`` >= ``Z_EDGE`` are
    reported, not flagged; the rest are labelled ``fixed_z`` or ``fitted_z``."""
    out: dict = {"flags": [], "flags_z_edge": []}
    for f in flags:
        system = f.split("@", 1)[0]
        if z_used[system] >= Z_EDGE:
            out["flags_z_edge"].append(f)
        else:
            out["flags"].append({"flag": f, "z": "fixed_z" if system in fixed else "fitted_z"})
    return out


def stage_2(cache: Path) -> dict:
    """PLAN_2.md: El Gordo stage 2b's rules 2a / 2b, frame tables, DJA reference, local noise,
    gate C2 shared, on this field; then deviation D2 / D3 on the flags."""
    import s3fs

    t, sha, nbytes = eg.fetch_table(DJA_PHOT)
    unit, grid = eg.flux_unit_jy(t, with_tot_corr=True), eg.grid_rotation_deg(t)
    res: dict = {"dja": {"sha256": sha, "bytes": nbytes, "rows": len(t)}, "rule_2a": unit}
    res["rule_2b"] = grid
    if not (unit["passes"] and grid["passes"]):
        res["outcome"] = "stopped at rule 2a / 2b; nothing read"
        return res
    mx.FIELDS[FIELD] = (IMAGES, MODEL, OBS, TILE)
    mx.IMAGES, mx.OUT, mx.S3 = OUT, FRAME, S3_2736
    fs = s3fs.S3FileSystem(anon=True, default_block_size=64 * 1024)
    cache.mkdir(parents=True, exist_ok=True)
    shapes = cache / Path(S3_2736.format(OBS, TILE, "cat.ecsv")).name
    if not shapes.exists():
        fs.get(S3_2736.format(OBS, TILE, "cat.ecsv")[5:], str(shapes))
    FRAME.mkdir(parents=True, exist_ok=True)
    frame = mx.run_field(FIELD, cache)  # rule 1
    (FRAME / "summary.json").write_text(
        json.dumps(_clean({FIELD: frame}), indent=1, allow_nan=False) + "\n"
    )
    w, d, pix = uref._read_sci(fs, S3_2736.format(OBS, TILE, "i2d.fits")[5:])
    rc = eg.to_cat_dja(t, unit["unit_jy"], grid["pixel_scale_arcsec"])
    inputs = {
        "reference": Path(DJA_PHOT).name,
        "reference_sha256": sha,
        "n_rows_kept": int(len(rc["flux"])),
        "flux": f"{eg.BAND}_flux_aper_{eg.APER_IX} x flux_auto / flux_aper_{eg.APER_IX} x "
        f"{eg.BAND}_tot_corr",
        "noise_image": S3_2736.format(OBS, TILE, "i2d.fits"),
        "noise_pixel_scale_arcsec": pix,
    }
    depth = cr.depth_mag(rc["flux"], rc["snr"])
    src_x, src_y = (
        np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(t["ra"], t["dec"], unit="deg"))
    )
    ip = Table.read(OUT / f"imageplane_{FIELD}_{MODEL}.ecsv")
    z_used = dict(zip(np.asarray(ip["system"], str), np.asarray(ip["z_used"], float), strict=False))
    fixed = json.loads((OUT / "summary.json").read_text())["stage1_fixed_z_systems"]
    res["frame"] = frame
    for contain, out in ((False, OUT / "without_contain"), (True, OUT)):  # rule 3 / rule 7
        out.mkdir(parents=True, exist_ok=True)
        r = uref._reference_run(
            rc, depth, inputs, out, contain, d, w, src_x, src_y, pix, FIELD, FRAME
        )
        res["frozen" if contain else "without_contain"] = r | split_z_edge(
            r["flags"], z_used, fixed
        )
    return res


def sheet(path: Path, classes=("absent", "ambiguous", "confused"), half: float = 3.0) -> int:
    """Post hoc contact sheet: F200W stamps of the stage-2 controls in ``classes`` (red circle =
    ``search_arcsec``, title = system, class and arc reference)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = Table.read(OUT / f"control_{FIELD}.ecsv")
    c = c[np.isin(np.asarray(c["class_local"], str), classes)]
    stamp, close = lc.image_stamper(S3_2736.format(OBS, TILE, "i2d.fits"), half)
    n = max(len(c), 1)
    fig, axes = plt.subplots(1, n, figsize=(2.6 * n, 2.9), squeeze=False)
    try:
        for ax, r in zip(axes[0], c, strict=False):
            sci, _, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
            med = np.nanmedian(sci)
            mad = np.nanmedian(np.abs(sci - med)) or 1.0
            ext = (xx.max(), xx.min(), yy.min(), yy.max())
            ax.imshow(
                sci,
                origin="lower",
                cmap="gray_r",
                vmin=med - 3 * mad,
                vmax=med + 30 * mad,
                extent=ext,
            )
            ax.add_patch(plt.Circle((0, 0), float(r["search_arcsec"]), fill=False, color="red"))
            ax.set_title(f"{r['system']} {r['class_local']} arc={bool(r['arc_ref'])}", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
    finally:
        close()
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    return len(c)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage", choices=("1", "2", "sheet"), default="1")
    ap.add_argument("--cache", type=Path, default=Path("data/cache/smacsx"))
    args = ap.parse_args(argv)
    if args.stage == "sheet":
        print(sheet(OUT / "sheet_2.png"), "controls")
        return 0
    if args.stage == "2":
        summary = json.loads((OUT / "summary.json").read_text())
        summary["stage2"] = _clean(stage_2(args.cache))
        (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
        s2 = summary["stage2"]
        print(json.dumps({k: s2.get(k) for k in ("rule_2a", "rule_2b", "outcome")}, indent=1))
        if "frozen" in s2:
            keys = ("gate_C2_shared", "flags", "flags_z_edge", "target_classes")
            print(json.dumps({k: s2["frozen"][k] for k in keys}, indent=1))
        return 0
    _, fixed = write_images()
    # fitted_z's grid starts at z_lens; take it from best.par (0.39), as lc.MODELS has none
    lc.MODELS[MODEL].setdefault("z_lens", lc.load_model(MODEL)[0].z_lens)
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    summary = {"stage1_fixed_z_systems": fixed}
    summary["stage1_gate_K"] = mi.run_field(FIELD, OUT, GRID_STEP)
    blob, total = eg.fetch_head(DJA_PHOT)
    summary["stage1_reference_inputs"] = stage2_inputs(blob, total) | {"url": DJA_PHOT}
    summary["rules"] = {
        "K_MAX_RMS_ARCSEC": mi.K_MAX_RMS_ARCSEC,
        "K_MIN_MATCHED": mi.K_MIN_MATCHED,
        "MATCH_ARCSEC": mi.MATCH_ARCSEC,
        "grid_step_arcsec": GRID_STEP,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    g = summary["stage1_gate_K"][MODEL]["gate_K"]
    print(json.dumps({"gate_K": g, "ref": summary["stage1_reference_inputs"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

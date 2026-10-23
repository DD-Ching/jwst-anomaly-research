"""RXJ2129: the reference-catalogue extra-image test, stage 1 (gate K and reference inputs).

Frozen plan: ``results/lens_consistency/rxj2129_djaref/PLAN_1.md``. The Dhanasingham+2026 image list
(arXiv:2604.11882, Appendix B table "Systems_RA_and_dec") is transcribed below; rule 1 runs gate K
of ``m0417_images`` on it with the pinned JWST-era Lenstool maps (``rxj2129-lenstool-maps``, no
frame offset); rule 2 reads only the FITS headers of the DJA v7.0 photometry (HTTP range) and
checks SMACS rule 2's columns. Writes ``summary.json``. ``PLAN_2.md`` (``--stage 2``): El Gordo
stage 2b's rules on GO 2767 F200W and the DJA v7.0 reference, with the system-3 knot group's
overlapping targets joined for flags (deviation D1).

    python scripts/rxj2129_djaref.py [--stage 2] [--cache data/cache/rxj2129x]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import a2744_uncoverref as uref
import arcs_rxj2129
import elgordo_djaref as eg
import m0417_canucsref as cr
import m0417_extra as mx
import m0417_images as mi
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from m0417_canucsref_posthoc import _clean
from smacs_djaref import stage2_inputs

from jwst_anomaly import lensmodel

FIELD = "rxj2129"
MODEL = "rxj2129-lenstool-maps"
OUT = Path("results/lens_consistency/rxj2129_djaref")
IMAGES = "dhanasingham26_images.ecsv"
DJA_PHOT = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/rxj2129-grizli-v7.0-fix_phot.fits"
GRID_STEP = 0.25  # m0417_images default (PLAN_1.md rule 1)
PAPER = "Dhanasingham et al. 2026, arXiv:2604.11882, Appendix B table Systems_RA_and_dec"

# (paper id, RA deg, Dec deg); observed, JWST NIRCam frame aligned to Gaia DR3 (paper Sect. 6)
_TABLE = {
    "1": ("0.6786", (("1a", 322.4149230, 0.0904156), ("1b", 322.4151893, 0.0889741),
                     ("1c", 322.4166530, 0.0867530))),
    "2": ("0.9160", (("2a", 322.4146395, 0.0923864), ("2b", 322.4162941, 0.0881023),
                     ("2c", 322.4165862, 0.0877771))),
    "3": ("1.5194", (("3a", 322.4159519, 0.0915111), ("3b", 322.4173436, 0.0906621),
                     ("3c", 322.4169561, 0.0903403), ("3d", 322.4185684, 0.0849278))),
    "4": ("1.5202", (("4a", 322.4155685, 0.0921535), ("4b", 322.4174895, 0.0901781),
                     ("4c", 322.4184261, 0.0853785))),
    "5": ("1.5210", (("5a", 322.4179753, 0.0932720), ("5b", 322.4201821, 0.0897722),
                     ("5c", 322.4204071, 0.0883142))),
    "6": ("3.0815", (("6a", 322.4137586, 0.0941984), ("6b", 322.4167395, 0.0877614),
                     ("6c", 322.4169923, 0.0873924))),
    "7": ("3.4270", (("7a", 322.4137266, 0.0920863), ("7b", 322.4144323, 0.0886347),
                     ("7c", 322.4175405, 0.0838683))),
    "G": ("9.5100", (("G1", 322.4202364, 0.0936441), ("G2", 322.4215532, 0.0916940),
                     ("G3", 322.4224106, 0.0889363))),
    "S3-1": ("1.5194", (("S3-1a", 322.4163448, 0.0912444), ("S3-1b", 322.4161198, 0.0911278),
                        ("S3-1c", 322.4163198, 0.0907945))),
    "S3-2": ("1.5194", (("S3-2a", 322.4169115, 0.0905028), ("S3-2b", 322.4167364, 0.0904193),
                        ("S3-2c", 322.4165032, 0.0905862))),
}  # fmt: skip
# stage 2 (PLAN_2.md; ASSUMPTIONs)
OBS, TILE = "o002", "t001"  # the GO 2767 F200W tile containing all 31 images (PLAN_2.md)
S3_2767 = "s3://stpubdata/jwst/public/jw02767/L3/t/{0}/jw02767-{0}_{1}_nircam_clear-f200w_{2}"
FRAME = OUT / "frame"
KNOT_GROUP = ("3", "S3-1", "S3-2")  # one source galaxy (PLAN_2.md deviation D1)
_KNOT_ORDER = (("recovered", "confused"), ("ambiguous",), ("absent",))


def image_table() -> Table:
    """PLAN_1.md inputs: the paper table with ids ``<family>.<n>`` and every row at its family z."""
    rows = []
    for fam, (z, imgs) in _TABLE.items():
        for n, (pid, ra, dec) in enumerate(imgs, 1):
            rows.append((f"{fam}.{n}", fam, pid, ra, dec, 1.0, float(z), True))
    t = Table(rows=rows, names=("image_id", "system", "paper_id", "ra", "dec", "a", "z", "z_spec"))
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{PAPER} (S3-1 / S3-2 at system 3's z and '--' rows at their family z, as the paper "
        "states; SN 2022riv row left out)"
    )
    return t


def merge_knot_targets(systems, ra, dec, search, classes, group=KNOT_GROUP) -> dict:
    """PLAN_2.md D1: join targets of different ``group`` families whose search circles overlap
    (transitively); a joined set takes the first member class in ``_KNOT_ORDER``. Returns the
    flags (``system@ra,dec``; a joined set is named by its first row) and the joined row sets."""
    n = len(systems)
    parent = list(range(n))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    cosd = np.cos(np.deg2rad(np.asarray(dec, float)))
    for i in range(n):
        for j in range(i + 1, n):
            if systems[i] == systems[j] or not {systems[i], systems[j]} <= set(group):
                continue
            dx = (ra[i] - ra[j]) * 3600.0 * 0.5 * (cosd[i] + cosd[j])
            sep = float(np.hypot(dx, (dec[i] - dec[j]) * 3600.0))
            if sep < search[i] + search[j]:
                parent[root(j)] = root(i)
    sets: dict[int, list[int]] = {}
    for i in range(n):
        sets.setdefault(root(i), []).append(i)
    flags, joined = [], []
    for g in sets.values():
        cls = [str(classes[k]) for k in g]
        merged = cls[0]
        if len(g) > 1:
            joined.append({"rows": g, "systems": [str(systems[k]) for k in g], "classes": cls})
            for want in _KNOT_ORDER:
                hit = [c for c in cls if c in want]
                if hit:
                    merged = hit[0]
                    break
            joined[-1]["class"] = merged
        if merged == "absent":
            k = g[0]
            flags.append(f"{systems[k]}@{float(ra[k]):.6f},{float(dec[k]):.6f}")
    return {"flags": flags, "knot_sets": joined}


def stage_2(cache: Path) -> dict:
    """PLAN_2.md: El Gordo stage 2b's rules 2a / 2b, frame tables, DJA reference, local noise,
    gate C2 shared, on this field; then deviation D1 on the flags."""
    import s3fs

    t, sha, nbytes = eg.fetch_table(DJA_PHOT)
    unit, grid = eg.flux_unit_jy(t, with_tot_corr=True), eg.grid_rotation_deg(t)
    res: dict = {"dja": {"sha256": sha, "bytes": nbytes, "rows": len(t)}, "rule_2a": unit}
    res["rule_2b"] = grid
    if not (unit["passes"] and grid["passes"]):
        res["outcome"] = "stopped at rule 2a / 2b; nothing read"
        return res
    arcs_rxj2129.stage_members(lensmodel.RXJ2129_LENSTOOL_MAPS)
    mx.FIELDS[FIELD] = (IMAGES, MODEL, OBS, TILE)
    mx.IMAGES, mx.OUT, mx.S3 = OUT, FRAME, S3_2767
    fs = s3fs.S3FileSystem(anon=True, default_block_size=64 * 1024)
    cache.mkdir(parents=True, exist_ok=True)
    shapes = cache / Path(S3_2767.format(OBS, TILE, "cat.ecsv")).name
    if not shapes.exists():
        fs.get(S3_2767.format(OBS, TILE, "cat.ecsv")[5:], str(shapes))
    FRAME.mkdir(parents=True, exist_ok=True)
    frame = mx.run_field(FIELD, cache)  # rule 1
    (FRAME / "summary.json").write_text(
        json.dumps(_clean({FIELD: frame}), indent=1, allow_nan=False) + "\n"
    )
    w, d, pix = uref._read_sci(fs, S3_2767.format(OBS, TILE, "i2d.fits")[5:])
    rc = eg.to_cat_dja(t, unit["unit_jy"], grid["pixel_scale_arcsec"])
    inputs = {
        "reference": Path(DJA_PHOT).name,
        "reference_sha256": sha,
        "n_rows_kept": int(len(rc["flux"])),
        "flux": f"{eg.BAND}_flux_aper_{eg.APER_IX} x flux_auto / flux_aper_{eg.APER_IX} x "
        f"{eg.BAND}_tot_corr",
        "noise_image": S3_2767.format(OBS, TILE, "i2d.fits"),
        "noise_pixel_scale_arcsec": pix,
    }
    depth = cr.depth_mag(rc["flux"], rc["snr"])
    src_x, src_y = (
        np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(t["ra"], t["dec"], unit="deg"))
    )
    res["frame"] = frame
    for contain, out in ((False, OUT / "without_contain"), (True, OUT)):  # rule 3 / rule 7
        out.mkdir(parents=True, exist_ok=True)
        r = uref._reference_run(
            rc, depth, inputs, out, contain, d, w, src_x, src_y, pix, FIELD, FRAME
        )
        pred = Table.read(out / f"predicted_{FIELD}.ecsv")
        k = merge_knot_targets(
            [str(x) for x in pred["system"]],
            np.asarray(pred["ra"], float),
            np.asarray(pred["dec"], float),
            np.asarray(pred["search_arcsec"], float),
            [str(x) for x in pred["class_local"]],
        )
        r["flags_rule_6_per_row"] = r.pop("flags")
        r["flags"] = k["flags"] if r["gate_C2_shared"]["passes"] else []
        r["knot_sets"] = k["knot_sets"]
        res["frozen" if contain else "without_contain"] = r
    return res


def stage_1() -> int:
    blob, total = eg.fetch_head(DJA_PHOT)  # rule 2 first: cheap, and fails before the solve
    OUT.mkdir(parents=True, exist_ok=True)
    t = image_table()
    t.write(OUT / IMAGES, format="ascii.ecsv", overwrite=True)
    arcs_rxj2129.stage_members(lensmodel.RXJ2129_LENSTOOL_MAPS)
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    summary: dict = {
        "images": {"n": len(t), "families": len(set(t["system"])), "source": PAPER},
        "stage1_gate_K": mi.run_field(FIELD, OUT, GRID_STEP),
        "stage1_reference_inputs": stage2_inputs(blob, total) | {"url": DJA_PHOT},
    }
    summary["rules"] = {
        "K_MAX_RMS_ARCSEC": mi.K_MAX_RMS_ARCSEC,
        "K_MIN_MATCHED": mi.K_MIN_MATCHED,
        "MATCH_ARCSEC": mi.MATCH_ARCSEC,
        "grid_step_arcsec": GRID_STEP,
    }
    ip = Table.read(OUT / f"imageplane_{FIELD}_{MODEL}.ecsv")
    d = np.asarray(ip["dtheta_arcsec"], float)
    summary["per_family_rms_arcsec"] = {
        fam: mi.rms_of(d[np.asarray(ip["system"], str) == fam])
        for fam in dict.fromkeys(np.asarray(ip["system"], str))
    }
    (OUT / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    g = summary["stage1_gate_K"][MODEL]["gate_K"]
    print(json.dumps({"gate_K": g, "ref": summary["stage1_reference_inputs"]}, indent=1))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage", choices=("1", "2"), default="1")
    ap.add_argument("--cache", type=Path, default=Path("data/cache/rxj2129x"))
    args = ap.parse_args(argv)
    if args.stage == "1":
        return stage_1()
    summary = json.loads((OUT / "summary.json").read_text())
    summary["stage2"] = _clean(stage_2(args.cache))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    s2 = summary["stage2"]
    print(json.dumps({k: s2.get(k) for k in ("rule_2a", "rule_2b", "outcome")}, indent=1))
    if "frozen" in s2:
        keys = ("gate_C2_shared", "flags", "knot_sets", "target_classes")
        print(json.dumps({k: s2["frozen"][k] for k in keys}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""MACS J1931 stage 3d: frame check F restricted to images with measurable F200W light at their
position (continuum criterion), validated on MACS0417 first; builder as stage V2 / 3b.

Frozen plan: ``results/lens_consistency/macs1931_stage3d/PLAN.md``.

    python scripts/m1931_stage3d.py validate          # rule 3d-2: F_c on MACS0417
    python scripts/m1931_stage3d.py f                 # rule 3d-3: F_c on MACS J1931
    python scripts/m1931_stage3d.py run --maps DIR    # rule 3d-4 (only after F_c passed)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import bullet_extra as be
import m0417_canucsref_posthoc as cp
import m0417_frame as mf
import m1423_canucsref as m14
import m1931_deblend as md
import m1931_stage3 as s3
import m1931_stage3c as s3c
import macs1931_gatek as mg
import numpy as np
from astropy.table import Table

OUT_VAL = Path("results/lens_consistency/macs0417_ownref_3d")
OUT = Path("results/lens_consistency/macs1931_stage3d")
CONT_AP_ARCSEC = 0.3  # ASSUMPTION (rule 3d-1), = the local-noise rule's aperture
CONT_MIN_SNR = 5.0  # ASSUMPTION (rule 3d-1), = m1423_canucsref.MIN_SIGMA
FC_MIN_N = 10  # ASSUMPTION (F_c)
FC_MIN_FRACTION = 0.5  # ASSUMPTION (F_c)
OFF_SECTION = "off_section"


def continuum_snr(data, w, ra, dec, src_x, src_y, pix: float, seed: int = m14.SEED) -> list:
    """Rule 3d-1: aperture sum / local noise at each position; NaN when the aperture leaves the
    section or holds a NaN pixel (``off_section``). One generator for all positions, in order."""
    rng = np.random.default_rng(seed)
    x, y = (np.asarray(v, float) for v in w.world_to_pixel_values(ra, dec))
    out = []
    for xi, yi in zip(x, y, strict=True):
        f = cp.aperture_sum(data, xi, yi, CONT_AP_ARCSEC / pix)
        if not np.isfinite(f):
            out.append(float("nan"))
            continue
        out.append(f / m14.local_noise(data, xi, yi, src_x, src_y, rng, pix))
    return out


def continuum_class(snr: float) -> str:
    """Rule 3d-1 per image: ``continuum``, ``faint`` or ``off_section``."""
    if not np.isfinite(snr):
        return OFF_SECTION
    return "continuum" if snr >= CONT_MIN_SNR else "faint"


def matched_each(ra, dec, cra, cdec) -> np.ndarray:
    """Per image: a catalogue row within the F matcher's radius (``F_MATCH_ARCSEC``)."""
    from astropy.coordinates import SkyCoord

    g = SkyCoord(np.asarray(ra, float), np.asarray(dec, float), unit="deg")
    c = SkyCoord(np.asarray(cra, float), np.asarray(cdec, float), unit="deg")
    return g.match_to_catalog_sky(c)[1].arcsec <= be.F_MATCH_ARCSEC


def frame_check_fc(ra, dec, snr, cra, cdec) -> dict:
    """F_c: the stage 3 matcher on the continuum-detected subset, every condition combined with
    "and": n_c >= FC_MIN_N, n_m >= FC_MIN_FRACTION * n_c, |median dx|, |median dy| <= 0.1",
    and the stage 3c chance control on the same subset."""
    cls = np.array([continuum_class(s) for s in snr])
    m = cls == "continuum"
    ra_c, dec_c = np.asarray(ra, float)[m], np.asarray(dec, float)[m]
    n_c = int(m.sum())
    f = be.frame_check_f(ra_c, dec_c, cra, cdec, 0) if n_c else {"n_matched": 0, "passes": False}
    f["passes_matcher_medians"] = bool(f["passes"])
    sh = s3c.shifted_counts(ra_c, dec_c, cra, cdec) if n_c else [0] * s3c.CHANCE_N_DIRECTIONS
    f["chance_shifted_counts"] = sh
    f["chance_mean"] = float(np.mean(sh))
    f["chance_passes"] = s3c.chance_passes(f["n_matched"], sh)
    f["n_continuum"] = n_c
    f["continuum_classes"] = {c: int((cls == c).sum()) for c in sorted(set(cls))}
    f["fc_rules"] = {
        "aperture_arcsec": CONT_AP_ARCSEC,
        "min_snr": CONT_MIN_SNR,
        "min_n": FC_MIN_N,
        "min_fraction": FC_MIN_FRACTION,
    }
    f["passes"] = bool(
        n_c >= FC_MIN_N
        and f["n_matched"] >= FC_MIN_FRACTION * n_c
        and f["passes_matcher_medians"]
        and f["chance_passes"]
    )
    return f


def use_fc(stored: dict):
    """Rule 3d-4 wiring: ``m1931_stage3.run``'s rule 5 outcome becomes the stored F_c outcome,
    after checking that ``run`` rebuilt the same catalogue (row count) as rule 3d-3."""
    want = int(stored["catalogue"]["n_rows_kept"])

    def extra_f(f: dict, images, cat) -> None:
        if len(cat["ra"]) != want:
            raise ValueError(f"catalogue has {len(cat['ra'])} rows, rule 3d-3 had {want}")
        f["passes_stage3_rule"] = f["passes"]
        f["frame_F_c"] = stored["frame_F_c"]
        f["passes"] = bool(stored["frame_F_c"]["passes"])

    return extra_f


def _fc_on(images: Table, sci, err, w, box, uri: str) -> tuple[dict, Table]:
    cat, meta, _, data, (sx, sy) = md.build(sci, err, w)
    snr = continuum_snr(data, w, images["ra"], images["dec"], sx, sy, meta["pixel_scale_arcsec"])
    f = frame_check_fc(images["ra"], images["dec"], snr, cat["ra"], cat["dec"])
    rows = Table(
        {
            "image_id": np.asarray(images["image_id"], str),
            "ra": np.asarray(images["ra"], float),
            "dec": np.asarray(images["dec"], float),
            "continuum_snr": snr,
            "continuum_class": [continuum_class(s) for s in snr],
            "matched_0p5": matched_each(images["ra"], images["dec"], cat["ra"], cat["dec"]),
        }
    )
    cat_info = {"image": uri, "section_y0_y1_x0_x1": [int(v) for v in box], **meta}
    return {"catalogue": cat_info, "frame_F_c": f}, rows


def validate() -> dict:
    """Rule 3d-2: F_c on the MACS0417 gold images (measured frame offset) on the V2 section."""
    off = json.loads((mf.OUT / "summary.json").read_text())["measurement"]["offset"]
    images = mf.shifted_loader(off["dra"], off["ddec"])("macs0417")
    uri, (sci, err, w, box) = md._section()
    s, rows = _fc_on(images, sci, err, w, box, uri)
    s["frame_offset_arcsec"] = [off["dra"], off["ddec"]]
    s["n_images"] = len(images)
    return s, rows


def field_f() -> dict:
    """Rule 3d-3: F_c on the 51 used MACS J1931 images on the 3b section."""
    import s3fs

    images = mg.load_images()
    fr = s3.OUT_3B / "frame"
    table = Table.read(fr / f"predicted_{s3.FIELD}.ecsv")
    ctrl = Table.read(fr / f"control_{s3.FIELD}.ecsv")
    pos = s3.section_positions(table, ctrl, images, True)
    fs = s3fs.S3FileSystem(anon=True, default_block_size=64 * 1024)
    uri = s3.S3_6882.format("i2d.fits")
    sci, err, w, box = md.read_section(fs, uri[5:], pos["ra"], pos["dec"])
    s, rows = _fc_on(images, sci, err, w, box, uri)
    s["n_images"] = len(images)
    return s, rows


def _write(out: Path, s: dict, rows: Table) -> None:
    out.mkdir(parents=True, exist_ok=True)
    s["plan"] = str(OUT / "PLAN.md")
    (out / "summary.json").write_text(json.dumps(cp._clean(s), indent=1, allow_nan=False) + "\n")
    rows.write(out / "images_fc.ecsv", overwrite=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("step", choices=["validate", "f", "run"])
    ap.add_argument("--maps", type=Path)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/m1931s3"))
    a = ap.parse_args(argv)
    if a.step == "validate":
        s, rows = validate()
        _write(OUT_VAL, s, rows)
        print(json.dumps(cp._clean(s["frame_F_c"]), indent=1))
        return 0
    val = json.loads((OUT_VAL / "summary.json").read_text())
    if not val["frame_F_c"]["passes"]:
        ap.error("rule 3d-2 validation failed on MACS0417: stage 3d stops (PLAN.md)")
    if a.step == "f":
        s, rows = field_f()
        _write(OUT, s, rows)
        print(json.dumps(cp._clean(s["frame_F_c"]), indent=1))
        return 0
    f = json.loads((OUT / "summary.json").read_text())
    if not f["frame_F_c"]["passes"]:
        ap.error("rule 3d-3 F_c failed on MACS J1931: no reading (PLAN.md)")
    if a.maps is None:
        ap.error("--maps is required for run")
    s = s3.run(a.maps, a.cache, OUT, True, md.BKG_BOX, use_fc(f))
    s["plan"] = str(OUT / "PLAN.md")
    s = cp._clean(s)
    (OUT / "run_summary.json").write_text(json.dumps(s, indent=1, allow_nan=False) + "\n")
    keys = ("rule2_gate_K_wls", "rule5_frame_F", "stopped")
    print(json.dumps({k: s[k] for k in keys if k in s}, indent=1))
    if "frozen" in s:
        print(json.dumps({k: s["frozen"][k] for k in ("gate_C2_shared", "flags")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

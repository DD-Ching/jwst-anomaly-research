"""Arc orientations in a second band: F444W shapes for the stored F200W ``arcs`` sources.

Lensing is achromatic, so a stretch set by the lens points the same way in F200W and F444W; a
blend of two differently coloured objects or a clumpy star-forming shape need not. For every
source of a stored ``lens_consistency.py arcs`` table (``results/arcs/<model>/arcs.ecsv``) the
nearest F444W pipeline-catalogue source within ``--match-arcsec`` supplies a second axis. The
F200W position and photo-z are kept; ``orientation_table`` re-classes the source with the F444W
axis against the same model. A source counts only if its F444W shape passes ``--min-e444`` and
``--min-snr444`` (ASSUMPTIONS: the F444W PSF is twice as wide, so ellipticities shrink).

Known-case check: background sources aligned in F200W must stay aligned in F444W; the fraction
that does is the test's own calibration. Outputs ``results/arcs_bands/<model>/bands.ecsv`` and
``summary.json`` (derived; the catalogues are observed).

    python scripts/arcs_bands.py --model abell370-canucs
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import lensmodel, paths, schema
from jwst_anomaly.photometry import fetch_catalog

_DIR = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("lens_consistency", _DIR / "lens_consistency.py")
lc = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("lens_consistency", lc)
_spec.loader.exec_module(lc)

S3_HTTPS = "https://stpubdata.s3.amazonaws.com/"

# F444W pipeline catalogue of the same observation as each stored F200W run, and its manifest.
FIELDS = {
    "abell2744-bergamini23": ("abell2744", "jw02561-o001_t003_nircam_clear-f444w_cat.ecsv"),
    "elgordo-caminha23": ("elgordo", "jw01176-o241_t012_nircam_clear-f444w_cat.ecsv"),
    "macs0416-canucs": ("macs0416", "jw01208-o004_t002_nircam_clear-f444w_cat.ecsv"),
    "macs1149-cats": ("macs1149", "jw01208-o008_t004_nircam_clear-f444w_cat.ecsv"),
    "abell370-canucs": ("abell370", "jw01208-o002_t001_nircam_clear-f444w_cat.ecsv"),
}


def manifest_sha(manifest: str, filename: str) -> str:
    t = Table.read(paths.repo_root() / "data" / "manifests" / f"{manifest}.ecsv")
    rows = t[t["productFilename"] == filename]
    if len(rows) != 1:
        raise SystemExit(f"error: {filename} not once in {manifest}.ecsv")
    return str(rows["sha256"][0])


def s3_url(filename: str) -> str:
    prog, obs = filename[2:7], filename.split("-")[1].split("_")[0]
    return f"{S3_HTTPS}jwst/public/jw{prog}/L3/t/{obs}/{filename}"


def match_second_band(stored: Table, shapes: Table, match_arcsec: float) -> Table:
    """``stored`` rows with the nearest ``shapes`` axis, ellipticity and S/N as ``*_444``."""
    c1 = SkyCoord(stored["ra"], stored["dec"], unit="deg")
    c2 = SkyCoord(shapes["ra"], shapes["dec"], unit="deg")
    idx, sep, _ = c1.match_to_catalog_sky(c2)
    ok = sep.to_value(u.arcsec) <= match_arcsec
    out = stored.copy(copy_data=True)
    for col, src in (("pa_444", "pa_obs"), ("e_444", "ellipticity"), ("snr_444", "snr")):
        out[col] = np.where(ok, np.asarray(shapes[src], float)[idx], np.nan)
    out["sep_444"] = np.where(ok, sep.to_value(u.arcsec), np.nan)
    out["dpa_bands"] = lensmodel.axis_offset_deg(out["pa_obs"], out["pa_444"])
    return out


def reclassify(model, rows: Table, args) -> Table:
    """``orientation_table`` with the F444W axis; returns ``class_444`` beside the F200W class."""
    keep = ["label", "ra", "dec", "z_phot", "z160", "z840"]
    src = Table({c: rows[c] for c in keep if c in rows.colnames})
    src["pa_obs"] = rows["pa_444"]
    t = lc.orientation_table(
        model, src, aligned_deg=args.aligned_deg, anti_deg=args.anti_deg, z_margin=args.z_margin
    )
    out = rows.copy(copy_data=True)
    out["class_444"] = t["orientation_class"]
    out["offset_min_444"] = t["offset_min"]
    return out


def summarize(t: Table, args) -> dict:
    test = t[t["strong_shear"] & (t["population"] == "background")]
    usable = np.isfinite(test["pa_444"]) & (test["e_444"] >= args.min_e444)
    usable &= test["snr_444"] >= args.min_snr444
    u_rows = test[usable]
    al = u_rows[u_rows["orientation_class"] == "aligned"]
    anti_200 = test[test["orientation_class"] == "anti"]
    sub = u_rows.copy()
    sub["orientation_class"] = sub["class_444"]
    return {
        "background_strong_shear": len(test),
        "usable_444": len(u_rows),
        "dpa_bands_median": round(float(np.median(u_rows["dpa_bands"])), 1)
        if len(u_rows)
        else None,
        "control_aligned_200": len(al),
        "control_stay_aligned_444": int(np.sum(al["class_444"] == "aligned")),
        "control_anti_444": int(np.sum(al["class_444"] == "anti")),
        "classes_444": lc.class_stats(sub, args.aligned_deg),
        "anti_expected_444": lc.anti_expected(sub, args.sigma_eps, args.min_e444, args.anti_deg),
        "anti_200": [
            {
                "label": int(r["label"]),
                "class_444": str(r["class_444"]) if np.isfinite(r["pa_444"]) else "no_match",
                "dpa_bands": None
                if not np.isfinite(r["dpa_bands"])
                else round(float(r["dpa_bands"]), 1),
                "e_444": None if not np.isfinite(r["e_444"]) else round(float(r["e_444"]), 2),
                "snr_444": None if not np.isfinite(r["snr_444"]) else round(float(r["snr_444"]), 1),
            }
            for r in anti_200
        ],
        "anti_444": [int(v) for v in sub["label"][(sub["orientation_class"] == "anti")]],
    }


def run(args) -> dict:
    stored = Table.read(args.arcs_dir / args.model / "arcs.ecsv")
    manifest, filename = FIELDS[args.model]
    cat = fetch_catalog(s3_url(filename), manifest_sha(manifest, filename))
    model, _, _ = lc.load_model(args.model)
    lc.apply_frame_offset(args.model, model)
    t = reclassify(model, match_second_band(stored, lc.load_shapes(cat), args.match_arcsec), args)
    t.meta.update(provenance=schema.Provenance.DERIVED.value, second_band=filename)
    summary = {
        "model": args.model,
        "second_band": filename,
        "assumptions": {
            k: getattr(args, k)
            for k in (
                "match_arcsec",
                "min_e444",
                "min_snr444",
                "aligned_deg",
                "anti_deg",
                "z_margin",
                "sigma_eps",
            )
        },
    } | summarize(t, args)
    out = args.out / args.model
    out.mkdir(parents=True, exist_ok=True)
    cols = [
        "label",
        "ra",
        "dec",
        "population",
        "strong_shear",
        "g_min",
        "orientation_class",
        "class_444",
        "pa_obs",
        "pa_444",
        "dpa_bands",
        "ellipticity",
        "e_444",
        "snr",
        "snr_444",
        "sep_444",
    ]
    t[cols].write(out / "bands.ecsv", overwrite=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", required=True, choices=sorted(FIELDS))
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--out", type=Path, default=paths.repo_root() / "results" / "arcs_bands")
    p.add_argument("--match-arcsec", type=float, default=0.3)
    p.add_argument("--min-e444", type=float, default=0.3)
    p.add_argument("--min-snr444", type=float, default=10.0)
    p.add_argument("--aligned-deg", type=float, default=30.0)
    p.add_argument("--anti-deg", type=float, default=60.0)
    p.add_argument("--z-margin", type=float, default=0.1)
    p.add_argument("--sigma-eps", type=float, default=0.26)
    args = p.parse_args(argv)
    print(json.dumps(run(args), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

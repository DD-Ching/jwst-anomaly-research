"""List Lenstool model potentials whose spectroscopic redshift is off the cluster plane (D-042).

A single-plane cluster model puts every member potential at the cluster redshift. A potential
whose galaxy has a spectroscopic redshift elsewhere is a foreground or background interloper; the
D-042 rule checks any extra image it produces against a model with that potential on its own plane
(``lens_consistency.py validate/images --plane NAME=Z``, D-046).

Each potential is matched to the nearest catalogue source with ``Z_SPEC > 0`` within
``--match-arcsec``. It is off-plane when ``|z_spec - z_cluster| / (1 + z_cluster) > --dz``.
Output: a table (``<data>/derived/offplane/<model>.ecsv``) and the ``--plane`` arguments.

Redshifts come from the field's CANUCS catalogue (``--field``; a z_spec is trusted when the
photo-z agrees), or, for models in ``MUSE``, from the authors' MUSE redshift catalogue on CDS
(trusted by its quality flag; only trusted rows are matched). These are matched in the model
frame: Caminha+23's positions are on the HST/RELICS frame of its model (ReadMe note 1), and
Richard+21's Abell 2744 and Abell 370 positions agree with the models' members to a median 0.08"
and 0.03" (measured 2026-10-14); Caminha+22's SMACS 0723 positions are moved by the measured
offset to Mahler+22's frame (``frame_offset_arcsec``). ASSUMPTION: a trusted z_spec within the
radius is the potential's own; a MUSE entry for a lensed knot on a member would pass, so every
off-plane row behind the cluster is checked against its photometry before it is interpreted.

Usage::

    python scripts/offplane_potentials.py --model macs0416-canucs --field macs0416
    python scripts/offplane_potentials.py --model abell2744-bergamini23
    python scripts/offplane_potentials.py --model abell370-canucs [--field abell370]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402
import orphan_pairs  # noqa: E402

from jwst_anomaly import paths  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

MATCH_ARCSEC = 0.5  # potential centre to catalogue source (ASSUMPTION)
DZ = 0.02  # |dz| / (1 + z_cluster) beyond which a galaxy is off the plane, ~6000 km/s (ASSUMPTION)
# |Z_ML - z_spec| / (1 + z_spec) within which the photo-z supports the z_spec (ASSUMPTION; a
# compiled z_spec can belong to a blended lensed image or a neighbour)
PHOTOZ_TOL = 0.15

_CDS = "https://cdsarc.cds.unistra.fr/ftp/J/A+A"
# MUSE redshift catalogues (fixed-width CDS files; 1-based byte ranges from each ReadMe). ``ok``:
# quality flags trusted (Richard+21 zconf 3 secure / 2 probable; Caminha+23 q_zsp 3 secure,
# 2 likely, 9 single line; 4 = star). Matched in the model frame (module docstring).
MUSE = {
    "abell2744-bergamini23": {
        "url": f"{_CDS}/646/A83/a2744.dat",
        "sha256": "89d4615fd73d3439c2aab3ca41c84cb1f26f393d375131363de0affa8441e985",
        "cols": {
            "SOURCE": (1, 5),
            "QUAL": (15, 15),
            "RA": (17, 25),
            "DEC": (27, 35),
            "Z_SPEC": (37, 46),
        },
        "ok": {2, 3},
        "ref": "Richard+21 MUSE (J/A+A/646/A83)",
    },
    "elgordo-caminha23": {
        "url": f"{_CDS}/678/A3/catalog.dat",
        "sha256": "c4ddb6e28415ae964ff8ee9828837cd87e5cbea4ff249e77aa78997247b036da",
        "cols": {
            "SOURCE": (1, 28),
            "RA": (30, 39),
            "DEC": (41, 51),
            "Z_SPEC": (53, 58),
            "QUAL": (61, 61),
        },
        "ok": {2, 3, 9},
        "ref": "Caminha+23 MUSE (J/A+A/678/A3)",
    },
    # Richard+21's Abell 370 positions agree with the CANUCS potentials to a median 0.03"
    # (measured 2026-10-14); the CANUCS DR1 z_spec compilation stays available with --field
    "abell370-canucs": {
        "url": f"{_CDS}/646/A83/a370.dat",
        "sha256": "05f524e56ca072b414b3dd1b4134f23e35b1eb99e844dba2d3a5d3935c97d0e0",
        "cols": {
            "SOURCE": (1, 5),
            "QUAL": (15, 15),
            "RA": (17, 25),
            "DEC": (27, 35),
            "Z_SPEC": (37, 46),
        },
        "ok": {2, 3},
        "ref": "Richard+21 MUSE (J/A+A/646/A83)",
    },
    "smacs0723-iclv2": {
        "url": f"{_CDS}/666/L9/muse.dat",
        "sha256": "81868f291a94da964b2a22252ea9fc1624f680d230ff3eedd79e95005277a113",
        "cols": {
            "SOURCE": (1, 2),
            "RA": (4, 15),
            "DEC": (17, 28),
            "Z_SPEC": (30, 35),
            "QUAL": (37, 37),
        },
        "ok": {2, 3, 9},
        "ref": "Caminha+22 MUSE (J/A+A/666/L9)",
        # (dRA cos dec, dDec) arcsec from Mahler+22's model frame to Caminha+22's positions:
        # median over the 24 trusted members within 1", scatter (MAD) 0.009" (2026-10-14)
        "frame_offset_arcsec": (-0.218, 0.262),
    },
}

# CANUCS models without a MUSE entry: the field whose catalogue carries their redshifts
MODEL_FIELD = {"macs0416-canucs": "macs0416"}


def read_muse(path: Path, spec: dict) -> Table:
    """The MUSE catalogue as RA, DEC, Z_SPEC, SOURCE, REF_SPECZ and Z_OK (trusted quality)."""
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    col = {k: [ln[a - 1 : b].strip() for ln in rows] for k, (a, b) in spec["cols"].items()}
    qual = np.array([int(q) for q in col["QUAL"]])
    ra, dec = np.array(col["RA"], float), np.array(col["DEC"], float)
    if "frame_offset_arcsec" in spec:  # catalogue frame -> model frame
        dx, dy = spec["frame_offset_arcsec"]
        sky = SkyCoord(ra * u.deg, dec * u.deg)
        moved = sky.spherical_offsets_by(-dx * u.arcsec, -dy * u.arcsec)
        ra, dec = moved.ra.deg, moved.dec.deg
    return Table(
        {
            "SOURCE": col["SOURCE"],
            "RA": ra,
            "DEC": dec,
            "Z_SPEC": np.array(col["Z_SPEC"], float),
            "REF_SPECZ": [spec["ref"]] * len(rows),
            "Z_OK": np.isin(qual, sorted(spec["ok"])),
        }
    )


def match_potentials(model, cat: Table, z_cluster: float, match_arcsec: float, dz: float) -> Table:
    """One row per potential: nearest ``Z_SPEC > 0`` source, separation, off-plane flag and
    ``z_trusted``: the source's photo-z (``Z_ML``) supports the z_spec, or, for a catalogue
    without ``Z_ML``, its quality flag (``Z_OK``) does."""
    keep = np.asarray(cat["Z_SPEC"], float) > 0
    if "Z_ML" not in cat.colnames:
        if "Z_OK" not in cat.colnames:
            raise ValueError("catalogue needs Z_ML (photo-z check) or Z_OK (quality flag)")
        # match trusted redshifts only, so a nearer low-quality entry cannot hide the member's own
        keep &= np.asarray(cat["Z_OK"], bool)
    spec = cat[keep]
    names = [c.name for c in model.components]
    ra, dec = model.to_sky([c.x for c in model.components], [c.y for c in model.components])
    pot = SkyCoord(ra * u.deg, dec * u.deg)
    src = SkyCoord(np.asarray(spec["RA"], float) * u.deg, np.asarray(spec["DEC"], float) * u.deg)
    idx, sep, _ = pot.match_to_catalog_sky(src)
    sep = sep.arcsec
    hit = sep <= match_arcsec
    zs = np.where(hit, np.asarray(spec["Z_SPEC"], float)[idx], np.nan)
    ids = np.asarray(spec["SOURCE"])
    has_zml = "Z_ML" in spec.colnames
    zp = np.where(hit, np.asarray(spec["Z_ML"], float)[idx], np.nan) if has_zml else None
    ref = np.where(hit, np.asarray(spec["REF_SPECZ"]).astype(str)[idx], "")
    out = Table(
        {
            "name": names,
            "ra": ra,
            "dec": dec,
            "v_disp": [c.v_disp for c in model.components],
            "r_cut": [c.r_cut for c in model.components],
            "z_model": [c.z_lens for c in model.components],
            "cat_id": np.where(hit, ids[idx], -1 if ids.dtype.kind in "iu" else ""),
            "sep_arcsec": sep,
            "z_spec": zs,
            "ref_specz": ref,
        }
    )
    with np.errstate(invalid="ignore"):
        out["offplane"] = np.abs(zs - z_cluster) / (1.0 + z_cluster) > dz
        if has_zml:
            out["z_phot"] = zp
            out["photoz_agrees"] = np.abs(zp - zs) / (1.0 + zs) <= PHOTOZ_TOL
            out["z_trusted"] = out["photoz_agrees"]
            out.meta["photoz_tol"] = PHOTOZ_TOL
        else:
            out["z_trusted"] = hit  # only trusted rows were matched
    out.meta.update(
        {
            "z_cluster": z_cluster,
            "match_arcsec": match_arcsec,
            "dz": dz,
            "label": "derived (positions from the model, redshifts observed)",
        }
    )
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", default="macs0416-canucs", choices=sorted(lc.MODELS))
    ap.add_argument(
        "--field",
        default=None,
        choices=sorted(orphan_pairs.FIELDS),
        help="CANUCS catalogue for the redshifts (default: the model's MUSE catalogue or field)",
    )
    ap.add_argument("--match-arcsec", type=float, default=MATCH_ARCSEC)
    ap.add_argument("--dz", type=float, default=DZ)
    args = ap.parse_args(argv)
    if lc.is_map_model(args.model):
        raise SystemExit(f"error: {args.model} is a map model (no potentials)")
    model, _, _ = lc.load_model(args.model)
    if args.field is None and args.model not in MUSE:
        args.field = MODEL_FIELD.get(args.model)
        if args.field is None:
            raise SystemExit(f"error: {args.model} has no MUSE catalogue; give --field")
    if args.field is None:
        spec = MUSE[args.model]
        cat = read_muse(fetch_catalog(spec["url"], spec["sha256"]), spec)
        z_cluster = float(np.median([c.z_lens for c in model.components]))
        frame = "model frame (MUSE positions on the model's HST/RELICS frame)"
        source = f"{spec['ref']} {spec['url']} sha256 {spec['sha256']}"
    else:
        lc.apply_frame_offset(args.model, model)  # model -> JWST frame of the catalogue (D-034)
        field = orphan_pairs.FIELDS[args.field]
        cat = Table.read(fetch_catalog(*field["catalogue"]), hdu=1)
        z_cluster = field["z_cluster"]
        frame = "JWST frame (CANUCS catalogue)"
        source = f"CANUCS DR1 {field['catalogue'][0]} sha256 {field['catalogue'][1]}"
    out = match_potentials(model, cat, z_cluster, args.match_arcsec, args.dz)
    out.meta.update({"frame": frame, "redshift_catalogue": source})
    dest = paths.data_root() / "derived" / "offplane" / f"{args.model}.ecsv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.write(dest, overwrite=True)
    n_spec = int(np.isfinite(out["z_spec"]).sum())
    off = out[out["offplane"]]
    print(
        f'{len(out)} potentials, {n_spec} with z_spec within {args.match_arcsec}", '
        f"{len(off)} off-plane ({int(off['z_trusted'].sum())} trusted)"
    )
    off.pprint(max_lines=-1, max_width=-1)
    trusted = off[off["z_trusted"]]
    print(" ".join(f"--plane {r['name']}={r['z_spec']:.4f}" for r in trusted))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

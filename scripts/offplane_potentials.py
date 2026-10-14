"""List Lenstool model potentials whose spectroscopic redshift is off the cluster plane (D-042).

A single-plane cluster model puts every member potential at the cluster redshift. A potential
whose galaxy has a spectroscopic redshift elsewhere is a foreground or background interloper; the
D-042 rule checks any extra image it produces against a model with that potential on its own plane
(``lens_consistency.py validate/images --plane NAME=Z``, D-046).

Each potential is matched to the nearest catalogue source with ``Z_SPEC > 0`` within
``--match-arcsec``. It is off-plane when ``|z_spec - z_cluster| / (1 + z_cluster) > --dz``.
Output: a table (``<data>/derived/offplane/<model>.ecsv``) and the ``--plane`` arguments.

Usage::

    python scripts/offplane_potentials.py --model macs0416-canucs --field macs0416
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


def match_potentials(model, cat: Table, z_cluster: float, match_arcsec: float, dz: float) -> Table:
    """One row per potential: nearest ``Z_SPEC > 0`` source, separation, off-plane flag and
    whether the source's photo-z (``Z_ML``) supports the z_spec."""
    spec = cat[np.asarray(cat["Z_SPEC"], float) > 0]
    names = [c.name for c in model.components]
    ra, dec = model.to_sky([c.x for c in model.components], [c.y for c in model.components])
    pot = SkyCoord(ra * u.deg, dec * u.deg)
    src = SkyCoord(np.asarray(spec["RA"], float) * u.deg, np.asarray(spec["DEC"], float) * u.deg)
    idx, sep, _ = pot.match_to_catalog_sky(src)
    sep = sep.arcsec
    hit = sep <= match_arcsec
    zs = np.where(hit, np.asarray(spec["Z_SPEC"], float)[idx], np.nan)
    zp = np.where(hit, np.asarray(spec["Z_ML"], float)[idx], np.nan)
    ref = np.where(hit, np.asarray(spec["REF_SPECZ"]).astype(str)[idx], "")
    out = Table(
        {
            "name": names,
            "ra": ra,
            "dec": dec,
            "v_disp": [c.v_disp for c in model.components],
            "r_cut": [c.r_cut for c in model.components],
            "z_model": [c.z_lens for c in model.components],
            "cat_id": np.where(hit, np.asarray(spec["SOURCE"])[idx], -1),
            "sep_arcsec": sep,
            "z_spec": zs,
            "ref_specz": ref,
            "z_phot": zp,
        }
    )
    with np.errstate(invalid="ignore"):
        out["offplane"] = np.abs(zs - z_cluster) / (1.0 + z_cluster) > dz
        out["photoz_agrees"] = np.abs(zp - zs) / (1.0 + zs) <= PHOTOZ_TOL
    out.meta.update(
        {
            "z_cluster": z_cluster,
            "match_arcsec": match_arcsec,
            "dz": dz,
            "photoz_tol": PHOTOZ_TOL,
            "label": "derived (positions from the model, redshifts observed)",
        }
    )
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", default="macs0416-canucs", choices=sorted(lc.MODELS))
    ap.add_argument("--field", default="macs0416", choices=sorted(orphan_pairs.FIELDS))
    ap.add_argument("--match-arcsec", type=float, default=MATCH_ARCSEC)
    ap.add_argument("--dz", type=float, default=DZ)
    args = ap.parse_args(argv)
    if lc.is_map_model(args.model):
        raise SystemExit(f"error: {args.model} is a map model (no potentials)")
    model, _, _ = lc.load_model(args.model)
    lc.apply_frame_offset(args.model, model)  # model frame -> JWST frame of the catalogue (D-034)
    field = orphan_pairs.FIELDS[args.field]
    cat = Table.read(fetch_catalog(*field["catalogue"]), hdu=1)
    out = match_potentials(model, cat, field["z_cluster"], args.match_arcsec, args.dz)
    dest = paths.data_root() / "derived" / "offplane" / f"{args.model}.ecsv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.write(dest, overwrite=True)
    n_spec = int(np.isfinite(out["z_spec"]).sum())
    off = out[out["offplane"]]
    print(
        f'{len(out)} potentials, {n_spec} with z_spec within {args.match_arcsec}", '
        f"{len(off)} off-plane ({int(off['photoz_agrees'].sum())} with photo-z support)"
    )
    off.pprint(max_lines=-1, max_width=-1)
    trusted = off[off["photoz_agrees"]]
    print(" ".join(f"--plane {r['name']}={r['z_spec']:.4f}" for r in trusted))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""`arcs` on SPT-CL J0615-5746: an independent test of the S1063 anti excess.

Plan (frozen before the run): ``results/arcs/spt0615-relics/PLAN.md``. Converts the RELICS BPZ
catalogue to a zout-like table (``z_phot`` = ``zb``, ``z160`` = ``zbmin``, ``z840`` = ``zbmax``,
stars with ``stel`` >= 0.8 dropped; the 95 % interval stands in for the 68 % one, an ASSUMPTION),
removes the BPZ-JWST astrometric offset when it exceeds 0.1" in either axis, then runs
``lens_consistency.py validate`` and ``arcs`` (D-024 defaults) with the RELICS GLAFIC v2 (primary)
and Lenstool v1 (second) maps on the program 4212 F200W catalogue and applies the frozen outcome
rule through ``arcs_rxj2129.decide``. Shapes are ``observed``, photo-z ``derived`` (RELICS),
predicted axes ``model_prediction``, classes, counts and P values ``derived``.

    python scripts/arcs_spt0615.py --catalog <jw04212-o001_t001 F200W _cat.ecsv> \\
        --bpz <hlsp_relics_hst_acs-wfc3ir_spt0615m57_multi_v1_cat.txt>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arcs_rxj2129 import decide  # noqa: E402

PRIMARY = "spt0615-relics-glafic"
SECOND = "spt0615-relics-lenstool"
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "spt0615-relics"
STAR_STEL = 0.8  # ASSUMPTION (plan): SExtractor stellarity at or above this is a star
PAIR_ARCSEC = 1.0  # offset-measurement pairs
SHIFT_MIN_ARCSEC = 0.1  # offset removed only above this, per axis (plan)
SNR_MIN = 10.0


def bpz_zout(bpz: Table) -> Table:
    """RELICS BPZ rows -> zout-like table (``ra``, ``dec``, ``z_phot``, ``z160``, ``z840``)."""
    keep = np.asarray(bpz["stel"], float) < STAR_STEL
    b = bpz[keep]
    return Table(
        {
            "id": np.asarray(b["id"]),
            "ra": np.asarray(b["RA"], float),
            "dec": np.asarray(b["Dec"], float),
            "z_phot": np.asarray(b["zb"], float),
            "z160": np.asarray(b["zbmin"], float),
            "z840": np.asarray(b["zbmax"], float),
        }
    )


def read_bpz(path: Path) -> Table:
    """Read the RELICS catalogue; column names come from its ``# N name description`` lines."""
    names = []
    for line in path.read_text().splitlines():
        if line.startswith("##"):
            continue
        if not line.startswith("#"):
            break
        parts = line[1:].split()
        if len(parts) >= 2 and parts[0].isdigit():
            names.append(parts[1])
    return ascii.read(path, names=names, comment="#", format="no_header")


def offset_arcsec(zout: Table, cat: Table) -> dict:
    """Median (JWST - BPZ) offset in arcsec (RA*cos dec, Dec) over pairs within 1"."""
    snr = np.asarray(cat["isophotal_flux"], float) / np.asarray(cat["isophotal_flux_err"], float)
    c = cat[np.isfinite(snr) & (snr >= SNR_MIN)]
    sky = c["sky_centroid"]
    c1 = SkyCoord(sky.ra, sky.dec)
    c2 = SkyCoord(zout["ra"], zout["dec"], unit="deg")
    idx, sep, _ = c1.match_to_catalog_sky(c2)
    ok = sep.arcsec <= PAIR_ARCSEC
    dra = (c1.ra.deg[ok] - c2.ra.deg[idx[ok]]) * np.cos(np.radians(c1.dec.deg[ok])) * 3600
    ddec = (c1.dec.deg[ok] - c2.dec.deg[idx[ok]]) * 3600
    return {
        "n_pairs": int(ok.sum()),
        "dra_arcsec": float(np.median(dra)),
        "ddec_arcsec": float(np.median(ddec)),
    }


def shift(zout: Table, off: dict) -> bool:
    """Apply the offset to the BPZ positions if either axis exceeds 0.1"; returns whether it did."""
    if not off.get("n_pairs", 1) or not np.isfinite([off["dra_arcsec"], off["ddec_arcsec"]]).all():
        raise SystemExit("error: no BPZ-JWST pairs within 1 arcsec; offset undefined")
    if max(abs(off["dra_arcsec"]), abs(off["ddec_arcsec"])) <= SHIFT_MIN_ARCSEC:
        return False
    zout["dec"] = zout["dec"] + off["ddec_arcsec"] / 3600
    zout["ra"] = zout["ra"] + off["dra_arcsec"] / 3600 / np.cos(np.radians(zout["dec"]))
    return True


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--bpz", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--zout", type=Path, help="where to write the converted photo-z table")
    args = ap.parse_args(argv)
    zout = bpz_zout(read_bpz(args.bpz))
    off = offset_arcsec(zout, Table.read(args.catalog))
    off["shift_applied"] = shift(zout, off)
    zpath = args.zout or args.bpz.with_name("spt0615_bpz_zout.fits")
    zout.write(zpath, overwrite=True)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "photoz_offset.json").write_text(
        json.dumps({"bpz": args.bpz.name, "n_bpz_galaxies": len(zout), **off}, indent=1) + "\n"
    )
    print(json.dumps(off))
    decide(PRIMARY, SECOND, args.out, args.catalog, zpath, "results/arcs/spt0615-relics/PLAN.md")
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())

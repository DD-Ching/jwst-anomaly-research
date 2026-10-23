"""RXJ2129: the reference-catalogue extra-image test, stage 1 (gate K and reference inputs).

Frozen plan: ``results/lens_consistency/rxj2129_djaref/PLAN_1.md``. The Dhanasingham+2026 image list
(arXiv:2604.11882, Appendix B table "Systems_RA_and_dec") is transcribed below; rule 1 runs gate K
of ``m0417_images`` on it with the pinned JWST-era Lenstool maps (``rxj2129-lenstool-maps``, no
frame offset); rule 2 reads only the FITS headers of the DJA v7.0 photometry (HTTP range) and
checks SMACS rule 2's columns. Writes ``summary.json``.

    python scripts/rxj2129_djaref.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import arcs_rxj2129
import elgordo_djaref as eg
import m0417_images as mi
import numpy as np
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


def main(argv: list[str] | None = None) -> int:
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    t = image_table()
    t.write(OUT / IMAGES, format="ascii.ecsv", overwrite=True)
    arcs_rxj2129.stage_members(lensmodel.RXJ2129_LENSTOOL_MAPS)
    mi.FIELDS[FIELD] = {"images": IMAGES, "models": (MODEL,)}
    summary: dict = {
        "images": {"n": len(t), "families": len(set(t["system"])), "source": PAPER},
        "stage1_gate_K": mi.run_field(FIELD, OUT, GRID_STEP),
    }
    blob, total = eg.fetch_head(DJA_PHOT)
    summary["stage1_reference_inputs"] = stage2_inputs(blob, total) | {"url": DJA_PHOT}
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
    (OUT / "summary.json").write_text(json.dumps(_clean(summary), indent=1) + "\n")
    g = summary["stage1_gate_K"][MODEL]["gate_K"]
    print(json.dumps({"gate_K": g, "ref": summary["stage1_reference_inputs"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

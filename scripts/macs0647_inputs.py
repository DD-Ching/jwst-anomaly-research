"""MACS0647: input check for the reference-catalogue extra-image route (no residual read).

Three inputs the route needs (CHANGELOG 2026-10-23 RXJ2129 handoff): a public JWST-era image list,
model maps, and a reference catalogue that deblends arcs from cluster light. The image list is the
Meena+2023 table (arXiv:2211.13334, ``lensed_systems.tex`` of the arXiv source), parsed by
``parse_lensed_systems`` and written to ``meena23_images.ecsv``; the Nishida+2025 GLAFIC
constraint list (arXiv:2505.00553 appendix) is parsed by ``parse_nishida25``
(``nishida25_images.ecsv``); the reference is the DJA v7.0 photometry, of which only the FITS
headers are read (HTTP range) and checked against SMACS rule 2's columns. Writes ``summary.json``.

    python scripts/macs0647_inputs.py --tex <lensed_systems.tex> --tex-n25 <ms.tex of 2505.00553>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import astropy.units as u
import elgordo_djaref as eg
import numpy as np
import smacs_djaref as sd
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/macs0647_inputs")
IMAGES = "meena23_images.ecsv"
IMAGES_N25 = "nishida25_images.ecsv"
PAPER_N25 = (
    "Nishida et al. 2025, arXiv:2505.00553, Appendix tables 'Table:MACS0647_1' and '_2' (ms.tex); "
    "rows with the JWST check mark"
)
PAPER = (
    "Meena et al. 2023, ApJL 944, L6 (arXiv:2211.13334), Table 'tab:lensed' (lensed_systems.tex)"
)
DJA_PHOT = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/macs0647-grizli-v7.0-fix_phot.fits"

_ROW = re.compile(r"^\s*([0-9]+\.[0-9]+)(\??)\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*([^&]*)&(.*)$")
_ZPHOT = re.compile(r"\$\s*([0-9.]+)_")
_DATA_ROW = re.compile(r"^\s*[0-9]+\.[0-9]+\S*\s*&")  # any uncommented row with an image id


def _unparsed(line: str) -> ValueError:
    return ValueError(f"table row not parsed: {line.strip()[:80]!r}")


def table_id(paper_id: str) -> tuple[str, str]:
    """(image_id, family) of a paper id: ``"4.2"`` -> ``("4.2", "4")``. System 6 (MACS0647-JD)
    lists three components with three images each (``6.<c><n>``); each component is its own
    family, ``"6.12"`` -> ``("6-1.2", "6-1")``, as in Hsiao et al. 2023's components."""
    fam, n = paper_id.split(".")
    if fam == "6" and len(n) == 2:
        fam, n = f"6-{n[0]}", n[1]
    return f"{fam}.{n}", fam


def parse_lensed_systems(tex: str) -> Table:
    """Rows of the Meena+2023 image table; TeX-commented rows (``%``) are left out, a trailing
    ``?`` on the id sets ``uncertain``, ``new_in_jwst`` marks the families whose first row's
    comment is "New system" (first identified in the JWST data), and ``z_phot`` is the EAZY median
    (NaN for ``--``). A data row or a redshift cell of any other form raises ``ValueError``."""
    rows, new = [], {}
    for line in tex.splitlines():
        if line.lstrip().startswith("%"):
            continue
        m = _ROW.match(line)
        if not m:
            if _DATA_ROW.match(line):
                raise _unparsed(line)
            continue
        pid, q, ra, dec, zcol, comment = m.groups()
        image_id, fam = table_id(pid)
        new.setdefault(fam.split("-")[0], "New system" in comment)
        z = _ZPHOT.search(zcol)
        if not z and zcol.strip() != "--":
            raise _unparsed(line)
        rows.append(
            (
                image_id,
                fam,
                pid + q,
                float(ra),
                float(dec),
                float(z.group(1)) if z else np.nan,
                bool(q),
            )
        )
    t = Table(
        rows=rows, names=("image_id", "system", "paper_id", "ra", "dec", "z_phot", "uncertain")
    )
    t["z_spec"] = np.zeros(len(t), bool)
    t["new_in_jwst"] = [new[f.split("-")[0]] for f in t["system"]]
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        f"{PAPER}; z_phot is the paper's EAZY median (no spectroscopic z in the table)"
    )
    return t


_ROW_N25 = re.compile(
    r"^\s*([0-9]+\.[0-9]+)\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*([^&]*)&\s*([^&]*)&"
)
_Z_N25 = re.compile(r"\$\s*([0-9.]+)\s*(?:\\pm\s*([0-9.]+))?\s*\$")


def parse_nishida25(tex: str) -> Table:
    """JWST rows of the Nishida+2025 MACS0647 tables. ``z_jwst`` is printed on a family's first row
    and carried to its other rows; ``z_jwst_pm`` is its printed ``\\pm`` (NaN when none is printed,
    the paper does not say what a value without one means; ``\\nodata`` gives NaN for both). A data
    row, or a family's first redshift cell, of any other form raises ``ValueError``."""
    rows, z, pm, last = [], np.nan, np.nan, None
    for line in tex.splitlines():
        if line.lstrip().startswith("%"):
            continue
        m = _ROW_N25.match(line)
        if not m:
            if _DATA_ROW.match(line):
                raise _unparsed(line)
            continue
        pid, ra, dec, jwst, zcol = m.groups()
        fam = pid.split(".")[0]
        if fam != last:  # the family's first printed row, whether or not it is a JWST row
            last = fam
            zm = _Z_N25.search(zcol)
            if not zm and zcol.strip() != "\\nodata":
                raise _unparsed(line)
            z = float(zm.group(1)) if zm else np.nan
            pm = float(zm.group(2)) if zm and zm.group(2) else np.nan
        if "checkmark" not in jwst:
            continue
        rows.append((pid, fam, pid, float(ra), float(dec), z, pm))
    t = Table(
        rows=rows, names=("image_id", "system", "paper_id", "ra", "dec", "z_jwst", "z_jwst_pm")
    )
    t.meta["provenance"] = "observed"
    t.meta["source"] = PAPER_N25
    return t


def table_sha256(t: Table) -> str:
    """sha256 of every column of every row (``str`` of each value, ``,`` and ``;`` separated)."""
    s = ";".join(",".join(str(v) for v in row) for row in t.iterrows())
    return hashlib.sha256(s.encode()).hexdigest()


def list_offsets(a: Table, b: Table, radius: float = 0.5) -> dict:
    """Nearest-neighbour match of ``b`` to ``a`` (one direction, not unique): counts within 0.3 /
    0.5 / 1.0 arcsec, the distinct ``a`` rows among them, and, over the pairs within ``radius``
    only, the median offset b - a (arcsec; RA times cos Dec)."""
    ca = SkyCoord(a["ra"], a["dec"], unit="deg")
    cb = SkyCoord(b["ra"], b["dec"], unit="deg")
    i, sep, _ = cb.match_to_catalog_sky(ca)
    s = sep.arcsec
    k = s < radius
    if not k.any():
        raise ValueError(f"no pair within {radius} arcsec")
    dra = ((cb.ra - ca[i].ra) * np.cos(cb.dec.rad)).to(u.arcsec).value[k]
    dde = (cb.dec - ca[i].dec).to(u.arcsec).value[k]
    return {
        "n_within": {str(r): int((s < r).sum()) for r in (0.3, 0.5, 1.0)},
        "n_distinct_a_within": {str(r): len(set(i[s < r])) for r in (0.3, 0.5, 1.0)},
        "radius_arcsec": radius,
        "median_dra_arcsec": float(np.median(dra)),
        "median_ddec_arcsec": float(np.median(dde)),
        "median_sep_arcsec": float(np.median(s[k])),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tex", type=Path, required=True)
    ap.add_argument("--tex-n25", type=Path, required=True)
    a = ap.parse_args(argv)
    raw25 = a.tex_n25.read_bytes()
    n25 = parse_nishida25(raw25.decode("utf-8", errors="replace"))
    raw = a.tex.read_bytes()
    t = parse_lensed_systems(raw.decode("utf-8", errors="replace"))
    for tab in (t, n25):
        if [lensmodel.image_family(i) for i in tab["image_id"]] != list(tab["system"]):
            raise ValueError("an image id does not parse back to its family")
    blob, total = eg.fetch_head(DJA_PHOT)  # network first, so a failure writes nothing
    fams = list(dict.fromkeys(t["system"]))
    summary = {
        "images": {
            "n": len(t),
            "families": len(fams),
            "uncertain": [str(p) for p in t["paper_id"][t["uncertain"]]],
            "families_new_in_jwst": list(dict.fromkeys(t["system"][t["new_in_jwst"]])),
            "families_with_one_image": [f for f in fams if np.sum(t["system"] == f) < 2],
            "n_z_phot_missing": int(np.isnan(t["z_phot"]).sum()),
            "tex_sha256": hashlib.sha256(raw).hexdigest(),
            "table_sha256": table_sha256(t),
            "source": PAPER,
        },
        "images_nishida25": {
            "n": len(n25),
            "families": len(set(n25["system"])),
            "families_z_with_pm": sorted(
                set(n25["system"][np.isfinite(n25["z_jwst_pm"])]), key=int
            ),
            "families_z_without_pm": sorted(
                set(n25["system"][np.isfinite(n25["z_jwst"]) & ~np.isfinite(n25["z_jwst_pm"])]),
                key=int,
            ),
            "families_z_nodata": sorted(set(n25["system"][~np.isfinite(n25["z_jwst"])]), key=int),
            "tex_sha256": hashlib.sha256(raw25).hexdigest(),
            "table_sha256": table_sha256(n25),
            "source": PAPER_N25,
        },
        "nishida25_vs_meena23": list_offsets(t, n25),
        "reference_inputs": sd.stage2_inputs(blob, total) | {"url": DJA_PHOT},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    t.write(OUT / IMAGES, format="ascii.ecsv", overwrite=True)
    n25.write(OUT / IMAGES_N25, format="ascii.ecsv", overwrite=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

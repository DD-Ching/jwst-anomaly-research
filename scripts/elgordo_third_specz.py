"""El Gordo: DJA NIRSpec redshifts at the Diego+2023 third-image candidates (D-034 addendum ay).

Frozen plan: ``results/lens_consistency/elgordo_third_specz/PLAN.md``. Rules 1-4 are the functions
below (rule 1 matching and rule 2 image z reuse ``macs0647_specz``). The DJA table is streamed, not
stored (cloud disk rule):

    curl -sS <TABLE url> | gunzip -c | python scripts/elgordo_third_specz.py /dev/stdin

"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import macs0647_specz as ms
import numpy as np
from astropy.table import Table

OUT = Path("results/lens_consistency/elgordo_third_specz")
IMAGES = Path("results/lens_consistency/elgordo_heldout/diego23_images.ecsv")
TABLE = ms.TABLE  # DJA v4.6, the bytes of addendum ai
PRIMARY = ("6c", "13c", "20c")
REPORTED = ("4c", "15c", "18c", "19c")
# Frozen rules (PLAN.md; ASSUMPTIONs)
Z_AGREE = ms.Z_AGREE  # 0.01 (1 + z)
Z_MIN_IMAGES = 3
Z_MIN_AGREE = 0.8


def agrees(z_image: float, z_family: float) -> bool:
    return bool(abs(z_image - z_family) / (1 + z_family) <= Z_AGREE)


def gate_z(pairs: list[tuple[float, float]]) -> dict:
    """Rule 3: (image z, family z) of the control images that have an image z."""
    n = len(pairs)
    frac = float(np.mean([agrees(a, b) for a, b in pairs])) if n else 0.0
    return {
        "n_images": n,
        "frac_agree": frac,
        "passes": bool(n >= Z_MIN_IMAGES and frac >= Z_MIN_AGREE),
    }


def reading(status: str, z_image: float | None, z_family: float, z_passes: bool) -> str:
    """Rule 4: a target's reading from its rule-2 status and z."""
    if not z_passes:
        return "not_decided"
    if status == "none":
        return "no_spectrum"
    if status == "ambiguous":
        return "ambiguous"
    return "confirmed" if agrees(z_image, z_family) else "z_differs"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("table_csv", type=Path)
    args = ap.parse_args(argv)
    img = Table.read(IMAGES, format="ascii.ecsv")
    img["image_id"] = np.asarray(img["image_id"], str)
    matches, n_rows = ms.match_rows(args.table_csv, img)
    if n_rows != TABLE["n_rows"]:  # a cut-off stream changes the matches
        raise SystemExit(f"read {n_rows} rows, expected {TABLE['n_rows']}: incomplete table")
    zfam = {str(r["image_id"]): float(r["z"]) for r in img}
    spec = {str(r["image_id"]): bool(r["z_spec"]) for r in img}
    s1 = {i: ms.image_z(matches.get(i, [])) for i in zfam}
    targets = PRIMARY + REPORTED
    controls = {
        i: (v[1], zfam[i])
        for i, v in s1.items()
        if v[0] == "spec_z" and spec[i] and i not in targets and _family(i) <= 23
    }
    gz = gate_z(list(controls.values()))
    summary = {
        "plan": str(OUT / "PLAN.md"),
        "table": TABLE,
        "rows_read": n_rows,
        "rules": {
            "MATCH_RADIUS_ARCSEC": ms.MATCH_RADIUS_ARCSEC,
            "SPEC_GRADE": ms.SPEC_GRADE,
            "Z_AGREE": Z_AGREE,
            "Z_MIN_IMAGES": Z_MIN_IMAGES,
            "Z_MIN_AGREE": Z_MIN_AGREE,
        },
        "matches": matches,
        "S1": {i: v[0] for i, v in s1.items() if i in matches},
        "controls": {
            i: {"z_image": a, "z_family": b, "agrees": agrees(a, b)}
            for i, (a, b) in controls.items()
        },
        "gate_Z": gz,
        "targets": {
            i: {
                "primary": i in PRIMARY,
                "status": s1[i][0],
                "z_image": s1[i][1],
                "z_family": zfam[i],
                "reading": reading(s1[i][0], s1[i][1], zfam[i], gz["passes"]),
            }
            for i in targets
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({"gate_Z": gz, "targets": summary["targets"]}, indent=1))
    return 0


def _family(image_id: str) -> int:
    """Leading integer of a Diego ID (``1.1a`` -> 1, ``20c`` -> 20)."""
    digits = ""
    for ch in image_id:
        if not ch.isdigit():
            break
        digits += ch
    return int(digits)


if __name__ == "__main__":
    raise SystemExit(main())

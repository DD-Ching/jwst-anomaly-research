"""Abell S1063 (RXC J2248.7-4431): the published multiple-image list of Chan et al. 2020
(ApJ 888, 35; VizieR J/ApJ/888/35 table1), in the image-list format of ``m0417_images.py``, so the
D-035 rms gate can be opened with spectroscopic system redshifts once the CATS v4.1 maps load.

Each image's z is its system's MUSE z_spec where any image of the system has one (z_spec True),
else the Chan model's geometric z ``Geoz`` (z_spec False). No lens-model reading is made here.

    python scripts/s1063_imagelist.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

OUT = Path("results/lens_consistency/abells1063_images")
CATALOG = "J/ApJ/888/35/table1"
SOURCE = (
    "Chan et al. 2020, ApJ 888, 35 (RXC J2248.7-4431 lensed images), VizieR J/ApJ/888/35 "
    "table1, streamed {date}; z of a system is its MUSE zspec where any image has one "
    "(z_spec True), else the Chan model Geoz (z_spec False)"
)


def _value(x) -> float:
    return float("nan") if np.ma.is_masked(x) else float(x)


def system_of(image_id: str) -> str:
    """'12.3' -> '12', '1a.1' -> '1a': the family is everything before the last dot."""
    return image_id.rsplit(".", 1)[0]


def build(rows: Table, date: str) -> tuple[Table, dict]:
    """Image list (id, ra, dec, z, z_spec, cl, dup_position) and a summary from the VizieR rows."""
    ids = [str(i).strip() for i in rows["ID"]]
    zspec = np.array([_value(z) for z in rows["zspec"]])
    geoz = np.array([_value(z) for z in rows["Geoz"]])
    systems = [system_of(i) for i in ids]
    sys_zspec: dict[str, float] = {}
    for s, z in zip(systems, zspec, strict=True):
        if np.isfinite(z):
            sys_zspec.setdefault(s, z)
    z_out, is_spec = [], []
    for s, g in zip(systems, geoz, strict=True):
        if s in sys_zspec:
            z_out.append(sys_zspec[s])
            is_spec.append(True)
        else:
            z_out.append(g)
            is_spec.append(False)
    coord = SkyCoord(rows["RAJ2000"], rows["DEJ2000"], unit=("hourangle", "deg"))
    # Rows sharing an exact position with another row (a transcription error in the source table;
    # 26.1 and 27.1-27.3 in J/ApJ/888/35): flagged, to be left out until checked against the paper
    pos = [
        (str(a).strip(), str(d).strip())
        for a, d in zip(rows["RAJ2000"], rows["DEJ2000"], strict=True)
    ]
    dup = [pos.count(q) > 1 for q in pos]
    out = Table(
        {
            "image_id": ids,
            "ra": coord.ra.deg,
            "dec": coord.dec.deg,
            "z": np.round(np.array(z_out, dtype=float), 4),
            "z_spec": is_spec,
            "cl": [int(c) for c in rows["CL"]],
            "dup_position": dup,
        }
    )
    out["ra"].unit = out["dec"].unit = "deg"
    out.meta = {"provenance": "observed", "source": SOURCE.format(date=date)}
    # Systems whose images disagree in zspec (a check on the system grouping)
    spread = {}
    for s in sorted(set(systems)):
        zs = zspec[[k for k, t in enumerate(systems) if t == s]]
        zs = zs[np.isfinite(zs)]
        if zs.size > 1:
            spread[s] = float(zs.max() - zs.min())
    # Model z vs spec z where both exist (derived; the Chan model may be fitted to these z)
    dz = [
        float(sys_zspec[s] - g)
        for s, g in zip(systems, geoz, strict=True)
        if s in sys_zspec and np.isfinite(g)
    ]
    summary = {
        "catalog": CATALOG,
        "retrieved": date,
        "n_images": len(out),
        "n_systems": len(set(systems)),
        "n_systems_zspec": len(sys_zspec),
        "n_images_system_zspec": int(np.sum(is_spec)),
        "dup_position_ids": [i for i, d in zip(ids, dup, strict=True) if d],
        "n_images_no_z": int(np.sum(~np.isfinite(np.array(z_out)))),
        "cl_counts": {str(c): int(np.sum(out["cl"] == c)) for c in sorted(set(out["cl"]))},
        "max_zspec_spread_within_system": max(spread.values()) if spread else None,
        "median_abs_zspec_minus_geoz": float(np.median(np.abs(dz))) if dz else None,
    }
    return out, summary


def main() -> None:
    from datetime import UTC, datetime

    from astroquery.vizier import Vizier

    rows = Vizier(row_limit=-1, columns=["**"]).get_catalogs(CATALOG)[0]
    date = datetime.now(UTC).strftime("%Y-%m-%d")
    out, summary = build(rows, date)
    OUT.mkdir(parents=True, exist_ok=True)
    out.write(OUT / "chan20_images.ecsv", overwrite=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

"""Is the RELICS Lenstool v1 amplitude offset a z_ref normalisation? (D-034 addendum bi)

Frozen plan: ``results/lens_consistency/relics_lenstool_zref/PLAN.md``. Tests H_U (one RELICS-wide
z_ref ≈ 2.09 fixed by the Abell 3192 factor) against the stored rule-7 amplitude factors of five
clusters, and H_C (a per-map normalisation keyword) on the published FITS headers, read by byte
range (header only).

    python scripts/relics_lenstool_zref.py            # headers + table -> summary.json
    python scripts/relics_lenstool_zref.py --offline  # H_U only, no network
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import lens_consistency as lc
import numpy as np
import slice_relics_macs0358 as m0358

OUT = Path("results/lens_consistency/relics_lenstool_zref")
Z_REF = 2.0946208510160416  # slice_relics_macs0358/posthoc.json amplitude_as_dls_ds_zref
TOLERANCE = 0.10  # PLAN.md H_U rule (ASSUMPTION)
HEADER_BYTES = 2880 * 10
#: keywords whose name or comment may carry a source z or a D_LS/D_S normalisation
NORM_PATTERN = re.compile(r"(^Z|ZS|ZSRC|ZSOURCE|DLS|DSL|DS|NORM|REDSH|SOURCE)", re.I)


def dls_ds(z_lens: float, z_s: float) -> float:
    from astropy.cosmology import FlatLambdaCDM

    cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3)
    return float(
        cosmo.angular_diameter_distance(z_lens, z_s) / cosmo.angular_diameter_distance(z_s)
    )


def stored_factors() -> dict[str, float]:
    """Median rule-7 factor per cluster from the stored runs (derived)."""
    base = Path("results/lens_consistency")
    scale = json.loads((base / "slice_relics_heldout/posthoc.json").read_text())["scale"]
    out = {c: float(np.median([f["scale"] for f in scale[c]["families"].values()])) for c in scale}
    rxj = json.loads((base / "slice_relics_rxj0152/summary.json").read_text())
    out["cl0152"] = float(rxj["models"]["cl0152-relics-lenstool"]["amplitude"]["median"])
    a3192 = json.loads((base / "slice_relics_macs0358/posthoc.json").read_text())
    out["a3192"] = float(
        a3192["amplitude_as_dls_ds_zref"]["a3192-relics-lenstool"]["median_factor"]
    )
    return out


def h_u(factors: dict[str, float]) -> dict:
    rows = {}
    for c, f in factors.items():
        z_l = lc.MODELS[f"{c}-relics-lenstool"]["z_lens"]
        pred = dls_ds(z_l, Z_REF)
        rows[c] = {
            "z_lens": z_l,
            "f_obs": f,
            "f_pred": pred,
            "diff": f - pred,
            "z_ref_implied": m0358.zref_for_factor(f, z_l),
            "consistent": abs(f - pred) <= TOLERANCE,
        }
    return {
        "z_ref": Z_REF,
        "tolerance": TOLERANCE,
        "rows": rows,
        "rejected": not all(r["consistent"] for r in rows.values()),
    }


def read_header(url: str):
    """Primary header by byte range (astropy Header)."""
    import requests
    from astropy.io import fits

    r = requests.get(url, headers={"Range": f"bytes=0-{HEADER_BYTES - 1}"}, timeout=60)
    r.raise_for_status()
    text = r.content.decode("ascii", "replace")
    end = re.search(r"END {77}", text)
    if end is None:
        raise ValueError(f"END card not within {HEADER_BYTES} bytes: {url}")
    return fits.Header.fromstring(text[: end.end()]), r.status_code


def h_c(clusters) -> dict:
    out = {}
    for c in clusters:
        files = lc.MODELS[f"{c}-relics-lenstool"]["files"]
        out[c] = {}
        for key in ("kappa_map", "alpha_x"):
            url = files[key][0]
            hdr, status = read_header(url)
            cards = [
                (cd.keyword, str(cd.value), cd.comment)
                for cd in hdr.cards
                if cd.keyword not in ("", "COMMENT", "HISTORY")
            ]
            notes = [str(cd.value) for cd in hdr.cards if cd.keyword in ("COMMENT", "HISTORY")]
            norm = [cd for cd in cards if NORM_PATTERN.search(cd[0]) or NORM_PATTERN.search(cd[2])]
            out[c][key] = {
                "url": url,
                "http": status,
                "cards": cards,
                "notes": notes,
                "normalisation_like": norm,
            }
    return out


def h_c_verdict(h_c_rows: dict) -> str:
    """PLAN.md H_C: "not_supported" when no map header has a normalisation-like card; otherwise
    "to_read", because the 0.05 match of an implied factor needs the card's meaning (by hand)."""
    found = any(m["normalisation_like"] for maps in h_c_rows.values() for m in maps.values())
    return "to_read" if found else "not_supported"


#: post hoc (not in PLAN.md): SMACS 0723 Lenstool v1 (addendum aw: |alpha| 1.14-1.56 x GLAFIC's)
SMACS_X = (
    "https://archive.stsci.edu/hlsps/relics/smacs0723m73/models/lenstool/"
    "hlsp_relics_model_model_smacs0723m73_lenstool_v1_x-arcsec-deflect.fits"
)


def posthoc_makers() -> dict:
    """Post hoc: map maker (MAPMKR_F / MAPMKR_L) and DATE of every RELICS Lenstool v1 deflection map
    registered in ``lens_consistency.MODELS``, plus SMACS 0723 (observed header values)."""
    urls = {
        name.removesuffix("-relics-lenstool"): spec["files"]["alpha_x"][0]
        for name, spec in lc.MODELS.items()
        if name.endswith("-relics-lenstool")
    }
    urls["smacs0723"] = SMACS_X
    out = {}
    for c, url in sorted(urls.items()):
        hdr, _ = read_header(url)
        out[c] = {
            "url": url,
            "map_maker": f"{hdr.get('MAPMKR_F', '')} {hdr.get('MAPMKR_L', '')}".strip(),
            "date": hdr.get("DATE"),
            "naxis1": hdr.get("NAXIS1"),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    factors = stored_factors()
    res = {"h_u": h_u(factors)}
    if not args.offline:  # --offline prints H_U only and keeps the stored summary
        res["h_c"] = h_c(sorted(factors))
        res["h_c_verdict"] = h_c_verdict(res["h_c"])
        res["posthoc_makers"] = posthoc_makers()
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "summary.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    for c, r in res["h_u"]["rows"].items():
        print(
            f"{c:9s} z_l {r['z_lens']:.3f} f_obs {r['f_obs']:.3f} f_pred {r['f_pred']:.3f} "
            f"diff {r['diff']:+.3f} z_ref {r['z_ref_implied']} ok {r['consistent']}"
        )
    print("H_U rejected:", res["h_u"]["rejected"])


if __name__ == "__main__":
    main()

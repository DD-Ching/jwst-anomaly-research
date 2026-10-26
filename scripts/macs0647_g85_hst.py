"""MACS0647 1433 g85 (post-hoc vetting): HST WFC3/IR F110W PID 17504 (2023-12, 2024-03) minus
CLASH 2011 F110W difference photometry at the JWST transient position.

Not part of the frozen plan (``results/transients/macs0647_1433/PLAN.md``); a post-hoc third-epoch
check. Reuses the byte-range reader, resampling, registration and aperture functions of
``scripts/sneos_hst.py`` / ``scripts/sneos_ratio.py`` (tested there). Fluxes are observed
(difference images); every number is derived.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sneos_hst as sh  # noqa: E402
import sneos_ratio as sr  # noqa: E402

RA, DEC = 101.99665, 70.250721  # arc_diff h1023 (results/transients/macs0647_1433/core/hits.ecsv)
TEMPLATE, TEMPLATE_CHECK = "ibffa2020", "ibffa9030"  # CLASH F110W 2011-10 / 2011-11
SCIENCE = ("ifaa27010", "ifaa29010", "ifaa26010", "ifaa25010", "ifaa30010", "ifaa28010")
OUT = Path(__file__).resolve().parents[1] / "results/transients/macs0647_1433/g85_hst.json"
C_AA = 2.99792458e18  # speed of light, Å/s


def ab_mag(flam_sum: float, photplam: float) -> float:
    """AB magnitude of an aperture sum in erg/s/cm²/Å (NaN when not positive)."""
    if not flam_sum > 0:
        return float("nan")
    fnu = flam_sum * photplam**2 / C_AA
    return float(-2.5 * np.log10(fnu) - 48.6)


def measure(obs: str, tmpl: sh.Drz) -> dict:
    m = sh.Drz(obs)
    sci, _, w = m.section(RA, DEC, sr.STAMP_HALF_ARCSEC)
    p = tuple(float(v) for v in w.world_to_pixel_values(RA, DEC))
    yy, xx = np.indices(sci.shape)
    mask = np.hypot(xx - p[0], yy - p[1]) <= sr.SN_MASK_ARCSEC / m.pix
    diff, info = sh.difference(sci, w, tmpl, RA, DEC, mask, m.pix)
    f = sh._single(diff, p, m.pix)
    s = sh.random_sigma(diff, p, p, m.pix)
    photplam = float(m.hdul["SCI"].header.get("PHOTPLAM", m.hdul[0].header.get("PHOTPLAM")))
    return {
        "used": bool(info["registered"]),  # unregistered rows are reported, never used
        "science": obs,
        "template": tmpl.obs,
        "date": m.date,
        "exptime": m.exptime,
        "flux_flam": f,
        "sigma_random_flam": s,
        "snr": f / s,
        "ab_mag_aperture": ab_mag(f, photplam),
        "ab_mag_5sigma_aperture": ab_mag(5 * s, photplam),
        **info,
    }


def main() -> int:
    tmpl, tchk = sh.Drz(TEMPLATE), sh.Drz(TEMPLATE_CHECK)
    rows = [measure(o, tmpl) for o in SCIENCE]
    rows.append(measure(SCIENCE[0], tchk))
    out = {
        "provenance": "derived (post hoc; not in PLAN.md)",
        "position": [RA, DEC],
        "aperture_arcsec": sh.APER_ARCSEC,
        "annulus_arcsec": list(sh.ANNULUS_ARCSEC),
        "rows": rows,
    }
    for r in rows:  # strict JSON: NaN -> null
        for k, v in r.items():
            if isinstance(v, float) and not np.isfinite(v):
                r[k] = None
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    for r in rows:
        print(
            r["science"], r["template"], r["date"], f"S/N {r['snr']:+.1f}",
            f"m {r['ab_mag_aperture']} lim5 {r['ab_mag_5sigma_aperture']:.2f}",
            f"reg {r['registered']} {r['shift_arcsec']:.3f}",
        )  # fmt: skip
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

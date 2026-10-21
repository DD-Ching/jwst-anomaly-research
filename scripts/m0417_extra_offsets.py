"""Post hoc diagnostic for ``m0417_extra.py`` (not a frozen rule): offset of each catalogued
image's 0.4"-recentred F200W peak from its published position, per field, and the reference-rule
quantities ``forced_check`` applies (peak / position flux, f(0.2")/f(0.4")).

    python scripts/m0417_extra_offsets.py
"""

from __future__ import annotations

import json

import lens_consistency as lc
import m0417_extra as mx
import numpy as np


def main() -> int:
    out = {}
    for field, (_, _, obs, tile) in mx.FIELDS.items():
        t = mx.load_images(field)
        stamp, close = lc.image_stamper(
            mx.S3.format(obs, tile, "i2d.fits"), lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.1
        )
        rows = []
        try:
            for r in t:
                s, e, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
                f, _ = lc.peak_flux(s, e, xx, yy)
                f0, _ = lc.aperture_snr(s, e, xx, yy, 0.0, 0.0)
                px, py = lc.peak_position(s, xx, yy)
                fb, _ = lc.aperture_snr(s, e, xx, yy, px, py, r_ap=0.4)
                rows.append(
                    {
                        "image_id": str(r["image_id"]),
                        "peak_dx_arcsec": float(px),
                        "peak_dy_arcsec": float(py),
                        "peak_over_position": float(f / f0) if f0 > 0 else float("nan"),
                        "concentration": float(f / fb) if fb > 0 else float("nan"),
                    }
                )
        finally:
            close()
        dx = np.array([r["peak_dx_arcsec"] for r in rows])
        dy = np.array([r["peak_dy_arcsec"] for r in rows])
        po = np.array([r["peak_over_position"] for r in rows])
        out[field] = {
            "provenance": "derived (post hoc)",
            "n_images": len(rows),
            "median_peak_dx_arcsec": float(np.nanmedian(dx)),
            "median_peak_dy_arcsec": float(np.nanmedian(dy)),
            "median_abs_offset_arcsec": float(np.nanmedian(np.hypot(dx, dy))),
            "fraction_peak_over_position_le_2": float(np.mean(po <= lc.MAX_PEAK_OVER_POSITION)),
            "rows": rows,
        }
        print(field, {k: v for k, v in out[field].items() if k != "rows"})
    (mx.OUT / "offsets_posthoc.json").write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

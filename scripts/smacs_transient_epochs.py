"""SMACS 0723 system 5 transient pair: small-aperture photometry in every public NIRCam epoch.

The two-epoch grid (D-027 addendum 2) found two point sources on the system 5 arc that are absent in
the 2022 ERO and present in VENUS o057 (2026-06-05). The position of image B is the target of JWST
DD 12774 and HST 18297 ("SMACS0723-GLSN-1", SN Atalanta). This measures both positions in the four
public NIRCam epochs (2022 ERO, 2024 GO-4043, 2026 VENUS, 2026 DD 12774) by S3 byte-range cutouts:
background-subtracted flux in an r = 0.15" aperture with a 0.4-0.7" annulus, no aperture correction,
so the numbers compare epochs, not total fluxes. Flux and the noise estimate are ``observed`` /
``derived``; the stamps go to the outputs directory (``smacs_transient_epochs.png``, untracked).

    python scripts/smacs_transient_epochs.py --out results/transients/smacs_grid/pair_epochs.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from jwst_anomaly import paths
from jwst_anomaly.cutouts import make_cutouts

PAIR = {"A": (110.825136, -73.452178), "B": (110.823747, -73.452975)}  # grid g03415, g03423
S3 = "s3://stpubdata/jwst/public/"
EPOCHS = (  # (UT date, program, i2d key without the band suffix, bands)
    (
        "2022-06-07",
        "2736 ERO",
        "jw02736/L3/t/o001/jw02736-o001_t001_nircam_clear-",
        ("f150w", "f444w"),
    ),
    ("2024-05-11", "4043", "jw04043/L3/t/o001/jw04043-o001_t001_nircam_clear-", ("f115w", "f444w")),
    (
        "2026-06-05",
        "6882 VENUS",
        "jw06882/L3/t/o057/jw06882-o057_t057_nircam_clear-",
        ("f150w", "f444w"),
    ),
    (
        "2026-07-06",
        "12774 DD",
        "jw12774/L3/t/o002/jw12774-o002_t032_nircam_clear-",
        ("f150w", "f444w"),
    ),
)
R_AP, R_IN, R_OUT = 0.15, 0.4, 0.7  # arcsec; ASSUMPTION: compact aperture against the arc light
MJYSR_ARCSEC2_TO_UJY = 1e12 / 4.25452e10  # 1 sr = 4.25452e10 arcsec^2


def aperture_flux(sci: np.ndarray, x0: float, y0: float, pix_arcsec: float) -> tuple[float, float]:
    """Background-subtracted aperture flux and a MAD noise estimate, in uJy, for a MJy/sr image."""
    yy, xx = np.mgrid[: sci.shape[0], : sci.shape[1]]
    rr = np.hypot(xx - x0, yy - y0) * pix_arcsec
    ap, an = rr < R_AP, (rr > R_IN) & (rr < R_OUT)
    bg = float(np.nanmedian(sci[an]))
    sig = 1.4826 * float(np.nanmedian(np.abs(sci[an] - bg)))
    k = pix_arcsec**2 * MJYSR_ARCSEC2_TO_UJY
    return float(np.nansum(sci[ap] - bg)) * k, sig * np.sqrt(ap.sum()) * k


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    targets = Table(
        {
            "source_uid": list(PAIR),
            "ra": [v[0] for v in PAIR.values()],
            "dec": [v[1] for v in PAIR.values()],
        }
    )
    cut_dir = paths.outputs_dir() / "smacs_transient_epochs"
    rows, stamps = [], []
    for date, prog, key, bands in EPOCHS:
        for band in bands:
            uri = f"{S3}{key}{band}_i2d.fits"
            ct = make_cutouts(
                uri, targets, size_arcsec=2.0, out_dir=cut_dir, storage_options={"anon": True}
            )
            for r in ct:
                with fits.open(r["path"]) as h:
                    hdu = h["SCI"] if "SCI" in h else h[0]
                    sci = np.asarray(hdu.data, dtype=float)
                    w = WCS(hdu.header)
                x0, y0 = w.world_to_pixel_values(*PAIR[r["source_uid"]])
                pix = float(
                    np.mean([abs(q.to_value("arcsec")) for q in w.proj_plane_pixel_scales()])
                )
                f, e = aperture_flux(sci, float(x0), float(y0), pix)
                rows.append(
                    {
                        "date": date,
                        "program": prog,
                        "band": band.upper(),
                        "image": r["source_uid"],
                        "flux_ujy_r015": round(f, 4),
                        "err_ujy": round(e, 4),
                        "quality_flag": str(r["quality_flag"]),
                    }
                )
                stamps.append((f"{date} {band.upper()} {r['source_uid']}", sci, x0, y0))
                print(json.dumps(rows[-1]))
    out = {
        "provenance": "observed (level-3 mosaic fluxes); err is a derived MAD estimate",
        "aperture": {
            "r_ap_arcsec": R_AP,
            "annulus_arcsec": [R_IN, R_OUT],
            "aperture_correction": None,
        },
        "pair": PAIR,
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n = len(stamps) // 2
        fig, axes = plt.subplots(2, n, figsize=(2 * n, 4.4))
        for i, (title, sci, x0, y0) in enumerate(stamps):
            a = axes[i % 2, i // 2]
            lo, hi = np.nanpercentile(sci, [5, 99.5])
            a.imshow(sci, origin="lower", vmin=lo, vmax=hi, cmap="gray")
            a.plot(x0, y0, "r+", ms=6)
            a.set_xticks([])
            a.set_yticks([])
            a.set_title(title, fontsize=6)
        fig.tight_layout()
        fig.savefig(cut_dir / "smacs_transient_epochs.png", dpi=120)
        print(cut_dir / "smacs_transient_epochs.png")
    except ImportError:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())

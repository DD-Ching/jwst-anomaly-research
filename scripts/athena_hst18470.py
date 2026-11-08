"""SN Athena reappearance check in public HST 18470 WFC3/IR F105W visits.

Plan and frozen rule: results/transients/athena_hst18470_v5/PLAN.md. Same photometry as
results/transients/confirm_r1/macs0417_hst18470/hst_check.py.txt (0.2" aperture, 0.6-1.0" annulus,
300 random-aperture robust errors) on each listed visit;
the rule compares the last listed visit to the first.

Usage: python scripts/athena_hst18470.py <data_dir> <out_dir> [visit ...]
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hit_lightcurve import aperture_photometry  # noqa: E402

POS = {"athena2_12781": (64.394375, -11.90758), "glafic_7p6": (64.394334, -11.907649)}
VISITS = ("ifu901010", "ifu902010", "ifu903010", "ifu904010", "ifu905010")
APPEAR_SIGMA = 5.0  # ASSUMPTION (PLAN.md)
NOT_SEEN_SIGMA = 3.0  # ASSUMPTION (PLAN.md)


def decide(delta: float, sigma_d: float, point_source_seen: bool) -> str:
    """Frozen rule of PLAN.md for one position: appeared / not seen / undetermined."""
    if not (np.isfinite(delta) and np.isfinite(sigma_d) and sigma_d > 0):
        return "undetermined"
    z = delta / sigma_d
    if z >= APPEAR_SIGMA and point_source_seen:
        return "appeared"
    if z < NOT_SEEN_SIGMA:
        return "not seen"
    return "undetermined"


def combine(decisions: list[str]) -> str:
    """Field-level outcome: appeared at either position wins; not seen needs both."""
    if "appeared" in decisions:
        return "appeared"
    if decisions and all(d == "not seen" for d in decisions):
        return "not seen"
    return "undetermined"


# Recorded after looking at stamps.png (PLAN.md: a point source must be visible for "appeared").
SEEN = {"athena2_12781": False, "glafic_7p6": False}


def evaluate(res: list[dict], seen: dict[str, bool]) -> dict:
    """Apply the frozen rule to the last listed visit against the first one."""
    first, last = res[0], res[-1]
    missing = [r["visit"] for r in (first, last) if "error" in r]
    if missing:
        return {"outcome": "undetermined", "error": f"visit not available: {missing}"}
    per = {}
    for name in POS:
        a, b = first[name], last[name]
        delta = b["flux_njy"] - a["flux_njy"]
        sd = float(np.hypot(a["err_njy"], b["err_njy"]))
        per[name] = {
            "delta_njy": round(delta, 1),
            "sigma_d_njy": round(sd, 1),
            "z": round(delta / sd, 2),
            "appeared_min_delta_njy": round(APPEAR_SIGMA * sd, 1),
            "not_seen_max_delta_njy": round(NOT_SEEN_SIGMA * sd, 1),
            "point_source_seen": seen[name],
            "decision": decide(delta, sd, seen[name]),
        }
    return {
        "reference_visit": first["visit"],
        "test_visit": last["visit"],
        "test_mjd": last["mjd"],
        "rule": f"appeared: z >= {APPEAR_SIGMA} and seen; not seen: z < {NOT_SEEN_SIGMA}",
        "positions": per,
        "outcome": combine([v["decision"] for v in per.values()]),
    }


def measure(data_dir: Path, visits: tuple[str, ...]):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astroquery.mast import Observations

    res, stamps = [], {}
    for v in visits:
        obs = Observations.query_criteria(obs_id=v)
        if not len(obs):
            res.append({"visit": v, "error": "not found"})
            continue
        p = Observations.get_product_list(obs)
        p = p[[str(x).endswith("_drz.fits") for x in p["productFilename"]]]
        dl = Observations.download_products(p[:1], download_dir=str(data_dir), cache=True)
        fn = dl["Local Path"][0]
        with fits.open(fn) as h:
            sci = h["SCI"].data.astype(float)
            hd, p0 = h["SCI"].header, h[0].header
            w = WCS(hd, h)
            photfnu = hd.get("PHOTFNU", p0.get("PHOTFNU"))
            sc = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600
            mjd, filt = p0.get("EXPSTART"), p0.get("FILTER")
        k = photfnu * 1e9  # e-/s -> nJy
        rec = {"visit": v, "filter": filt, "mjd": round(mjd, 2), "pixscale": round(sc, 4)}
        for name, (ra, dec) in POS.items():
            x, y = (float(t) for t in w.world_to_pixel_values(ra, dec))
            f = aperture_photometry(sci, x, y, 0.2 / sc, (0.6 / sc, 1.0 / sc)) * k
            rng = np.random.default_rng(0)
            rnd = []
            for _ in range(300):
                xx, yy = x + rng.uniform(-60, 60), y + rng.uniform(-60, 60)
                if np.hypot(xx - x, yy - y) * sc < 1.5:
                    continue
                q = aperture_photometry(sci, xx, yy, 0.2 / sc, (0.6 / sc, 1.0 / sc))
                if np.isfinite(q):
                    rnd.append(q * k)
            e = 1.4826 * np.median(np.abs(np.array(rnd) - np.median(rnd)))
            rec[name] = {"flux_njy": round(float(f), 1), "err_njy": round(float(e), 1)}
            iy, ix = int(round(y)), int(round(x))
            if 20 <= iy < sci.shape[0] - 20 and 20 <= ix < sci.shape[1] - 20:
                stamps[(v, name)] = sci[iy - 20 : iy + 21, ix - 20 : ix + 21] * k
            else:
                print(f"{v} {name}: stamp off the image edge", flush=True)
        res.append(rec)
        print(rec, flush=True)
    return res, stamps


def main() -> None:
    warnings.filterwarnings("ignore")
    data_dir, out = Path(sys.argv[1]), Path(sys.argv[2])
    visits = tuple(sys.argv[3:]) or VISITS
    out.mkdir(parents=True, exist_ok=True)
    res, stamps = measure(data_dir, visits)
    (out / "forced.json").write_text(json.dumps(res, indent=1) + "\n")
    summary = evaluate(res, SEEN)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, len(visits), figsize=(2.4 * len(visits), 5), squeeze=False)
    for j, v in enumerate(visits):
        for i, name in enumerate(POS):
            a = ax[i, j]
            s = stamps.get((v, name))
            if s is not None and s.size:
                ref = stamps.get((visits[0], name))
                lo, hi = np.nanpercentile(ref if ref is not None else s, (5, 99.5))
                a.imshow(s, origin="lower", cmap="gray", vmin=lo, vmax=hi)  # shared scale per row
                a.plot(20, 20, "r+")
            a.set_title(f"{v} {name}", fontsize=7)
            a.set_xticks([])
            a.set_yticks([])
    plt.tight_layout()
    plt.savefig(out / "stamps.png", dpi=90)
    print("wrote", out)


if __name__ == "__main__":
    main()

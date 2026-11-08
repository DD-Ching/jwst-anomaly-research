"""Does the RELICS GLAFIC v3 model of RXC J0600 treat row 581 as a cluster member?

Galaxy-scale κ peaks as a membership proxy, validated on NED spectroscopic redshifts
(``results/arcs/kappa-members-581/PLAN.md``, frozen before any κ value was read). The κ maps and the
NED table are read in the scratchpad and deleted; only ``summary.json`` is written.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import tempfile
from pathlib import Path

import numpy as np
from arcs_photoz_contam import NED_TAP, row_position
from arcs_shear_error import RELICS
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.wcs import WCS

OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "kappa-members-581"
MODELS = RELICS + "rxc0600m20/models/"
MAPS = {
    "glafic_v3": MODELS + "glafic/v3/hlsp_relics_model_model_rxc0600m20_glafic_v3_kappa.fits",
    "lenstool_v1": MODELS + "lenstool/hlsp_relics_model_model_rxc0600m20_lenstool_v1_kappa.fits",
}
DECIDING = "glafic_v3"
RUN = "rxc0600-relics/rxc0600-relics-glafic"
ROW = 581
CENTRE = (90.0340, -20.1358)
NED_RADIUS_DEG = 2.0 / 60.0
# Frozen rules (PLAN.md; ASSUMPTIONs)
Z_CL = 0.43
DZ = 0.03
SPEC_FLAGS = {"SLS", "S1L", "SST", "IST"}
TYPES = {"G", "G_Lens"}
BOX = 2.0  # arcsec, full width
PEAK_R = 0.4  # arcsec
PEAK_RATIO = 1.5
EDGE = 1.0  # arcsec
N_MIN = 8
RECALL_MIN = 0.70
BG_MAX = 0.10


def z_class(z: float) -> str:
    lim = DZ * (1 + Z_CL)
    if abs(z - Z_CL) <= lim:
        return "member"
    return "background" if z > Z_CL + lim else "foreground"


def spectroscopic(zflag: str, phytype: str) -> bool:
    return str(zflag).strip() in SPEC_FLAGS and str(phytype).strip() in TYPES


def peak_at(kappa: np.ndarray, x: float, y: float, pix: float) -> dict | None:
    """Peak statistic at 0-based pixel (x, y); None when within EDGE of the map edge."""
    ny, nx = kappa.shape
    e = EDGE / pix
    if x < e or y < e or x > nx - 1 - e or y > ny - 1 - e:
        return None
    h = BOX / 2 / pix
    x0, x1 = int(np.floor(x - h)), int(np.ceil(x + h)) + 1
    y0, y1 = int(np.floor(y - h)), int(np.ceil(y + h)) + 1
    sub = kappa[max(y0, 0) : min(y1, ny), max(x0, 0) : min(x1, nx)]
    iy, ix = np.unravel_index(np.nanargmax(sub), sub.shape)
    d = float(np.hypot(ix + max(x0, 0) - x, iy + max(y0, 0) - y)) * pix
    ratio = float(sub[iy, ix] / np.nanmedian(sub))
    return {
        "max": float(sub[iy, ix]),
        "dist_arcsec": round(d, 3),
        "ratio": round(ratio, 3),
        "peak": bool(d <= PEAK_R and ratio >= PEAK_RATIO),
    }


def gate(member_peaks: list[bool], bg_peaks: list[bool]) -> dict:
    nm, nb = len(member_peaks), len(bg_peaks)
    rec = float(np.mean(member_peaks)) if nm else float("nan")
    fb = float(np.mean(bg_peaks)) if nb else float("nan")
    ok = nm >= N_MIN and nb >= N_MIN and rec >= RECALL_MIN and fb <= BG_MAX
    return {"n_members": nm, "n_background": nb, "recall": rec, "bg_rate": fb, "passes": bool(ok)}


def label(gate_passes: bool, peak: bool | None) -> str:
    if not gate_passes or peak is None:
        return "undetermined"
    return "GLAFIC member (model input)" if peak else "GLAFIC background (model input)"


def fetch(url: str, dest: Path) -> str:
    import requests

    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        h = hashlib.sha256()
        with dest.open("wb") as f:
            for chunk in r.iter_content(1 << 20):
                h.update(chunk)
                f.write(chunk)
    return h.hexdigest()


def ned_table() -> Table:
    import requests

    q = (
        "SELECT prefname, ra, dec, z, zflag, prefphytype FROM NEDTAP.objdir WHERE "
        "CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {CENTRE[0]}, {CENTRE[1]}, {NED_RADIUS_DEG})) = 1 "
        "AND z IS NOT NULL"
    )
    p = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}
    r = requests.get(NED_TAP + "/sync", params=p, timeout=120)
    r.raise_for_status()
    return ascii.read(io.StringIO(r.text), format="csv")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    ned = ned_table()
    spec = [r for r in ned if spectroscopic(r["zflag"], r["prefphytype"])]
    pos = row_position(RUN, ROW)
    out: dict = {
        "params": {
            "z_cl": Z_CL,
            "dz": DZ,
            "spec_flags": sorted(SPEC_FLAGS),
            "types": sorted(TYPES),
            "box_arcsec": BOX,
            "peak_r_arcsec": PEAK_R,
            "peak_ratio": PEAK_RATIO,
            "edge_arcsec": EDGE,
            "n_min": N_MIN,
            "recall_min": RECALL_MIN,
            "bg_max": BG_MAX,
        },
        "ned": {"n_rows": len(ned), "n_spec": len(spec)},
        "row": {"label": ROW, "ra": pos["ra"], "dec": pos["dec"]},
        "maps": {},
    }
    with tempfile.TemporaryDirectory() as tmp:
        for name, url in MAPS.items():
            p = Path(tmp) / Path(url).name
            sha = fetch(url, p)
            with fits.open(p) as hd:
                kap = np.asarray(hd[0].data, float)
                w = WCS(hd[0].header)
            p.unlink()
            pix = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600.0
            per = []
            for r in spec:
                x, y = w.world_to_pixel_values(float(r["ra"]), float(r["dec"]))
                st = peak_at(kap, float(x), float(y), pix)
                if st is None:
                    continue
                per.append(
                    {
                        "name": str(r["prefname"]),
                        "z": float(r["z"]),
                        "class": z_class(float(r["z"])),
                        **st,
                    }
                )
            mem = [d["peak"] for d in per if d["class"] == "member"]
            bg = [d["peak"] for d in per if d["class"] == "background"]
            x, y = w.world_to_pixel_values(pos["ra"], pos["dec"])
            st581 = peak_at(kap, float(x), float(y), pix)
            g = gate(mem, bg)
            out["maps"][name] = {
                "url": url,
                "sha256": sha,
                "pixel_arcsec": round(pix, 5),
                "gate": g,
                "row_581": st581,
                "label": label(g["passes"], None if st581 is None else st581["peak"]),
                "galaxies": per,
            }
    out["outcome"] = out["maps"][DECIDING]["label"]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps(
            {k: {kk: v[kk] for kk in ("gate", "row_581", "label")} for k, v in out["maps"].items()},
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

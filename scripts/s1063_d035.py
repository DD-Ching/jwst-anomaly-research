"""Abell S1063 D-035 revisit: image-plane rms of the CATS v4.1 deflection map on the Chan et al.
2020 image list with its published system redshifts (frozen plan:
``results/lens_consistency/abells1063_d035/PLAN.md``).

The CATS files live on ``archive.stsci.edu``, whose chain ends at emSign Root CA - G1. Cloud
sessions' trust store lacks that root; Mozilla's (certifi) has it, so pass a bundle that holds
both, e.g. ``SSL_CERT_FILE=<proxy bundle + certifi bundle>``. Verification stays on.

    python scripts/s1063_d035.py
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel, schema  # noqa: E402
from jwst_anomaly.paths import repo_root  # noqa: E402

MODEL = "abells1063-cats"
LIST = repo_root() / "results/lens_consistency/abells1063_images/chan20_images.ecsv"
OUT = repo_root() / "results/lens_consistency/abells1063_d035"


@dataclass(frozen=True)
class Params:
    quoted_rms_arcsec: float = 0.48  # CATS v4.1 S1063 quoted rms (D-035 table)
    factor: float = lc.IMAGE_LIST_RMS_FACTOR  # ASSUMPTION (D-035): whole-list gate 1.5 x quoted
    system_max_rms: float = lc.GATE_MAX_RMS_ARCSEC  # ASSUMPTION: per-system opening, 0.5"
    grid_step: float = 0.25  # solver grid step, arcsec (validate's default)
    frame_warn_arcsec: float = 0.1  # report a median Chan - CATS offset above this
    match_arcsec: float = 1.0  # ASSUMPTION: Chan - CATS arcs.txt nearest-match radius


def chan_images(path: Path = LIST) -> Table:
    """The Chan list in ``load_lenstool_images`` form (``a`` = 1: no position error is used);
    ``dup_position`` rows are left out (VizieR defect, see the list's RESULT.md)."""
    t = Table.read(path, format="ascii.ecsv")
    t = t[~np.asarray(t["dup_position"], bool)]
    ids = [str(i) for i in t["image_id"]]
    out = Table(
        {
            "image_id": ids,
            "system": [lensmodel.image_family(i) for i in ids],
            "ra": np.asarray(t["ra"], float),
            "dec": np.asarray(t["dec"], float),
            "a": np.ones(len(t)),
            "z": np.asarray(t["z"], float),
            "z_spec": np.asarray(t["z_spec"], bool),
        }
    )
    src = path.relative_to(repo_root()) if path.is_relative_to(repo_root()) else path
    out.meta.update(provenance=schema.Provenance.OBSERVED.value, source=str(src))
    return out


def whole_list_gate(dtheta, quoted: float, factor: float, shared=None) -> dict:
    """R1: passes when every image is solved (finite dtheta) AND rms(dtheta) <= factor x quoted.
    ``n_shared_match`` (images whose nearest prediction is shared, i.e. a predicted image is
    missing or merged) is reported only; the frozen rule does not use it."""
    d = np.asarray(dtheta, float)
    ok = np.isfinite(d)
    rms = float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else float("nan")
    limit = factor * quoted
    return {
        "n_images": int(len(d)),
        "n_solved": int(ok.sum()),
        "n_shared_match": None if shared is None else int(np.sum(shared)),
        "rms_arcsec": rms,
        "limit_arcsec": limit,
        "passes": bool(ok.all() and len(d) > 0 and rms <= limit),
    }


def outcome(gate_s: dict, opened_s: list[str]) -> str:
    """Frozen outcome: the whole S list, else the S systems passing R2, else map only."""
    if gate_s["passes"]:
        return "open: all z_spec systems"
    if opened_s:
        return f"open: z_spec systems {', '.join(opened_s)}"
    return "map only"


def match_cats(chan: Table, cats: Table, max_sep: float):
    """Nearest CATS ``arcs.txt`` image to each Chan image: ``(idx, sep_arcsec, matched)``, where
    ``matched`` keeps sep <= ``max_sep`` and, for a CATS image nearest to several Chan images,
    only the closest one (one-to-one)."""
    a = SkyCoord(chan["ra"], chan["dec"], unit="deg")
    b = SkyCoord(cats["ra"], cats["dec"], unit="deg")
    idx, sep, _ = a.match_to_catalog_sky(b)
    sep = sep.arcsec
    m = sep <= max_sep
    for k in np.unique(idx[m]):
        both = np.where(m & (idx == k))[0]
        m[both[both != both[np.argmin(sep[both])]]] = False
    return idx, sep, m


def frame_check(chan: Table, cats: Table, max_sep: float) -> dict:
    """Median Chan - CATS offset over one-to-one matches within ``max_sep`` arcsec."""
    idx, sep, m = match_cats(chan, cats, max_sep)
    a = SkyCoord(chan["ra"][m], chan["dec"][m], unit="deg")
    b = SkyCoord(cats["ra"][idx[m]], cats["dec"][idx[m]], unit="deg")
    dra, ddec = a.spherical_offsets_to(b)
    return {
        "n_chan": len(chan),
        "n_matched": int(m.sum()),
        "median_sep_arcsec": float(np.median(sep[m])) if m.any() else None,
        "median_dra_arcsec": float(np.median(dra.arcsec)) if m.any() else None,
        "median_ddec_arcsec": float(np.median(ddec.arcsec)) if m.any() else None,
    }


def main() -> int:
    p = Params()
    model, files, _ = lc.load_model(MODEL)
    files.update(lc.model_files(MODEL, ("arcs.dat",)))
    cats = lensmodel.load_lenstool_images(files["arcs.dat"])
    chan = chan_images()
    summary: dict = {"model": MODEL, "model_sha256": model.sha256, "params": asdict(p)}
    summary["frame_check"] = frame_check(chan, cats, p.match_arcsec)
    fc = summary["frame_check"]
    summary["frame_offset_flag"] = bool(
        fc["median_sep_arcsec"] is not None and fc["median_sep_arcsec"] > p.frame_warn_arcsec
    )
    bt = lensmodel.backtrace_images(model, chan, {})
    half = lc.grid_half_width(model, chan)
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, p.grid_step), half, p.grid_step
    )
    ip = lensmodel.imageplane_residuals(model, grid, bt)
    ip["z_spec"] = chan["z_spec"]
    ip.meta["provenance"] = schema.Provenance.MODEL_PREDICTION.value
    zs = np.asarray(ip["z_spec"], bool)
    sh = np.asarray(ip["shared_match"], bool)
    gates = {
        "S": whole_list_gate(ip["dtheta_arcsec"][zs], p.quoted_rms_arcsec, p.factor, sh[zs]),
        "A": whole_list_gate(ip["dtheta_arcsec"], p.quoted_rms_arcsec, p.factor, sh),
    }
    opened = lc.systems_passing_rms(ip[zs], p.system_max_rms)
    per_system = {}
    systems = np.asarray(ip["system"], str)
    for s in dict.fromkeys(systems):
        sel = systems == s
        d = np.asarray(ip["dtheta_arcsec"][sel], float)
        per_system[s] = {
            "z_spec": bool(zs[sel][0]),
            "n": int(len(d)),
            "n_predicted": int(ip["n_predicted"][sel][0]),
            "rms_arcsec": float(np.sqrt(np.mean(d**2))) if np.isfinite(d).all() else None,
        }
    summary.update(
        gates=gates,
        r2_open_systems_S=opened,
        outcome=outcome(gates["S"], opened),
        per_system=per_system,
        grid={"half_width_arcsec": half, "step_arcsec": p.grid_step},
    )
    OUT.mkdir(parents=True, exist_ok=True)
    ip.write(OUT / "imageplane.ecsv", format="ascii.ecsv", overwrite=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("frame_check", "gates", "outcome")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

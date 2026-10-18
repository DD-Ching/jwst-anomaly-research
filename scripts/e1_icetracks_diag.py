"""Post hoc diagnostic for E1 + IceTracks-DR2 (not part of the frozen plan): the DR2-ICECAT control
channel's long-lag counts against four null variants, to locate the deficit seen in
``results/e1_icetracks/counts.ecsv``.

Variants: ``main`` (as the run: ICECAT by scramble_jit, DR2 by the uptime jitter); ``in_dr2``
(ICECAT alerts after the last DR2 track dropped); ``uptime`` (ICECAT also by the DR2 uptime jitter,
since its alerts come from the same detector); ``both``; ``dr2_only`` (ICECAT fixed, DR2 jittered)
and ``icecat_only`` (DR2 fixed, ICECAT by scramble_jit), to see which side carries the deficit;
``in_uptime`` (only ICECAT alerts inside a DR2 good run, both sides by the uptime jitter).

  python scripts/e1_icetracks_diag.py [--n 1000]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_icetracks as it  # noqa: E402
import e1_swift as sw  # noqa: E402
import nf_ghost_dr2 as dr2  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly import neutrino_frontier as nf  # noqa: E402

OUT = ROOT / "results" / "e1_icetracks" / "control_null_variants.json"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    a = ap.parse_args(argv)
    raw, up = dr2.fetch()
    ev = sw.load(sw.fetch(False))
    d = en.Sample.from_table(it.dr2_events(raw))
    ic_all = en.Sample.from_table(ev["ICECAT"])
    keep = ic_all.mjd <= d.mjd.max()
    ic_in = en.Sample(
        "ICECAT", ic_all.mjd[keep], ic_all.ra[keep], ic_all.dec[keep], ic_all.sigma[keep],
        ic_all.year[keep],
    )  # fmt: skip
    nlag = len(it.P.lag_edges) - 1
    out = {
        "provenance": "derived (post hoc diagnostic)",
        "n_scrambles": a.n,
        "labels": en.lag_labels(it.P),
        "icecat_all": len(ic_all.mjd),
        "icecat_in_dr2": int(keep.sum()),
    }
    inu = nf.in_uptime(ic_all.mjd, up[:, 0], up[:, 1])
    ic_up = en.Sample(
        "ICECAT", ic_all.mjd[inu], ic_all.ra[inu], ic_all.dec[inu], ic_all.sigma[inu],
        ic_all.year[inu],
    )  # fmt: skip
    out["icecat_in_uptime"] = int(inu.sum())
    out["icecat_north_in_uptime"] = int((inu & (ic_all.dec > it.DEC_MIN)).sum())
    out["icecat_north"] = int((ic_all.dec > it.DEC_MIN).sum())
    print({k: v for k, v in out.items() if k.startswith("icecat")}, flush=True)
    variants = (
        ("in_uptime", ic_up, True),
        ("main", ic_all, False),
        ("in_dr2", ic_in, False),
        ("uptime", ic_all, True),
        ("both", ic_in, True),
        ("dr2_only", ic_all, None),
        ("icecat_only", ic_all, "fix"),
    )
    for name, ic, upt in variants:
        obs = en.count_channel(d, ic, it.P, same=False)[:nlag, :2]
        rng = np.random.default_rng(8_000_000)
        null = []
        for _ in range(a.n):
            dj = (
                d
                if upt == "fix"
                else nf.jitter_uptime(d, up[:, 0], up[:, 1], it.P.jitter_s / en.DAY, rng)
            )
            if upt is None:
                ij = ic
            elif upt is True:
                ij = nf.jitter_uptime(ic, up[:, 0], up[:, 1], it.P.jitter_s / en.DAY, rng)
            else:
                ij = en.scramble_jit(ic, rng, it.P)
            null.append(en.count_channel(dj, ij, it.P, same=False)[:nlag, :2])
        null = np.array(null, float)
        z = (obs - null.mean(0)) / np.maximum(null.std(0), 1e-9)
        out[name] = {
            "obs_same_wide": obs.tolist(),
            "null_mean": np.round(null.mean(0), 2).tolist(),
            "z": np.round(z, 2).tolist(),
        }
        print(name, "wide z", np.round(z[:, 1], 2).tolist(), flush=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

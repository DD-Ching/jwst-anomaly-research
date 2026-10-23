"""Re-fetch the inputs of the STScI PSF check (results/w2_j2244_stpsf/sources.json) through the
download ledger and verify their sha256, so the run can be repeated from the recorded sources:

    NET_LEDGER_CAFILE=<proxy bundle + certifi roots> python scripts/w2_j2244_stpsf_fetch.py OUT_DIR
    python scripts/w2_j2244_stpsf.py OUT_DIR RESULT_DIR      (isolated env with photutils)

Inputs: 12 FLT exposures (S3), the STDPSF model (www.stsci.edu) and the Gaia DR3 rows (ESA TAP; a
byte mismatch is accepted when the parsed rows equal the recorded ones).
"""

import hashlib
import json
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import net_ledger as nl  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/w2_j2244_stpsf"
CAP = 1_000_000_000  # owner limit for this unit, bytes (all downloads)


def main(out: str) -> None:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    src = json.load(open(RES / "sources.json"))
    ledger = nl.Ledger(RES / "net_ledger.jsonl", CAP)
    items = [(f"{r}_flt.fits", v["url"], v["bytes"], v["sha256"]) for r, v in src["flt"].items()]
    p = src["psf_model"]
    items.append(("STDPSF_WFC3IR_F160W.fits", p["url"], p["bytes"], p["sha256"]))
    for name, url, size, sha in items:
        b = nl.fetch(url, ledger, est_bytes=size)
        ok = hashlib.sha256(b).hexdigest() == sha
        (out / name).write_bytes(b)
        print(name, len(b), "sha256 ok" if ok else "SHA256 MISMATCH")
        if not ok:
            raise SystemExit(f"{name}: sha256 differs from sources.json")
    g = src["gaia_dr3_tap"]
    url = (
        g["url"]
        + "?"
        + urllib.parse.urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "QUERY": g["query"]}
        )
    )
    b = nl.fetch(url, ledger, est_bytes=g["bytes"])
    if hashlib.sha256(b).hexdigest() != g["sha256"]:
        j = json.loads(b)
        names = [c["name"] for c in j["metadata"]]
        rows = [dict(zip(names, r, strict=True)) for r in j["data"]]
        if rows != g["rows"]:
            raise SystemExit("Gaia rows differ from sources.json")
        print("gaia3.json: bytes differ, rows equal")
    else:
        print("gaia3.json sha256 ok")
    (out / "gaia3.json").write_bytes(b)
    print("ledger used", ledger.used())


if __name__ == "__main__":
    main(sys.argv[1])

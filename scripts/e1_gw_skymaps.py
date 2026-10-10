"""E1: reduce the public GWTC sky maps to one nside-32 (NESTED) probability map per event
(D-074 next step; observed sky maps -> derived low-resolution maps).

The four Zenodo tarballs (GWTC-2.1, GWTC-3, GWTC-4.1, GWTC-5.0 parameter-estimation releases) are
streamed member by member and never stored (cloud disk rule); only the reduced maps are written
under the data root, with a manifest of the tarball digests. Flat (``PROB``, RING or NESTED) and
multi-order (``UNIQ`` + ``PROBDENSITY``) maps are both read. One map per event: the first member
in the preference order ``PREFER`` (ASSUMPTION: any of the release's waveform maps is adequate at
1.8 deg pixels).

Usage: python scripts/e1_gw_skymaps.py [--out DIR] [--tar KEY ...]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import pickle
import re
import tarfile
import time
from pathlib import Path

import numpy as np

NSIDE = 32
ORDER = 5  # log2(NSIDE)
NPIX = 12 * NSIDE**2

#: Zenodo record -> tarball (checked 2026-10-09).
TARS = {
    "GWTC-2.1": ("6513631", "IGWN-GWTC2p1-v2-PESkyMaps.tar.gz"),
    "GWTC-3": ("8177023", "IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz"),
    "GWTC-4.1": ("20275769", "IGWN-GWTC4p1-18965dda8_5-Archived_Skymaps.tar.gz"),
    "GWTC-5.0": ("20348005", "IGWN-GWTC5p0-29ebe06b7_25-Archived_Skymaps.tar.gz"),
}
#: Preferred map per event, first match wins (ASSUMPTION).
PREFER = ("Mixed", "IMRPhenomXPHM", "SEOBNRv4PHM", "")
MANIFEST = Path(__file__).resolve().parents[1] / "data" / "manifests" / "e1_gw_skymaps.ecsv"
EVENT_RE = re.compile(r"(GW\d{6}_\d{6})")


def to_nested32(table, header) -> np.ndarray:
    """A sky-map FITS table -> probability per nside-32 NESTED pixel (sums to 1)."""
    from astropy_healpix.core import ring_to_nested

    cols = {c.upper(): c for c in table.dtype.names}
    out = np.zeros(NPIX)
    if "UNIQ" in cols:
        uniq = np.asarray(table[cols["UNIQ"]], dtype=np.int64)
        dens = np.asarray(table[cols["PROBDENSITY"]], dtype=float)
        from astropy_healpix import uniq_to_level_ipix

        order, ipix = (np.asarray(x, dtype=np.int64) for x in uniq_to_level_ipix(uniq))
        prob = dens * 4 * np.pi / (12 * 4.0**order)
        hi = order >= ORDER
        np.add.at(out, ipix[hi] >> (2 * (order[hi] - ORDER)), prob[hi])
        for o in np.unique(order[~hi]):  # coarser pixels: spread over their children
            m = order == o
            k = 4 ** (ORDER - o)
            child = (ipix[m][:, None] * k + np.arange(k)[None, :]).ravel()
            np.add.at(out, child, np.repeat(prob[m] / k, k))
    else:
        prob = np.asarray(table[cols["PROB"]], dtype=float).ravel()
        nside = int(np.sqrt(prob.size / 12))
        if str(header.get("ORDERING", "RING")).upper().startswith("RING"):
            nest = np.empty_like(prob)
            nest[ring_to_nested(np.arange(prob.size), nside)] = prob
            prob = nest
        f = (nside // NSIDE) ** 2
        out = prob.reshape(NPIX, f).sum(1)
    out = np.clip(out, 0, None)
    return out / out.sum()


def maps_digest(names, prob) -> str:
    """Content digest of the reduced maps (the npz file itself carries zip timestamps)."""
    h = hashlib.sha256("\n".join(map(str, names)).encode())
    h.update(np.ascontiguousarray(prob, dtype=np.float32).tobytes())
    return h.hexdigest()


def _rank(name: str) -> int:
    for i, p in enumerate(PREFER):
        if p in name:
            return i
    return len(PREFER)


#: ASSUMPTION: range-request size and connections per tarball (Zenodo: ~0.4 MB/s per connection).
CHUNK = 8 << 20
CONNECTIONS = 12


def fetch_bytes(url: str, connections: int = CONNECTIONS) -> bytes:
    """Whole file into memory by parallel HTTP range requests (never written to disk)."""
    from concurrent.futures import ThreadPoolExecutor

    import requests

    size = None
    for attempt in range(6):  # size from a 1-byte range request; the proxy drops some requests
        try:
            r = requests.get(url, headers={"Range": "bytes=0-0"}, timeout=60, stream=True)
            if r.status_code == 206:
                size = int(r.headers["Content-Range"].split("/")[1])
                break
            r.close()  # a server that ignores Range would send the whole body
        except (requests.RequestException, KeyError, ValueError):
            time.sleep(2**attempt)
    if size is None:
        raise OSError(f"no size for {url}")
    buf = bytearray(size)

    def get(lo: int) -> None:
        hi = min(lo + CHUNK, size) - 1
        for attempt in range(6):
            try:
                r = requests.get(url, headers={"Range": f"bytes={lo}-{hi}"}, timeout=120)
                if r.status_code == 206 and len(r.content) == hi - lo + 1:
                    buf[lo : hi + 1] = r.content
                    return
            except requests.RequestException:
                pass
            time.sleep(2**attempt)
        raise OSError(f"range {lo}-{hi} of {url} failed")

    with ThreadPoolExecutor(connections) as ex:
        list(ex.map(get, range(0, size, CHUNK)))
    return buf  # no bytes() copy: peak memory stays one archive


def reduce_tar(key: str) -> dict:
    """Fetch one tarball into memory; return {event: (rank, member, map)}, its sha256 and size."""
    from astropy.io import fits

    rec, fname = TARS[key]
    url = f"https://zenodo.org/records/{rec}/files/{fname}?download=1"
    t0 = time.time()
    data = fetch_bytes(url)
    print(f"{key}: {len(data) / 1e6:.0f} MB in {time.time() - t0:.0f} s", flush=True)
    best: dict[str, tuple[int, str, np.ndarray]] = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
        for m in tf:
            if not m.isfile() or not m.name.endswith((".fits", ".fits.gz", ".fits.fz")):
                continue
            ev = EVENT_RE.search(m.name)
            if ev is None:
                continue
            rank = _rank(Path(m.name).name)
            if ev[1] in best and best[ev[1]][0] <= rank:
                continue
            raw = tf.extractfile(m).read()
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            try:
                with fits.open(io.BytesIO(raw)) as h:
                    hdu = h[1]
                    best[ev[1]] = (rank, m.name, to_nested32(hdu.data, hdu.header))
            except Exception as e:  # noqa: BLE001 - a bad member is logged, not fatal
                print(f"{key}: skip {m.name}: {e}", flush=True)
    sha = hashlib.sha256(data).hexdigest()
    return {
        "key": key,
        "maps": best,
        "sha256": sha,
        "bytes": len(data),
        "url": url,
        "retrieved": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t0)),
    }


def durable_copy(out: Path) -> bool:
    """S-1: fetch the published reduced maps (``derived_store``) and accept them only when their
    content digest equals the one pinned in ``MANIFEST``. False means: reduce from the source."""
    from astropy.table import Table

    from jwst_anomaly import derived_store

    npz = out / "gw_skymaps_nside32.npz"
    if not MANIFEST.exists():
        return False
    want = Table.read(MANIFEST, format="ascii.ecsv").meta.get("maps_sha256")

    def ok() -> bool:
        if want is None or not npz.exists():
            return False
        z = np.load(npz)
        return maps_digest(z["name"], z["prob"]) == want

    # a local reduction is accepted by its content digest (npz bytes differ between runs)
    return ok() or (derived_store.fetch(npz.name, npz) is not None and ok())


def main(argv=None) -> None:
    from astropy.table import Table

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=None, help="default: <data root>/e1_events")
    ap.add_argument("--tar", nargs="*", default=list(TARS))
    ap.add_argument(
        "--refresh", action="store_true", help="skip the durable copy; reduce from the tarballs"
    )
    a = ap.parse_args(argv)
    if a.out is None:
        from jwst_anomaly.paths import data_root

        a.out = data_root() / "e1_events"
    a.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if not a.refresh and set(a.tar) == set(TARS) and durable_copy(a.out):
        print(f"durable copy verified in {time.time() - t0:.0f} s; nothing to reduce")
        return
    # one tarball per call; for speed run one process per --tar KEY, then a plain run merges the
    # caches. Each reduced tarball is cached so a failed fetch never loses the others.
    res = []
    for key in a.tar:
        cache = a.out / f"gw_skymaps_{key}.pkl"
        if cache.exists():
            res.append(pickle.loads(cache.read_bytes()))
            continue
        r = reduce_tar(key)
        cache.write_bytes(pickle.dumps(r))
        res.append(r)
    if set(a.tar) != set(TARS):
        print(f"cached {a.tar}; run without --tar to write the merged maps and the manifest")
        return
    names, members, cats, maps = [], [], [], []
    man = Table(
        names=("catalog", "url", "sha256", "bytes", "n_events", "retrieved"),
        dtype=(str, str, str, int, int, str),
    )
    for r in res:
        man.add_row(
            (
                r["key"],
                r["url"],
                r["sha256"],
                r["bytes"],
                len(r["maps"]),
                r.get("retrieved", "2026-10-10 (cache without a timestamp)"),
            )
        )
        for ev, (_, mem, p) in sorted(r["maps"].items()):
            names.append(ev)
            members.append(mem)
            cats.append(r["key"])
            maps.append(p.astype(np.float32))
    npz = a.out / "gw_skymaps_nside32.npz"
    np.savez_compressed(
        npz,
        name=np.array(names),
        member=np.array(members),
        catalog=np.array(cats),
        prob=np.array(maps),
    )
    man.meta["provenance"] = (
        "observed (GWTC PE sky maps), streamed; reduced to nside 32 NESTED: derived"
    )
    man.meta["reason"] = (
        "E1 GW channels need sky directions (D-074 addendum 4); the maps exist only in these "
        "whole-release tarballs, which are streamed into memory and never stored"
    )
    man.meta["preference"] = list(PREFER)
    man.meta["maps_sha256"] = maps_digest(np.array(names), np.array(maps))
    man.write(MANIFEST, overwrite=True)
    print(man)
    print(f"{len(names)} maps in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

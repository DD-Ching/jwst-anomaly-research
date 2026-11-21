"""VizieR server selection: CDS first, then the CfA mirror (D-006 addendum 2026-11-21).

CDS (``vizier.cds.unistra.fr``) has answered HTTP 404 with a generic "Page unavailable" page since
2026-11-10 (#819), while the CfA mirror (``vizier.cfa.harvard.edu``, astroquery's documented
alternative server) serves the same catalogues. Callers get an astroquery ``Vizier`` bound to the
first server that answers a one-row probe. The mirror is only a fallback: CDS stays the reference
whenever it answers, because the mirror may lag it.

``JWST_ANOMALY_VIZIER_SERVER`` pins a server and skips the probe (tests, reproducing a run).
"""

from __future__ import annotations

import logging
import os
import threading

log = logging.getLogger(__name__)

CDS_VIZIER = "vizier.cds.unistra.fr"
CFA_VIZIER = "vizier.cfa.harvard.edu"
VIZIER_SERVERS = (CDS_VIZIER, CFA_VIZIER)
ENV_VAR = "JWST_ANOMALY_VIZIER_SERVER"
# One Hipparcos row: small, stable since 1997, present on every VizieR mirror.
PROBE_QUERY = "/viz-bin/asu-tsv?-source=I/239/hip_main&-out=HIP&-out.max=1&HIP=1"
PROBE_TIMEOUT_S = 15.0  # ASSUMPTION: a live server answers the one-row probe well inside 15 s
PROBE_ATTEMPTS = 2  # ASSUMPTION: the CfA mirror occasionally fails one request (seen 2026-11-21)


def probe_ok(text: str) -> bool:
    """True when an asu-tsv probe body carries the HIP 1 data row (not an HTML error page)."""
    if "<html" in text[:500].lower():
        return False
    return any(line.strip() == "1" for line in text.splitlines() if not line.startswith("#"))


def _probe(server: str, timeout: float) -> bool:
    import requests

    try:
        r = requests.get(f"https://{server}{PROBE_QUERY}", timeout=timeout)
    except requests.RequestException as exc:
        log.info("VizieR probe %s failed: %s", server, exc)
        return False
    ok = r.status_code == 200 and probe_ok(r.text)
    if not ok:
        log.info("VizieR probe %s: HTTP %s, no data row", server, r.status_code)
    return ok


def _bare_host(value: str) -> str:
    """``https://host/`` or ``host/`` -> ``host`` (astroquery adds the scheme itself)."""
    v = value.strip()
    for scheme in ("https://", "http://"):
        if v.lower().startswith(scheme):
            v = v[len(scheme) :]
    return v.split("/", 1)[0]


_LOCK = threading.Lock()
_SELECTED: list[str] = []  # the probed server for this process (at most one entry)


def vizier_server() -> str:
    """``$JWST_ANOMALY_VIZIER_SERVER`` if set, else the first answering ``VIZIER_SERVERS`` entry.

    The probe runs once per process under a lock (threads wait for it); :func:`forget` drops the
    choice so the next call probes again. Raises ``ConnectionError`` when no server answers; that
    failure is not cached.
    """
    pinned = _bare_host(os.environ.get(ENV_VAR, ""))
    if pinned:
        return pinned
    with _LOCK:
        if not _SELECTED:
            _SELECTED.append(_first_answering(PROBE_TIMEOUT_S))
        return _SELECTED[0]


def forget() -> None:
    """Drop the probed server (call after a query on it failed, so a retry can switch hosts)."""
    with _LOCK:
        _SELECTED.clear()


def _first_answering(timeout: float) -> str:
    for server in VIZIER_SERVERS:
        if any(_probe(server, timeout) for _ in range(PROBE_ATTEMPTS)):
            if server != CDS_VIZIER:
                log.warning("CDS VizieR unreachable; using the mirror %s", server)
            return server
    raise ConnectionError(f"no VizieR server answered: {', '.join(VIZIER_SERVERS)}")


def vizier(**kwargs):
    """``astroquery.vizier.Vizier(**kwargs)`` bound to :func:`vizier_server`."""
    from astroquery.vizier import Vizier

    return Vizier(vizier_server=vizier_server(), **kwargs)


def tag(source: str) -> str:
    """``source`` with the serving VizieR host appended, for ``meta["source"]`` provenance."""
    return f"{source} [VizieR host {vizier_server()}]"


def asu_tsv_url() -> str:
    """The asu-tsv endpoint on the selected server (for scripts that call it directly)."""
    return f"https://{vizier_server()}/viz-bin/asu-tsv"

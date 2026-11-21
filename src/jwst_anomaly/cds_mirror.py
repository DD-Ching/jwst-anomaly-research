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
from functools import lru_cache

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


def vizier_server(timeout: float = PROBE_TIMEOUT_S) -> str:
    """``$JWST_ANOMALY_VIZIER_SERVER`` if set, else the first answering ``VIZIER_SERVERS`` entry.

    The probe result is cached for the process. Raises ``ConnectionError`` when none answers.
    """
    pinned = os.environ.get(ENV_VAR, "").strip()
    return pinned or _first_answering(timeout)


@lru_cache(maxsize=1)
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


def asu_tsv_url() -> str:
    """The asu-tsv endpoint on the selected server (for scripts that call it directly)."""
    return f"https://{vizier_server()}/viz-bin/asu-tsv"

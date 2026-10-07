#!/usr/bin/env python
"""Check that the http(s) links in Markdown files resolve.

Standard library only, so it runs before the project environment exists. Usage::

    python scripts/check_links.py SOURCES.md docs/landscape.md DECISIONS.md

Each unique URL gets a HEAD request with redirects followed. If HEAD does not return 2xx (many
servers answer 403/405, loop, or drop the connection), a GET is tried instead. A final 2xx, or a
3xx that is not an unfollowable redirect (e.g. 300, 304), counts as OK. Anything else (4xx, 5xx,
redirect loops, redirects without a usable Location, network errors) is a failure, and the exit
code is 1. Requests to one host are limited (default 2 at a time) and 429/503 responses are retried
after ``Retry-After``. Not checked: URLs in fenced code blocks, single-line code spans and HTML
comments; localhost/loopback URLs; templated URLs such as ``https://host/<name>`` or
``https://host/{id}``.
"""

from __future__ import annotations

import argparse
import http.client
import http.cookiejar
import ipaddress
import re
import socket
import ssl
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

USER_AGENT = (
    "Mozilla/5.0 (compatible; jwst-anomaly-check-links/1.0; "
    "+https://github.com/DD-Ching/jwst-anomaly-research)"
)
DEFAULT_TIMEOUT = 20.0
DEFAULT_WORKERS = 6
DEFAULT_PER_HOST = 2
MAX_WAIT = 30.0
RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
REDIRECT_STATUSES = frozenset({301, 302, 303, 307, 308})

# CommonMark: a backtick fence's info string cannot contain backticks (so ```x``` is inline code).
_FENCED_CODE = re.compile(
    r"^[ \t]{0,3}(?P<f>`{3,})[^`\n]*\n.*?^[ \t]{0,3}(?P=f)`*[ \t]*$"
    r"|^[ \t]{0,3}(?P<t>~{3,})[^\n]*\n.*?^[ \t]{0,3}(?P=t)~*[ \t]*$",
    re.M | re.S,
)
_INLINE_CODE = re.compile(r"(`+)(?!`).+?(?<!`)\1(?!`)")  # single line: a stray ` can't eat text
_HTML_COMMENT = re.compile(r"<!--.*?-->", re.S)
# Stop at whitespace and characters that delimit URLs in Markdown/HTML. Parentheses are allowed;
# unbalanced trailing ones are trimmed afterwards (Markdown link syntax).
_URL = re.compile(r"""https?://[^\s<>"'`|\[\]{}\\]+""", re.I)
_TRAILING_PUNCT = ".,;:!?*_~"
# What follows a URL match decides whether it is a template: "<name>" / "{id}" placeholders are,
# HTML tags ("<br>", "</td>") and kramdown attributes ("{.class}", "{: ...}") are not.
_PLACEHOLDER = re.compile(r"<(?P<tag>[A-Za-z][\w-]*)>|\{[A-Za-z_][\w-]*\}")
_HTML_TAGS = frozenset(
    "a abbr b br code del details div em hr i img ins kbd li ol p pre s small span strong sub "
    "summary sup table tbody td th thead tr u ul".split()
)
_LOCAL_HOSTS = frozenset({"localhost", "localhost.localdomain"})

_SSL_CONTEXT = ssl.create_default_context()
_host_lock = threading.Lock()
_host_slots: dict[str, threading.BoundedSemaphore] = {}


@dataclass(frozen=True)
class Result:
    url: str
    ok: bool
    status: int | None = None
    method: str = ""
    final_url: str = ""
    error: str = ""
    retry_after: float | None = None


class Response(NamedTuple):
    status: int
    final_url: str
    retry_after: float | None = None
    redirect_error: str = ""  # set when a 3xx could not be followed (loop, no Location, ...)


def _strip_code(text: str) -> str:
    text = _FENCED_CODE.sub(" ", text)
    text = _HTML_COMMENT.sub(" ", text)
    return _INLINE_CODE.sub(" ", text)


def _trim(url: str) -> str:
    while url:
        if url[-1] in _TRAILING_PUNCT:
            url = url[:-1]
        elif url[-1] == ")" and url.count(")") > url.count("("):
            url = url[:-1]
        else:
            break
    return url


def _is_template(text: str, end: int) -> bool:
    m = _PLACEHOLDER.match(text, end)
    if not m:
        return False
    return m.group("tag") is None or m.group("tag").lower() not in _HTML_TAGS


def extract_urls(text: str) -> list[str]:
    """Return the unique http(s) URLs in Markdown ``text``, in order of first appearance."""
    text = _strip_code(text)
    seen: dict[str, None] = {}
    for m in _URL.finditer(text):
        if _is_template(text, m.end()):
            continue  # e.g. https://archive.stsci.edu/hlsp/<name>: a pattern, not a link
        url = _trim(m.group(0))
        if urllib.parse.urlsplit(url).hostname:
            seen.setdefault(url, None)
    return list(seen)


def is_local(url: str) -> bool:
    """True for localhost, loopback and unspecified addresses (never checked)."""
    host = (urllib.parse.urlsplit(url).hostname or "").lower()
    if host in _LOCAL_HOSTS or host.endswith(".localhost"):
        return True
    try:
        addr = ipaddress.ip_address(host)
    except ValueError:
        return False
    return addr.is_loopback or addr.is_unspecified


class _KeepMethodRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Follow redirects without turning HEAD into GET (urllib does that before Python 3.13)."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is not None and req.get_method() == "HEAD":
            new.method = "HEAD"
        return new


def _quote(url: str) -> str:
    # Percent-encode non-ASCII characters and spaces; keep reserved characters as written.
    return urllib.parse.quote(url, safe=":/?#[]@!$&'()*+,;=%~")


def _retry_after(headers) -> float | None:
    value = (headers or {}).get("Retry-After", "")
    try:
        return max(0.0, float(value))
    except (TypeError, ValueError):
        return None  # HTTP-date form or absent: fall back to our own backoff


def _request(url: str, method: str, timeout: float) -> Response:
    """Send one request, following redirects. HTTP error statuses are returned, not raised.

    Network-level failures (DNS, TLS, timeouts, resets) raise.
    """
    opener = urllib.request.build_opener(
        urllib.request.HTTPSHandler(context=_SSL_CONTEXT),
        _KeepMethodRedirectHandler(),
        urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()),
    )
    req = urllib.request.Request(
        _quote(url),
        method=method,
        headers={"User-Agent": USER_AGENT, "Accept": "*/*", "Accept-Language": "en"},
    )
    try:
        with opener.open(req, timeout=timeout) as resp:
            return Response(resp.status, resp.geturl())
    except urllib.error.HTTPError as exc:
        exc.close()
        # urllib follows 301/302/303/307/308 itself, so one surfacing here was not followable.
        redirect_error = str(exc.reason) if exc.code in REDIRECT_STATUSES else ""
        return Response(exc.code, exc.geturl() or url, _retry_after(exc.headers), redirect_error)


def _transient(exc: BaseException) -> bool:
    """Whether a network error might succeed on retry (timeouts, resets) or never will."""
    reason = getattr(exc, "reason", None)
    for err in (exc, reason):
        if isinstance(err, (socket.gaierror, ssl.SSLCertVerificationError, ValueError)):
            return False
    return True


def _attempt(url: str, method: str, timeout: float, retries: int) -> Result:
    retries = max(0, retries)
    for attempt in range(retries + 1):
        try:
            resp = _request(url, method, timeout)
        except (urllib.error.URLError, http.client.HTTPException, OSError, ValueError) as exc:
            reason = getattr(exc, "reason", None) or exc
            res = Result(url, False, None, method, "", f"{type(exc).__name__}: {reason}")
            if attempt < retries and _transient(exc):
                time.sleep(min(2.0 * (attempt + 1), MAX_WAIT))
                continue
            return res
        ok = 200 <= resp.status < 400 and not resp.redirect_error
        res = Result(
            url, ok, resp.status, method, resp.final_url, resp.redirect_error, resp.retry_after
        )
        if resp.status in RETRY_STATUSES and attempt < retries:
            time.sleep(min(resp.retry_after or 2.0 * (attempt + 1), MAX_WAIT))
            continue
        return res
    return res


def _host_slot(url: str, limit: int) -> threading.BoundedSemaphore:
    host = (urllib.parse.urlsplit(url).hostname or "").lower()
    with _host_lock:
        if host not in _host_slots:
            _host_slots[host] = threading.BoundedSemaphore(max(1, limit))
        return _host_slots[host]


def check_url(
    url: str,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = 1,
    per_host: int = DEFAULT_PER_HOST,
) -> Result:
    """HEAD first; unless that gives 2xx, fall back to GET (servers that reject HEAD)."""
    with _host_slot(url, per_host):
        head = _attempt(url, "HEAD", timeout, retries=0)
        if head.ok and head.status is not None and head.status < 300:
            return head
        if head.status in RETRY_STATUSES:  # rate-limited: wait before asking again
            time.sleep(min(head.retry_after or 2.0, MAX_WAIT))
        return _attempt(url, "GET", timeout, retries)


def collect(paths: list[Path]) -> dict[str, list[str]]:
    """Map each unique URL to the files that contain it."""
    found: dict[str, list[str]] = {}
    for path in paths:
        for url in extract_urls(path.read_text(encoding="utf-8")):
            found.setdefault(url, []).append(str(path))
    return found


def _format(res: Result) -> str:
    code = str(res.status) if res.status is not None else "ERR"
    line = f"{'OK' if res.ok else 'FAIL':4} {code:>3} {res.method:4} {res.url}"
    if res.final_url and res.final_url.rstrip("/") != _quote(res.url).rstrip("/"):
        line += f" -> {res.final_url}"
    if res.error:
        line += f"  [{res.error}]"
    return line


def _non_negative_int(value: str) -> int:
    n = int(value)
    if n < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return n


def _positive_int(value: str) -> int:
    n = int(value)
    if n < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return n


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("files", nargs="+", type=Path, help="Markdown files to scan")
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT, help="seconds/request")
    parser.add_argument("--workers", type=_positive_int, default=DEFAULT_WORKERS)
    parser.add_argument("--per-host", type=_positive_int, default=DEFAULT_PER_HOST)
    parser.add_argument("--retries", type=_non_negative_int, default=1, help="GET retries")
    parser.add_argument(
        "--exclude", action="append", default=[], metavar="REGEX", help="skip matching URLs"
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="print failures only")
    args = parser.parse_args(argv)

    missing = [str(p) for p in args.files if not p.is_file()]
    if missing:
        print(f"error: not a file: {', '.join(missing)}", file=sys.stderr)
        return 2

    found = collect(args.files)
    excludes = [re.compile(p) for p in args.exclude]
    skipped: list[str] = []
    todo: list[str] = []
    for url in found:
        local_or_excluded = is_local(url) or any(p.search(url) for p in excludes)
        (skipped if local_or_excluded else todo).append(url)

    def run(url: str) -> Result:
        return check_url(url, args.timeout, args.retries, args.per_host)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(run, todo))

    failures = [r for r in results if not r.ok]
    for res in sorted(results, key=lambda r: (r.ok, r.url)):
        if res.ok and args.quiet:
            continue
        print(_format(res))
        if not res.ok:
            print(f"          in: {', '.join(found[res.url])}")
    for url in skipped:
        if not args.quiet:
            print(f"SKIP     {url}")
    n_ok = len(results) - len(failures)
    print(
        f"\n{len(found)} unique URLs in {len(args.files)} files: {n_ok} ok, "
        f"{len(failures)} failed, {len(skipped)} skipped."
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

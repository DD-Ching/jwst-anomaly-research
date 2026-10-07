"""Offline tests for scripts/check_links.py (no internet; one test uses a loopback server)."""

from __future__ import annotations

import importlib.util
import socket
import sys
import threading
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_links.py"
_spec = importlib.util.spec_from_file_location("check_links", _SCRIPT)
check_links = importlib.util.module_from_spec(_spec)
sys.modules["check_links"] = check_links  # dataclasses need the module registered
_spec.loader.exec_module(check_links)


# --- URL extraction ---------------------------------------------------------------------------


def test_extract_markdown_forms():
    text = (
        "See [docs](https://example.org/docs) and <https://example.org/auto>.\n"
        "Bare https://example.org/bare, then (https://example.org/paren).\n"
        '[t](https://example.org/titled "Title") | https://example.org/cell |\n'
        "**https://example.org/bold** ~~https://example.org/struck~~ _https://example.org/em_\n"
    )
    assert check_links.extract_urls(text) == [
        "https://example.org/docs",
        "https://example.org/auto",
        "https://example.org/bare",
        "https://example.org/paren",
        "https://example.org/titled",
        "https://example.org/cell",
        "https://example.org/bold",
        "https://example.org/struck",
        "https://example.org/em",
    ]


def test_extract_keeps_balanced_parentheses_and_queries():
    text = "[w](https://en.wikipedia.org/wiki/Lens_(optics)) http://example.org/q?a=1&b=2#frag."
    assert check_links.extract_urls(text) == [
        "https://en.wikipedia.org/wiki/Lens_(optics)",
        "http://example.org/q?a=1&b=2#frag",
    ]


def test_extract_skips_code_comments_and_templates():
    text = (
        "```python\nurl = 'https://example.org/in-fence'\n```\n"
        "~~~\nhttps://example.org/in-tilde-fence\n~~~\n"
        "Inline `https://example.org/in-code` and <!-- https://example.org/comment -->.\n"
        "Pattern https://archive.stsci.edu/hlsp/<name> and https://host.org/{id}/x.\n"
        "Not a URL: s3://stpubdata/jwst and ftp://example.org/x; https:// alone.\n"
        "Kept: https://example.org/kept\n"
    )
    assert check_links.extract_urls(text) == ["https://example.org/kept"]


def test_extract_html_tags_and_attributes_are_not_templates():
    text = (
        "| https://a.org/x<br>https://b.org/y |\n"
        "https://c.org/z</td> https://d.org/w<sup>1</sup> [l](https://e.org/v){.cls}\n"
    )
    assert check_links.extract_urls(text) == [
        "https://a.org/x",
        "https://b.org/y",
        "https://c.org/z",
        "https://d.org/w",
        "https://e.org/v",
    ]


def test_inline_triple_backticks_are_not_a_fence():
    text = (
        "```inline``` then https://a.org/kept\nmore https://a.org/kept2\n"
        "```\nhttps://a.org/in-fence\n```\nafter https://a.org/after\n"
    )
    assert check_links.extract_urls(text) == [
        "https://a.org/kept",
        "https://a.org/kept2",
        "https://a.org/after",
    ]


def test_extract_dedupes_in_order():
    text = "https://b.org/x https://a.org/y https://b.org/x."
    assert check_links.extract_urls(text) == ["https://b.org/x", "https://a.org/y"]


@pytest.mark.parametrize(
    "url, local",
    [
        ("http://localhost:8080/x", True),
        ("http://LOCALHOST/x", True),
        ("http://api.localhost/x", True),
        ("http://127.0.0.1:5000/", True),
        ("http://[::1]/", True),
        ("http://0.0.0.0/", True),
        ("https://archive.stsci.edu/", False),
        ("https://10.0.0.1/", False),
    ],
)
def test_is_local(url, local):
    assert check_links.is_local(url) is local


# --- HTTP logic with a monkeypatched transport ------------------------------------------------


@pytest.fixture
def sleeps(monkeypatch):
    calls = []
    monkeypatch.setattr(check_links.time, "sleep", calls.append)
    return calls


def _fake_transport(monkeypatch, responses):
    """responses: {(url, method): status | Response | exception}; returns the call log."""
    calls = []

    def fake_request(url, method, timeout):
        calls.append((url, method))
        outcome = responses[(url, method)]
        if isinstance(outcome, BaseException):
            raise outcome
        if isinstance(outcome, check_links.Response):
            return outcome
        return check_links.Response(outcome, url)

    monkeypatch.setattr(check_links, "_request", fake_request)
    return calls


def test_head_success_needs_no_get(monkeypatch, sleeps):
    calls = _fake_transport(monkeypatch, {("https://a.org", "HEAD"): 200})
    res = check_links.check_url("https://a.org")
    assert res.ok and res.status == 200 and res.method == "HEAD"
    assert calls == [("https://a.org", "HEAD")]


@pytest.mark.parametrize("head", [405, 403, 404, ConnectionResetError("reset")])
def test_head_rejected_falls_back_to_get(monkeypatch, sleeps, head):
    _fake_transport(monkeypatch, {("u", "HEAD"): head, ("u", "GET"): 200})
    res = check_links.check_url("u")
    assert res.ok and res.status == 200 and res.method == "GET"


def test_non_redirect_3xx_counts_as_ok(monkeypatch, sleeps):
    calls = _fake_transport(monkeypatch, {("u", "HEAD"): 304, ("u", "GET"): 300})
    res = check_links.check_url("u")
    assert res.ok and res.status == 300
    assert calls == [("u", "HEAD"), ("u", "GET")]  # a HEAD 3xx is double-checked with GET


def test_unfollowable_redirect_fails(monkeypatch, sleeps):
    loop = check_links.Response(302, "u", None, "redirect loop")
    _fake_transport(monkeypatch, {("u", "HEAD"): loop, ("u", "GET"): loop})
    res = check_links.check_url("u")
    assert not res.ok and res.status == 302 and "loop" in res.error


def test_not_found_fails(monkeypatch, sleeps):
    _fake_transport(monkeypatch, {("u", "HEAD"): 404, ("u", "GET"): 404})
    res = check_links.check_url("u")
    assert not res.ok and res.status == 404


def test_transient_network_error_is_retried(monkeypatch, sleeps):
    err = urllib.error.URLError(TimeoutError("timed out"))
    calls = _fake_transport(monkeypatch, {("u", "HEAD"): err, ("u", "GET"): err})
    res = check_links.check_url("u", retries=1)
    assert not res.ok and res.status is None and "timed out" in res.error
    assert calls == [("u", "HEAD"), ("u", "GET"), ("u", "GET")]


@pytest.mark.parametrize(
    "err",
    [
        urllib.error.URLError(socket.gaierror(11001, "getaddrinfo failed")),
        ValueError("unknown url type"),
    ],
)
def test_permanent_network_error_is_not_retried(monkeypatch, sleeps, err):
    calls = _fake_transport(monkeypatch, {("u", "HEAD"): err, ("u", "GET"): err})
    assert not check_links.check_url("u", retries=3).ok
    assert calls == [("u", "HEAD"), ("u", "GET")]
    assert sleeps == []


def test_rate_limit_waits_for_retry_after(monkeypatch, sleeps):
    seq = iter([check_links.Response(429, "u", retry_after=7.0), check_links.Response(200, "u")])

    def fake_request(url, method, timeout):
        if method == "HEAD":
            return check_links.Response(429, url, retry_after=3.0)
        return next(seq)

    monkeypatch.setattr(check_links, "_request", fake_request)
    assert check_links.check_url("u", retries=1).ok
    assert sleeps == [3.0, 7.0]


def test_negative_retries_are_clamped(monkeypatch, sleeps):
    _fake_transport(monkeypatch, {("u", "HEAD"): 405, ("u", "GET"): 503})
    res = check_links.check_url("u", retries=-1)
    assert not res.ok and res.status == 503


def test_main_reports_failures_and_skips_local(monkeypatch, tmp_path, capsys):
    a = tmp_path / "a.md"
    b = tmp_path / "b.md"
    a.write_text("[ok](https://ok.org) https://bad.org/x http://localhost:8000/\n", "utf-8")
    b.write_text("again https://ok.org and https://skip.me/y\n", "utf-8")
    checked = []

    def fake_check(url, timeout, retries, per_host):
        checked.append(url)
        ok = url == "https://ok.org"
        return check_links.Result(url, ok, 200 if ok else 404, "GET", url)

    monkeypatch.setattr(check_links, "check_url", fake_check)
    code = check_links.main([str(a), str(b), "--exclude", r"skip\.me"])
    out = capsys.readouterr().out
    assert code == 1
    assert sorted(checked) == ["https://bad.org/x", "https://ok.org"]
    assert "FAIL 404 GET  https://bad.org/x" in out
    assert "SKIP     http://localhost:8000/" in out
    assert "4 unique URLs in 2 files: 1 ok, 1 failed, 2 skipped." in out


def test_main_all_ok_missing_file_and_bad_args(monkeypatch, tmp_path):
    md = tmp_path / "ok.md"
    md.write_text("https://ok.org\n", "utf-8")
    monkeypatch.setattr(
        check_links,
        "check_url",
        lambda u, t, r, p: check_links.Result(u, True, 200, "HEAD", u),
    )
    assert check_links.main([str(md), "-q"]) == 0
    assert check_links.main([str(tmp_path / "missing.md")]) == 2
    with pytest.raises(SystemExit):
        check_links.main([str(md), "--retries", "-1"])


# --- Real urllib path against a loopback server (still offline) -------------------------------


class _Handler(BaseHTTPRequestHandler):
    def _reply(self, head_only):
        self.server.seen.append((self.command, self.path))  # before replying: no race
        if self.path == "/nohead" and head_only:
            self.send_response(405)
        elif self.path == "/redirect":
            self.send_response(301)
            self.send_header("Location", "/ok")
        elif self.path == "/loop":
            self.send_response(302)
            self.send_header("Location", "/loop")
        elif self.path == "/noloc":
            self.send_response(302)
        elif self.path in ("/ok", "/nohead"):
            self.send_response(200)
        else:
            self.send_response(404)
        self.send_header("Content-Length", "0")
        self.end_headers()

    def do_HEAD(self):  # noqa: N802 (http.server API)
        self._reply(head_only=True)

    def do_GET(self):  # noqa: N802
        self._reply(head_only=False)

    def log_message(self, *args):
        pass


@pytest.fixture
def server(monkeypatch):
    # Bypass any environment/registry proxy so loopback requests stay local.
    monkeypatch.setattr(urllib.request, "getproxies", dict)
    srv = ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    srv.seen = []
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    yield srv
    srv.shutdown()
    srv.server_close()


def test_against_loopback_server(server):
    base = f"http://127.0.0.1:{server.server_address[1]}"
    ok = check_links.check_url(f"{base}/ok", timeout=5)
    assert ok.ok and ok.method == "HEAD"

    nohead = check_links.check_url(f"{base}/nohead", timeout=5)
    assert nohead.ok and nohead.method == "GET" and nohead.status == 200

    server.seen.clear()
    redirect = check_links.check_url(f"{base}/redirect", timeout=5)
    assert redirect.ok and redirect.final_url == f"{base}/ok"
    assert server.seen == [("HEAD", "/redirect"), ("HEAD", "/ok")]  # HEAD kept across redirect

    for path in ("/loop", "/noloc"):
        res = check_links.check_url(f"{base}{path}", timeout=5, retries=0)
        assert not res.ok and res.status == 302 and res.error, path

    missing = check_links.check_url(f"{base}/missing", timeout=5, retries=0)
    assert not missing.ok and missing.status == 404

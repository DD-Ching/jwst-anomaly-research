"""Offline tests of scripts/net_ledger.py (local HTTP server, no external network)."""

import sys
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import net_ledger as nl  # noqa: E402


@pytest.fixture
def server(tmp_path, monkeypatch):
    monkeypatch.setenv("NO_PROXY", "127.0.0.1,localhost")
    monkeypatch.setenv("no_proxy", "127.0.0.1,localhost")
    (tmp_path / "www").mkdir()
    (tmp_path / "www" / "blob.bin").write_bytes(bytes(range(256)) * 400)  # 102400 B
    handler = partial(SimpleHTTPRequestHandler, directory=str(tmp_path / "www"))
    handler.log_message = lambda *a, **k: None
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    yield f"http://127.0.0.1:{httpd.server_port}/blob.bin"
    httpd.shutdown()


def test_reserve_refuses_over_cap_and_persists(tmp_path):
    p = tmp_path / "l.jsonl"
    lg = nl.Ledger(p, 1000)
    lg.reserve(600, "a")
    with pytest.raises(nl.BudgetExceeded):
        lg.reserve(500, "b")
    # a "restart": a new object on the same file still sees the open reservation
    assert nl.Ledger(p, 1000).used() == 600


def test_concurrent_reservations_are_all_recorded(tmp_path):
    lg = nl.Ledger(tmp_path / "l.jsonl", 10**6)

    def work():
        for _ in range(20):
            nl.Ledger(lg.path, lg.cap).reserve(100, "t")

    ts = [threading.Thread(target=work) for _ in range(8)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    assert lg.used() == 8 * 20 * 100


def test_done_replaces_reservation_by_actual(tmp_path):
    lg = nl.Ledger(tmp_path / "l.jsonl", 1000)
    rid = lg.reserve(600, "a")
    lg.done(rid, 250, "ok")
    lg.external(100, "wheel")
    assert lg.used() == 350


def test_fetch_counts_body_header_and_tls(tmp_path, server):
    lg = nl.Ledger(tmp_path / "l.jsonl", 10**6)
    body = nl.fetch(server, lg, est_bytes=102400)
    assert len(body) == 102400
    used = lg.used()
    assert 102400 + nl.TLS_ALLOWANCE < used < 102400 + nl.TLS_ALLOWANCE + 1000


def test_fetch_refused_before_request(tmp_path, server):
    lg = nl.Ledger(tmp_path / "l.jsonl", 50_000)
    with pytest.raises(nl.BudgetExceeded):
        nl.fetch(server, lg, est_bytes=102400)
    assert lg.used() == 0


def test_fetch_cut_mid_stream_on_underestimate(tmp_path, server):
    cap = 60_000
    lg = nl.Ledger(tmp_path / "l.jsonl", cap)
    with pytest.raises(nl.BudgetExceeded):
        nl.fetch(server, lg, est_bytes=1000)
    assert lg.used() <= cap

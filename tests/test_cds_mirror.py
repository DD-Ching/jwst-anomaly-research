"""VizieR server selection (D-006 addendum 2026-11-21)."""

import pytest

from jwst_anomaly import cds_mirror as cm

TSV = "#\n#Name: I/239/hip_main\n#Column\tHIP\nHIP\n\n------\n     1\n\n"
HTML = "<!doctype html>\n<html lang='en'><title>Page unavailable — CDS</title>1\n</html>"


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    monkeypatch.delenv(cm.ENV_VAR, raising=False)
    cm.forget()
    yield
    cm.forget()


def test_probe_ok_accepts_the_data_row_and_rejects_error_pages():
    assert cm.probe_ok(TSV)
    assert not cm.probe_ok(HTML)
    assert not cm.probe_ok("#\nHIP\n------\n\n")  # header only, no row


def test_cds_first_then_mirror_then_error(monkeypatch):
    calls = []

    def fake(answering):
        def probe(server, timeout):
            calls.append(server)
            return server in answering

        return probe

    monkeypatch.setattr(cm, "_probe", fake({cm.CDS_VIZIER, cm.CFA_VIZIER}))
    assert cm.vizier_server() == cm.CDS_VIZIER and calls == [cm.CDS_VIZIER]

    cm.forget()
    calls.clear()
    monkeypatch.setattr(cm, "_probe", fake({cm.CFA_VIZIER}))
    assert cm.vizier_server() == cm.CFA_VIZIER
    assert calls == [cm.CDS_VIZIER] * cm.PROBE_ATTEMPTS + [cm.CFA_VIZIER]

    cm.forget()
    calls.clear()
    monkeypatch.setattr(cm, "_probe", fake({cm.IDIA_VIZIER}))
    assert cm.vizier_server() == cm.IDIA_VIZIER
    assert calls == [cm.CDS_VIZIER, cm.CDS_VIZIER, cm.CFA_VIZIER, cm.CFA_VIZIER, cm.IDIA_VIZIER]

    cm.forget()
    monkeypatch.setattr(cm, "_probe", fake(set()))
    with pytest.raises(cm.NoVizierServer):
        cm.vizier_server()


def test_env_pin_skips_the_probe_even_after_caching(monkeypatch):
    monkeypatch.setattr(cm, "_probe", lambda s, t: s == cm.CDS_VIZIER)
    assert cm.vizier_server() == cm.CDS_VIZIER
    monkeypatch.setenv(cm.ENV_VAR, "vizier.example.org")
    assert cm.vizier_server() == "vizier.example.org"
    assert cm.asu_tsv_url() == "https://vizier.example.org/viz-bin/asu-tsv"
    assert cm.tag("VizieR X") == "VizieR X [VizieR host vizier.example.org]"
    monkeypatch.setenv(cm.ENV_VAR, " https://vizier.example.org/viz-bin/ ")
    assert cm.vizier_server() == "vizier.example.org"


def test_tag_never_probes(monkeypatch):
    def boom(server, timeout):
        raise AssertionError("tag() must not probe")

    monkeypatch.setattr(cm, "_probe", boom)
    assert cm.tag("X") == "X [VizieR host none queried in this process]"
    monkeypatch.setattr(cm, "_probe", lambda s, t: s == cm.CFA_VIZIER)
    cm.vizier_server()
    assert cm.tag("X") == f"X [VizieR host {cm.CFA_VIZIER}]"


def test_choice_is_cached_until_forget(monkeypatch):
    up = {cm.CDS_VIZIER}
    monkeypatch.setattr(cm, "_probe", lambda s, t: s in up)
    assert cm.vizier_server() == cm.CDS_VIZIER
    up = {cm.CFA_VIZIER}  # CDS goes down mid-run
    assert cm.vizier_server() == cm.CDS_VIZIER  # cached
    cm.forget()
    assert cm.vizier_server() == cm.CFA_VIZIER


def test_vizier_binds_the_selected_server(monkeypatch):
    monkeypatch.setenv(cm.ENV_VAR, cm.CFA_VIZIER)
    v = cm.vizier(columns=["HIP"], row_limit=1)
    assert v.VIZIER_SERVER == cm.CFA_VIZIER and v.ROW_LIMIT == 1


@pytest.mark.network
def test_live_some_server_answers():
    assert cm.vizier_server() in cm.VIZIER_SERVERS

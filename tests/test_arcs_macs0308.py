"""Tests for scripts/arcs_macs0308.py: model registration and wiring (no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_models_registered_as_planned():
    lc = _load("lens_consistency")
    am = _load("arcs_macs0308")
    for name in (am.PRIMARY, am.SECOND):
        spec = lc.MODELS[name]
        assert spec["kind"] == "maps"
        assert spec["z_lens"] == 0.356  # Acebron et al. 2018, Table 1
        assert set(spec["files"]) == {"alpha_x", "alpha_y", "kappa_map", "mag_z2"}
        assert all(len(sha) == 64 for _, sha in spec["files"].values())
    assert "glafic/v2" in lc.MODELS[am.PRIMARY]["files"]["alpha_x"][0]
    assert "zitrin-ltm-gauss" in lc.MODELS[am.SECOND]["files"]["alpha_x"][0]


def test_main_wires_conversion_and_decide(tmp_path, monkeypatch):
    am = _load("arcs_macs0308")
    calls = {}
    monkeypatch.setattr(am, "read_bpz", lambda p: "bpz")
    monkeypatch.setattr(am, "bpz_zout", lambda t: _Z())
    monkeypatch.setattr(am.Table, "read", lambda p: "cat")
    monkeypatch.setattr(
        am, "offset_arcsec", lambda z, c: {"n_pairs": 5, "dra_arcsec": 0.3, "ddec_arcsec": 0.0}
    )
    monkeypatch.setattr(am, "shift", lambda z, off: True)
    monkeypatch.setattr(am, "decide", lambda *a: calls.setdefault("decide", a))
    bpz = tmp_path / "bpz.txt"
    assert am.main(["--catalog", "c.ecsv", "--bpz", str(bpz), "--out", str(tmp_path / "o")]) == 0
    primary, second, out, _, zpath, plan = calls["decide"]
    assert (primary, second) == ("macs0308-relics-glafic", "macs0308-relics-zitrin")
    assert zpath == tmp_path / "macs0308_bpz_zout.fits"
    assert plan == "results/arcs/macs0308-relics/PLAN.md"
    assert '"shift_applied": true' in (out / "photoz_offset.json").read_text()


class _Z(list):
    def write(self, path, overwrite=False):
        Path(path).write_text("z")

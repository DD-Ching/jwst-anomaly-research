"""Offline end-to-end tests of the runner with contract-valid fake stages.

``FakeStages``/``install_fakes``/``write_config`` are also used by ``test_cli.py``.
"""

from __future__ import annotations

import json
import pathlib
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.table import Table

from jwst_anomaly import (
    acquire,
    catalog,
    crossmatch,
    cutouts,
    features,
    paths,
    pipeline,
    provenance,
    quality,
    query,
    rank,
    schema,
    viz,
)
from jwst_anomaly.candidates import CandidateStore

N_SOURCES = 12

TEST_CONFIG = {
    "name": "fake_sample",
    "archive": {
        "collection": "JWST",
        "calib_level": 3,
        "data_rights": "PUBLIC",
        "product_subgroups": ["CAT", "I2D"],
    },
    "samples": [
        {
            "id": "field_a",
            "role": "science",
            "description": "fake cluster field",
            "proposal_id": "2736",
            "instrument_name": "NIRCAM/IMAGE",
            "ref_band": "F200W",
            "obs_ids": [
                "jw02736-o001_t001_nircam_clear-f150w",
                "jw02736-o001_t001_nircam_clear-f200w",
            ],
        },
        {
            "id": "field_b",
            "role": "control",
            "description": "fake control field",
            "proposal_id": "1345",
            "instrument_name": "NIRCAM/IMAGE",
            "ref_band": "F200W",
            "obs_ids": ["jw01345-o001_t021_nircam_clear-f200w"],
        },
    ],
    "cloud": {
        "s3_bucket": "stpubdata",
        "l3_key_pattern": "jwst/public/jw{proposal:05d}/L3/t/o{observation:03d}/{obs_id}_{suffix}",
    },
    "stages": {
        "catalog": {"merge_radius_arcsec": 0.1},
        "quality": {"grid_arcsec": 1.0, "min_ranked": 2},
        "rank": {"methods": ["robust_z", "isolation_forest"], "random_state": 0},
        "cutouts": {"top_k": 3, "size_arcsec": 3.0, "bands": ["F200W"]},
        "crossmatch": {"top_k": 4, "radius_arcsec": 1.0, "services": ["simbad"]},
    },
    "outputs": {"candidates_top_k": 5},
}


def _meta(table: Table, prov: schema.Provenance, source: str) -> Table:
    table.meta.update(provenance=prov.value, source=source)
    return table


class FakeStages:
    """Deterministic stand-ins honouring the stage signatures and ``schema`` contracts."""

    def __init__(self) -> None:
        self.calls: dict[str, list] = defaultdict(list)

    def query_observations(self, **criteria):
        self.calls["query_observations"].append(criteria)
        obs_ids = list(criteria.get("obs_id") or [])
        n = len(obs_ids)
        t = Table(
            {
                "obsid": np.arange(n, dtype=int) + 1000,
                "obs_id": obs_ids,
                "proposal_id": [criteria.get("proposal_id", "")] * n,
                "instrument_name": [criteria.get("instrument_name", "")] * n,
                "filters": [pipeline.band_from_name(o) for o in obs_ids],
                "calib_level": [3] * n,
                "t_exptime": [1000.0] * n,
                "s_ra": [110.8] * n,
                "s_dec": [-73.45] * n,
                "dataRights": ["PUBLIC"] * n,
            }
        )
        return _meta(t, schema.Provenance.OBSERVED, "fake MAST query")

    def list_products(self, observations, subgroups=("CAT", "I2D"), calib_level=3):
        self.calls["list_products"].append((len(observations), tuple(subgroups), calib_level))
        rows = []
        for obs in observations:
            for sub, suffix in (("CAT", "cat.ecsv"), ("I2D", "i2d.fits")):
                if sub in subgroups:
                    name = f"{obs['obs_id']}_{suffix}"
                    rows.append(
                        (obs["obsid"], obs["obs_id"], name, sub, f"mast:JWST/product/{name}")
                    )
        t = Table(
            rows=rows,
            names=("obsID", "obs_id", "productFilename", "productSubGroupDescription", "dataURI"),
        )
        t["size"] = [
            3_000_000 if s == "CAT" else 1_800_000_000 for s in t["productSubGroupDescription"]
        ]
        t["calib_level"] = calib_level
        return _meta(t, schema.Provenance.OBSERVED, "fake MAST product list")

    def fetch_products(self, products, data_root=None, manifest_path=None):
        self.calls["fetch_products"].append(list(products["productSubGroupDescription"]))
        t = Table(
            {
                "dataURI": list(products["dataURI"]),
                "productFilename": list(products["productFilename"]),
                "local_path": [f"cache/{n}" for n in products["productFilename"]],
                "size": list(products["size"]),
                "sha256": ["0" * 64] * len(products),
                "retrieved_utc": ["2026-10-07T00:00:00Z"] * len(products),
                "pipeline_version": ["1.0"] * len(products),
            }
        )
        return _meta(t, schema.Provenance.OBSERVED, "fake manifest")

    def load_pipeline_catalog(self, path):
        self.calls["load_pipeline_catalog"].append(Path(path))
        band = pipeline.band_from_name(Path(path).name)
        i = np.arange(N_SOURCES)
        mag = 25.0 + 0.1 * i
        mag[3] = 19.0  # one obvious outlier
        t = Table(
            {
                "label": i + 1,
                "ra": 110.8 + 1e-4 * i,
                "dec": -73.45 + 1e-4 * i,
                "aper50_abmag": mag,
            }
        )
        t["ra"].unit = t["dec"].unit = "deg"
        t.meta["band"] = band
        return _meta(t, schema.Provenance.OBSERVED, f"fake catalog {Path(path).name}")

    def merge_bands(self, catalogs, ref_band, radius_arcsec=0.1):
        self.calls["merge_bands"].append((sorted(catalogs), ref_band, radius_arcsec))
        ref = catalogs[ref_band]
        t = Table(
            {
                "source_uid": [f"{ref_band.lower()}-{lab:04d}" for lab in ref["label"]],
                "ra": ref["ra"],
                "dec": ref["dec"],
                "ref_band": [ref_band] * len(ref),
                "n_bands": [len(catalogs)] * len(ref),
            }
        )
        for band, cat in catalogs.items():
            t[schema.band_column(band, "aper50_abmag")] = cat["aper50_abmag"]
        return _meta(t, schema.Provenance.DERIVED, "fake merge")

    def build_features(self, sources):
        ref = sources["ref_band"][0]
        t = Table({"source_uid": sources["source_uid"]})
        t["mag_ref"] = sources[schema.band_column(ref, "aper50_abmag")]
        t.meta["feature_spec"] = {"mag_ref": "aper50 AB magnitude in the reference band"}
        return _meta(t, schema.Provenance.DERIVED, "fake features")

    def score_anomalies(self, features, methods=("robust_z",), random_state=0):
        self.calls["score_anomalies"].append((tuple(methods), random_state))
        mag = np.asarray(features["mag_ref"], float)
        base = np.abs(mag - np.median(mag))
        t = Table({"source_uid": features["source_uid"], "score": base})
        for j, m in enumerate(methods):
            t[f"score_{m}"] = base * (1 + 0.1 * j)
        order = np.argsort(-base, kind="stable")
        rnk = np.empty(len(t), int)
        rnk[order] = np.arange(1, len(t) + 1)
        t["rank"] = rnk
        t["top_features"] = ["mag_ref"] * len(t)
        return _meta(t, schema.Provenance.MODEL_PREDICTION, "fake scores")

    def make_cutouts(self, image_uri, targets, size_arcsec=3.0, out_dir=None, **kw):
        self.calls["make_cutouts"].append((image_uri, len(targets), size_arcsec))
        self.calls["make_cutouts_kw"].append({"out_dir": str(out_dir), **kw})
        band = pipeline.band_from_name(image_uri)
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        paths_ = []
        for uid in targets["source_uid"]:
            p = out_dir / f"{uid}_{band}.fits"
            p.write_text("fake cutout")
            paths_.append(str(p))
        n = len(targets)
        t = Table(
            {
                "source_uid": targets["source_uid"],
                "band": [band] * n,
                "path": paths_,
                "frac_nan": [0.25] + [0.0] * (n - 1),
                "on_edge": [True] + [False] * (n - 1),
                "quality_flag": ["edge"] + ["ok"] * (n - 1),
            }
        )
        return _meta(t, schema.Provenance.DERIVED, f"fake cutouts of {image_uri}")

    def crossmatch(self, targets, radius_arcsec=1.0, services=("simbad",)):
        self.calls["crossmatch"].append((len(targets), radius_arcsec, tuple(services)))
        n = len(targets)
        t = Table(
            {
                "source_uid": targets["source_uid"],
                "n_matches": [1] + [0] * (n - 1),
                "is_known_object": [True] + [False] * (n - 1),
                "is_star": [True] + [False] * (n - 1),
                "best_match_id": ["Fake Star 1"] + [""] * (n - 1),
                "best_match_type": ["*"] + [""] * (n - 1),
                "best_match_service": ["simbad"] + [""] * (n - 1),
                "best_match_sep_arcsec": [0.12] + [np.nan] * (n - 1),
            }
        )
        t.meta["services"] = list(services)
        t.meta["radius_arcsec"] = radius_arcsec
        return _meta(t, schema.Provenance.OBSERVED, "fake crossmatch")

    def sample_weight_map(self, image_uri, step=None, **kw):
        # No WHT map: quality.assess_sources then runs only its PSF-sharpness test.
        self.calls["sample_weight_map"].append((image_uri, step))
        return None

    def query_matches(self, targets, radius_arcsec=1.0, services=("gaia", "simbad"), **kw):
        # Offline: no catalog knows these sources (classify then finds no stars).
        self.calls["query_matches"].append((len(targets), radius_arcsec, tuple(services)))
        return _match_table([], services, radius_arcsec)

    def contact_sheet(self, cutouts_table, out_png, *, ranks=None, title=None, **kw):
        self.calls["contact_sheet"].append((len(cutouts_table), dict(ranks or {}), title))
        out_png = Path(out_png)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        out_png.write_bytes(b"fake png")
        return out_png


_CLASSIFY = "classify.classify_sources"


def install_fakes(monkeypatch: pytest.MonkeyPatch, fakes: FakeStages | None = None) -> FakeStages:
    """Replace every stage function with ``fakes`` (module attributes, as the runner calls them)."""
    fakes = fakes or FakeStages()
    for module, name in (
        (query, "query_observations"),
        (query, "list_products"),
        (acquire, "fetch_products"),
        (catalog, "load_pipeline_catalog"),
        (catalog, "merge_bands"),
        (features, "build_features"),
        (rank, "score_anomalies"),
        (cutouts, "make_cutouts"),
        (cutouts, "sample_weight_map"),
        (crossmatch, "query_matches"),
        (crossmatch, "crossmatch"),
        (viz, "contact_sheet"),
    ):
        monkeypatch.setattr(module, name, getattr(fakes, name))
    return fakes


def write_config(tmp_path: Path, **overrides) -> Path:
    config = json.loads(json.dumps(TEST_CONFIG))
    for key, value in overrides.items():
        if value is None:
            config.pop(key, None)
        else:
            config[key] = value
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


@pytest.fixture
def env(tmp_path, monkeypatch):
    """Isolated data/output roots; returns the outputs directory."""
    monkeypatch.setenv(paths.DATA_ENV, str(tmp_path / "data"))
    monkeypatch.setenv(paths.OUTPUTS_ENV, str(tmp_path / "outputs"))
    return tmp_path / "outputs"


def _store(outputs: Path) -> CandidateStore:
    return CandidateStore(pipeline.default_db_path(outputs), create=False)


def test_run_end_to_end_with_fakes(tmp_path, env, monkeypatch):
    fakes = install_fakes(monkeypatch)
    config_path = write_config(tmp_path)
    run_id = pipeline.run(config_path)

    assert provenance.RUN_ID_PATTERN.match(run_id)
    run_dir = env / "runs" / run_id
    with _store(env) as store:
        run = store.get_run(run_id)
        assert run["status"] == "completed"
        assert run["config_name"] == "fake_sample"
        assert run["context"]["config"]["sha256"] == provenance.file_sha256(config_path)
        assert run["summary"]["n_candidates"] == 10
        cands = store.list_candidates(run_id=run_id)
        a = store.list_candidates(run_id=run_id, sample_id="field_a")

    assert len(cands) == 10 and len(a) == 5
    assert [c["rank"] for c in a] == [1, 2, 3, 4, 5]
    top = a[0]
    assert top["source_uid"] == "f200w-0004"  # the injected outlier
    assert top["status"] == "new"
    assert top["ra"] == pytest.approx(110.8003)
    assert set(top["method_scores"]) == {"robust_z", "isolation_forest"}
    assert top["flags"] == {"F200W": {"quality_flag": "edge", "on_edge": True, "frac_nan": 0.25}}
    assert top["cutouts"]["F200W"] == "field_a/cutouts/f200w-0004_F200W.fits"
    assert (run_dir / top["cutouts"]["F200W"]).is_file()
    assert top["xmatch"]["best_match_id"] == "Fake Star 1"
    assert a[1]["xmatch"]["n_matches"] == 0 and a[1]["xmatch"]["best_match_sep_arcsec"] is None
    assert a[3]["flags"] == {} and a[3]["xmatch"] is not None  # rank 4: crossmatch only
    assert a[4]["xmatch"] is None  # rank 5: beyond both top_k

    # Stage parameters come from the config; only catalogs are fetched.
    assert fakes.calls["fetch_products"] == [["CAT", "CAT"], ["CAT"]]
    assert fakes.calls["score_anomalies"][0] == (("robust_z", "isolation_forest"), 0)
    assert fakes.calls["merge_bands"][0] == (["F150W", "F200W"], "F200W", 0.1)
    assert fakes.calls["make_cutouts"][0] == (
        "s3://stpubdata/jwst/public/jw02736/L3/t/o001/jw02736-o001_t001_nircam_clear-f200w_i2d.fits",
        3,
        3.0,
    )
    assert fakes.calls["crossmatch"][0] == (4, 1.0, ("simbad",))
    criteria = fakes.calls["query_observations"][0]
    assert criteria == {
        "obs_collection": "JWST",
        "proposal_id": "2736",
        "instrument_name": "NIRCAM/IMAGE",
        "calib_level": 3,
        "obs_id": TEST_CONFIG["samples"][0]["obs_ids"],
    }
    assert fakes.calls["load_pipeline_catalog"][0] == (
        env.parent / "data" / "cache" / "jw02736-o001_t001_nircam_clear-f150w_cat.ecsv"
    )

    # Persisted run files.
    assert (run_dir / "config.yaml").read_bytes() == config_path.read_bytes()
    context = json.loads((run_dir / "run_context.json").read_text(encoding="utf-8"))
    assert context["run_id"] == run_id
    record = json.loads((run_dir / "run_record.json").read_text(encoding="utf-8"))
    assert record["status"] == "completed"
    assert {s["status"] for s in record["stages"] if s["stage"] != _CLASSIFY} == {"ok"}
    for name in ("observations", "products", "manifest", "sources", "features", "scores"):
        assert (run_dir / "field_a" / f"{name}.ecsv").is_file(), name
    scores = Table.read(run_dir / "field_a" / "scores.ecsv")
    assert scores.meta["provenance"] == "model_prediction"
    targets = Table.read(run_dir / "field_a" / "targets.ecsv")
    assert list(targets["source_uid"][:2]) == [a[0]["source_uid"], a[1]["source_uid"]]

    report = (run_dir / "report.md").read_text(encoding="utf-8")
    assert f"# Anomaly-ranking run `{run_id}`" in report
    assert "not evidence of new physics" in report
    assert "**Status:** completed" in report
    assert "#### Top 5 candidates" in report
    assert "[F200W](field_a/cutouts/f200w-0004_F200W.fits)" in report
    assert "F200W: edge, edge, NaN 25%" in report
    assert "- Contact sheet: [field_a/contact_sheet.png](field_a/contact_sheet.png)" in report
    assert (run_dir / "field_a" / "contact_sheet.png").is_file()
    sheet_calls = fakes.calls["contact_sheet"]
    assert sheet_calls and sheet_calls[0][1]["f200w-0004"] == 1  # ranks passed through
    assert 'simbad: Fake Star 1 (*, 0.12", star, known; 1 matches)' in report
    assert "no counterpart within 1.0 arcsec" in report
    assert "## Limitations" in report
    assert "`mag_ref`: aper50 AB magnitude in the reference band" in report


def test_reference_config_runs_with_fakes(env, monkeypatch):
    """The real reference config parses and drives the full chain (with fake stages)."""
    install_fakes(monkeypatch)
    config_path = paths.repo_root() / "configs" / "reference_sample.yaml"
    config = pipeline.load_config(config_path)
    run_id = pipeline.run(config_path)
    with _store(env) as store:
        cands = store.list_candidates(run_id=run_id)
    assert {c["sample_id"] for c in cands} == {s["id"] for s in config["samples"]}
    report = (env / "runs" / run_id / "report.md").read_text(encoding="utf-8")
    # MIRI sample has no F200W image: cutouts fall back to its ref band.
    assert "cutout bands ['F200W'] not in this sample; used ref_band F770W" in report


def not_implemented(*args, **kwargs):
    """Stand-in for a skeleton stub, independent of whether the real stage has landed."""
    raise NotImplementedError("bootstrap unit 1: archive access")


def test_required_stage_not_implemented_fails_cleanly(tmp_path, env, monkeypatch):
    """A required stage that is still a stub fails the run with a clear, recorded message."""
    install_fakes(monkeypatch)
    monkeypatch.setattr(query, "query_observations", not_implemented)
    with pytest.raises(pipeline.PipelineError) as info:
        pipeline.run(write_config(tmp_path))
    err = info.value
    assert "query.query_observations" in str(err)
    assert "stage not implemented yet" in str(err)
    assert "field_a" in str(err)
    with _store(env) as store:
        run = store.get_run(err.run_id)
        assert run["status"] == "failed"
        assert "stage not implemented" in run["error"]
        assert store.list_candidates(run_id=err.run_id) == []
    report = err.report_path.read_text(encoding="utf-8")
    assert "**Status:** failed" in report
    assert "stage not implemented yet" in report


def test_optional_stage_failure_is_recorded_and_run_continues(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)

    def broken(*args, **kwargs):
        raise RuntimeError("S3 unreachable")

    monkeypatch.setattr(cutouts, "make_cutouts", broken)
    monkeypatch.setattr(crossmatch, "crossmatch", broken)
    run_id = pipeline.run(write_config(tmp_path))
    with _store(env) as store:
        run = store.get_run(run_id)
        cands = store.list_candidates(run_id=run_id)
    assert run["status"] == "completed"
    failed = [s for s in run["summary"]["stages"] if s["status"] == "failed"]
    assert {s["stage"] for s in failed} == {"cutouts.make_cutouts[F200W]", "crossmatch.crossmatch"}
    assert all(not s["required"] for s in failed)
    assert all(c["flags"] == {} and c["xmatch"] is None for c in cands)
    report = (env / "runs" / run_id / "report.md").read_text(encoding="utf-8")
    assert "RuntimeError: S3 unreachable" in report
    assert "| not run | not run | none |" in report


def test_optional_stages_skipped_when_not_configured(tmp_path, env, monkeypatch):
    fakes = install_fakes(monkeypatch)
    stages = {k: v for k, v in TEST_CONFIG["stages"].items() if k not in ("cutouts", "crossmatch")}
    run_id = pipeline.run(write_config(tmp_path, stages=stages))
    assert "make_cutouts" not in fakes.calls and "crossmatch" not in fakes.calls
    with _store(env) as store:
        run = store.get_run(run_id)
        cands = store.list_candidates(run_id=run_id, sample_id="field_a")
    skipped = {s["stage"] for s in run["summary"]["stages"] if s["status"] == "skipped"}
    assert skipped == {"cutouts.make_cutouts", "crossmatch.crossmatch", _CLASSIFY}
    assert len(cands) == 5  # outputs.candidates_top_k


def test_disabled_stage_is_skipped_with_reason(tmp_path, env, monkeypatch):
    fakes = install_fakes(monkeypatch)
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages["crossmatch"]["enabled"] = False
    run_id = pipeline.run(write_config(tmp_path, stages=stages), samples=["field_b"])
    assert "crossmatch" not in fakes.calls
    with _store(env) as store:
        stages_run = store.get_run(run_id)["summary"]["stages"]
    assert {(s["stage"], s["message"]) for s in stages_run if s["status"] == "skipped"} == {
        ("crossmatch.crossmatch", "disabled in config"),
        (_CLASSIFY, "not configured"),
    }


def test_interrupt_is_recorded(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)

    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(features, "build_features", interrupted)
    with pytest.raises(KeyboardInterrupt):
        pipeline.run(write_config(tmp_path))
    with _store(env) as store:
        (run,) = store.list_runs()
    assert run["status"] == "failed" and run["error"] == "interrupted"
    assert (env / "runs" / run["run_id"] / "report.md").is_file()


def test_duplicate_scored_uid_aborts(tmp_path, env, monkeypatch):
    fakes = install_fakes(monkeypatch)

    def dup_scores(features, methods=("robust_z",), random_state=0):
        t = fakes.score_anomalies(features, methods, random_state)
        t["source_uid"][1] = t["source_uid"][0]
        return t

    monkeypatch.setattr(rank, "score_anomalies", dup_scores)
    with pytest.raises(pipeline.PipelineError, match="duplicate source_uid 'f200w-0001'"):
        pipeline.run(write_config(tmp_path))
    with _store(env) as store:
        assert store.list_runs()[0]["status"] == "failed"


def test_list_valued_top_features_and_cutout_merge_failure(tmp_path, env, monkeypatch):
    """Odd but valid stage outputs must not break storage; optional merge errors stay optional."""
    fakes = install_fakes(monkeypatch)

    def list_top(features, methods=("robust_z",), random_state=0):
        t = fakes.score_anomalies(features, methods, random_state)
        t["top_features"] = [["mag_ref", "color"]] * len(t)
        return t

    def mixed_paths(image_uri, targets, size_arcsec=3.0, out_dir=None):
        t = fakes.make_cutouts(image_uri, targets, size_arcsec, out_dir)
        if "f150w" in image_uri:  # valid alone, but cannot be stacked with the F200W table
            t["path"] = np.arange(len(t), dtype=float)
        return t

    monkeypatch.setattr(rank, "score_anomalies", list_top)
    monkeypatch.setattr(cutouts, "make_cutouts", mixed_paths)
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages["cutouts"]["bands"] = ["F150W", "F200W"]
    run_id = pipeline.run(write_config(tmp_path, stages=stages), samples=["field_a"])
    with _store(env) as store:
        run = store.get_run(run_id)
        top = store.list_candidates(run_id=run_id)[0]
    assert run["status"] == "completed"
    assert top["top_features"] == "mag_ref, color"
    assert top["flags"] == {} and top["xmatch"] is not None
    (merge,) = [s for s in run["summary"]["stages"] if s["stage"] == "cutouts (combine bands)"]
    assert merge["status"] == "failed" and not merge["required"]


def test_relative_outputs_dir_gives_run_relative_cutout_links(tmp_path, monkeypatch):
    install_fakes(monkeypatch)
    monkeypatch.setenv(paths.DATA_ENV, str(tmp_path / "data"))
    monkeypatch.chdir(tmp_path)
    run_id = pipeline.run(write_config(tmp_path), outputs_dir="out", samples=["field_a"])
    with CandidateStore(tmp_path / "out" / "candidates.sqlite", create=False) as store:
        top = store.list_candidates(run_id=run_id)[0]
    rel = top["cutouts"]["F200W"]
    assert rel == "field_a/cutouts/f200w-0004_F200W.fits"
    assert (tmp_path / "out" / "runs" / run_id / rel).is_file()


def test_contract_violation_aborts(tmp_path, env, monkeypatch):
    fakes = install_fakes(monkeypatch)

    def bad_scores(features, methods=("robust_z",), random_state=0):
        t = fakes.score_anomalies(features, methods, random_state)
        t.remove_column("rank")
        return t

    monkeypatch.setattr(rank, "score_anomalies", bad_scores)
    with pytest.raises(pipeline.PipelineError, match=r"rank\.score_anomalies.*missing required"):
        pipeline.run(write_config(tmp_path))


def test_missing_provenance_aborts(tmp_path, env, monkeypatch):
    fakes = install_fakes(monkeypatch)

    def no_prov(sources):
        t = fakes.build_features(sources)
        del t.meta["provenance"]
        return t

    monkeypatch.setattr(features, "build_features", no_prov)
    with pytest.raises(pipeline.PipelineError, match="provenance"):
        pipeline.run(write_config(tmp_path))


def test_ref_band_missing_aborts(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)
    samples = json.loads(json.dumps(TEST_CONFIG["samples"]))
    samples[0]["ref_band"] = "F444W"
    with pytest.raises(pipeline.PipelineError, match="ref_band F444W not among"):
        pipeline.run(write_config(tmp_path, samples=samples))


@pytest.mark.parametrize(
    "case, note",
    [
        ("proprietary", "dropped 1 rows with dataRights != PUBLIC"),
        (
            "missing",
            "configured obs_ids not returned by MAST: ['jw02736-o001_t001_nircam_clear-f150w']",
        ),
    ],
)
def test_observation_filtering(tmp_path, env, monkeypatch, case, note):
    """Unconfigured obs_ids are dropped; non-public rows and missing obs_ids are reported."""
    fakes = install_fakes(monkeypatch)
    real_query = fakes.query_observations

    def noisy_query(**criteria):
        t = real_query(**criteria)  # rows: f150w, f200w
        t.add_row(t[1])
        t["obs_id"] = [*t["obs_id"][:2], "jw09999-o001_t001_nircam_clear-f200w"]  # unconfigured
        if case == "proprietary":
            t["dataRights"] = ["EXCLUSIVE_ACCESS", "PUBLIC", "PUBLIC"]
        else:
            t.remove_row(0)
        return t

    monkeypatch.setattr(query, "query_observations", noisy_query)
    run_id = pipeline.run(write_config(tmp_path), samples=["field_a"])
    with _store(env) as store:
        run = store.get_run(run_id)
    notes = [s["message"] for s in run["summary"]["stages"] if s["stage"] == "filter observations"]
    assert notes == [note]
    sample = run["summary"]["samples"][0]
    assert sample["bands"] == ["F200W"] and sample["n_observations"] == 1


def test_sample_selection(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)
    config_path = write_config(tmp_path)
    run_id = pipeline.run(config_path, samples=["field_b"])
    with _store(env) as store:
        assert {c["sample_id"] for c in store.list_candidates(run_id=run_id)} == {"field_b"}
    with pytest.raises(pipeline.ConfigError, match="unknown sample"):
        pipeline.run(config_path, samples=["nope"])


def test_parquet_tables_and_explicit_paths(tmp_path, monkeypatch):
    install_fakes(monkeypatch)
    monkeypatch.setenv(paths.DATA_ENV, str(tmp_path / "data"))
    out, db = tmp_path / "elsewhere", tmp_path / "db" / "store.sqlite"
    config_path = write_config(tmp_path, outputs={"table_format": "parquet"})
    run_id = pipeline.run(config_path, outputs_dir=out, db_path=db)
    scores = Table.read(out / "runs" / run_id / "field_a" / "scores.parquet")
    assert scores.meta["provenance"] == "model_prediction"
    with CandidateStore(db, create=False) as store:
        # default candidates_top_k = max(cutouts.top_k, crossmatch.top_k) = 4
        assert len(store.list_candidates(run_id=run_id, sample_id="field_a")) == 4


@pytest.mark.parametrize(
    "mutate, match",
    [
        (lambda c: c.pop("name"), "'name'"),
        (lambda c: c.update(samples=[]), "non-empty list"),
        (lambda c: c["samples"].append(dict(c["samples"][0])), "duplicate sample id"),
        (lambda c: c["samples"][0].pop("ref_band"), "ref_band"),
        (lambda c: c["samples"][1].update(obs_ids=[], proposal_id=None), "obs_ids"),
        (lambda c: c.update(stages=[1]), "'stages' must be a mapping"),
        (lambda c: c.update(outputs={"table_format": "csv"}), "table_format"),
        (lambda c: c["samples"][0].update(id="a/b"), "may only use"),
        (lambda c: c["stages"]["cutouts"].update(top_k=0), "positive integer"),
        (lambda c: c["stages"]["cutouts"].update(enabled="false"), "true or false"),
        (lambda c: c["stages"].update(cutouts=True), "stages.cutouts must be a mapping"),
        (lambda c: c.update(outputs={"candidates_top_k": -1}), "positive integer"),
    ],
)
def test_load_config_rejects(tmp_path, mutate, match):
    config = json.loads(json.dumps(TEST_CONFIG))
    mutate(config)
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(pipeline.ConfigError, match=match):
        pipeline.load_config(path)


def test_load_config_missing_file(tmp_path):
    with pytest.raises(pipeline.ConfigError, match="not found"):
        pipeline.load_config(tmp_path / "nope.yaml")


def test_l3_image_uri():
    cloud = TEST_CONFIG["cloud"]
    assert pipeline.l3_image_uri(cloud, "jw01345-o001_t021_nircam_clear-f200w") == (
        "s3://stpubdata/jwst/public/jw01345/L3/t/o001/jw01345-o001_t021_nircam_clear-f200w_i2d.fits"
    )
    with pytest.raises(pipeline.ConfigError):
        pipeline.l3_image_uri(None, "jw01345-o001_t021_nircam_clear-f200w")
    with pytest.raises(ValueError, match="obs_id"):
        pipeline.l3_image_uri(cloud, "not-an-obs-id")


@pytest.mark.parametrize(
    "name, band",
    [
        ("jw02736-o001_t001_nircam_clear-f200w_cat.ecsv", "F200W"),
        ("jw02736-o002_t001_miri_f1800w", "F1800W"),
        ("jw01345-o001_t021_nircam_clear-f150w2_cat.ecsv", "F150W2"),
        ("jw01234-o001_t001_nircam_f444w-f470n", "F470N"),
        ("no band here", None),
    ],
)
def test_band_from_name(name, band):
    assert pipeline.band_from_name(name) == band


def test_catalog_meta_band_takes_precedence():
    t = Table({"label": [1]})
    t.meta["band"] = "f356w"
    assert pipeline._band_of(t, "jw02736-o001_t001_nircam_clear-f200w") == "F356W"
    t.meta["band"] = "CLEAR"  # not a filter name: fall back to the obs_id
    assert pipeline._band_of(t, "jw02736-o001_t001_nircam_clear-f200w") == "F200W"


def test_xmatch_cell_does_not_invent_a_non_detection():
    assert pipeline._xmatch_cell(None) == "not run"
    assert pipeline._xmatch_cell({"n_matches": None}) == "result missing"
    assert pipeline._xmatch_cell({"n_matches": 0, "radius_arcsec": 1.0}) == (
        "no counterpart within 1.0 arcsec"
    )


def test_quality_gate_excludes_flagged_sources_from_ranking(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)

    def gate(sources, weight_map, *, ref_band=None, **kw):
        uids = [str(u) for u in sources["source_uid"]]
        t = Table(
            {
                "source_uid": uids,
                "rel_weight": [1.0] * len(uids),
                "edge_dist_arcsec": [5.0] * len(uids),
                "sharper_than_psf": [False] * len(uids),
                "quality_ok": [u != uids[0] for u in uids],
                "quality_reason": ["low_weight" if u == uids[0] else "" for u in uids],
            }
        )
        t.meta.update(provenance="derived", source="fake gate")
        return t

    monkeypatch.setattr(quality, "assess_sources", gate)
    run_id = pipeline.run(write_config(tmp_path))
    run_dir = env / "runs" / run_id
    sources = Table.read(run_dir / "field_a" / "sources.ecsv")
    scores = Table.read(run_dir / "field_a" / "scores.ecsv")
    assert len(scores) == len(sources) - 1
    assert str(sources["source_uid"][0]) not in {str(u) for u in scores["source_uid"]}
    assert (run_dir / "field_a" / "quality.ecsv").is_file()
    report = (run_dir / "report.md").read_text(encoding="utf-8")
    n = len(sources)
    assert f"Quality gate (derived, D-011): {n - 1} of {n} passed; flags: low_weight 1" in report


def test_quality_gate_failure_ranks_everything(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)

    def broken(image_uri, step=None, **kw):
        raise OSError("S3 unreachable")

    monkeypatch.setattr(cutouts, "sample_weight_map", broken)
    run_id = pipeline.run(write_config(tmp_path))
    run_dir = env / "runs" / run_id
    sources = Table.read(run_dir / "field_a" / "sources.ecsv")
    scores = Table.read(run_dir / "field_a" / "scores.ecsv")
    assert len(scores) == len(sources)
    report = (run_dir / "report.md").read_text(encoding="utf-8")
    assert "quality gate failed: all sources ranked (ungated)" in report
    assert "Quality gate (derived, D-011): not applied (all sources ranked)" in report


def test_quality_gate_with_too_few_survivors_ranks_everything(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)

    def strict_gate(sources, weight_map, *, ref_band=None, **kw):
        uids = [str(u) for u in sources["source_uid"]]
        t = Table(
            {
                "source_uid": uids,
                "rel_weight": [0.1] * len(uids),
                "edge_dist_arcsec": [5.0] * len(uids),
                "sharper_than_psf": [False] * len(uids),
                "quality_ok": [i == 0 for i in range(len(uids))],
                "quality_reason": ["" if i == 0 else "low_weight" for i in range(len(uids))],
            }
        )
        t.meta.update(provenance="derived", source="fake gate")
        return t

    monkeypatch.setattr(quality, "assess_sources", strict_gate)
    run_id = pipeline.run(write_config(tmp_path))
    run_dir = env / "runs" / run_id
    sources = Table.read(run_dir / "field_a" / "sources.ecsv")
    scores = Table.read(run_dir / "field_a" / "scores.ecsv")
    assert len(scores) == len(sources)  # a gate that empties the sample is not applied
    report = (run_dir / "report.md").read_text(encoding="utf-8")
    assert "quality gate unusable (only 1 of" in report


def test_quality_gate_skipped_when_not_configured(tmp_path, env, monkeypatch):
    fakes = install_fakes(monkeypatch)
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages.pop("quality")
    run_id = pipeline.run(write_config(tmp_path, stages=stages))
    assert fakes.calls["sample_weight_map"] == []
    report = (env / "runs" / run_id / "report.md").read_text(encoding="utf-8")
    assert "Quality gate (derived, D-011): not applied (all sources ranked)" in report


def _match_table(rows, services=("gaia", "simbad"), radius_arcsec=0.5, failed=()):
    t = Table(
        rows=rows,
        names=("source_uid", "service", "match_id", "match_type", "sep_arcsec", "is_star"),
        dtype=(str, str, str, str, float, bool),
    )
    t.meta.update(
        provenance="observed",
        source="fake",
        services_requested=list(services),
        services_ok=[s for s in services if s not in failed],
        services_failed=list(failed),
        radius_arcsec=radius_arcsec,
    )
    return t


def _fake_matches(star_uids):
    def query_matches(targets, radius_arcsec=1.0, services=("gaia", "simbad"), **kw):
        rows = [
            (str(u), "gaia", f"Gaia DR3 {i}", "astrometric_star", 0.05, True)
            for i, u in enumerate(targets["source_uid"])
            if str(u) in star_uids
        ]
        return _match_table(rows, services, radius_arcsec)

    return query_matches


def test_stars_are_ranked_as_their_own_stratum(tmp_path, env, monkeypatch):
    fakes = install_fakes(monkeypatch)
    probe = pipeline.run(write_config(tmp_path))  # learn the fake source uids
    sources = Table.read(env / "runs" / probe / "field_a" / "sources.ecsv")
    star_uids = {str(u) for u in sources["source_uid"][:6]}
    monkeypatch.setattr(crossmatch, "query_matches", _fake_matches(star_uids))
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages["classify"] = {"star_top_k": 3, "min_stars": 2}
    run_id = pipeline.run(write_config(tmp_path, stages=stages), samples=["field_a"])
    run_dir = env / "runs" / run_id
    galaxy_scores = Table.read(run_dir / "field_a" / "scores.ecsv")
    star_scores = Table.read(run_dir / "field_a-stars" / "scores.ecsv")
    assert {str(u) for u in star_scores["source_uid"]} <= star_uids
    assert not ({str(u) for u in galaxy_scores["source_uid"]} & star_uids)
    assert (run_dir / "field_a" / "populations.ecsv").is_file()
    with _store(env) as store:
        star_cands = store.list_candidates(run_id=run_id, sample_id="field_a-stars")
    assert 0 < len(star_cands) <= 3
    report = (run_dir / "report.md").read_text(encoding="utf-8")
    assert "### field_a-stars" in report
    assert "ranked among themselves (D-012)" in report
    star_cut_dirs = [
        c["out_dir"] for c in fakes.calls["make_cutouts_kw"] if "field_a-stars" in c["out_dir"]
    ]
    assert star_cut_dirs  # the star stratum ran its own cutouts
    assert "6 stars ranked separately as field_a-stars (D-012)" in report
    assert "applied in field_a; this stratum holds its" in report


def test_classify_failure_keeps_one_ranking(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)

    def failed(targets, radius_arcsec=1.0, services=("gaia",), **kw):
        return _match_table([], services, radius_arcsec, failed=("gaia",))

    monkeypatch.setattr(crossmatch, "query_matches", failed)
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages["classify"] = {}
    run_id = pipeline.run(write_config(tmp_path, stages=stages), samples=["field_a"])
    run_dir = env / "runs" / run_id
    assert not (run_dir / "field_a-stars").exists()
    report = (run_dir / "report.md").read_text(encoding="utf-8")
    assert "star/galaxy separation failed: one ranking" in report


def test_too_few_stars_are_not_ranked(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)
    probe = pipeline.run(write_config(tmp_path))
    sources = Table.read(env / "runs" / probe / "field_a" / "sources.ecsv")
    monkeypatch.setattr(crossmatch, "query_matches", _fake_matches({str(sources["source_uid"][0])}))
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages["classify"] = {"min_stars": 5}
    run_id = pipeline.run(write_config(tmp_path, stages=stages), samples=["field_a"])
    report = (env / "runs" / run_id / "report.md").read_text(encoding="utf-8")
    assert "star/galaxy strata not used (1 stars" in report
    scores = Table.read(env / "runs" / run_id / "field_a" / "scores.ecsv")
    assert str(sources["source_uid"][0]) in {str(u) for u in scores["source_uid"]}  # kept
    assert not (env / "runs" / run_id / "field_a-stars").exists()


def test_mostly_stars_falls_back_to_one_ranking(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)
    probe = pipeline.run(write_config(tmp_path))
    sources = Table.read(env / "runs" / probe / "field_a" / "sources.ecsv")
    all_but_one = {str(u) for u in sources["source_uid"][1:]}
    monkeypatch.setattr(crossmatch, "query_matches", _fake_matches(all_but_one))
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages["classify"] = {"min_stars": 2}
    run_id = pipeline.run(write_config(tmp_path, stages=stages), samples=["field_a"])
    with _store(env) as store:
        assert store.get_run(run_id)["status"] == "completed"
    report = (env / "runs" / run_id / "report.md").read_text(encoding="utf-8")
    assert "1 others; min_stars 2): one ranking" in report


@pytest.mark.parametrize(
    "bad, expect",
    [
        ({"star_top_k": -1}, "star_top_k >= 1"),
        ({"star_top_k": "ten"}, "invalid literal"),
        ({"min_stars": 1}, "min_stars >= 2"),
    ],
)
def test_invalid_classify_params_fail_only_the_stage(tmp_path, env, monkeypatch, bad, expect):
    install_fakes(monkeypatch)
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages["classify"] = bad
    run_id = pipeline.run(write_config(tmp_path, stages=stages), samples=["field_a"])
    with _store(env) as store:
        run = store.get_run(run_id)
    assert run["status"] == "completed"
    failed = [s for s in run["summary"]["stages"] if s["stage"] == _CLASSIFY]
    assert failed and failed[0]["status"] == "failed" and expect in failed[0]["message"]


def test_stratum_id_collision_fails_classify(tmp_path, env, monkeypatch):
    install_fakes(monkeypatch)
    config = json.loads(json.dumps(TEST_CONFIG))
    clash = json.loads(json.dumps(config["samples"][0]))
    clash["id"] = f"{config['samples'][0]['id']}-stars"
    config["samples"].append(clash)
    config["stages"]["classify"] = {}
    path = tmp_path / "clash.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    run_id = pipeline.run(path, samples=[config["samples"][0]["id"]])
    report = (env / "runs" / run_id / "report.md").read_text(encoding="utf-8")
    assert "star/galaxy separation failed: one ranking" in report


def test_cutouts_reuse_the_gate_reference_weight_for_every_stratum(tmp_path, env, monkeypatch):
    from types import SimpleNamespace

    fakes = install_fakes(monkeypatch)
    probe = pipeline.run(write_config(tmp_path))
    sources = Table.read(env / "runs" / probe / "field_a" / "sources.ecsv")
    monkeypatch.setattr(
        cutouts, "sample_weight_map", lambda uri, **kw: SimpleNamespace(reference_weight=2.5)
    )

    def gate(sources, weight_map, *, ref_band=None, **kw):
        uids = [str(u) for u in sources["source_uid"]]
        t = Table(
            {
                "source_uid": uids,
                "rel_weight": [1.0] * len(uids),
                "edge_dist_arcsec": [5.0] * len(uids),
                "sharper_than_psf": [False] * len(uids),
                "quality_ok": [True] * len(uids),
                "quality_reason": [""] * len(uids),
            }
        )
        t.meta.update(provenance="derived", source="fake gate")
        return t

    monkeypatch.setattr(quality, "assess_sources", gate)
    monkeypatch.setattr(
        crossmatch, "query_matches", _fake_matches({str(u) for u in sources["source_uid"][:6]})
    )
    stages = json.loads(json.dumps(TEST_CONFIG["stages"]))
    stages["classify"] = {"min_stars": 2}
    fakes.calls["make_cutouts_kw"].clear()
    pipeline.run(write_config(tmp_path, stages=stages), samples=["field_a"])
    kws = fakes.calls["make_cutouts_kw"]
    assert {"field_a", "field_a-stars"} <= {pathlib.Path(k["out_dir"]).parent.name for k in kws}
    assert all(k.get("weight_ref") == 2.5 for k in kws)


def test_rejected_gate_does_not_leak_its_reference_weight(tmp_path, env, monkeypatch):
    from types import SimpleNamespace

    fakes = install_fakes(monkeypatch)
    monkeypatch.setattr(
        cutouts, "sample_weight_map", lambda uri, **kw: SimpleNamespace(reference_weight=2.5)
    )

    def reject_all_but_one(sources, weight_map, *, ref_band=None, **kw):
        uids = [str(u) for u in sources["source_uid"]]
        t = Table(
            {
                "source_uid": uids,
                "rel_weight": [0.1] * len(uids),
                "edge_dist_arcsec": [5.0] * len(uids),
                "sharper_than_psf": [False] * len(uids),
                "quality_ok": [i == 0 for i in range(len(uids))],
                "quality_reason": ["" if i == 0 else "low_weight" for i in range(len(uids))],
            }
        )
        t.meta.update(provenance="derived", source="fake gate")
        return t

    monkeypatch.setattr(quality, "assess_sources", reject_all_but_one)
    pipeline.run(write_config(tmp_path), samples=["field_a"])
    kws = fakes.calls["make_cutouts_kw"]
    assert kws and all("weight_ref" not in k for k in kws)

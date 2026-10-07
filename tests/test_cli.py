import json

import pytest
from test_pipeline import install_fakes, not_implemented, write_config  # shared fake stages

from jwst_anomaly import cli, paths, query


@pytest.fixture
def outputs(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.DATA_ENV, str(tmp_path / "data"))
    monkeypatch.setenv(paths.OUTPUTS_ENV, str(tmp_path / "outputs"))
    monkeypatch.setenv(cli.AUTHOR_ENV, "cli-test")
    return tmp_path / "outputs"


@pytest.fixture
def run_id(tmp_path, outputs, monkeypatch, capsys):
    install_fakes(monkeypatch)
    assert cli.main(["run", "--config", str(write_config(tmp_path))]) == 0
    out = capsys.readouterr().out
    assert "completed: 10 candidates" in out
    return out.split()[1]


def test_run_and_list(run_id, outputs, capsys):
    assert (outputs / "runs" / run_id / "report.md").is_file()
    assert cli.main(["candidates", "list", "--run", run_id, "--sample", "field_a"]) == 0
    out = capsys.readouterr().out
    lines = out.strip().splitlines()
    assert lines[0].split() == ["RUN_ID", "SAMPLE_ID", "RANK", "SOURCE_UID", "SCORE", "STATUS"]
    assert len(lines) == 6 and "f200w-0004" in lines[1] and lines[1].endswith("new")

    assert cli.main(["candidates", "list", "--json", "--status", "new", "--limit", "3"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert len(rows) == 3 and rows[0]["run_id"] == run_id

    assert cli.main(["runs", "list"]) == 0
    assert run_id in capsys.readouterr().out


def test_vetting_workflow(run_id, capsys):
    uid = ["f200w-0004", "--sample", "field_a"]
    # Ambiguous: the same uid exists in both samples of the fake run.
    assert cli.main(["candidates", "set-status", "f200w-0004", "triaged", "--note", "x"]) == 1
    assert "ambiguous" in capsys.readouterr().err

    assert cli.main(["candidates", "set-status", *uid, "unexplained", "--note", "x"]) == 1
    assert "invalid transition new -> unexplained" in capsys.readouterr().err

    assert cli.main(["candidates", "set-status", *uid, "triaged", "--note", "first look"]) == 0
    assert "status new -> triaged" in capsys.readouterr().out

    assert cli.main(["candidates", "set-status", *uid, "artifact", "--note", "spike"]) == 1
    assert "vetting note" in capsys.readouterr().err

    vet = ["--test", "diffraction_spike", "--outcome", "fail", "--provenance", "observed"]
    assert cli.main(["candidates", "add-vetting", *uid, *vet, "--evidence", "cutout"]) == 0
    assert "diffraction_spike -> fail" in capsys.readouterr().out
    assert cli.main(["candidates", "set-status", *uid, "artifact", "--note", "spike"]) == 0
    capsys.readouterr()

    assert cli.main(["candidates", "show", *uid]) == 0
    cand = json.loads(capsys.readouterr().out)
    assert cand["status"] == "artifact"
    assert cand["vetting_notes"][0]["author"] == "cli-test"
    assert [h["new_status"] for h in cand["status_history"]] == ["triaged", "artifact"]

    assert cli.main(["candidates", "list", "--status", "artifact"]) == 0
    assert "f200w-0004" in capsys.readouterr().out


def test_run_with_stub_stage_fails_cleanly(tmp_path, outputs, monkeypatch, capsys):
    install_fakes(monkeypatch)
    monkeypatch.setattr(query, "query_observations", not_implemented)
    assert cli.main(["run", "--config", str(write_config(tmp_path))]) == 1
    err = capsys.readouterr().err
    assert "required stage query.query_observations failed" in err
    assert "stage not implemented yet" in err
    assert "Traceback" not in err
    assert "report:" in err


def test_errors_and_usage(tmp_path, outputs, capsys):
    assert cli.main(["candidates", "list"]) == 1
    assert "no candidate store" in capsys.readouterr().err
    assert cli.main(["run", "--config", str(tmp_path / "missing.yaml")]) == 1
    assert "config file not found" in capsys.readouterr().err
    bogus = tmp_path / "not-a-db.sqlite"
    bogus.write_text("name: not a database\n" * 100)
    assert cli.main(["candidates", "list", "--db", str(bogus)]) == 1
    assert "error:" in capsys.readouterr().err
    assert cli.main([]) == 2
    assert cli.main(["candidates", "set-status", "x", "bogus", "--note", "n"]) == 2
    assert cli.main(["--version"]) == 0
    assert "jwst-anomaly" in capsys.readouterr().out


def test_mutating_commands_need_run_when_uid_in_several_runs(run_id, tmp_path, capsys):
    assert cli.main(["run", "--config", str(write_config(tmp_path))]) == 0
    second = capsys.readouterr().out.split()[1]
    args = ["candidates", "set-status", "f200w-0004", "triaged", "--sample", "field_a"]
    assert cli.main([*args, "--note", "x"]) == 1
    assert "occurs in 2 runs" in capsys.readouterr().err
    assert cli.main([*args, "--note", "x", "--run", second]) == 0
    assert f"run {second}" in capsys.readouterr().out

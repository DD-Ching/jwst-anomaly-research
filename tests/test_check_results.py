"""Offline tests of scripts/check_results.py (docs/efficiency-rules.md)."""

import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import check_results as cr  # noqa: E402


def test_size_rule():
    files = [
        "results/a/rows.json",
        "results/a/rows.json.gz",
        "results/a/net_ledger.jsonl",
        "results/a/fig.png",
        "results/a/old.json",
        "src/x.json",
    ]
    sizes = {
        "results/a/rows.json": 600_000,
        "results/a/rows.json.gz": 900_000,
        "results/a/net_ledger.jsonl": 10,
        "results/a/fig.png": 900_000,
        "results/a/old.json": 700_000,
        "src/x.json": 900_000,
    }
    bad = cr.size_violations(files, sizes, {"results/a/old.json"})
    assert len(bad) == 2
    assert bad[0].startswith("results/a/rows.json:") and "ledger" in bad[1]


def test_review_flags(tmp_path):
    (tmp_path / "out.json").write_text(json.dumps({"label": "derived", "n": 1732, "x": 0.0409}))
    with gzip.open(tmp_path / "ph.json.gz", "wt") as fh:
        fh.write(json.dumps({"label": "derived (post hoc)", "lost": 177}))
    (tmp_path / "PLAN.md").write_text("threshold 0.25\n")
    (tmp_path / "RESULT.md").write_text(
        "\n".join(
            [
                "- n is 1,732 and x 0.041; threshold 0.25.",  # all sourced (rounding)
                "- lost 177 injections.",  # post hoc number without label
                "- lost 177 injections (post hoc).",  # labelled: fine
                "- a share of 8.9 points.",  # not stored
                "- a share of 8.9 points (not stored).",  # labelled: fine
                "- this happens mainly because of gaps.",  # causal, no number
                "- confirmed real dimming.",
                "- **Deviations:** none.",
                "  - a note",
            ]
        )
        + "\n"
    )
    f = cr.review(tmp_path, [])
    joined = "\n".join(f)
    assert sum("item 2" in x for x in f) == 1 and "[177]" in joined
    assert sum("item 1" in x for x in f) == 1 and "[8.9]" in joined
    assert sum("item 3" in x for x in f) == 1
    assert any("confirmed real" in x for x in f)
    assert any("item 6" in x for x in f)
    assert "1,732" not in joined and "0.041" not in joined and "0.25" not in joined


def test_md_sections_newest_entries(tmp_path):
    c = tmp_path / "CHANGELOG.md"
    c.write_text("# Changelog\n\nintro\n\n## new\n- a 12.5\n\n## old\n- b 99.9\n")
    assert cr.md_sections(c) == ["## new\n- a 12.5\n"]
    d = tmp_path / "DECISIONS.md"
    d.write_text("# D\n\n## D-001 one\n- x\n\n## D-002 two\n- y\n")
    assert cr.md_sections(d) == ["## D-002 two\n- y\n"]


def test_hash_skip_keeps_plain_numbers():
    toks = [t for t, _ in cr.numbers_in_line("commit 038fdd0 read 1234567 bytes")]
    assert toks == ["1234567"]

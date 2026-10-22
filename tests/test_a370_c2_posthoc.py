"""Post hoc shared-pair counting of gate C2 controls (scripts/a370_c2_posthoc.py)."""

import importlib.util
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_SCRIPTS))
_spec = importlib.util.spec_from_file_location("a370_c2_posthoc", _SCRIPTS / "a370_c2_posthoc.py")
ph = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ph)


def test_shared_pair_across_systems_counts_once_and_absent_wins():
    systems = ["2", "102", "202", "4", "4"]
    pairs = [("r1", "c1")] * 3 + [("r2", "c2")] * 2
    classes = ["absent", "recovered", "absent", "recovered", "absent"]
    merged, groups = ph.merge_shared_pairs(systems, pairs, classes)
    # systems 2 / 102 / 202 merge into one absent control; one system's own images are not merged
    assert groups == [[0, 1, 2]]
    assert merged == ["absent", "recovered", "absent"]


def test_no_shared_pairs_leaves_classes_unchanged():
    merged, groups = ph.merge_shared_pairs(
        ["1", "2"], [("a", "b"), ("c", "d")], ["recovered", "absent"]
    )
    assert groups == []
    assert merged == ["recovered", "absent"]

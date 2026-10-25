"""Checks behind docs/efficiency-rules.md (owner rule 2026-10-26). Standard library only.

python scripts/check_results.py --sizes
    Hard check (pre-commit hook): a tracked plain-text output under results/ (.json, .jsonl, .csv,
    .ecsv, .tsv, .txt) larger than 500 kB, or a tracked raw download ledger
    (net_ledger*.jsonl), fails unless it is listed in scripts/check_results_grandfathered.txt
    (files that predate the rule). Write such outputs as .gz from the start.

python scripts/check_results.py --review results/<dir> [--md CHANGELOG.md DECISIONS.md ...]
    Advisory pre-review report for docs/pre-review-checklist.md items 1, 2, 3, 5 and 6. Numbers in
    <dir>/RESULT.md (and the newest CHANGELOG entry / newest DECISIONS record if given) that are
    not found in a committed data file of <dir> (.json, .jsonl, .csv, .ecsv, and their .gz) or in
    <dir>/PLAN.md, and are not on a line saying "not stored"; numbers found only in files labelled
    post hoc but on a line without "post hoc"; causal phrases on lines without a number; risky
    phrasings; a "Deviations: none" heading followed by notes. Exit 0 unless --strict and findings.
    The report is a reading aid: every flagged line still needs a person's judgement.
"""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRANDFATHER = ROOT / "scripts" / "check_results_grandfathered.txt"
TEXT_SUFFIXES = {".json", ".jsonl", ".csv", ".ecsv", ".tsv", ".txt"}
MAX_PLAIN = 500_000  # bytes
DATA_SUFFIXES = (".json", ".jsonl", ".csv", ".ecsv", ".json.gz", ".jsonl.gz", ".csv.gz", ".ecsv.gz")

NUM = re.compile(r"(?<![\w.#§/-])(\d{1,3}(?:,\d{3})+|\d+(?:\.\d+)?)(?![\w/])")
SKIP_CONTEXT = re.compile(
    r"\d{4}-\d{2}-\d{2}|\d{1,2}:\d{2}(:\d{2})?|#\d+|§\s*\d+|D-\d+|\bQ\d+\b|\bv\d+(\.\d+)*\b|"
    r"\b(?=[0-9a-f]*[a-f])[0-9a-f]{7,40}\b|arXiv:\S+|DOI\s*\S+|10\.\d{4,}/\S+|KIC\s*\d+|\bkepid\b"
)
CAUSAL = re.compile(
    r"\bbecause\b|\bdue to\b|\bdriven by\b|\bcaused by\b|\bmainly\b|\bmostly\b|\bcomes from\b|"
    r"\bresult of\b|因為|主要來自|導致",
    re.I,
)
RISKY = [
    (
        re.compile(r"confirmed real|確認真實", re.I),
        re.compile(r"not vetted|unvetted|未核查", re.I),
        "item 5: 'confirmed real' on an object that may not be vetted",
    ),
    (
        re.compile(r"achromatic", re.I),
        re.compile(r"not tested|not applicable|未測", re.I),
        "item 5: single band must read 'achromatic: not tested'",
    ),
    (
        re.compile(r"false[- ]positive rate|false[- ]pass rate|誤報率", re.I),
        re.compile(r"\bnot\b|conditional|條件|不是", re.I),
        "item 5: a conditional false-pass rate written as a sky false-positive rate",
    ),
    (
        re.compile(r"PDCSAP.*independent|independent.*PDCSAP", re.I),
        re.compile(r"\bnot\b|never|不能", re.I),
        "item 5: SAP and PDCSAP are not independent",
    ),
]


# --------------------------------------------------------------------------------------- --sizes


def tracked_files() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True)
    return [line for line in out.stdout.splitlines() if line]


def grandfathered() -> set[str]:
    if not GRANDFATHER.exists():
        return set()
    lines = GRANDFATHER.read_text().splitlines()
    return {s.strip() for s in lines if s.strip() and not s.startswith("#")}


def size_violations(files, sizes, allowed) -> list[str]:
    """Pure rule (tested): which tracked files break the output-size or raw-ledger rule."""
    bad = []
    for f in files:
        p = Path(f)
        if f in allowed:
            continue
        if re.fullmatch(r"net_ledger.*\.jsonl", p.name):
            bad.append(f"{f}: raw download ledger is tracked; commit a .jsonl.gz snapshot instead")
        elif f.startswith("results/") and p.suffix in TEXT_SUFFIXES and sizes[f] > MAX_PLAIN:
            bad.append(f"{f}: {sizes[f]} B plain text > {MAX_PLAIN} B; write it as .gz")
    return bad


def check_sizes() -> int:
    files = [f for f in tracked_files() if (ROOT / f).exists()]
    sizes = {f: (ROOT / f).stat().st_size for f in files}
    bad = size_violations(files, sizes, grandfathered())
    for b in bad:
        print(b)
    return 1 if bad else 0


# -------------------------------------------------------------------------------------- --review


def _load(p: Path) -> str:
    raw = gzip.open(p, "rt").read() if p.suffix == ".gz" else p.read_text(errors="replace")
    return raw


def _values(text: str) -> list[float]:
    vals = []
    for m in NUM.finditer(text):
        try:
            vals.append(float(m.group(1).replace(",", "")))
        except ValueError:
            pass
    return vals


def _is_post_hoc(text: str) -> bool:
    head = text[:2000].lower()
    return "post hoc" in head or "posthoc" in head


def data_numbers(d: Path) -> tuple[set[float], set[float]]:
    """All numbers in committed data files of d (recursively) and PLAN.md; second set: numbers
    found only in files labelled post hoc."""
    regular, post = set(), set()
    for p in sorted(d.rglob("*")):
        name = p.name
        if not p.is_file() or name.startswith("net_ledger"):
            continue
        if name.endswith(DATA_SUFFIXES) or name == "PLAN.md":
            text = _load(p)
            target = post if _is_post_hoc(text) and name != "PLAN.md" else regular
            target.update(_values(text))
    return regular, post - regular


def _decimals(tok: str) -> int:
    return len(tok.split(".")[1]) if "." in tok else 0


def _matches(x: float, dec: int, pool: set[float]) -> bool:
    tol = 0.5 * 10 ** (-dec) + 1e-12
    return any(abs(v - x) <= tol for v in pool) or any(abs(100 * v - x) <= tol for v in pool)


def numbers_in_line(line: str) -> list[tuple[str, float]]:
    clean = SKIP_CONTEXT.sub(" ", line)
    out = []
    for m in NUM.finditer(clean):
        tok = m.group(1)
        v = float(tok.replace(",", ""))
        if "." not in tok and v < 10:
            continue  # small integers (counts like 3, cell indices) are too common to check
        out.append((tok, v))
    return out


def md_sections(path: Path) -> list[str]:
    """RESULT.md whole; CHANGELOG.md newest entry; DECISIONS.md newest record; else whole file."""
    text = path.read_text()
    if path.name == "CHANGELOG.md":
        parts = re.split(r"\n(?=## )", text)
        return [next(p for p in parts if p.startswith("## "))]
    if path.name == "DECISIONS.md":
        parts = re.split(r"\n(?=## D-)", text)
        return [parts[-1]]
    return [text]


def review(d: Path, extra: list[Path]) -> list[str]:
    regular, post_only = data_numbers(d)
    findings = []
    for path in [d / "RESULT.md", *extra]:
        if not path.exists():
            findings.append(f"{path}: missing")
            continue
        for sec in md_sections(path):
            lines = sec.splitlines()
            for i, line in enumerate(lines):
                low = line.lower()
                where = f"{path.name}: '{line.strip()[:90]}'"
                for tok, v in numbers_in_line(line):
                    dec = _decimals(tok)
                    if _matches(v, dec, regular):
                        continue
                    if _matches(v, dec, post_only):
                        if "post hoc" not in low:
                            msg = f"item 2 [{tok}] only in a post hoc file, no 'post hoc' label"
                            findings.append(f"{msg}: {where}")
                        continue
                    if "not stored" not in low and "not counted" not in low:
                        findings.append(
                            f"item 1 [{tok}] not found in committed data or PLAN: {where}"
                        )
                if CAUSAL.search(line) and not re.search(r"\d", line) and "not counted" not in low:
                    findings.append(f"item 3 causal wording without a number: {where}")
                for pat, ok, msg in RISKY:
                    if pat.search(line) and not ok.search(line):
                        findings.append(f"{msg}: {where}")
                if re.search(r"Deviations:?\**:?\s*none\.?\s*$", line, re.I):
                    nxt = lines[i + 1] if i + 1 < len(lines) else ""
                    if nxt.startswith("  -"):
                        findings.append(f"item 6 'Deviations: none' followed by notes: {where}")
    return findings


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", action="store_true")
    ap.add_argument("--review", type=Path)
    ap.add_argument("--md", type=Path, nargs="*", default=[])
    ap.add_argument("--strict", action="store_true")
    a = ap.parse_args(argv)
    if a.sizes:
        return check_sizes()
    if a.review:
        f = review(a.review, a.md)
        for x in f:
            print(x)
        print(f"{len(f)} advisory finding(s); check each against docs/pre-review-checklist.md")
        return 1 if (a.strict and f) else 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())

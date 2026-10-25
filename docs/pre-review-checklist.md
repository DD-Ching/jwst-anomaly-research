# Pre-review checklist (owner rule 2026-10-26)

Work through this list **before** sending a result to the independent review. The review and the
targeted re-check use the same list as their baseline. The list changes no method, threshold, review
count or quality requirement; it only moves the wording checks before the review.

**Advisory run** (items 1, 2, 3, 5, 6):
`python scripts/check_results.py --review results/<dir> --md CHANGELOG.md DECISIONS.md`.
Every flagged line still needs a person's judgement, and an unflagged line is not proof of correctness.

1. **Every number has a committed source.** Each number in RESULT, CHANGELOG, DECISIONS and TASKS
   appears in a committed JSON / ECSV / CSV file (gzip included) or in the frozen PLAN. Otherwise the
   line says "not stored". Derived arithmetic (a percentage, a sum, a magnitude) counts as not stored
   unless the derived value itself is in a file. [script: partly; TASKS manual]
2. **Post hoc is labelled everywhere,** including the summary at the top and the CHANGELOG entry.
   [script: numbers found only in files labelled post hoc]
3. **Causal statements carry a checkable number.** This covers "because …", "mainly from …",
   "driven by …". Without a number, rewrite the line as a description or label it "not counted".
   [script: flags causal words on lines without a number; manual judgement]
4. **By-construction rejections are stated as such.** If a control family is rejected because of how
   it was built (its model can never pass a test), the write-up says so. It is not reported as a
   measured capability. [manual]
5. **No over-strong labels:**
   - nothing unvetted is "confirmed real";
   - a conditional false-pass rate is not a sky false-positive rate;
   - single band means "achromatic: not tested", never "passed";
   - SAP and PDCSAP (or any two products of one observation) are not independent confirmations.

   [script: phrase flags; manual]
6. **Deviations heading:** if a "Deviations: none" heading has sub-bullets, write
   "Deviations: none. Notes:". [script]

Also from `scripts/CLAUDE.md` "Before writing results": frozen rules as tested functions, narrow
wording, labelled deviations, lint before committing.

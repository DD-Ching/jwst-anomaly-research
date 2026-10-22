# Abell 370 gate C2 failure: post hoc diagnosis (result)

Plan: `PLAN.md` (committed first). Code: `scripts/a370_c2_posthoc.py`; test: `tests/test_a370_c2_posthoc.py`.
Numbers: `posthoc.json`. **Everything here is post hoc.** The frozen outcome of `../RESULT.md` stands: C2 fails,
no flag is read. Table diagnosis only: no image was looked at in this step, nothing was re-measured.

- **Reproduction:** `gate_c2` on the stored control tables gives the frozen arc-referenced counts on all five
  CANUCS-reference fields (Abell 370 19 / 25, MACS0416 58 / 59, MACS0417 10 / 10, MACS1149 5 / 5, MACS1423 2 / 4).
- **Abell 370 absent arc-referenced controls (6):**
  - 2 (2.4+2.5), 102 (102.3+102.4) and 202 (202.4+202.5) are three BUFFALO systems that use one CANUCS reference
    source (2101809) and one counterpart source (2101806); observed / predicted 0.135–0.173 each. The same
    segment pair is counted as three controls.
  - 3.1 (ratio 0.281; own segment a3 2.32″ vs reference 0.78″), 15.3 (ratio 0.018; own a3 0.17″ vs reference 1.86″)
    and 42.2 (ratio 0.301): three different source pairs. Cause not tested.
- **Counting each shared (reference, counterpart) pair once across systems** (absent if any member is absent):
  Abell 370 arc-referenced 19 / 23 = 0.826, all 35 / 40 → would pass (threshold 0.80, ASSUMPTION; one more absent
  control would fail it). MACS0416 49 / 50 (passes either way); MACS0417, MACS1149 and MACS1423 have no such
  pairs and are unchanged. So the counting changes no other field's outcome and does not explain MACS1423.
- **Reading (narrow):** half of the Abell 370 failure is one segment pair counted three times; the other half is
  three unrelated pairs. The counting rule was found on Abell 370, so it cannot be validated there.
- **Proposal (not applied):** freeze "one control per shared (reference, counterpart) pair, absent if any member is
  absent" before the next field this route runs on. Whether Abell 370's flags may be read under a post hoc pass is
  the owner's decision; this run reads none.

# MACS J1931 stage 3b: section covers every image of frame check F (deviation, frozen before its run)

**Deviation, post failure.** Stage 3 (`../macs1931_counter/PLAN.md`, "Stage 3: frozen rules") stopped at frame check F
(18 / 51 < 20; `RESULT.md`). Code review then found a plan defect: rule 4 bounds the F200W section by the control and
predicted positions (spec-z systems only), while rule 5 matches all 51 used images, and 3 of them (13.1, 13.2, 13.3,
photo-z) lie outside that section, so they could not match by construction (not stored before this plan: the
check was a session computation; 3b's `posthoc.json` stores `inside_section`). With 18 matched, F could flip if all three matched
(21 ≥ 20), so stage 3's F reading is not conclusive. This deviation was written after seeing 18 / 51; it changes only
the section and nothing that was tuned on the outcome.

## The one change
- **3b-1:** rule 4's section bounds the control and predicted positions **and all 51 used images**, plus 15″, clipped to
  the mosaic (`m1931_stage3.section_positions(..., all_images=True)`). Everything else, rules 1–10, is unchanged:
  builder parameters, F (0.5″, ≥ 20, |median| ≤ 0.1″), chain, footprint rule (now applied before local noise; it had
  not run), gate C2, flags.
- Output `../macs1931_stage3b/`. Run: `python scripts/m1931_stage3.py --maps DIR --plan-3b`, then `--posthoc --plan-3b`.
- Reading: if F fails again, no stage 3 reading (as before). If F passes, stage 3's rules 6–10 decide, and the result
  names this deviation. A widened section changes the background map and detections near the old edges, so the 3b
  catalogue is not identical to stage 3's.

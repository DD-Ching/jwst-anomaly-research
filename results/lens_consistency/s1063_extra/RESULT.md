# Abell S1063 WSLAP+ JWST 2025: extra-image test, WLS source, DJA v7 reference (2026-10-24)

Plan frozen in fdedcac with the functions and tests (`PLAN.md`), before any DJA row, F200W pixel or prediction was
read. Script `scripts/s1063_extra.py` (tests `tests/test_s1063_extra.py`); numbers from `summary.json`; tables
`predicted_abells1063.ecsv` / `control_abells1063.ecsv` (rules 5–7), `frame/` (rule 3), `without_contain/` (rule 8);
`imageplane_abells1063_abells1063-wslap.ecsv` (rule 2). Provenance: positions observed, predictions
model_prediction, classes derived.

**Deviation (labelled in the plan):** this stage gates on the WLS source after the mean-source gate K failed
(`../s1063_wslap/RESULT.md`); the WLS statistic was declared before that run (addendum k / Bullet precedent).

| Rule | Found | Status |
|---|---|---|
| 1 F: Chan+2020 vs DJA v7 frame | 22 / 36 matched ≤ 0.5″ (≥ 18), median dx +0.008″, dy +0.068″ (≤ 0.1″) | **passes** |
| 2 K-WLS repeated | rms 0.5805″, 36 / 36 within 3″; equals `../s1063_wslap/summary.json` | **passes** |
| 4a unit | zero point 23.90 mag → µJy (1,609 rows) | **passes** |
| 4b grid | rotation 0.0002°, determinant < 0, 0.0400″ / px | **passes** |
| 6 gate C2 (shared, local noise, containing-source rule) | 19 / 20 (0.95); arc-referenced 8 / 8 | **passes** |
| 7 flags | 0 | **no flag** |

**Reading.** The DJA v7 reference route (El Gordo stage 2b form) passes gate C2 on Abell S1063, but **none of the 9
predicted targets is testable**: 6 `inconsistent_reference` (systems 12 and 19: their catalogued images disagree in
μ-corrected flux beyond the frozen 3× spread), 1 `no_reference` (system 10, `candidate`), 1 `undetectable` (system 13,
`faint`), 1 `undetectable_local` (system 52). So this is no flag **and no reach**: the test does not constrain the
WSLAP+ counter-image predictions of S1063 under these rules. Image classes in the frame tables: 37 observed, 3
candidate, 5 missing, 1 faint, 8 demagnified (`rule3_frame_tables`); system 10.3 is unpredicted under the WLS source.

**Reported, not rules.** Without the containing-source rule C2 fails (15 / 20, arc 6 / 8) and there is no flag; without
the local-noise rule C2 passes (22 / 23, arc 10 / 10) with no flag; the pipeline gate C of rule 3 passes (8 / 8).

**Contact sheet** (`sheet.png`, display only; `scripts/s1063_extra_sheet.py`, 6″ GO 3293 F200W stamps, search circle in
red). Looked at before writing: the system 10 and 12 circles hold or touch bright cluster members, and system 13's
stamp is a noisy region crossed by a diagonal streak. Nothing in the stamps changes a class; the sheet is not used for
any verdict.

**Limitations.** The DJA v7 table has 1,621 rows (1,587 kept; reference depth 27.39 mag, `frozen.inputs`). Which
exposures it was built from was not checked (HTTP Last-Modified 2023-09-26, not stored; GO 3293 F200W is dated 2024-09-28
in its header, not stored); whether a deeper catalogue would make systems 12 and 19 consistent is untested. The reference fluxes are not bCG-subtracted. z is fixed at z_spec;
the 36 images are WSLAP+'s own constraints (not held out).

**Data.** Two α members re-fetched by range reads (sha256 in `summary.json`, equal to the gate-K run), held in the
session scratchpad and deleted after the run; DJA table (4,389,120 B) in memory only; GO 3293 SCI read once by a range
stream; the pipeline `_cat.ecsv` in the session cache, deleted by the script after rule 3.
After `/code-review` the code was revised (the model is passed to `m0417_extra.run_field` as a `model=` argument
instead of through the `install_model` hook named in PLAN rule 3, with the same behaviour: no registry lookup, no frame
offset; cache deletion; contact-sheet edge clipping) and the run repeated: `summary.json` and every table identical
(diff of the two runs not stored). Deviation in wording only.

## Next
- S1063 extra-image route closed as no reach under these rules. Reopen with a deeper JWST-era S1063 photometry
  catalogue (e.g. one built on GO 3293) that makes systems 12 / 19 consistent, or a released WSLAP+ image list.

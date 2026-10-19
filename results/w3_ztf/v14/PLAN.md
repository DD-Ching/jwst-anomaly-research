# W3 ZTF v14 — finish the #370 pending checks; then one exposure-independent area with the frozen v13 rule

Plan fixed 2026-10-19, before any new cutout and before any data of the new area are read. Owner authorization: the
round after #370 (4 h, ≤ 3 GB). No new discriminator, and the threshold stays fixed.

## 1. #370 pending items (≤ 40 new cutouts)
- **Terms:**
  - 96 **flagged runs** (pre-screen) in 6 groups;
  - **0 final survivors** (none passed vetting);
  - the 5 groups checked by image in #370 were a **diagnostic sample** of flagged groups, not candidates.
- **Inventory** (`pending_inventory_script.py.txt` → `pending_inventory.json`):
  - The 1.85 GB partition is not re-read. Run exposures are rebuilt from IRSA TAP and DR24 light-curve API
    detections, matched by expid.
  - An ND is an in-footprint, infobits-0 exposure of a member quadrant with no member detection.
  - The frozen good-exposure rule is not re-applied (approximation). Over all 96 runs, the rebuilt NDs (691) are
    close to the stored dim epochs (735); per run they agree exactly for 2 of 96 runs.
- **Selection** (`pending_select_script.py.txt` → `pending_selection.json`, frozen before any cutout):
  - group 724208200014557 (never checked) with v6 sampling;
  - for the other groups, pending NDs stratified by (run, field, CCD, quadrant). Strata on a detector quadrant not
    checked in #370 come first. One ND per stratum, ≤ 2 strata per group, plus 1 detected reference;
  - official mask cutouts (`_mskimg.fits`) for the selected NDs.
  - 26 cutouts in total (17 science + 9 mask). Everything else stays **pending, not excluded**.
- **Method:** the v6/v13 image-check core (`../v13/image_script.py.txt`, from the IBE search onward) is reused
  unchanged on this selection. The v13 column check is applied to the new cutouts. Mask bits are read at the source
  position and named from the mask header.
- **Wording:** a dark column without an official mask flag or native detector position stays "compatible with a data
  problem". The column relation is **not** added as a global veto.
- **Priority:** if any item shows a reproducible, image-confirmed dimming without an ordinary explanation, the rest of
  the budget goes to it and section 2 is not started.

## 2. One exposure-independent area: Npix 8170 (chosen in #368 by metadata, never read)
- **Re-check before any read** (from `../v12/obs_conditions.json`, `../v12/budget.json`):
  - worst-cell necessary fraction 0.6615 ≥ 0.5;
  - product ratio 1.05 ≥ 0.8;
  - fields 634, 635, 681 and 1679, which share no field with the development data (2382: 723, 724, 1722) or any
    earlier sample;
  - light-curve partition 2,207.9 MB.
  - If any check fails, the round stops without another area.
- **Budget:** 2,207.9 + objects (parent 2042, ≤ 100 MB) + TAP + part 1 (≈ 15 MB) + cutouts (≤ 10 MB) ≈ 2.35 GB, under
  3 GB.
- **Sample** (`select()`, cuts unchanged; dedupe against `../v12/dedupe.json` plus the v12/v13 2382 groups):
  - ≤ 5,500 groups (seed 7);
  - calibration = the lowest 500 by crc32("v14:" + key), search = the rest (≤ 5,000);
  - control stars: the comparison range without any sample member.
  - Groups bind all objectids of a source.
  - Brightness and sampling percentiles are reported for both sets; search-group sampling is reported only after the
    search.

## 3. Validation with the frozen v13 rule (`--edge time`, EDGE_FRAC 0.25, EDGE_NSIG 2; other thresholds unchanged)
- **Injections:** 10 cells × 200, seeds 5100 + cell (offset 4000), all injections in the denominator.
- **Controls:**
  - eclipse, PSPL and common mode × 200;
  - **fade × 1,000** (new option `--n-fade`; same fade parameter distribution; seed 5201);
  - the fade rate is reported with a Wilson 95 % interval.
- **Pass, all required:**
  - every cell ≥ 0.5;
  - every control class ≤ 10 % (point estimate, the existing criterion);
  - ≤ 5 real calibration survivors.
- **Pass:** one search of the frozen search group follows directly.
- **Fail:** record and stop. No rule change, no new seed and no other area.
- Simulated control rates are not real-sky false-alarm rates, and passing this control family does not cover all
  ordinary variables.

## 4. Image check after the search (only for final survivors)
- **Order:** vetted survivors by key. **Limits:** ≤ 5 groups, ≤ 64 cutouts.
- **Steps:**
  1. exposure, footprint, pairing and same-exposure comparisons;
  2. before / during / after the run;
  3. ordinary explanations versus W3 only after confirmation.
- Pending items stay pending. A missing record is unknown, never zero flux.

## 5. Stop
- At least 30 min are kept for review and records. Small inputs (samples, inventory, seeds, per-event outputs) are
  kept in the repo, and the large caches are deleted.

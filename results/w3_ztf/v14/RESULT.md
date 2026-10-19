# W3 ZTF v14 — #370 pending checks done (official masks); Npix 8170 passes validation; 5,000 groups searched; 1 vetted survivor, not a real dimming

- **Plan:** `PLAN.md`, committed in 2ed8929 with the frozen part-1 selection, before any new cutout.
- **8170 sample:** committed before any light curve (`8170/sample.json.gz`). The uncompressed file is 1.3 MB, over the
  repo's 1 MB rule, so it is git-ignored and rebuilt from the `.gz`; ra/dec are rounded to 1e-7 deg.
- **Labels:** TAP, DR24, IBE images and masks `observed`; counts `derived`; injections `simulated`.
- **Scope:**
  - the frozen v13 W3 rule, catalogue stage, the 10 calibrated cells;
  - no abundance limit;
  - no statement about wormholes in general.

## 1. #370 pending items
- **Terms (#370):**
  - 96 flagged runs in 6 groups (pre-screen);
  - 0 final survivors;
  - the image check of 5 groups was a **diagnostic sample of flagged groups**, not candidates.
- **Inventory** (`pending_inventory.json`; TAP + DR24 light-curve API, expid matching; no partition re-read):
  - 691 rebuilt NDs in the 96 run windows against 735 stored dim epochs. Per run they agree exactly for 2 of 96,
    because the good-exposure rule is not re-applied (approximation).
  - 26 of the NDs had been cut out in #370.
- **Frozen selection** (`pending_selection.json`):
  - the never-checked group 724208200014557 (6 science cutouts);
  - one ND per detector-quadrant stratum, up to 2 per group, plus 1 reference for 4 other groups (11 science
    cutouts). 724208200008944 had no pending stratum;
  - 9 official mask cutouts (`_mskimg`).
  - 26 cutouts in total (cap 40). **656 rebuilt NDs stay pending, not excluded.**
- **Results** (`pending_image.json`, `pending_image.png`, `mask_check.json`):
  - **3 new auto-classified "dimming" NDs** (723205400006654: 1; 724208200007532: 2). The official mask at the source
    position has value 516 = bit 2 **LOW RESPONSIVITY** + bit 9 **DEAD (UNRESPONSIVE)**, bit names from the mask
    header. The other masked ND of 723205400006654 (class 5) also carries 516 at the centre.
    - So for these exposures the official pipeline flags the pixels under the source as dead or low-response. That
      supports, beyond the #370 image-only dark-column relation, a detector origin for these NDs.
    - The 17 #370 NDs were not re-masked (cap).
    - The relation is **not** made a global veto.
  - **724208200014557** (6th group): 2 NDs are undetermined (class 5). Their masks show bit 7 (pixel spike) and bit 1
    (a SExtractor detection inside the aperture), and the only post-run reference is a poor image. **Not
    image-confirmed as dimming**; it stays undetermined.
  - **724208400014500** (86 runs): 2 more NDs show the source present (class 3, ratio 9.3 and 6.0), with mask bits
    1 and 8 (**SATURATED**). Read as catalogue rejections of a saturated source, not dimming.
  - **724208200006955:** 1 ND undetermined (mask bits 1 and 7).
- No item needed priority follow-up, so section 2 went ahead.

## 2. Npix 8170 (exposure-independent of the development data)
- **Re-check passed** (`sample_info.json` → recheck): worst-cell necessary fraction 0.6615, product ratio 1.05.
  Fields 634, 635, 681 and 1679 share none with the development fields (723, 724, 1722) or earlier samples.
- **Sample:**
  - 5,500 groups, 500 calibration and 5,000 search (crc32 rank);
  - 20,048 control stars, none a sample member, with no objectid shared between calibration and search;
  - refmag r (5/50/95 %): calibration 14.88 / 16.54 / 17.41; search 14.84 / 16.51 / 17.42.
  - r-band series epochs (5/50/95 %; `coverage.json`, computed after the search): calibration 1,036 / 1,849 / 3,526;
    search 1,017 / 1,849 / 3,526.
- **Data:** 49,509 light-curve rows (`8170/fetch_info.json`); TAP 61,483 rows.

## 3. Validation with the frozen v13 rule (`8170/decision.json`): **passed**

| cell (ρ, blend, t_E) | recovered / 200 |
|---|---|
| 0.01, 1, 10 | 136 |
| 0.01, 1, 30 | 161 |
| 0.01, 1, 100 | 156 |
| 0.1, 1, 10 | 123 |
| 0.1, 1, 30 | 152 |
| 0.01, 0.7, 10 | 141 |
| 0.01, 0.7, 30 | 158 |
| 0.01, 0.7, 100 | 127 |
| 0.1, 0.7, 10 | 114 (lowest, 0.57) |
| 0.1, 0.7, 30 | 149 |

- **Controls passing vetting:**
  - **fade 50 / 1,000 = 5.0 %** (Wilson 95 % interval 3.8–6.5 %);
  - eclipse 0 / 200, PSPL 0 / 200 and common mode 0 / 200 (each 0–1.9 %).
- **Real calibration light curves:** 999 series, 0 flagged runs.
- **Fade control across sets, under the same rule:** 17 / 200 (8.5 %) in v13 and 50 / 1,000 (5.0 %) here. The
  intervals overlap (v13 Wilson 5.4–13.2 %, computed in the session, not stored). The estimate is of the same order
  but not tightly stable.
- Simulated control rates are not real-sky false-alarm rates.

## 4. Search (`8170/search_holdout_v14_time.json`)
- **One search of the 5,000 frozen groups:** 9,995 series, 8 flagged runs in 4 groups, **1 vetted survivor**:
  635216100004733 (g band, 3 epochs, MJD 58976.48–58978.41).
- The other 7 runs fail the one-off test (and one also fails the achromatic test).

## 5. Image check of the survivor (`8170/image_check.json`, `image_check.png`; 7 cutouts, 0.56 MB)
- **Before the run:** 2 detected references.
- **The 3 run NDs:**
  - 122247916 and 122339453: maglimit 19.25 and 19.30, below the image check's 19.5 quality cut (class 2, poor
    image);
  - 122440822: the source is present at 0.88 × the reference flux (class 3, "catalogue missing it").
- **After the run:** detections 7–9 min after the third ND (122441409, 122441457).
- **The contact sheet shows the star in all three NDs.**
- **Reading:**
  - the source (g ≈ 18.9 AB from the 100.6 µJy reference) lies 0.3–0.7 mag above those exposures' limits;
  - the catalogue dropped it in shallow images;
  - this is **not a real dimming**, and no ordinary-versus-W3 comparison applies.
- Post hoc note, not a rule change: the frozen pre-screen counts any in-footprint ND as dim, whatever the source
  brightness relative to the exposure limit. That lets near-limit g-band sources pass.

## 6. Cost (`cost.json`)
- **Downloads ≤ 2,344.6 MB in total:**
  - 8170 partition 2,207.9 MB;
  - objects ≤ 93.1 MB;
  - TAP 13.2 MB;
  - part 1 TAP and LC API about 9.4 MB per inventory run (3 runs);
  - cutouts 1.4 + 0.2 + 0.6 MB.

## Deviations and notes
- **Part 1 inventory:** rebuilt without the good-exposure rule (approximation, stated above).
- **Sample storage:** gzipped (see the header).
- **Image script:** restricted to vetted survivors (PLAN section 4); that is the only change from `../v13/image_script.py.txt`.

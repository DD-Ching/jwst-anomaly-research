# W3 ZTF v3: counting fix and what a catalogue non-detection means — result

Plan `PLAN.md` (committed before any v3 number). Development half only (already seen, not blind); holdout unread.
Thresholds unchanged. Every number is `derived` from the files named; image work is exploratory.

## A. Counting fix (`../../../scripts/w3_ztf.py` v3; tests `tests/test_w3_ztf.py`)
- Exposures now come from IRSA's official metadata (`ztf.ztf_current_meta_sci`, 38,026 exposures of the 26
  quadrants; ~95 s by TAP, not stored); light-curve times are matched to each exposure's heliocentric mid-time within
  0.001 d (25,289 of 21,334,085 detections unmatched, 0.12 %; `v2_vs_v3.json`). Each comparison star counts once
  per exposure, the target's own group is excluded, and only comparison stars with detections before and after the
  exposure are in the denominator. Regression tests fail on the v2 code (duplicates, target as its own comparison,
  out-of-span stars).
- **Changes v2 → v3** (`v2_vs_v3.json`):
  - non-detections per development r-band star: median 0.04 % of good exposures (90th percentile 4.3 %, 99th 48 %),
    so v1's "3.4 %" was a binning artefact; the tail is real catalogue gaps (B below);
  - real survivors 9 → 8 (826210200000603 dropped; the other 8 unchanged, still in the two windows MJD 58657–58661
    and 59106–59109); flagged runs 1,619 → 1,520;
  - W3 recovery (blend 1): ρ 0.01: 0.46 / 0.605 / 0.69 / 0.655 for t_E 3 / 10 / 30 / 100 d (v2 0.455 / 0.56 / 0.71 /
    0.615); ρ 0.1: 0.41 / 0.585 / 0.525 / **0.30** (v2 0.335 / 0.565 / 0.62 / 0.24). Controls passing vetting
    unchanged: eclipses 0, fades 3 / 200, PSPL 0, common mode 0 (by construction).
  - The frozen pass conditions still fail (ρ 0.1, t_E 100 d: 0.30 < 0.5; 8 survivors > 5).

## B. What the catalogue non-detections are (`image_selection.json`, `forced_photometry.json`, `stamps.png`)
- Three v2 survivors fixed in the plan: 848201300006405 and 826215400000661 (window MJD 58655–58661) and
  826215100008147 (window 59106–59109); 23 exposures chosen by the plan rule before any flux; IRSA IBE positional
  search plus 120″ science-image cutouts, 1.34 MB downloaded in all (IBE answered 503 or reset several times; retries
  succeeded). Forced photometry at the fixed sky position (2″ aperture, 4–6″ annulus, MAGZP, local random-aperture
  error), with 1–5 field stars in the same cutouts as a common-change check (their median ratio to the pre-window
  level is 0.82–1.18 where measured).
- **The 11 catalogue non-detections:**
  - **class 1, no coverage or no valid observation: 6** — for 848201300006405 all 4 window exposures do not contain
    the position at all (IBE position search; the field pointing put the star off the quadrant); one exposure of
    826215400000661 has infobits 67108864 (flux there 1.07 × its median, but invalid by the rule); one of
    826215100008147 has non-finite pixels in the aperture (chip edge);
  - **class 3, source present in the image but missing from the catalogue: 5** — 826215400000661 at 0.56, 0.98 and
    1.15 × its pre-window flux (the 0.56 exposure has 3.9″ seeing, which a 2″ aperture undercounts), and
    826215100008147 at 0.98 and 1.02 × (in those exposures the star fell in the neighbouring quadrant q4 of the same
    CCD, so the catalogue entry used for this objectid had no detection);
  - **class 4, a real drop: 0.** No verified dimming.
- So the shared non-detection windows are coverage and catalogue effects (pointing offsets moving stars off or
  across quadrants, a flagged exposure), not astrophysical changes, for the three sources looked at.

## Answers to the two questions
1. **What the light-curve tables are still usable for (after the fix):** partial dims of detected stars that stay
   above the limit (W3 with a blend, or other dips ≥ 1 mag), at the measured stage efficiencies above (ρ = 0.01:
   0.61–0.69 for t_E 10–100 d; t_E 3 d: 0.41–0.46; large ρ t_E slow edges: 0.21–0.30). A non-detection alone is
   **not** usable as a disappearance: in the cases checked it meant no coverage, an invalid exposure, or a catalogue
   miss with the star still there.
2. **Can a few images fill the gaps reliably and cheaply?** Yes for confirming or rejecting flagged runs: cutouts are
   ~30–300 kB each, the classification separated all 23 exposures without ambiguity, and a whole check of 23 exposures
   cost 1.34 MB and a few minutes (IBE's intermittent 503s need retries). The remaining coverage gap can be closed
   from metadata alone: `ztf_current_meta_sci` carries each quadrant's corner coordinates (ra1–ra4, dec1–dec4), so
   "position inside the exposure" can be tested for every star without images, and objectids of the same star in
   neighbouring quadrants can be merged.

## If a new independent validation is run (not started)
- Scale: a new development set of 1,000 stars in a different HATS partition (fresh, unseen), with three additions
  frozen first: per-star footprint containment from the corner coordinates, merging of the same star's objectids
  across quadrants, and IBE forced photometry of every non-detection exposure of each vetted survivor (expected
  ≲ 10 survivors × ≤ 8 exposures). Measured costs: one light-curve partition 57–61 s, exposure metadata ~95 s per
  four fields, development search ~1 min, 4,000 injections ~3 min (session observations), cutouts ~60 kB average
  (`forced_photometry.json`: 1.34 MB for 23).
- Remaining gaps that this would not close: full disappearance is calibrated only through the image check (the
  catalogue has no upper limits); t_E ≈ 3 d events are recovered at ~0.4; slow-edged W3 (ρ t_E ≳ 10 d) at ≤ 0.30
  because the fade veto removes them — a model-versus-ordinary-fade overlap, not a bug. W3 remains a restricted
  negative-mass hypothesis; nothing here is a wormhole candidate.

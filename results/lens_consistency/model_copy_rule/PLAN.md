# D-041 model-copy rule as a function (plan, frozen before the run)

**Question.** Can D-041's model-copy rule be computed instead of read by eye, so that `forced_check` flags
(`absent` / `confused`) that are the model's own copy of a catalogued, unpredicted image are labelled
automatically on every field?

**Rule (ASSUMPTION thresholds, from D-041).** A predicted row with `forced_class` in {`absent`, `confused`} is a
`model_copy` of catalogued image I when I belongs to the same system, I is in the run's
`unpredicted_catalogued_images`, and the separation is ≤ 1.75 × the model's image-plane rms. Nearest such I wins.
D-041 also requires I to be detected (S/N ≥ 5); that needs photometry that the stored tables lack, so the function
takes an optional detected set and records whether the detection condition was applied.

**Known case (gate G).** MACS0717 `macs0717-cats` (stored `results/lens_consistency/macs0717-cats/`), rms 3.21″
(D-041 Evidence). The hand classification (docs/fields/macs0717.md) counts 28 `confused` copy rows (systems 1 ×3,
3 ×2, 7, 12, 17, 24 ×2, 32 ×2, 34, 45 ×2, 50, 56, 62 ×2, 67, 71 ×2, 72, 75, 76, 79, 80 ×2; the 45 / 79 / 80
three are copies of undetected images) and 3 `absent` copy rows (systems 5, 6, 34). G passes when the
positional rule (no detection condition) flags exactly these 31 rows by per-system and per-class counts among
the `absent` + `confused` rows. Otherwise the rule is not adopted as a class and the mismatches are listed.

**Ordinary explanations this targets.** Model position error (rms 3.2″) plus a neighbour inside the search circle.
**What changes the plan.** A pass makes the rule a reported column (`copy_of`) for future `images --forced-image`
runs; a fail keeps the per-row reading by eye.

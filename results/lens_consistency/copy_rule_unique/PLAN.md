# D-042 second copy rule as a function (plan, frozen before the run)

**Question.** D-041's positional copy rule (1.75 × rms, nearest image) fails gate G on MACS0717
(`../model_copy_rule/RESULT.md`). D-042 used a stricter, one-to-one second rule by hand on MACS0416. Can that rule be
computed, so `absent` / `confused` rows that are the model's copy of an unpredicted catalogued image are labelled
automatically?

**Rule (ASSUMPTION thresholds, D-042).** Within one system, take the rows with `forced_class` in {`absent`,
`confused`} and the system's catalogued images that are unpredicted (no prediction matched to them). A row and an
image are a pair when their separation is ≤ 6.5 × the model's image-plane rms. A row is a `copy_of` image I when the
row has exactly one such image (I), I has exactly one such row, and the model's magnification at I (same source
redshift) has the row's sign (same parity). Everything else gets no label.

**Known case (gate U).** MACS0416 `macs0416-cats` (rms 0.76″ as D-042 records it, so the radius is 4.94″) `images` re-run with the recorded settings (docs/fields/macs0416.md:
F277W forced image by S3 byte range, frame offset +0.208″, −0.025″). D-042's hand reading: copy rows 1, 15, 122 ×2
(confused) and not 45 (opposite parity, μ −9 at 1.89″ from 45.2). U passes when the rule labels exactly these four
rows among the `absent` + `confused` rows (per-system counts) and does not label 45. Before the gate, the re-run must
reproduce the recorded forced counts (recovered 3, confused 13, absent 2); otherwise the run stops as `route
failed`.

**Then (diagnostic, only if U passes).** Apply the rule to the stored MACS0717 table (hand list in
`../model_copy_rule/`), MACS1149 and the two Abell 2744 extra-image tables and list agreements and differences with
the recorded hand verdicts. No verdict changes from this step.

**Ordinary explanation this targets.** Model position error plus a neighbour inside the search circle.
**What changes the plan.** A pass with agreement on the other fields makes the rule a candidate `copy_of` column
for `images --forced-image` runs; a fail keeps the by-eye reading and closes this route.

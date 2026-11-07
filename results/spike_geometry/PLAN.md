# Spike-geometry check for single-band elongated detections (TASKS Then 6, D-014 "Revisit if")

Frozen 2026-11-07 before any measurement on the target sample.

**Hypothesis.** Single-band detections that pass the D-014 gate only through a coincident DJA object, and that are
highly elongated next to a bright source, are diffraction-spike fragments (CEERS `f200w_2842`, run
`20261007T202419Z-65728e80`: F200W only, ellipticity 0.95, bright source 4.6″ away, DJA object 0.04″ off).

**Ordinary explanations of a positive flag.** A real edge-on or lensed galaxy that happens to point at a bright
neighbour (chance alignment); a tidal feature of the neighbour.

**Rule (ASSUMPTIONs, all frozen here).** A source is `spike_geom` when all hold:
1. detected in one band only in the merged pipeline catalog (`n_bands == 1`);
2. ellipticity ≥ 0.8;
3. a brighter source (≥ 3 mag brighter in the same band) lies within 10″;
4. the source's major axis points at that source: |Δθ| ≤ 10° between the major-axis angle and the bearing to it.

Radial alignment to a compact bright source is what a spike does; a lensed arc is tangential to its deflector, so
condition 4 should spare arcs.

**Validation before trusting it.**
- Known positive: CEERS `f200w_2842` must be flagged.
- Known negatives: catalogued SMACS lensed images / arcs (MJR2023 and ICLv2 multiple images) in the gated sample:
  ≤ 1 flagged, else the rule is rejected.
- Chance rate: rotate each source's major axis by 90° (a null that keeps positions and neighbours); the flagged
  fraction under rotation estimates chance alignment. The rule is useful only if real-sample flags exceed 3× the
  rotated count.

**What changes the plan.** If the positive is missed, or arcs are flagged, the rule is not adopted and the spike
line stays at D-014's confirmation gate. If it passes, it is reported as a quality flag (ranks unchanged) and a
wiring proposal, not wired in this PR.

## Addendum B (frozen 2026-11-07 ~08:20Z, after the rule above failed and after viewing the positive's cutout)
The 16″ cutout shows `f200w_2842` on a straight line crossing the whole cutout in F200W (fainter in F356W) that does
not pass through the 19.3 mag source 4.6″ away. New hypothesis: single-band streak detections lie on long linear
features (off-field spikes or trails), which a trail finder sees regardless of any neighbour.

**Tool (reuse, no new code path):** `scripts/w1_streak_mask.py` functions with the W1 settings unchanged
(`stream_binned`, `normalise`, `find_trails`, status-2 trails, `VETO_ARCSEC` 1″; results/w1_streak_veto/PLAN.md),
run on the CEERS t021 and SMACS F200W `_i2d` (streamed, not stored). Arm 2: the same plus the row/column-segment
stripe mask at `--seg-len` 64 (results/w1_row_veto/PLAN.md; ASSUMPTION).

**Adopt an arm as a proposed quality flag `on_streak` only if all hold (computed by `streak_verdict`):**
1. the positive `f200w_2842` lies within 1″ of the mask;
2. ≤ 1 of the matched ICLv2 SMACS images lies within 1″ of the mask;
3. in each field with ≥ 1 masked single-band source, the in-mask fraction of single-band sources is ≥ 3× that of
   sources detected in ≥ 3 bands (enrichment; a mask that hits everything equally is not a streak flag).

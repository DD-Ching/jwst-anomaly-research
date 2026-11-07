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

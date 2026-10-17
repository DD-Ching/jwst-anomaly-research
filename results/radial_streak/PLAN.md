# Radial (W1) screen with the calibrated streak veto: plan (2026-10-17, written before any run)

**Question.** Does the adopted W1 streak veto (`results/w1_streak_veto`, `w1_streak_mask.VETO_ARCSEC` = 1″ around
`findsat_mrt` status-2 trails) remove the streak-segment "arcs" that build `radial` centres, and does any `radial`
result change on the twelve W1 fields?

**Hypothesis under test (ordinary).** Linear streaks cut into several elongated catalogue segments whose major axes
are collinear, so their lines meet near the streak and inflate `max_lines` (seen by eye at MACS0647's 6-line centre,
`results/radial_more/RESULT.md`). Exotic reading (W1 radially stretched, demagnified images around a dark centre) is a
hypothesis only.

**Method.** `scripts/radial_more.py --streak-masks results/w1_streak_veto/masks`: shape rows within 1″ of a masked
pixel are dropped before `exotic_screens.cmd_radial`; everything else as in `results/radial_more/PLAN.md` (CLI
defaults, Gaia spike veto, spec `max_radius`). The null draws use the same vetoed table.

**Known case first (gate K).** MACS0647's strongest centre (6 lines, p 0.075): at least 4 of its 5 streak segments
must be vetoed and its line count must drop. If K fails, the veto is not used for `radial` and the field runs are
reported as a diagnostic only.

**Then.** All twelve W1 fields with the veto (the five of `results/radial_more` and the seven D-049 fields that have
masks). Flag rule unchanged: a centre with p_random < 0.05 (ASSUMPTION). Any flag goes to `/vet-candidate`; models
that fail the `arcs` control (RXJ2129, MACS0647) give no result.

**What would change the plan.** A flag in a control-passing field that the veto newly creates (fewer rows also
lowers the null) → vet it before any reading; no change anywhere → the radial streak item closes.

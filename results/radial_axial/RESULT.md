# Radial (W1) screen: centre-level axial (collinear-streak) rule — result (2026-10-17)

Plan: `PLAN.md` (committed before any run). Code: `exotic_screens.axial_concentration`, gate script
`scripts/radial_axial.py`. Table: `gate_k_macs0647.json`. Labels: shapes `observed`; centres, azimuths and A
`derived`; every threshold an ASSUMPTION.

## Inputs reproduced
`scripts/w1_stage_inputs.py --fields macs0647`, then `scripts/radial_more.py --fields macs0647` (no streak masks):
`max_lines` 6, `p_random_max` 0.075, 2 centres, as in `results/radial_more` (`gate_k_macs0647.json`).

## Gate K (known case): **failed**
The 6-line MACS0647 centre (model frame (2.5, −18.5)″) has axial concentration A = 0.731, below the frozen 0.9
(`gate_k_pass: false`). So the rule is **not adopted** and gate I (injections) was not run, as the plan says.

Why (sheet `sheet_macs0647_top.png`, looked at): the faint diagonal streak passes through the centre and carries
five of the six contributing arcs; one arc lies off it. The grid centre sits up to the 1″ line tolerance off the
streak, so streak segments a few arcsec away take azimuths spread over 48°–67° (mod 180), and the one on the other
side 98°: the azimuth about a grid point does not isolate a streak. The off-streak arc lowers A further.

## Post-hoc diagnostic (not a rule)
`post_hoc_max_on_one_line` (most contributors within 0.5″ of one straight line through two of them): **5 of 6** at
the MACS0647 centre, 2 of 3 at the second centre. Collinearity of the arcs themselves, independent of the grid
point, does pick out the streak in this one case. Since this case has now been looked at, it cannot serve as a blind
gate for that rule.

## Reading
No change to any `radial` conclusion (`results/radial_more`, `results/radial_streak`, D-049); nothing for
`/vet-candidate`. **Failed approach:** azimuthal axial concentration about the grid centre as a streak rule.

## Next (handoff)
A collinearity rule (contributors on one line within a tolerance count once) needs its own frozen plan, a different
known streak case (or injected streaks) as the gate, and W1 injections: one W1 source's two images are collinear with
the lens, so injected lenses with few lensed sources could lose recoveries to it.

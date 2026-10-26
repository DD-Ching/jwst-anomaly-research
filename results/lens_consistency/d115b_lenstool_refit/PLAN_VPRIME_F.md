# D115b Lenstool re-fit: stages V′ and F (plan, frozen 2026-10-26 11:27Z, before any run of this round)

Continues `PLAN.md` / `RESULT.md` (stage V failed as frozen on `best.par`'s 6-decimal reference). Tool
`scripts/d115b_lenstool_refit.py` (`--stage vprime-f`), rules tested in `tests/test_d115b_lenstool_refit.py`, both
committed before the run. Lenstool 8.15.5 (conda-forge), installed as in D-109.

## Stage V′ (a labelled deviation from stage V)
**Hypothesis:** with `input.par`'s full-precision reference (`64.03821668520617 -24.06750122335985`) in place of
`best.par`'s, and nothing else changed from the stage V run file, Lenstool reproduces the header χ²_pos.
The post hoc stage V value (344.42) is known, so this gate is a confirmation by the committed tool, not a blind test.
**Gate V′ (frozen):** `gate_v(chimul)`: |χ²_pos − 344.298308| ≤ 1 % (ASSUMPTION) of 344.298308.
Fail → stop; no stage F number is read.

## Stage F: χ² profile in member 8515's σ, nothing else changed (only if V′ passes)
Run file per factor f ∈ {1.00, 0.90, 0.85, 0.80, 0.75, 0.70}: the V′ file with `potentiel 8515`'s `v_disp`
set to f × 120.600977 km/s and `image 1 arcs.dat` added to `runmode` (Lenstool's `image.all`, as in D-109).
Lenstool's image-plane χ² uses each family's back-traced barycentre as the source; no parameter or source is fitted,
so this is Lenstool's analogue of the port's "nothing refitted" σ × 0.7 cost (D115b record, +12.6 with only D115
refitted), not the free re-fit.

Per f (all `model_prediction`): χ²_pos(f) (`chimul`), Δχ²(f) = χ²_pos(f) − χ²_pos(1.00), and
flag(f) = any `image.all` image of D115b.1 or of D115b.2 within 0.7″ of the #458 flag (64.048551, −24.064825).

**Rule F (frozen, function `rule_f`):**
- `undetermined` if flag(1.00) is false, or flag(f) is true at every f.
- Otherwise f_L* = the largest f with flag(f) false; `cheap` if Δχ²(f_L*) ≤ 6 (ASSUMPTION, system-51 precedent, as in
  the port scans), else `costly`.
- `cheap` → the extra image goes away in Lenstool's own χ² at an accepted cost; the record's verdict question becomes
  whether a re-fit keeps that (a proposal for `/vet-candidate`, not a verdict change here).
- `costly` → agrees in direction with the port scans (+12.28 at f* 0.75 with all sources fitted); the free re-fit
  stays the only step.

## Reported, not gated
- Stage F2 cost probe: Lenstool `inverse 3` (input.par's rate 0.015) with 8515's `v_disp` free in [60, 200] km/s and
  potentials O1–O4 freed with `input.par`'s `limit 1`–`limit 4` ranges, started from `best.par` with the V′
  reference; wall cap 10 min. Only its run time and whether it finished are recorded; no χ² from an unfinished run is
  read. It sizes the F2 plan.

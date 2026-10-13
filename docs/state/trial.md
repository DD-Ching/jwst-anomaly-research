# High-risk exploration trial (owner-approved 2026-10-12 → 2026-11-09; D-083)

**Goal:** let research-motivated high-risk units deepen across many runs and produce testable predictions,
reliable methods or real candidates. PR count, test count and compute time are **not** progress measures.
Every run reads this file during dispatch (research-cycle step 1). **Rules version: trial-2026-10-13.1** (adds the
every-run tally to trial-2026-10-12.1). **Every run**, routine searches included, appends a tally line below. A run
that works on a unit also records this rules version and its routine prompt version in the unit's session log; the
prompt version comes from its own prompt text, or is "none" for an interactive session or an older prompt. CLAIM
comments carrying `rules: trial-2026-10-12.1` remain valid.

## Rules
**T1. Weekly split** (a trial setting, checked per ISO week, not per run). Shares are of run time:

| Share | Use |
|---|---|
| 45 % | Routine searches |
| 30 % | High-risk units (at most **two** in depth) |
| 15 % | Verification (independent checks, vetting, held-out data) |
| 10 % | Upkeep (owner reviews, CI, base merges, this index) |

- A run reads the week's tally below.
- When an authorized unit has an open stage with budget left and no blocker, and its category is under its share,
  the run takes that unit before a routine-search unit.
- Unused shares are not filled with invented work.

**T2. WIP counting** (replaces "all needs-human PRs block new work"):
- **Awaiting merge:** non-draft agent PRs **without** `needs-human`. The cap of 3 applies to these.
- **Awaiting owner review:** non-draft `needs-human` PRs. They never block an authorized unit's continuation,
  base merges, fixes, or routine-search draft claims.
- **Owner-review backlog:** with ≥ 10 such PRs, open **at most one new needs-human PR per ISO week** from
  high-risk or theory work, and only from a unit stage that passed its gate.
- **Not limited by the backlog rule:**
  - a candidate, meaning a screen hit that survived `/vet-candidate`'s ordinary-explanation tests;
  - a routine-search result that needs owner judgement (for example a chain survivor, as in #145). Holding it back
    could hide a real signal.
  - Both are opened as needs-human PRs when they arise.
- Unvetted screen hits are not candidates. They stay on their draft claim PR until vetted.
- Draft claim PRs and unit PRs are in-flight work and count toward neither.
- Human merge rights for scientific results are unchanged.

**T3. Units.** Each authorized unit has a state file in `docs/state/units/<id>.md`. The file holds: stages with
pre-written gates, the budget per stage and in total, budget used, reusable artifacts (paths and commits), the
next concrete action, and stop conditions.
- A run continues from the file. It never re-downloads, rebuilds or re-reads what the file marks as verified.
- A stage that passes its gate continues automatically within budget, with no owner question.
- "Not computed yet" or "needs a new tool" is not a stop: first check whether the budget covers it.
- Time already spent is never a reason to continue past a failed gate.
- A failed gate pauses the unit and frees its slot. A replacement comes from a screened (S0-passed) direction;
  if none qualifies, the slot stays empty.
- Paused lines reopen only with a new method or input listed in `docs/state/pending.md` or the unit file.

**T4. Routine searches** (W3 MOA-II fields and similar): the "≥ 10 % change in a limit cell" figure is a
**reference only**.
- Judge a batch of fields on: cumulative limit gain, cost per field, new parameter-space coverage, and
  dependence on the injection model.
- A single field below 10 % is not a stop by itself.

**T5. Labels.** Results are labelled observed / derived / simulated / figure-read / assumption / hypothesis.
- Simulations and figure readings are never observational evidence.
- Data already used to choose a model or tune a method is never called held-out. Without true held-out data, the
  analysis is labelled exploratory.

**T6. Notifications** go to the owner only for: a stage completed, an important result, a failure, or a decision
needed. Unchanged CI, subscriptions, waiting for review and routine syncs are not reported.

## Units
| Unit | State file | Branch / PR | Status |
|---|---|---|---|
| X1 far-side SV hot spots (Sgr A*) | `docs/state/units/X1.md` | `claude/x1-s3` (#179; S2 in #177) | authorized; F0, S1 done 2026-10-12; S2 done 2026-10-13 (formal pass, a = 2.1 M, flux contrast only; fragile); S3 done 2026-10-13 (exploratory); next S4 |
| X2 | — | — | empty; filled only by an S0-passed direction that differs from X1 in mechanism, observable or data |

## Weekly tally (hours of run time; **every** run appends one line, routine searches included)
T1 needs every category, or the shares cannot be computed. A run that only works on a routine-search unit still
adds its line, with a note naming the unit.
| ISO week | Routine search | High-risk | Verification | Upkeep | Note |
|---|---|---|---|---|---|
| 2026-W42 | — | 1.5 (X1 F0) | — | 1.0 (trial set-up) | set-up session, not a routine run |
| 2026-W42 | — | 0.5 (X1 S1) | — | 0.2 (orientation, PR) | cloud routine 01CJcx3Y |
| 2026-W42 | ≈ 0.6 (gb14 pre-screen chunks 17–45) | — | — | — | cloud routine, 00:09–00:46Z; back-filled from the gb14 handoff by an interactive session |
| 2026-W42 | — | 0.8 (X1 S2, S3) | — | 0.15 (orientation, takeover, PR) | cloud routine 01ELiBXy |

## Outcome ledger (four columns only)
| Method improvements | Conditional theory results | Real candidates | Independent verifications |
|---|---|---|---|
| M2 stretch grid + free per-band scale (#167, provisional); W3 chain fixes (D-068) | Goal rounds 1–10 (#152–#166, provisional); far-side image bound vs GRAVITY summary (#172, provisional); X1 S2: no astrometric far-side discriminant at 9 µas, a fragile flux-contrast one for a = 2.1 M (#177); X1 S3: far side Δχ² +28 to +37 on the 2018-07-22 flare (a ≤ 7 M, exploratory, #179) | none (n0153 inconclusive, #61) | none of a candidate (SkyLLH D-080 is a benchmark reproduction) |

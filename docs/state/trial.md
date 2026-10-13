# High-risk exploration trial (owner-approved 2026-10-12 → 2026-11-09; D-083)

**Goal:** let research-motivated high-risk units deepen across many runs and produce testable predictions,
reliable methods or real candidates. PR count, test count and compute time are **not** progress measures.
Every run reads this file during dispatch (research-cycle step 1). **Rules version: trial-2026-10-13.2** (13.1 added the
every-run tally; 13.2 lets T1 outrank stale routine claims and fixes the 'ordinary model' wording). **Every run**, routine searches included, appends a tally line below. A run
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
- When an authorized unit has an open stage with budget left, no blocker and no live claim, and its category is
  under its share, the run takes that unit before a routine-search unit. This holds even when a stale
  routine-search claim exists (trial-2026-10-13.2). Each run completes and saves one step of the unit's current
  stage: a tool, theory or data step.
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

**T7. S0 screens for new units.** A screen takes ≤ 2 h of high-risk time and compares at most three
directions. Each must differ from earlier units in mechanism, observable or data. For each, answer:
1. the concrete unknown, and how it bears on channel existence, control or propagation;
2. how the prediction differs from **explicit ordinary models over their plausible parameter ranges** (not
   a proof that no ordinary explanation could exist);
3. the minimum verification: which data, tools or derivations it needs, which already exist, and which fit in
   the budget;
4. which result continues the unit and which stops it.

A direction does not need an existing anomaly, and "not computed yet" is not a failure. A new theory name alone
is not a reason to pick one. If one qualifies, write its state file with stage budgets and start it under the
existing authorization. If none qualifies, record the gaps and stay idle.

## Units
| Unit | State file | Branch / PR | Status |
|---|---|---|---|
| X1 far-side SV hot spots (Sgr A*) | `docs/state/units/X1.md` | `claude/x1-s4b` (#181) | **closed 2026-10-13** (gates passed as written; 3.7 / 20 h; tool cross-checked; position comparison exploratory; the S2 flux-contrast gate observable untested by data; no candidate). Conditional results: no astrometric discriminant at 9 µas; a fragile flux-contrast one at a = 2.1 M; S3's Δχ² +23 to +37 needs the traced centroid (model-free only +6.8, a ≤ 3 M). Tracer angles and weights verified by second routes. Reopen only with the inputs in its Closure |
| X2 SV-throat echoes vs published post-merger limits | `docs/state/units/X2.md` | `claude/x2-s1-s2` | S1–S3 done for one event (GW150914; phenomenological model P1–P4; 5.3 / 16 h). η_max 0.70–1.24 in the applicable a-range (a − 2M ≈ 2.2e-4 to 7.1e-4 M); no constraint at η ~ 1. S1 check 4 (a − 2M = 1e-3) and S3 spacing (+4.1 %) fail the pre-written 1 % (reflection delay; timing sensitivity recorded). **Paused for owner review**; next if resumed: S2b multi-event (≤ 2 h). Missing input: SV spin-2 equations, the source excitation η |

## Weekly tally (hours of run time; **every** run appends one line, routine searches included)
T1 needs every category, or the shares cannot be computed. A run that only works on a routine-search unit still
adds its line, with a note naming the unit.
| ISO week | Routine search | High-risk | Verification | Upkeep | Note |
|---|---|---|---|---|---|
| 2026-W42 | — | 1.5 (X1 F0) | — | 1.0 (trial set-up) | set-up session, not a routine run |
| 2026-W42 | — | 0.5 (X1 S1) | — | 0.2 (orientation, PR) | cloud routine 01CJcx3Y |
| 2026-W42 | ≈ 0.6 (gb14 pre-screen chunks 17–45) | — | — | — | cloud routine, 00:09–00:46Z; back-filled from the gb14 handoff by an interactive session |
| 2026-W42 | ≈ 0.6 (gb14 pre-screen 46–50, fits) | — | — | 0.1 (orientation, base merge) | cloud routine 01RHJ74R, 01:08–01:50Z |
| 2026-W42 | ≈ 0.6 (gb14 injection chunks 1–4 of 6) | — | — | 0.1 (orientation, claim) | cloud routine 01BiV92m, 02:08–02:52Z |
| 2026-W42 | — | 0.8 (X1 S2, S3) | — | 0.15 (orientation, takeover, PR) | cloud routine 01ELiBXy |
| 2026-W42 | ≈ 0.55 (gb14 injections 5–6, re-vet, limit, combine) | — | — | 0.15 (orientation, base merge, PR) | cloud routine 01B15J5c, 04:08–04:50Z |
| 2026-W42 | — | — | 0.5 (X1 S4: model-free disk bound, einsteinpy check) | 0.15 (orientation, claim, PR) | cloud routine 01BSu1H7, 05:08–05:50Z; S4 is an independent check, so it counts as verification |
| 2026-W42 | — | — | 0.3 (X1 S4b: image weights, second route) | 0.15 (orientation, claim, X1 closure, PR) | cloud routine 01QLViGu, 06:08–06:30Z |
| 2026-W42 | ≈ 0.35 (Abell 2744 `fluxratios`, D-032 addendum) | — | — | 0.15 (orientation, claim, PR) | cloud routine 01EaLT96, 07:08–07:35Z |
| 2026-W42 | ≈ 0.25 (`fluxratios --mu-samples`, Abell 2744 + El Gordo) | — | — | 0.15 (orientation, claim, PR) | cloud routine 01WdD16D, 08:09–08:45Z |
| 2026-W42 | — | — | ≈ 0.3 (halo-subtracted photometry of Abell 2744 flux flags) | 0.15 (orientation, claim, PR) | cloud routine 01FwFCPF, 09:08–09:45Z |
| 2026-W42 | — | — | ≈ 0.25 (halo_phot final review: injection through the chain, pool, errors) | 0.1 (takeover, PR) | cloud routine 0117Bja4, 10:08–10:45Z |
| 2026-W42 | ≈ 0.45 (halo_phot on all images of three fields; El Gordo 1c vetted) | — | — | 0.15 (orientation, claim, PR) | cloud routine 01CTJqMS, 11:08–11:50Z |
| 2026-W42 | — | 1.2 (S0 screen after X1) | — | 0.5 (X1 scope wording, prompt 2026-10-13.1, rules 13.2) | interactive session 01AzKMDv |
| 2026-W42 | — | — | ≈ 0.5 (halo_phot elongation gate, star calibration, three-field re-run) | 0.15 (orientation, claim, PR) | cloud routine 01Tan1i4, 12:10–12:58Z |
| 2026-W42 | ≈ 0.35 (`fluxratios --halo-phot`, Abell 2744; A2744-QSO1 flag vetted) | — | — | 0.15 (orientation, claim, PR) | cloud routine 01LzuV9A, 13:09–13:45Z |
| 2026-W42 | — | 1.4 (X2 S0 addendum: GrayHawk-v2, echo code, B critique, delay check) | — | 0.2 (state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 2.3 (X2 gate fix, S1, S2 one event) | 0.4 (X2 S3 time-domain check) | 0.2 (state files, PR) | interactive session 01AzKMDv |

## Outcome ledger (four columns only)
| Method improvements | Conditional theory results | Real candidates | Independent verifications |
|---|---|---|---|
| M2 stretch grid + free per-band scale (#167, provisional); W3 chain fixes (D-068) | Goal rounds 1–10 (#152–#166, provisional); far-side image bound vs GRAVITY summary (#172, provisional); X1 S2: no astrometric far-side discriminant at 9 µas, a fragile flux-contrast one for a = 2.1 M (#177); X1 S3: far side Δχ² +23 to +37 on the 2018-07-22 flare (a ≤ 7 M, exploratory, #179); X1 S4: model-free bound gives only +6.8 (a ≤ 3 M) and ≤ 0 (a ≥ 5 M), so S3 rests on the traced centroid (#180) | none (n0153 inconclusive, #61) | none of a candidate (SkyLLH D-080 is a benchmark reproduction); X1 tracer swept angles and image weights checked by second routes (#180, #181); X2 test-field scattering: Numerov vs DOP853, time domain vs frequency domain (this PR) |

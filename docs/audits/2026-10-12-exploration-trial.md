# Operations audit and a four-week trial for high-risk exploration (2026-10-12)

> **Adopted 2026-10-12 with owner corrections (D-083).** The operative rules are in `docs/state/trial.md`
> (version trial-2026-10-12.1). They supersede R1–R5 below where the two differ. Corrections:
> - §2 item 1 overstated the WIP effect;
> - the 10 % figure in R5 is a reference only (T4);
> - the split is checked per week, not per run;
> - X1's stages were re-budgeted after the F0 feasibility check (`docs/state/units/X1.md`).

Owner brief 2026-10-12: diagnose how the long-running process actually works, then propose how to fit high-risk
exploration into it without rebuilding the agent system. **This is a plan. Nothing here changes the schedule,
the routine prompt, settings or any merge.** The owner decides what to adopt.

## 1. What is running, and how sure we are

| Item | Status | Evidence |
|---|---|---|
| Hourly cloud routine `trig_01PNAmgcfqef8CvhPAY8ggbP` (cron `7 * * * *`) is enabled | **Inferred, not verified.** This session has no tool to read routine settings or run logs. | W3 claim PRs open at minute :08–:11 of the hour: #173 16:09, #171 07:30, #170 05:35, #168 22:10, #160 09:11, #151 03:09, #148 00:08, #147 23:09 (UTC). The newest, #173 (gb10), opened 16:09 UTC today. 27 distinct session IDs committed to main since 2026-10-11. |
| The routine produces merged results | Observed | 15 commits on main since 2026-10-10, nearly all W3 MOA-II fields (gb1–gb9, 0 survivors each). They were self-merged under the merge policy. |
| The routine reads the latest conclusions | **Partly.** | The prompt reads TASKS.md "Now", the DECISIONS headings and the newest CHANGELOG entry on main, plus `head -40` of CHANGELOG on each open agent PR branch. 14 non-draft `needs-human` PRs (#150–#172) carry conclusions that are **not on main**. |
| The routine's mission text is current | **No.** | `docs/cloud-routine-prompt.md` still states the D-023 cluster-lensing mission. Owner decisions since then (any dataset, 2026-10-08; System A/B; the goal rounds; S1; exploration units) appear only through TASKS.md. |

Concrete redo risk found:
- TASKS.md on main still lists S1 "Next: … TTE for short and faint bursts, an s ≠ 1 chain".
- Both were done and paused in #167 and #169, which are unmerged.
- A routine that orients on main could repeat them.

## 2. Rules that are missing or work against the goal

1. **WIP cap side effect.**
   - 14 needs-human PRs exceed the cap of 3 permanently, so the routine may only open draft claim PRs that it can
     self-merge.
   - In practice that means routine field sweeps (W3).
   - *Corrected after owner review:* high-risk work is not blocked outright. Draft claim PRs are exempt, and
     research that touches no open PR file stays allowed (recorded in an issue). The real effect is that such
     work cannot end in a new PR while the cap is exceeded, so it piles up or stops, and nothing lets an
     authorized unit continue across runs.
2. **No home on main for unmerged conclusions.** Paused or settled findings live only on PR branches, which feeds
   the redo risk above.
3. **No stage gates for multi-run units.** Units are sized to one ~40-minute run. Nothing lets a screened unit
   continue across runs with a budget, or makes it stop when a gate fails.
4. **No outcome ledger.** Progress is visible only as PRs and CHANGELOG entries. Method improvements, conditional
   theory, candidates and independent checks are not counted separately.
5. **The idle rule is implicit.** Rule 13 of the prompt ("stop instead of inventing low-value work") exists, but
   field sweeps always qualify as work. Nothing asks whether the next field still changes a limit materially:
   no field moves any cell by more than 19 %.
6. **Reopen rule.** Paused lines (S1 energy, the goal rounds, wh_farside) record their reopen inputs only in
   their PRs.

## 3. Proposed rules (to adopt via owner-approved edits)

**R1. Two WIP counters.**
- Agent PRs awaiting owner **review** (`needs-human`) are counted separately from agent PRs awaiting **merge**.
- The cap of 3 applies to the second counter.
- A separate owner-backlog cap: with ≥ 10 open needs-human PRs, at most **one new needs-human PR per week**, and
  only from a unit that passed a stage gate.
- Human merge rights are unchanged.

**R2. Pending-conclusions index on main.**
- One agent-mergeable file, `docs/state/pending.md`, lists each open needs-human PR in one line: unit, status
  (paused / settled / awaiting decision), reopen inputs, and a link.
- It claims no science beyond the PR's own words and touches no guarded file, so it is self-mergeable under the
  current policy.
- TASKS.md "Next" lines point to it instead of repeating stale items.

**R3. Stage gates for high-risk units.** Each unit has the stages below, and each stage states: the unknowns,
the budget cap, the pass result and the stop result. A unit may span runs and weeks within its cap.

| Stage | Purpose | Budget |
|---|---|---|
| S0 screen | Duplicate check against DECISIONS and open PRs, one model, one observable, an ordinary control | ≤ 2 h |
| S1 tool | Reuse check, then reproduce a known published case | ≤ 6 h |
| S2 discriminant | Prediction against an explicit ordinary model, with achievable precision | ≤ 8 h |
| S3 data | Smallest public data comparison, with validation data held out | ≤ 8 h |
| S4 independent check | A second code, an analytic limit or an independent data set | ≤ 4 h |

- "Not yet computed" or "needs a new tool" is **not** a stop. It triggers S1 if the work fits the cap.
- Time already spent is **never** a reason to continue past a failed gate.
- Reopening a paused line needs a new method or input named in its reopen list.

**R4. Outcome ledger** (`docs/state/ledger.md`): four columns only.

| Method improvements | Conditional theory results | Real candidates | Independent verifications |
|---|---|---|---|

- PR counts, test counts and compute hours are excluded.
- Seed from the record, for example:
  - method improvements: the M2 stretch grid (#167), the W3 chain fixes (D-068);
  - conditional theory: goal rounds 1–10, wh_farside;
  - candidates: none confirmed (n0153 inconclusive, #61);
  - independent verifications: none of a candidate. SkyLLH (D-080) is a benchmark reproduction.

**R5. Idle rule.**
- *(Adopted as a reference value only; see trial rule T4.)* A routine-search unit runs only if its expected effect on a limit cell is ≥ 10 % or it covers new parameter
  space.
- A high-risk unit runs only if it has an open stage with budget left.
- Otherwise the run stops after dispatch: no claim PR, no documentation-only edits.

## 4. Resource split (inside the current subscription; no paid compute)

| Share of routine run time | Use | Trade-off |
|---|---|---|
| **45 %** | Routine searches (W3 remaining large fields gb13/gb14, W5/W1 when new data) | W3 limits tighten more slowly; marginal fields first get the ≥ 10 % test |
| **30 %** | High-risk exploration: **at most two units in depth** | Fewer field PRs per week |
| **15 %** | Verification: independent checks, vetting, held-out data | Not optional, even when no candidate exists |
| **10 %** | Upkeep: owner reviews, CI, base merges, the pending index | Bounded so it cannot crowd out research |

Unused shares are **not** reassigned to new invented work. The run idles (R5).

## 5. Two high-risk directions, differing in mechanism and data

| Unit | Mechanism | Observable | Data | State |
|---|---|---|---|---|
| **X1 Far-side light at Sgr A*** (#172) | Light crossing a Simpson–Visser throat | Image-weighted centroid, light-curve peaks, Doppler weights of far-side hot spots | GRAVITY flare publications, EHT | S0 done (#172). Next: S1 tool |
| **X2 (slot, chosen by an S0 screen in week 2)** | Must differ from X1, e.g. a propagation effect over cosmological paths rather than a horizon-scale geometry | Time + spectrum + position jointly (not a single summary number) | Must be public and not already settled (check S1/S2/D1/E1) | Open; left empty if no candidate passes S0 |

## 6. First multi-week unit: X1, ray-traced far-side signatures

- **Question:** does any observable separate far-side SV hot spots from an explicit ordinary model (GRAVITY's
  same-side orbiting hot spot plus a measured static component) at precision GRAVITY reaches?
- **Total cap:** 20 h over four weeks.

| Week | Stage | Work and unknowns | Gate to continue | Stop if |
|---|---|---|---|---|
| 1 | S1 tool (≤ 6 h) | Reuse check: Gyoto (used for SV-like wormholes by Vincent and colleagues; licence and SV support to verify), else a 2D geodesic integrator for the SV metric (spherical symmetry keeps it small). Reproduce a published Schwarzschild hot-spot centroid loop and the far-side confinement of 2408.07350 | Known cases reproduced within 5 % | No tool reproduces a known case within the cap |
| 2 | S2 discriminant (≤ 8 h) | For i ≈ 155°, R = 9M and a ∈ {2.1, 3, 5, 7, 9}M: far-side vs ordinary centroid loop, flux modulation, peak count, image weights. Unknowns: static-flux offset, polarisation transport (only a flagged approximation) | At least one observable differs by more than GRAVITY's published per-epoch precision for a stated ordinary model | Every difference is below the precision or degenerate with the ordinary model's free parameters |
| 3 | S3 data (≤ 4 h) | Compare with the published per-flare light curves and astrometry: figures and any CDS tables of GRAVITY 2018 and 2023. Keep one flare unseen until the method is frozen | A prediction is frozen before the held-out flare is opened | Data insufficient: record the minimum reduced data needed |
| 4 | S4 check (≤ 2 h) | Analytic limits (b_c, the Kepler period) against the ray tracer; write-up | Agreement | Disagreement unexplained: the result is not reported as a finding |

The result counts in the ledger as a **conditional theory result** or a **method improvement**. It counts as a
candidate only if a held-out comparison favours the far-side model over the stated ordinary model and survives S4.

## 7. Plan versus running process

- **Adopted in the repository** (D-083):
  - `docs/state/trial.md`, `docs/state/pending.md` and `docs/state/units/X1.md`;
  - the research-cycle skill, CLAUDE.md, TASKS.md and `docs/cloud-routine-prompt.md`.
- **Not done from a session:** updating the routine's own prompt text. No tool reads or edits routines here. The
  routine reads the repository files above on every run, so the dispatch rules apply from its next run. The prompt
  copy (version 2026-10-12.1) must be pasted into the routine by the owner.
- **Acceptance** requires all three:
  - a routine run's CLAIM comment shows the version it used;
  - a run starts X1 S1;
  - a later run continues from X1's state file.

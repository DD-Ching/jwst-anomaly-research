# Exploration trial review (D-083, 2026-10-12 → 2026-11-09)

Written by a cloud routine run on the trial's last day (2026-11-09 03:13Z), as D-083 "Revisit if" asks.
**The decision (continue, adjust or stop) is the owner's.** This review changes no rule; until the owner rules,
runs keep applying `docs/state/trial.md` rules version trial-2026-10-13.2 (ASSUMPTION: the trial rules do not
lapse on their own; D-083 asks for a review, not an automatic stop).

## 1. Run time by category (derived from the weekly tally)

Source: every `| 2026-W… |` row of `docs/state/trial.md` (565 rows, 2026-W42 to the 02:08Z run of 2026-W46).
The hours are the runs' own estimates (the first number in each cell); they are not measured wall-clock time.

| ISO week | Rows | Routine | High-risk | Verification | Upkeep | Total h | Shares (R / H / V / U) |
|---|---|---|---|---|---|---|---|
| 2026-W42 | 177 | 36.6 | 27.2 | 23.9 | 25.1 | 112.7 | 32 / 24 / 21 / 22 % |
| 2026-W43 | 162 | 14.5 | 4.2 | 41.7 | 20.7 | 81.1 | 18 / 5 / 51 / 26 % |
| 2026-W44 | 118 | 1.9 | 0.0 | 7.2 | 7.6 | 16.7 | 11 / 0 / 43 / 46 % |
| 2026-W45 | 105 | 0.0 | 0.0 | 7.7 | 8.5 | 16.2 | 0 / 0 / 48 / 52 % |
| 2026-W46 (part) | 3 | 0.0 | 0.0 | 0.1 | 0.3 | 0.4 | — |
| **All** | 565 | 52.9 | 31.4 | 80.6 | 62.2 | 227.0 | 23 / 14 / 35 / 27 % |

Target split (T1): 45 / 30 / 15 / 10 %. Observed: the target was approached only in W42. From W44 on, high-risk
time was zero and routine-search time near zero; recorded work fell from about 113 h to about 16 h per week, and
most of it was verification of the `arcs` line plus upkeep. The low totals are consistent with T1's rule that unused
shares are not filled with invented work.

## 2. Units (from `docs/state/units/`, the Units table of `docs/state/trial.md` and TASKS.md)

| Unit | Status at trial end | Budget used | What blocks it |
|---|---|---|---|
| X1 far-side SV hot spots (Sgr A*) | closed 2026-10-13, gates passed as written; conditional theory results only | 3.7 / 20 h | reopen inputs in `X1.md` "Closure" |
| X2 SV-throat echoes | paused for owner review after S6a; writeup `X2_s6a_summary.md` | ≈ 13.9 / 16 h | owner's call on an expert review; polar SV equations, source excitation η |
| X3 W3 in Gaia epoch photometry | paused 2026-10-15; G1 held-out null, Γ₉₅ < 1 × 10⁻⁵ /star/yr (conditional) | ≈ 5 / 16 h | ESA DR4 epoch-astrometry documentation (2026-12-02) |
| X4 wave-optics GW lensing by negative mass | paused after S1k / exact-score round (2026-10-22); R1 undetermined | ≈ 10.8 / 18 h per the Units table (later S1g–S1k were owner-requested rounds) | owner: a precession parametrisation, GW data-analysis advice or the not-posted expert query |

S0 screens run: 5 (`S0-2026-10-13` … `S0-2026-10-16`). They authorized X2, X3 and X4. The two owner-check screens
(2026-10-15b, 2026-10-16) stopped without a unit.

## 3. Outcomes (the four ledger columns)

- **Method improvements:** M2 stretch grid (#167, provisional); W3 chain fixes (D-068); X3 frozen Gaia W3 screen with
  measured blind spots; W3 dimming coverage terms; macro-shear tool (D-106); Kepler W3 pilot (D-107); lens-model
  route checks (Sharon / Keeton / Lenstool / GLAFIC maps, Lenstool port D-109); the `arcs` anti line's frozen gates
  (D-024 addenda 8–20).
- **Conditional theory results:** X1 S2–S4, X2 S4–S6a, X4 S1–S1k (all simulated or derived; none is observational
  evidence, T5).
- **Real candidates:** none. Screen hits that stayed `inconclusive` after vetting: Sunrise n0153 (#61), Abell 370
  1565 (#199), Abell 2744 g09725 (#553), MACS0647 g85 (#556); Abell 2744 700.1 stayed `undetermined` (D-036). The
  `arcs` single-component anti excess did not survive vetting (D-024 addendum 20, P 0.440); one row, MACS0308 2157,
  stays open and needs spectroscopy or a JWST grism.
- **Independent verifications:** SkyLLH TXS 0506+056 benchmark (D-080, a reproduction, not a candidate check); X1
  tracer second routes; X2 Numerov vs DOP853 and the 1309.2448 reproduction; X4 PyCBC score reproduction (3.1e-6).

The outcome ledger table in `docs/state/trial.md` holds one aggregated row written early in the trial; the later
outcomes above are in CHANGELOG.md and DECISIONS.md but were never folded into it.

## 4. Did the trial's mechanisms work?

| Mechanism | Evidence | Assessment |
|---|---|---|
| Split WIP counters (T2) | needs-human PRs stopped blocking self-mergeable work; PR numbers ran from #174 to #806; awaiting merge at trial end: 1 (#751) | worked |
| Owner-review backlog rule | 27 open `needs-human` PRs at trial end (22 non-draft, 5 draft), oldest #61 (2026-10-11); 10 are goal rounds #152–#166 and 2 are draft physics rounds (#592, #619) | the backlog grew instead of shrinking; it is now the binding constraint |
| Per-unit state files (T3) | every run that touched a unit resumed from its file; no verified step was redone | worked |
| Weekly split (T1) | met only in W42; high-risk 0 % from W44 | not met; all four units paused on owner or external inputs, and T1 rightly did not invent work |
| `pending.md` | provisional conclusions of unmerged PRs stayed visible; no "done; do not redo" item was redone | worked |
| Idle rule | W46 runs stopped in ≈ 0.05 h each when nothing was unblocked | worked; each idle hourly run still costs ≈ 0.05 h of upkeep |

## 5. Recommendation (for the owner's decision)

**Adjust**, not stop:
1. Keep T2 (split counters), T3 (unit files), T5 (labels) and `pending.md` as standing rules; they worked.
2. Replace T1's fixed 45 / 30 / 15 / 10 split with "report the shares; no target", since the split was not
   reachable while all units waited on the owner.
3. Clear the owner-review backlog first. Lowest-effort route: close the 10 goal-round PRs #152–#166 as
   "recorded, not merged" (their conclusions stay in `pending.md`), then review the 6 candidate-labelled PRs
   (#61, #145, #199, #237, #553, #556; #145 and #237 are already explained as red noise and SN Atalanta).
4. Decide whether the S0 authorization (T7) continues. If yes, the next run can screen up to three new high-risk
   directions; if no, high-risk work waits for a new owner brief.
5. Reduce idle cost: with the queue blocked on dated inputs (CDS for #751, DR4 on 2026-12-02, Euclid DR1, HST
   18470 visits), consider a lower routine frequency until one of them arrives.

Stop conditions not met: no unit overran its budget without a gate, and no rule was broken.

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
| X2 SV-throat echoes vs published post-merger limits | `docs/state/units/X2.md` | `claude/x2-s6a-5d` | S1–S3 done for one event (GW150914; phenomenological model P1–P4). S4: matter-consistent axial equations (Yang & Lan; NED + phantom) give η_max 1.52 at the 0.45 s cell, outside the proxy bracket; S4b derived the strain/flux conversion. S5: the background is radially unstable (e-folding 9.9 M, 139 e-folds per echo spacing), so echo searches on it stop. S6 screen: rotation is not an established stabiliser. S6a: in 5D equal-J rotating Ellis wormholes the two unstable modes continue past their merger as an oscillating growing pair (verified at J/A 0.121 and 0.142; odd parity screened only; near-extremal unresolved; writeup drafted, no expert review yet; about 13.9 / 16 h). η_max 0.70–1.24 in the applicable a-range (a − 2M ≈ 2.2e-4 to 7.1e-4 M); η ~ 1 disfavoured under the scalar and EM proxies only, not proxy-independently. S1 check 4 (a − 2M = 1e-3) and S3 spacing (+4.1 %) fail the pre-written 1 % (reflection delay; timing sensitivity recorded). **Paused for owner review**; next if resumed: S2b multi-event (≤ 2 h). Missing input: polar SV perturbation equations, the source excitation η |
| X3 W3 in unselected Gaia epoch photometry (GAPS → DR4) | `docs/state/units/X3.md` | `claude/s0-2026-10-14` | **Paused 2026-10-15** (reopen when the DR4 epoch-astrometry fields and precision are confirmed, from documents or the release). Authorized by S0-2026-10-14 (A: PTA Shapiro sign test and B: FRB exotic lensing not qualified). G0 passed (post-hoc astrometric veto v1). **G1 done:** full screen calibrated on dev (W3 eff 0.32, PSPL 0, dev flags 0 / 1,995); held-out GAPS run once: 395,764 exposed, 196 flags, vetted, **no candidate**; Γ₉₅ < 1 × 10⁻⁵ /star/yr for photocentre-stable W3 (conditional; flags excluded by rules, RUWE figure a model estimate). **G2 prerequisite:** astrometry mostly repeats the ordinary blend effect; the lens image shift is measurable only for G ≲ 17, θ_E ≳ 3 mas, t_E ≈ 100 d (assumed σ); DR4 epoch-astrometry content unknown. Next: wait for DR4 documentation. ~5 / 16 h |
| X4 wave-optics GW lensing by negative-mass objects | `docs/state/units/X4.md` | `claude/s0-2026-10-15` | Authorized by S0-2026-10-15 (B: VASCO and C: quantum-processor wormhole not qualified). S0 minimal step: closed-form repulsive-lens F(w) (low-pass), separable from a positive point mass at observed SNR ρ ≈ 3.1–16.5 for M_Lz ≥ 300 M☉, y ≤ 1, intrinsic parameters fixed (strongly demagnified cases). **S1 done:** 18 / 36 simulated cases keep a residual after free aligned-spin binary (q ≤ 15) and positive-lens fits (frozen refinement; IMRPhenomD agrees); extreme q absorbs part. **S1b done:** R2, R3 no longer pass the screening threshold against precessing, higher-mode binaries (ρ_req > 30). **S1c done (R1):** residual kept against validated aligned-spin fits (ρ_req 5.8, bound-conditional, provisional); undetermined with XPHM (best ρ_req 9.2). **S1d:** the SciPy DE reference search fails its injection checks; R1 still undetermined. **S1e:** joint search: U reliable, P 3/4; tool still fails; R1 undetermined. **S1f:** validated U search finds no counter-example (FF 0.9135, ρ_req 7.4). No pilot proposed. ~10.8 / 18 h |

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
| 2026-W42 | ≈ 0.45 (Sunrise image positions in the Scofield+2025 maps) | — | — | 0.15 (orientation, claim, PR) | cloud routine 01XzRkq7, 14:09–14:50Z |
| 2026-W42 | ≈ 0.4 (Sunrise deflection unit, system 2 redshift scan) | — | — | 0.15 (orientation, #192 base merge, claim, PR) | cloud routine 016ZY9Gt, 15:08–15:50Z |
| 2026-W42 | ≈ 0.2 (Sunrise system 2 in the four RELICS models) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01UitTtD, 16:09–16:35Z |
| 2026-W42 | — | 0.5 (X2 S4: literature check, axial NED comparison) | 0.2 (reproduction of Yang–Lan QNMs) | 0.1 (state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.45 (`arcs` on Abell 2744 and El Gordo, shape-noise null, chain spread) | — | — | 0.1 (orientation, claim, PR) | cloud routine 0111axHH, 17:08–17:50Z |
| 2026-W42 | ≈ 0.2 (`arcs` on MACS1149 and MACS0416, cutouts of every anti) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01MAK8Jd, 18:08–18:20Z (#197) |
| 2026-W42 | ≈ 0.25 (`arcs` on Abell 370; five-field shape-noise sensitivity) | — | — | 0.05 (claim, PR) | cloud routine 01MAK8Jd, 18:20–18:45Z (#198) |
| 2026-W42 | ≈ 0.1 (`arcs` shape S/N cut on the stored five-field tables) | — | — | 0.1 (PR, CI) | cloud routine 01SirLB6, 21:18–21:40Z |
| 2026-W42 | — | 0.6 (X2 S4b: strain/flux derivation, single-cell reruns) | 0.2 (strain-integral check, coupling-off diagnostic) | 0.2 (S4c literature, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.3 (`arcs_bands.py`: F444W vs F200W arc axes, five fields) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01VwkNGL, 22:08–22:45Z |
| 2026-W42 | ≈ 0.25 (`arcs_curvature.py`: arc curvature vs −α, five fields) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01GW8SqA, 23:08–23:24Z (#203) |
| 2026-W42 | ≈ 0.2 (`arcs_curvature_model.py`: forward-modelled curvature, five fields) | — | — | 0.05 (PR) | cloud routine 01GW8SqA, 23:24–23:45Z |
| 2026-W42 | — | 1.4 (X2 S5: radial equations with NED, growing-mode search) | 0.2 (Ellis–Bronnikov reference, convergence) | 0.1 (state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 0.7 (X2 S6: rotation-stability literature and equation check) | — | 0.1 (S5 scope fix, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.5 (`arcs_counter.py`: counter-images of uncatalogued arcs, five fields, control) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01QyGL8z, 00:08–00:48Z |
| 2026-W42 | — | — | ≈ 0.4 (vet Abell 370 3881: CATS system 7, second-model images; exclusion fix, five-field reruns and controls) | 0.15 (orientation, takeover, base merge, PR) | cloud routine 01AQVssd, 01:08–01:50Z |
| 2026-W42 | ≈ 0.25 (SMACS system 17 halo-subtracted counter-image) | — | ≈ 0.15 (system 17 vetting: z scan, RELICS models, literature) | 0.1 (orientation, claim, PR) | cloud routine 014mrZNP, 02:08–02:45Z |
| 2026-W42 | ≈ 0.3 (SMACS systems 11, 16, 26: sibling references in three models, two bands; cutouts) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01LDvpj4, 03:07–03:40Z |
| 2026-W42 | — | 1.8 (X2 S6a: 5D derivation, background, mode tracking) | 0.9 (static/merger gates, shooting cross-check, constraint residuals, convergence) | 0.6 (review and fixes, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | — | 0.3 (X2 S6a writeup: prior-work search, shooting audit) | 0.5 (summary, reproduction package, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.3 (`--halo-filter-arcsec` on Abell 2744, MACS1149, MACS0416; cutouts of 16 confused rows) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01WHQem9, 04:08–04:40Z |
| 2026-W42 | ≈ 0.4 (`member_sub.py`: Sérsic member subtraction, three fields, controls) | — | — | 0.1 (orientation, claim, PR) | cloud routine 015QCMes, 05:08–05:50Z |
| 2026-W42 | ≈ 0.5 (`parity_test.py`: observational parity, synthetics, Abell 2744 run; not validated on known systems) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01GMu4HZ, 06:08–06:50Z |
| 2026-W42 | — | — | ≈ 0.3 (parity acceptance test: segmentation, resolution gate, start convergence; failed) | 0.15 (orientation, claim, PR) | cloud routine 01NCPwCw, 07:08–07:45Z |
| 2026-W42 | — | — | ≈ 1.2 (PSF-commuting parity fit, empirical PSF, acceptance rerun; failed, tool paused) | 0.15 (orientation, claim, PR) | cloud routine 01VnLVRw, 08:09–09:45Z |
| 2026-W42 | ≈ 0.4 (MACS0416 `halo_phot` + `fluxratios --halo-phot`, CANUCS model; cutouts) | — | — | 0.1 (orientation, claim, PR) | cloud routine 019kpPJS, 09:08–09:45Z |
| 2026-W42 | — | 1.4 (S0-2026-10-14 screen: PTA reach estimate, literature/data check; X3 G0 harness and gate runs) | 1.0 (X3 G0 positive controls and false-alarm calibration on two dev samples) | 0.3 (orientation, claim, DataLink fix, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.35 (`sunrise_counter.py`: Sunrise counter-images in four RELICS models; control; cutouts) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01Dfw2ou, 11:08–11:45Z |
| 2026-W42 | ≈ 0.35 (UNCOVER v2.0 maps: route, validation on both image lists; gate fails) | — | — | 0.15 (orientation, claim, PR) | cloud routine 01Mkg8y2, 12:09–12:55Z |
| 2026-W42 | ≈ 0.45 (UNCOVER v2.0 gate fix; 4.2c second-model photometry and injections) | — | — | 0.1 (orientation, takeover, base merge, PR) | cloud routine 012BA5rZ, 13:08–13:50Z |
| 2026-W42 | ≈ 0.4 (D-041 cross-model rule on every Abell 2744 forced row, UNCOVER v2.0; cutouts) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01DJ1D1o, 14:09–14:45Z |
| 2026-W42 | ≈ 0.3 (Abell 2744 sibling reference fluxes for both-model extras, `a2744_xmodel_refs.py`) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01GcHRcx, 15:08–15:40Z |
| 2026-W42 | ≈ 0.25 (`--plane` CLI; MACS0416 8757 two-plane χ² check) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01XYMGpL, 16:08–16:45Z |
| 2026-W42 | ≈ 0.35 (`offplane_potentials.py`: MACS0416 CANUCS off-plane potentials, O5 σ scan, `images --plane`) | — | — | 0.1 (orientation, claim, PR) | cloud routine 014CmXZZ, 17:08–17:45Z |
| 2026-W42 | ≈ 0.35 (off-plane potentials with MUSE z_spec: Abell 2744, El Gordo z 0.63 group σ scan, `images --plane`) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01FdC23M, 18:08–18:45Z |
| 2026-W42 | ≈ 0.35 (off-plane potentials: Abell 370 MUSE vs CANUCS z_spec, 180 σ controls; SMACS frame offset, `images --plane`) | — | — | 0.1 (orientation, claim, PR) | cloud routine 0181RxsQ, 19:08–19:45Z |
| 2026-W42 | — | — | ≈ 0.4 (vet five Abell 2744 follow-ups: published identities, Bergamini+23 counter-images, 264 knot colours; cutouts) | 0.1 (orientation, claim, PR) | cloud routine 0183gAYL, 20:08–20:50Z |
| 2026-W42 | — | — | ≈ 0.6 (vet El Gordo 2828/1601 and Sunrise 1756/1907/1578: spikes, DJA, persistence, proper motion, RELICS and Caminha+23 counter-images; cutouts) | 0.1 (orientation, claim, PR) | cloud routine 01EwxqDZ, 21:08–21:50Z |
| 2026-W42 | — | 0.3 (X3 G1: bulk-file route, sample counts) | 2.2 (X3 G1: full-pipeline calibration, references, selection bound, held-out run and vetting) | 0.3 (state files, review, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 0.2 (X3 G2 prerequisite: DR4 content check) | 0.5 (X3 G2 prerequisite: photocentre model, identifiability cases) | 0.2 (G1 corrections, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.3 (`m1149_refs.py`: MACS1149 `no_reference` rows, sibling references, Sharon v4cor; stamps) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01CLbKJM, 22:08–22:55Z |
| 2026-W42 | — | 1.5 (S0-2026-10-15: three directions; GW repulsive-lens closed form and distinguishability) | — | 0.2 (state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.3 (`m1149_sys16.py`: MACS1149 system 16 z scan in two models, far-image photometry; stamps) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01Pnu5BB, 23:08–23:45Z |
| 2026-W42 | ≈ 0.3 (`transient_grid.py`: Sunrise o010/o120 forced-photometry grid on every catalogued position, injections; contact sheets) | — | — | 0.1 (orientation, PR) | cloud routine 01QaBwbk, 00:08–00:45Z |
| 2026-W42 | — | 3.3 (X4 S1: IMR fitting tool, 36 cases, frozen refinement, wide-q and waveform cross-checks) | — | 0.3 (review fixes, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 1.4 (X4 S1b: multi-start fits, fit-error injections, IMRPhenomXPHM fits and continuations) | — | 0.2 (state files, review, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 1.9 (X4 S1c: unified fitter, lens grid, validation, R1 case fits in two families) | — | 0.1 (plan, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 1.0 (X4 S1d: reuse check, DE fitter, frozen validation) | — | 0.3 (plan, state files, review, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 1.4 (X4 S1e: diagnostics, joint DE, frozen validation) | — | 0.4 (record corrections, plan, state files, review, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 0.45 (X4 S1f: unlensed R1 fits, independent verification) | — | 0.15 (gate record, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 0.5 (S0-2026-10-15b Kapitza screen: literature check, derivation, thin-shell toy) | — | 0.2 (record, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | 0.3 (S0-2026-10-15b re-check: background existence, withdrawn claims) | — | 0.1 (record, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.2 (MACS0416 o211/o212 grid; contact sheets) | — | ≈ 0.2 (validation on Yan+2023 transients, grid fixes, persistence test) | 0.1 (orientation, claim, PR) | cloud routine 01AGqJkf, 03:08–03:45Z |
| 2026-W42 | ≈ 0.1 (MACS0416 o211/o213 grid; contact sheets, persistence test) | — | — | 0.05 (PR) | cloud routine 01AGqJkf, 03:30–03:45Z |
| 2026-W42 | ≈ 0.35 (`arc_diff.py`: MACS0416 o211/o212 arc differencing, validation on Yan+2023, three fix iterations; contact sheets) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01JZas31, 04:08–04:45Z |
| 2026-W42 | ≈ 0.4 (`arc_diff.py` pixel kernels: five PSF-match variants on o211/o212 and o211/o213, validation on Yan+2023; contact sheets) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01VwokZ2, 05:08–05:50Z |
| 2026-W42 | — | — | ≈ 0.4 (vet MACS0416 h108: five epoch pairs incl. CANUCS o004, `arc_diff_vet.py`, validation on Yan+2023; contact sheets) | 0.1 (orientation, claim, PR) | cloud routine 01VE5srt, 06:08–06:50Z |
| 2026-W42 | ≈ 0.4 (MACS0416 arc tiling: 16 / 28 tiles o211/o212; contact sheets) | — | ≈ 0.15 (`arc_diff_morph.py` validation on Yan+2023) | 0.1 (orientation, claim, PR) | cloud routine 012PG6nH, 07:08–07:48Z |
| 2026-W42 | ≈ 0.55 (MACS0416 arc tiling: remaining o211/o212 tiles with `--maps`; stalled by thread oversubscription, restarted) | — | — | 0.1 (orientation, claim) | cloud routine 016oUw9Y, 08:08–08:50Z |
| 2026-W42 | ≈ 0.35 (Abell 2744 DD 2756 grid, pair choice, `i2d_region.py`; contact sheets) | — | ≈ 0.1 (persistence test, UNCOVER third-epoch stamps) | 0.1 (orientation, claim, PR) | cloud routine 01HeG8Gq, 09:09–09:45Z |
| 2026-W42 | ≈ 0.45 (MACS0416 arc tiling: 9 more o211/o212 tiles; dipole verdict; core-peak check; contact sheets) | — | — | 0.15 (orientation, base merge, claim) | cloud routine 01QcfDQb, 10:09–10:50Z |
| 2026-W42 | ≈ 0.2 (Abell 2744 UNCOVER o001 vs 2756 o002 grid; contact sheets) | — | ≈ 0.15 (persistence test on 5 visits, dipole aperture checks, g11 core check) | 0.1 (orientation, claim, PR) | cloud routine 01JkMNZG, 11:08–11:45Z |
| 2026-W42 | — | — | ≈ 0.4 (MACS0416 t+2+2 h002/h003 persistence, h063 with CANUCS o004; `_cal` and difference contact sheets) | 0.15 (orientation, base merge, claim, PR) | cloud routine 01UyLcH7, 12:09–12:45Z |
| 2026-W42 | — | — | ≈ 0.5 (MACS0416 morphology veto: whole-chain injections on 11 tiles, control, rosette lobe-ratio fix, t−1−4 ring passes) | 0.1 (orientation, claim, PR) | cloud routine 01UfNGSL, 14:08–14:50Z |
| 2026-W42 | ≈ 0.3 (MACS0416 t−1−4 and 6 unrecorded tiles re-run with maps; contact sheets) | — | ≈ 0.2 (ring and sign_flip verdicts: injections at 27.5–29.5, Yan+2023 re-validation) | 0.1 (orientation, claim, PR) | cloud routine 0148wS49, 15:08–15:45Z |
| 2026-W42 | — | — | ≈ 0.4 (MACS0416 9 open arc_diff passes vs c004/o213; 7 ordinary, h044/h030 open) | 0.1 (orientation, claim, PR) | cloud routine 019UYvw9, 16:08–16:45Z |
| 2026-W42 | — | — | ≈ 0.4 (MACS0416 h044/h030 knot test, o211→o212 west box) | 0.1 (orientation, claim, PR) | cloud routine 01CGPk1D, 17:08–17:45Z |
| 2026-W42 | — | — | ≈ 0.45 (MACS0416 h044/h030 knot test on o211→c004, o211→o213 and the o212→o213 control; `arc_diff_knots.py`) | 0.1 (orientation, claim, PR) | cloud routine 018oYfhm, 18:08–18:50Z |
| 2026-W42 | — | — | ≈ 0.4 (MACS0416 h044/h030 per-dither `_cal` photometry and centroids, o211/o212/o213/c004; `cal_dither_phot.py`) | 0.15 (orientation, claim, review, PR) | cloud routine 01GzDwPS, 19:08–19:50Z |
| 2026-W42 | ≈ 0.6 (MACS0416 o211→o213 arc tiling, 23 / 28 tiles, morphology, contact sheets, persistence) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01G17H3Y, 20:08–20:50Z |
| 2026-W42 | ≈ 0.35 (MACS0416 o211→o213 last 5 tiles, h002 persistence, h014 per-dither photometry) | — | — | 0.1 (orientation, claim, review, PR) | cloud routine 01JkQLMW, 21:08–21:40Z |
| 2026-W42 | ≈ 0.25 (Sunrise o010/o120, o010/o052 and SMACS/VENUS grids re-run under origin recentring and mixed-sign hits; contact sheets) | — | — | 0.1 (orientation, claim, PR) | cloud routine 018nALFM, 22:08–22:35Z |
| 2026-W42 | ≈ 0.3 (Abell 2744 UNCOVER vs ALT F356W grid, PLATE-23a check, ALT o001 control, contact sheets) | — | ≈ 0.1 (g1/g2 multi-band and centroid checks) | 0.1 (orientation, claim, PR) | cloud routine 01LF3BFs, 23:08–23:35Z |
| 2026-W42 | — | — | — | 0.3 (process: merge=union check, operations notes, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | — | — | 0.1 (X3 DR4 reopen check: official DR4 documents only) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.7 (Abell 2744 UNCOVER vs ALT F356W `arc_diff` tiling, 33 tiles, population and ALT o001 tests, 76 survivors looked at) | — | — | 0.1 (orientation, claim) | cloud routine 01A8yfJB, 00:08–01:05Z |
| 2026-W42 | — | — | ≈ 0.4 (vet the 6 tiling priority hits: `hit_lightcurve.py` forced photometry in 46–64 covering mosaics, 53-source control calibration; stamps looked at) | 0.15 (orientation, takeover, PR) | cloud routine 01P4Y7VH, 01:08–01:50Z |
| 2026-W42 | — | 0.4 (S0-2026-10-16: unstable-wormhole message window, literature and code check) | — | 0.1 (record, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | — | ≈ 0.25 (Abell 2744 h182 arc-pixel controls: `pick-arc`, 78 controls in 64 mosaics; stamps looked at) | 0.1 (orientation, claim, PR) | cloud routine 01VqBbgX, 02:08–02:40Z |
| 2026-W42 | — | — | ≈ 0.4 (Abell 2744 `arc_diff` whole-chain injections, 4 tiles × 2 pairs; peer-step reconstruction; vetoed-star sheet looked at) | 0.1 (orientation, claim, PR) | cloud routine 01CgfDzc, 03:07–03:35Z |
| 2026-W42 | ≈ 0.4 (Abell 2744 `arc_diff` extended-light cut: injections, real passes, sheets; not adopted) | — | — | 0.1 (orientation, claim, PR) | cloud routine 012QTtJg, 04:08–04:50Z |
| 2026-W42 | ≈ 0.3 (Abell 2744 `arc_diff` stripe/spike/offset-core vetoes: injections, real passes, sheets; not adopted) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01KUpQtW, 05:08–05:35Z |
| 2026-W42 | — | — | ≈ 0.3 (PEARLS 1176 o221 registration vs UNCOVER: catalogue offsets, held-out check, image WCS check, h182 re-calibration) | 0.1 (plan, claim, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | — | — | 0.7 (S6a archival brief, PDF, README entry) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.8 (Sunrise arc `arc_diff`: o052 pair, sheet; o120 pair stopped) | — | — | 0.2 (orientation, claim, PR) | cloud routine 01FjRPtm, 06:08–07:25Z |
| 2026-W42 | — | 0.4 (Abell 2744 orphan-pair dark-deflector screen: frozen rules, alignment, injections, top-15 sheet) | — | 0.1 (orientation, claim, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | — | — | 0.3 (process: "Before writing results" rules; orphan vetting paused) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.15 (`arcs` on Abell S1063) | — | ≈ 0.2 (S1063 anti: F444W, faint split, ICL-subtracted shapes, three team models; sheets looked at) | 0.1 (orientation, claim) | cloud routine 01P8KV8A, 07:08–07:35Z |
| 2026-W42 | — | — | ≈ 0.3 (S1063 1316 / 2184: measurement-noise injections, HST axes, cutouts) | — | cloud routine trig01PNAm, 08:12–08:45Z |
| 2026-W42 | — | — | 0.6 (Abell 2744 orphan per-pair vetting: E2 controls, VizieR E5/E6, #7 check, classification) | 0.1 (state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | ≈ 0.3 (SMACS 0723 `arc_diff`: two pairs, Atalanta control, h061 re-box, h000 light curve; sheets looked at) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01De8CTn, 09:08–09:40Z |
| 2026-W42 | ≈ 0.3 (Sunrise arc_diff o010→o120 LW: frozen rules, run, morphology, distances, sheet) | — | — | 0.1 (claim, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | — | ≈ 0.4 (S1063 1557 deblend with own-flux gate, third-component and clump-orientation checks; sheet looked at) | 0.1 (orientation, takeover, merge, PR) | cloud routine trig01PNAm-1210, 12:10–12:50Z |
| 2026-W42 | — | — | ≈ 0.3 (S1063 clump-pair census of 24 strong-shear rows: plan, own-flux gates, clump-aware sample test; sheet looked at) | 0.1 (orientation, takeover, PR) | cloud routine trig01PNAm-1316, 13:17–13:45Z |
| 2026-W42 | — | — | ≈ 0.4 (full-chain known-case check: arc mask, SMACS / MACS0416 / Sunrise known cases, SMACS in-mask injections; mask figures looked at) | 0.1 (claim, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | — | ≈ 0.4 (review-list chain: position as label, known cases, whole-box injections, burden on 5 stored pairs, sheets looked at) | 0.1 (claim, state files, PR) | interactive session 01AzKMDv |
| 2026-W42 | — | — | ≈ 0.6 (S1063 1316/1557 public-spectroscopy feasibility: VizieR, ESO ObsCore, SODA cutouts, wording corrections) | 0.1 (PR) | interactive session 01P8KV8A |
| 2026-W42 | — | — | ≈ 0.3 (no-injection control at vetoed SMACS injection positions; sheet looked at) | 0.1 (orientation, claim, state files, PR) | cloud routine 01LuctU7, 14:08–14:40Z |
| 2026-W42 | ≈ 0.4 (`arcs` on MACS0417: RELICS v1 and CANUCS v1 maps, cutouts looked at) | — | — | 0.1 (orientation, claim, PR) | cloud routine 012exdqt, 16:08–16:45Z |
| 2026-W42 | ≈ 0.15 (`arcs` on MACS0647, CLASH maps, two epochs; no result) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01KYMo38, 17:08–17:25Z |
| 2026-W42 | ≈ 0.2 (`arcs` on MACS1423, CANUCS DR1 maps; cutouts looked at) | — | ≈ 0.2 (MACS1423 anti: ICL-subtracted shapes, measurement-noise injections) | 0.1 (orientation, claim, PR) | cloud routine 01Brfkw8, 18:08–18:45Z |
| 2026-W42 | — | — | 0.3 (MACS1423 `arcs` vetting step 2: clump-pair census, P 0.259, explained; #282) | 0.1 (orientation, claim, state files) | cloud routine 01HysHAx, 19:08–19:40Z |
| 2026-W42 | ≈ 0.3 (`arcs` on RXJ2129, Lenstool and GLAFIC maps by zip range reads; no result) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01Lvrs8W, 20:08–20:40Z |
| 2026-W42 | ≈ 0.35 (`arcs` on RXJ2129 with the DJA catalogue: PA convention gate, run, sheet; no result) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01ExftNC, 21:08–21:45Z |
| 2026-W42 | ≈ 1.0 (El Gordo arc_diff one-box pilot: screen, full-chain injections, 4 entries vetted) | — | ≈ 0.4 (spike, on-spike and bright-source controls) | 0.2 (orientation, claim, state files) | interactive session 01P8KV8A |
| 2026-W42 | — | — | ≈ 0.15 (`arcs` PA convention check: pa_sep vs pa_obs on stored census rows) | 0.1 (orientation, claim, PR) | cloud routine 01UdRNer, 22:08–22:35Z |
| 2026-W42 | ≈ 0.4 (El Gordo boxes 2–3 `arc_diff` screen and injections) | — | ≈ 0.9 (`arc_diff` fading calibration with pre-matching injections; registration fix; review-driven geometry fix; held-out reruns; #286 re-screen) | 0.3 (plan addenda, speed-up comparison, state files, review, PR) | cloud routine 01P8KV8A, 23:32–01:30Z (#288) |
| 2026-W42 | — | — | ≈ 0.25 (MACS1423 clump-pair census re-run with the fixed #273 code; sheet looked at) | 0.1 (orientation, takeover, PR) | cloud routine 01WSfhqY, 00:08–00:45Z |
| 2026-W42 | ≈ 0.4 (`arcs` on WHL0137: plan, GLAFIC model pin, run, sheet) | — | — | 0.1 (orientation, claim, PR) | cloud routine 016CxjVy, 01:08–01:45Z |
| 2026-W42 | — | — | ≈ 0.4 (kernel_loo: kernel-level and full-chain B2 equivalence checks) | ≈ 0.4 (shared-design stamp kernels, D-094; plan, state files, review, PR) | interactive session 01P8KV8A, 02:23–03:30Z (#291) |
| 2026-W42 | — | — | — | 0.3 (MOA-II manifest rows and retrieval dates; orientation, PR) | cloud routine 01RSWD6E, 02:08–02:45Z |
| 2026-W42 | — | — | ≈ 0.2 (SMACS ERO → VENUS `arc_diff` tile-median registration re-screen; sheet looked at) | 0.1 (orientation, claim, PR) | cloud routine 017PWozQ, 03:08–03:30Z |
| 2026-W42 | — | — | ≈ 0.25 (SFFT vs arc_diff on B1/B2 injections) | ≈ 0.25 (SFFT integration check, VENUS pair enumeration, plan, state files, PR) | interactive session 01P8KV8A, 04:02–04:45Z (#293) |
| 2026-W42 | ≈ 0.3 (W1 shear screen on Abell S1063, injections, five-field re-combination; peak sheet looked at) | — | — | 0.1 (orientation, claim, PR) | cloud routine 018YnfpM, 04:08–04:40Z |
| 2026-W42 | ≈ 0.3 (W1 shear aperture-coverage cut: plan, code, five-field c = 0 / 0.5 injections; not adopted) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01H2i3uj, 05:07–05:35Z |
| 2026-W42 | ≈ 0.5 (VENUS batch: 20 boxes, per-box calibration, 5 reviews) | — | — | ≈ 0.3 (#293 wording fixes, plan, box selection, state files, PR) | interactive session 01P8KV8A, 06:01–07:00Z (#296) |
| 2026-W42 | ≈ 0.4 (W1 Abell 2744 B-mode floor diagnostic: sheet looked at, contributors, rotation null, bright-neighbour screen) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01QRPt8P, 06:08–06:40Z |
| 2026-W42 | ≈ 0.4 (W1 Abell 2744 PSF anisotropy from stars: plan, `psf_chi` option, screens, injections later withdrawn after review; map looked at) | — | ≈ 0.1 (PLCK G165 h002 record correction: NS_0481 3.77″, literature subagent; #299) | 0.1 (orientation, claim, PR) | cloud routine 019HJb2R, 07:08–07:45Z |
| 2026-W42 | ≈ 0.35 (W1 STPSF model PSF vs stars and the Abell 2744 floor: plan, PSF grid, area-matched stars, screens; PSF figure looked at) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01TXkzEW, 08:08–08:40Z |
| 2026-W42 | ≈ 0.2 (W1 smoothed STPSF vs star χ: plan, five variants, S3 fails) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01VZp3Mr, 09:08–09:30Z |
| 2026-W42 | ≈ 0.4 (W1 shear on SMACS 0723 + El Gordo: plan, screens, injections, A2744 reproduction, sheet looked at) | — | — | 0.15 (orientation, claim, state files, PR) | cloud routine fb53b841, 10:08–10:55Z |
| 2026-W42 | ≈ 0.3 (W1 shear on MACS0417, MACS1423, WHL0137: plan, screens, injections, ten-field combination; sheets looked at) | — | ≈ 0.1 (MACS1423 hit: frozen streak/star test) | 0.1 (orientation, claim, state files, PR) | cloud routine 01Td5iwv, 11:08–11:35Z |
| 2026-W42 | ≈ 0.2 (W1 shear on RXJ2129 and MACS0647: plan, screens, injections, eleven-field combination; sheet looked at) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01E49N6v, 12:08–12:40Z |
| 2026-W42 | ≈ 0.4 (W1 streak veto: plan, findsat_mrt masks on 12 fields, veto screens and injections, eleven-field limit; sheets looked at) | — | — | 0.1 (orientation, claim, input staging, state files, PR) | cloud routine 01H5UuYw, 13:08–13:45Z |
| 2026-W42 | ≈ 0.3 (radial screen on the five W1 shear-only fields: plan, S1063 reproduction, runs, sheet looked at) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 015p2WtK, 15:08–15:45Z |
| 2026-W42 | ≈ 0.3 (radial screen with the W1 streak veto: plan, mask regeneration, gate K, twelve-field diagnostic; sheet looked at) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01CgPQ37, 16:08–16:50Z |
| 2026-W42 | ≈ 0.2 (radial centre-level axial rule: plan, MACS0647 reproduction, gate K failed, post-hoc collinearity; sheet looked at) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01Mk3R1x, 17:08–17:30Z |
| 2026-W42 | ≈ 0.2 (W1 shear stacked at radial centres: plan, four fields, two post hoc vet steps) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 0129VzsK, 18:08–18:28Z |
| 2026-W42 | ≈ 0.1 (W1 radial-centre stack with a second MACS0417 model) | — | — | 0.05 (state files, PR) | cloud routine 0129VzsK, 18:28–18:45Z |
| 2026-W42 | ≈ 0.3 (W1 injections at the `radial` centres through the stack; frozen rule degenerate, post hoc noise-aware limit) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01VzSwtX, 19:08–19:28Z |
| 2026-W42 | ≈ 0.1 (W1 stack injections: mass points 5 × 10¹¹ and 1 × 10¹² M☉) | — | — | 0.05 (plan, state files, PR) | cloud routine 01VzSwtX, 19:28–19:40Z |
| 2026-W42 | ≈ 0.3 (W5 counts at `radial` centres: plan, script, Gaia via CfA VizieR mirror, planned run, 2″-neighbour diagnostic) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01U6Y4Kd, 20:08–20:50Z |
| 2026-W42 | ≈ 0.25 (W1 shear on MACS0717: ASTRODEEP-FF photo-z adapter, streak mask, screen, injections, S_max sheet) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01HW32mP, 21:08–21:35Z |
| 2026-W42 | — | — | — | 0.3 (#319 review follow-ups: column-value hash pin, masked-flag tests, handoff) | cloud routine, 22:08–22:30Z |
| 2026-W42 | ≈ 0.1 (radial screen on MACS0717 with ASTRODEEP-FF photo-z: plan, run, sheet) | — | — | 0.05 (claim, state files, PR) | cloud routine, 22:27–22:45Z |
| 2026-W42 | ≈ 0.45 (W1 volume density over lens redshift: plan, `--lens-z`, gate Z, 55-job grid, volume integral) | — | — | 0.1 (orientation, claim, staging, state files, PR) | cloud routine 01MxAyxM, 23:08–23:55Z |
| 2026-W42 | ≈ 0.35 (W1 shear in the five CANUCS NCF: plan, gate R, streak masks, 25-job grid, flag vetting, volume combine) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01JAwoGc, 00:08–00:50Z (2026-10-18, ISO W42) |
| 2026-W42 | ≈ 0.25 (W1 faint-streak veto: plan, `--mrt-threshold`, two variants on two NCF, K1n failed; sheets looked at) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01DA23en, 01:08–01:40Z (2026-10-18) |
| 2026-W42 | ≈ 0.4 (W1 row/column-segment veto: plan, `segment_stripes`, NCF staging, two variants × five NCF masks and screens, cost check failed; sheets looked at) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 0188maj8, 02:08–02:50Z (2026-10-18) |
| 2026-W42 | ≈ 0.4 (W1 shear in GOODS-N: plan, footprint choice, staging, gate P, streak masks, 36-job grid, flag vetting and sheet, volume combine) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 017BFTDU, 03:08–03:45Z (2026-10-18) |
| 2026-W42 | ≈ 0.6 (W1 blank-field R image simulations, `results/w1_r_metacal/`) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01SnbR8J, 05:08–05:50Z (2026-10-18) |
| 2026-W42 | ≈ 0.75 (W1 shear in CEERS: plan, footprint script, staging, gate P, streak masks, 10 screens, two flags vetted with sheet, 48-job grid, volume combine) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01ASLQCx, 04:08–05:00Z (2026-10-18) |
| 2026-W42 | ≈ 0.4 (W1 shear in PRIMER: frozen plan, mosaic tiling, catalogue pins) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01RnsypG, 06:08–06:45Z |
| 2026-W42 | ≈ 0.5 (W1 shear in PRIMER run 4: last 23 injection jobs, R sensitivity, forty-field limit, RESULT) | — | ≈ 0.05 (x3y0 visit check) | 0.15 (orientation, takeover, base merge, PR) | cloud routine 018J4abj, 09:08–09:50Z |
| 2026-W42 | — | — | ≈ 0.35 (W1 PRIMER x3y0 `_cal` exposure check, `w1_cal_trail.py`; bright-star artefact, post hoc) | 0.1 (orientation, claim, PR) | cloud routine 0145mJQ8, 10:08–10:45Z |
| 2026-W42 | ≈ 0.2 (W1 lower-mass filter pilot: plan, staging, 6 injection jobs, D-103) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01Jc3ZxK, 11:08–11:35Z |
| 2026-W42 | — | — | ≈ 0.25 (W1 real-image metacal R: plan, tool, simulated gate, six GOODS-N tiles, sheet looked at) | 0.1 (orientation, claim, PR) | cloud routine 01SS72dG, 12:08–12:45Z |
| 2026-W42 | — | — | ≈ 0.2 (W1 metacal R v2: frozen plan, mask/floor/unique/tile gate, simulated gate, 22 tiles, sheet looked at) | 0.1 (orientation, claim, PR) | cloud routine 01E7j7xV, 13:08–13:35Z |
| 2026-W42 | — | — | ≈ 0.3 (W1 metacal R v3: frozen plan, ePSF tool, gate E failed, post hoc simulated + real diagnostic, sheet looked at) | 0.1 (orientation, claim, PR) | cloud routine 01NtqjLe, 14:08–14:45Z |
| 2026-W42 | — | — | — | 0.3 (orientation; queue check: every open line waits on the owner or is optional; `pending.md` sync for #199, #237, #279, #273; PR) | cloud routine 01RKeHhJ, 15:08–15:40Z |
| 2026-W42 | ≈ 0.25 (S2 DR10 z < 22: plan, `--survey dr10`, gate Q, fetch, three fits; control gate failed) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01JBqic5, 16:08–16:40Z |
| 2026-W42 | ≈ 0.3 (`arcs` on RXC J0600: map survey, plan, two RELICS models, sheets; does not replicate) | — | — | 0.1 (orientation, claim, PR) | cloud routine 013yHurq, 17:08–17:40Z |
| 2026-W42 | ≈ 0.3 (`arcs` on SPT0615: map/photo-z field survey incl. LAM atlas, plan, RELICS BPZ route, two RELICS models, sheets; no result) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01SmgqdA, 18:08–18:40Z |
| 2026-W42 | ≈ 0.25 (`arcs` on MACS0308: plan, two RELICS models, RELICS BPZ, sheets; does not replicate) | — | — | 0.1 (orientation, claim, PR) | cloud routine 018BWoMp, 19:08–19:35Z |
| 2026-W42 | ≈ 0.25 (Abell 370 halo fluxes: per-system rms gate, `halo_phot`, `fluxratios`, stamps) | — | — | 0.1 (orientation, claim, review, PR) | cloud routine 01StJeNF, 20:08–20:45Z |
| 2026-W42 | ≈ 0.3 (E1 + Swift BAT: plan, driver, 2 × 10⁴ scrambles, injections; control-channel diagnostics) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 016CcRnf, 21:08–21:45Z |
| 2026-W42 | ≈ 0.3 (E1 + Einstein Probe: GCN-circular node table, plan, driver, 2 × 10⁴ scrambles, injections; joint-event diagnostic) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01UFgoqS, 22:08–22:45Z |
| 2026-W42 | ≈ 0.2 (E1 + IceTracks-DR2: plan, driver, 2 × 10⁴ scrambles, injections) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01MS6Fgg, 23:08–23:40Z |
| 2026-W43 | — | — | ≈ 0.2 (E1 DR2–ICECAT control deficit: frozen plan and rule, three null variants; explained) | 0.1 (orientation, claim, state files, PR) | cloud routine 01EJiGRJ, 00:08–00:35Z |
| 2026-W43 | — | — | ≈ 0.25 (E1 BAT–GBM control deficit: frozen plan and rule; frozen rule undetermined) | 0.1 (orientation, claim, state files, PR) | cloud routine 01CSzjJs, 01:08–01:40Z |
| 2026-W43 | ≈ 0.2 (SMACS AT 2024jfz counter-images and Fermat delays in ICLv2; frozen testability rule; post hoc MAST epochs) | — | — | 0.15 (orientation, claim, state files, PR) | cloud routine 01CPQ1eo, 02:08–02:40Z |
| 2026-W43 | ≈ 0.25 (AT 2024jfz counter-image in NIRISS 4368: frozen plan and rule, injection control; post hoc z_s scan) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01VEAQr3, 03:08–03:40Z |
| 2026-W43 | ≈ 0.1 (AT 2024jfz host redshift: frozen plan and rule, DJA + MUSE match; untested) | — | — | 0.15 (orientation, #359 base merge, claim, state files, PR) | cloud routine 015gUbYY, 04:08–04:45Z |
| 2026-W43 | ≈ 0.15 (AT 2024jfz z_s 1.55–1.70 PSF-matched difference: frozen plan and rule, injection, control; untestable) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 015u2xwD, 05:09–05:35Z |
| 2026-W43 | — | — | ≈ 0.25 (review-list h039 / h005: frozen plan and rule, forced photometry, matched controls, stamps looked at; inconclusive, post hoc within control scatter) | 0.1 (orientation, claim, state files, PR) | cloud routine 01QCVQcw, 06:16–06:45Z |
| 2026-W43 | — | — | — | 0.3 (orientation; queue check: the two open-looking lensing items, the halo-phot elongation gate and the off-arc veto control, were already done; TASKS synced; PR) | cloud routine 01Hea1QD, 07:12–07:45Z |
| 2026-W43 | — | — | — | 0.25 (orientation; the in-flight W3 ZTF v11 branch was merged as #365 by its own session during this run; queue still blocked on owner decisions per the 07:12 handoff; no unit started) | cloud routine 017cBF6W, 08:09–08:35Z |
| 2026-W43 | ≈ 0.25 (W3 dimming Abell 2744 forced stage: frozen plan, 300 controls, calibrated, null, limits) | — | — | 0.15 (orientation, claim, state files, PR) | cloud routine 0125PxUc, 09:09–09:45Z |
| 2026-W43 | — | — | ≈ 0.3 (NEXUS W3 survivors: frozen plan, `_cal` exposure tool, chromatic veto (post hoc), 5 survivors + controls vetted, sheets looked at; 0 W3 events) | 0.1 (orientation, claim, state files, PR) | cloud routine 01Goi7d2, 11:12–11:50Z |
| 2026-W43 | — | — | — | 0.6 (orientation; brought W3 ZTF v13 #370 up to date with main, full offline tests, CI, merge per policy; queue otherwise blocked on owner decisions; no unit started) | cloud routine 01Mfyque, 12:09–12:50Z |
| 2026-W43 | ≈ 0.4 (W3 dimming MACS0416 regeneration: frozen plan and gate, S3 catalogue route, forced + inject, joint three-field `combine`) | — | — | 0.15 (orientation, claim, state files, PR) | cloud routine 01Ep8g9c, 13:12–13:55Z |
| 2026-W43 | ≈ 0.15 (W3 dimming: chromatic veto folded into `forced`, unit test) | — | — | 0.3 (orientation; issue #374 network-test failures diagnosed as upstream MAST/NED connection errors; state files, PR) | cloud routine 01CCHGKi, 14:12–14:50Z |
| 2026-W43 | ≈ 0.15 (W3 dimming: D-039 persistence-swap scoping, frozen rule, workload + veto-off what-if on three fields, reproduction exact) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 01Eubzuj, 15:12–15:30Z |
| 2026-W43 | ≈ 0.45 (W3 dimming: D-039 persistence run, frozen plan + tested rules; NEXUS 17 targets, partial controls on MACS0416 / Abell 2744) | — | — | 0.2 (orientation, merged #376 per policy, claim, state files, PR) | cloud routine 01WrZoCP, 16:08–16:50Z |
| 2026-W43 | ≈ 0.2 (W3 dimming: NEXUS `o014` forced check on the 12 `on_sky` targets, frozen plan + tested rule, sheet looked at; 0 survivors, mosaic-edge / depth) | — | — | 0.1 (orientation, claim, state files, PR) | cloud routine 019eUf3y, 17:09–17:35Z |
| 2026-W43 | ≈ 0.7 (W3 dimming: D-039 persistence run, MACS0416 e2–3 + Abell 2744 e1: all 20 targets `on_sky`, one epoch each; controls 194 / 194) | — | — | 0.15 (orientation, claim, state files, PR) | cloud routine 01GSZUHb, 18:08–19:10Z |
| 2026-W43 | — | — | ≈ 0.2 (W3 dimming injection edge coverage: frozen plan, rule and gate, `WHT` grids on 15 images, 4-variant what-if, sheet looked at; gate V failed as written, post hoc diagnostic) | 0.1 (orientation, claim, state files, PR) | cloud routine 01RmqBCN, 19:09–19:40Z |
| 2026-W43 | — | — | — | 0.5 (#380 W3 ZTF v14 brought up to date with `main` and merged under the policy; #374 live re-check of the three failing network tests) | cloud routine 01QjKQDx, 20:08–20:50Z |
| 2026-W43 | ≈ 0.4 (W3 dimming: D-039 persistence-test controls, frozen addendum; Abell 2744 e2, NEXUS e0/e1/e3/e4; NEXUS fails rule 1, exposure-count dependence post hoc) | — | — | 0.15 (orientation, #374 NED re-check, claim, state files, PR) | cloud routine 01CCiYbc, 23:08–23:50Z |
| 2026-W43 | — | — | — | 0.5 (#385 W3 ZTF v16 brought up to date with `main` (clean merge, offline suite green), `agent` label added; #374 NED re-check: still connection reset) | cloud routine 01EYNJfo, 00:08–00:45Z |
| 2026-W43 | — | — | ≈ 0.15 (W3 dimming: frozen `forced` + depth / footprint test on the 20 MACS0416 / Abell 2744 `on_sky` targets; 0 survivors, sheets looked at) | 0.1 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 01L4CY8s, 03:08–03:30Z |

## Outcome ledger (four columns only)
| Method improvements | Conditional theory results | Real candidates | Independent verifications |
|---|---|---|---|
| M2 stretch grid + free per-band scale (#167, provisional); W3 chain fixes (D-068); X3 G1: frozen W3 screen on unselected Gaia light curves with measured blind spots (RUWE veto of blended events, colour test, short t_E) and a held-out GAPS null, Γ₉₅ < 1 × 10⁻⁵ /star/yr conditional (#228) | Goal rounds 1–10 (#152–#166, provisional); far-side image bound vs GRAVITY summary (#172, provisional); X1 S2: no astrometric far-side discriminant at 9 µas, a fragile flux-contrast one for a = 2.1 M (#177); X1 S3: far side Δχ² +23 to +37 on the 2018-07-22 flare (a ≤ 7 M, exploratory, #179); X1 S4: model-free bound gives only +6.8 (a ≤ 3 M) and ≤ 0 (a ≥ 5 M), so S3 rests on the traced centroid (#180); X2 S4: with the NED + phantom axial equations, GW150914 η_max = 1.52 at a − 2M = 3.6e-4, outside the test-field bracket (#194); X2 S5: the S4 background (SV, phantom + magnetic NED, a − 2M = 3.6e-4) is radially unstable, e-folding 9.9 M (#204); X2 S6a: 5D equal-J rotating Ellis wormholes remain unstable just past the mode merger (oscillating pair verified at J/A 0.121, 0.142; even parity; #208) | none (n0153 inconclusive, #61) | none of a candidate (SkyLLH D-080 is a benchmark reproduction); X1 tracer swept angles and image weights checked by second routes (#180, #181); X2 test-field scattering: Numerov vs DOP853, time domain vs frequency domain (#192); X2 S4 two-channel tool reproduces Yang–Lan's a = 2.20 M QNMs (#194); X2 S6a: re-derived 5D equations reproduce 1309.2448's static eigenvalue and merger point; dense spectrum vs Evans shooting agree to ≤ 1e-7 (#208) |
| 2026-W43 | ≈ 0.25 (W3 dimming edge coverage: gate V2, stamp-mask sensitivity) | — | — | 0.1 (orientation, #374 re-check, claim, PR) | cloud routine 01FdEN2q, 21:09–21:35Z |
| 2026-W43 | ≈ 0.15 (W3 dimming: weight-map coverage term in `screen` / `inject`, gates R0 / R1, corrected headline) | — | — | 0.1 (orientation, claim, state files, review, PR) | cloud routine 01AUvG96, 22:08–22:40Z |
| 2026-W43 | — | — | ≈ 0.25 (W3 dimming: exposure-count `on_sky` efficiency, frozen plan, NEXUS e5–e7 controls held out; ≤ 4 class holds, e2 stopped) | 0.2 (orientation, #374 NED re-check, claim, state files, review, PR) | cloud routine 012vL5CG, 01:08–01:50Z |
| 2026-W43 | ≈ 0.25 (W3 dimming: D-039 persistence swap in `inject`, frozen plan, gates G0 / G1, joint re-`combine`, swap-cleared flags checked) | — | — | 0.1 (orientation, claim, state files, review, PR) | cloud routine 01EEgGY3, 04:08–04:45Z |
| 2026-W43 | — | — | ≈ 0.2 (W3 dimming: NEXUS `undetected` control stamps, measurable-pair efficiency; decision R2, NaN-tolerant aperture next) | 0.1 (orientation, claim, state files, review, PR) | cloud routine 01XJja1E, 05:08–05:40Z |
| 2026-W43 | ≈ 0.25 (W3 dimming NaN-tolerant aperture: k = 2 on NEXUS e1, e3–e7; G1 fails at 2 px, D2) | — | ≈ 0.1 (G0 identity, G1 fill accuracy on 42 known apertures) | 0.15 (orientation, #374 check, claim, review, PR) | cloud routine 01FXV1y5, 06:11–06:45Z |
| 2026-W43 | — | — | ≈ 0.35 (W3 dimming k = 1 aperture: frozen plan, G0 on e0 / e2, first NEXUS e2 controls at k = 0 and 1; R3, > 4 class holds) | 0.1 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 01MqRna1, 07:13–07:50Z |
| 2026-W43 | — | — | ≈ 0.25 (W3 dimming: MAST check for new NEXUS epochs; o015 Tile-6 frozen k = 1 check, R3, no monitored pair in its footprint) | 0.15 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 019GernwYABxMegDvtii19QF, 08:08–08:45Z |
| 2026-W43 | — | — | ≈ 0.1 (W3 dimming: held-out NEXUS `inconclusive` pairs split by cause, frozen plan; 12 / 12 C1 first-dither, post hoc `on_sky` time check) | 0.1 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 017FBUD4, 09:09–09:40Z |
| 2026-W43 | — | — | ≈ 0.15 (W3 ZTF fade control across 25 areas: frozen plan, homogeneity + inventory predictors (outcome C); post hoc sky-region clustering) | 0.15 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 01P3Lu4A, 10:09–10:40Z |
| 2026-W43 | ≈ 0.15 (W3 macro-shear: tool, gates G1–G3′, five cases, maps looked at) | — | ≈ 0.1 (validation on the isolated lens and the Chang-Refsdal caustic) | 0.15 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 01LNxfWW, 12:09–12:45Z |
| 2026-W43 | — | — | ≈ 0.2 (W3 macro-shear follow-up /code-review of merged #403: whole-window saddle minima, umbra chords, bin-area fix, re-run, maps looked at) | 0.1 (orientation, #374 status read, PR) | cloud routine 0168d1Q4, 13:09–13:40Z |
| 2026-W43 | ≈ 0.25 (`arcs` on RXC J0018+16: RELICS Lenstool maps, DJA v7.5 photo-z, frozen plan, run, sheet looked at) | — | — | 0.1 (orientation, claim, state files, review, PR) | cloud routine 014y1cGz, 14:08–14:45Z |
| 2026-W43 | — | — | ≈ 0.25 (`arcs` pooled anti tally over six replication fields: frozen plan, σ_ε scan, S1063 vs pooled rate; exploratory) | 0.1 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 01Cg2SZG, 15:15–15:50Z |
| 2026-W43 | — | — | ≈ 0.2 (`arcs` null calibrated on weak-shear rows, #407: frozen plan, σ_ε and isotropic-fraction fits, propagation; no leverage) | 0.15 (orientation, #374 status read, claim, state files, review, PR, this tally line) | cloud routine 01UexdJ1, 16:09–16:35Z |
| 2026-W43 | — | — | ≈ 0.2 (#409 W5 PSF radial: fourth /code-review, NaN/cache/PA fixes, EAS re-run reproduces summary.json; σ and coverage findings checked, not defects) | 0.1 (orientation, #374 status read, needs-human label + owner question on #409, this tally line) | cloud routine 01VV1Vfi, 18:09–18:35Z |
| 2026-W43 | — | — | ≈ 0.15 (MACS1149 system 16 published spec-z check: frozen plan, Treu+16, DJA NIRSpec v4.4; not decided) | 0.15 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 01RJnu6v, 19:08–19:45Z |
| 2026-W43 | — | — | ≈ 0.2 (Abell 2744 700.1 z scan in Bergamini+23: frozen plan, known-case check 39 / 40, 121 z, sheet looked at; outcome B, model-limited) | 0.15 (orientation, DJA photo-z index recheck: no new cluster, claim, state files, review, PR) | cloud routine 01Az7Cqh, 20:08–20:45Z |
| 2026-W43 | — | — | ≈ 0.4 (Abell 2744 700.1 in UNCOVER v2.0: frozen gate fails, tool not validated; post hoc 700.1 = UNCOVER 78, z_spec 2.3151, a/b merged in both models; undetermined, /vet-candidate next) | 0.1 (orientation, #374 status read, this tally line) | cloud routine 017ECfaD, 21:09–21:55Z |
| 2026-W43 | — | — | ≈ 0.4 (Abell 2744 700.1 /vet-candidate: frozen plan, 7-band pair colours on 20 known pairs (gate fails, 10 / 20), perturbers, sheet looked at; post hoc Bergamini+23 components) | 0.1 (orientation, claim, state files, review, PR) | cloud routine 012djKdn, 22:09–22:45Z |
| 2026-W43 | — | — | ≈ 0.3 (MACS1149 Sharon v4cor pinned: frozen plan, G1 maps pass, G2 fails as written on system 9 off the map; CATS system 13 needs its fitted z) | 0.1 (orientation, claim, state files, review, PR) | cloud routine 01NLjt8m, 23:09–23:45Z |
| 2026-W43 | — | ≈ 3.1 (X4 S1g: frozen budget top-up, 6 / 6 P checks pass, R1 seeds disagree, undetermined; 6.8 CPU-h on 3 cores, run in parallel) | ≈ 1.2 (S1 CHIME Cat 2 30-event validation: frozen method does not transfer; W2 HST Moffat chain: negative control S/N 13, validation fails; two /code-review rounds and re-runs) | ≈ 1.0 (orientation, #153/#374 checks, coordination, state files, review, PR) | interactive owner parallel round 01BBF8ia, 00:05–03:40Z (`docs/plans/2026-10-21-parallel-round.md`) |
| 2026-W43 | — | — | ≈ 0.15 (MACS1149 D-041 cross-model with the pinned Sharon: frozen plan, gate 38 / 39, 16 rows; one frozen flag, system 16 at its fitted z, D-037's case) | 0.1 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 01D7FXqw, 00:08–00:40Z |
| 2026-W43 | — | — | ≈ 0.25 (MACS0717 Sharon v4cor + Keeton v4 pinned: frozen plan, Table 7 known cases, G1 / G2 per model; Keeton passes, Sharon G2 fails as written on 6.2 off the map) | 0.15 (orientation, #374 status read, claim, state files, review, PR) | cloud routine 017CvAPJ, 01:08–01:45Z |
| 2026-W43 | — | — | ≈ 0.3 (MACS0717 D-041 cross-model with the pinned Keeton / Sharon: frozen plan, CATS images re-run, gate fails at CATS-fitted z 0.55 / 0.59; post hoc arcs.txt-z systems 0 to vet) | 0.15 (orientation, #374 status read, claim, TLS route, state files, review, PR) | cloud routine 01C96NzD, 02:08–02:50Z |
| 2026-W43 | ≈ 0.35 (MACS0417 / MACS1423 published images vs CANUCS v1 maps: image-list search, frozen plan, gate K fails (one central image each), post hoc leave-one-out, RELICS map check) | — | — | 0.15 (orientation, #374 status read, claim, TLS route, state files, PR) | cloud routine 0169xSUp, 03:08–03:45Z |
| 2026-W43 | — | — | — | 0.15 (orientation, #374 status read, queue triage: no unblocked step fits a slot; CANUCS `extra_predicted` forced photometry not started, since gate K failed as frozen and a central-image-free gate would be post hoc) | cloud routine 01VuUwtu, 04:08–04:35Z |
| 2026-W43 | ≈ 0.45 (W3 ZTF v25: frozen sibling-prediction plan, objects + samples, 4 areas validated, 2346 searched, 0 survivors; 4 / 4 as predicted, 2 carry over) | — | — | 0.15 (orientation, #374 status read, claim, ASAS-SN probe refused by session policy, state files, review, PR) | cloud routine 01VXHtxQ, 05:09–05:50Z |
| 2026-W43 | ≈ 0.2 (W3 ZTF v25 carry-over: 2347 and 4990 decided, 2347 searched, 0 survivors; sibling tally 6 / 6, H-sib adopted as queue rule) | — | — | 0.15 (orientation, claim, state files, review, PR) | cloud routine 01QcXKc4, 06:12–06:40Z |

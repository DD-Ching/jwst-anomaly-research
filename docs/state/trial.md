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
| X4 wave-optics GW lensing by negative-mass objects | `docs/state/units/X4.md` | `claude/s0-2026-10-15` | Authorized by S0-2026-10-15 (B: VASCO and C: quantum-processor wormhole not qualified). S0 minimal step: closed-form repulsive-lens F(w) (low-pass), separable from a positive point mass at observed SNR ρ ≈ 3.1–16.5 for M_Lz ≥ 300 M☉, y ≤ 1, intrinsic parameters fixed (strongly demagnified cases). **S1 done:** 18 / 36 simulated cases keep a residual after free aligned-spin binary (q ≤ 15) and positive-lens fits (frozen refinement; IMRPhenomD agrees); extreme q absorbs part. **S1b done:** R2, R3 no longer pass the screening threshold against precessing, higher-mode binaries (ρ_req > 30). **S1c done (R1):** residual kept against validated aligned-spin fits (ρ_req 5.8, bound-conditional, provisional); undetermined with XPHM (best ρ_req 9.2). **S1d:** the SciPy DE reference search fails its injection checks; R1 still undetermined. No pilot proposed. ~8.4 / 18 h |

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
| 2026-W42 | ≈ 0.2 (MACS0416 o211/o212 grid; contact sheets) | — | ≈ 0.2 (validation on Yan+2023 transients, grid fixes, persistence test) | 0.1 (orientation, claim, PR) | cloud routine 01AGqJkf, 03:08–03:45Z |
| 2026-W42 | ≈ 0.1 (MACS0416 o211/o213 grid; contact sheets, persistence test) | — | — | 0.05 (PR) | cloud routine 01AGqJkf, 03:30–03:45Z |
| 2026-W42 | ≈ 0.35 (`arc_diff.py`: MACS0416 o211/o212 arc differencing, validation on Yan+2023, three fix iterations; contact sheets) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01JZas31, 04:08–04:45Z |
| 2026-W42 | ≈ 0.4 (`arc_diff.py` pixel kernels: five PSF-match variants on o211/o212 and o211/o213, validation on Yan+2023; contact sheets) | — | — | 0.1 (orientation, claim, PR) | cloud routine 01VwokZ2, 05:08–05:50Z |
| 2026-W42 | — | — | ≈ 0.4 (vet MACS0416 h108: five epoch pairs incl. CANUCS o004, `arc_diff_vet.py`, validation on Yan+2023; contact sheets) | 0.1 (orientation, claim, PR) | cloud routine 01VE5srt, 06:08–06:50Z |
| 2026-W42 | ≈ 0.35 (Abell 2744 DD 2756 grid, pair choice, `i2d_region.py`; contact sheets) | — | ≈ 0.1 (persistence test, UNCOVER third-epoch stamps) | 0.1 (orientation, claim, PR) | cloud routine 01HeG8Gq, 09:09–09:45Z |

## Outcome ledger (four columns only)
| Method improvements | Conditional theory results | Real candidates | Independent verifications |
|---|---|---|---|
| M2 stretch grid + free per-band scale (#167, provisional); W3 chain fixes (D-068); X3 G1: frozen W3 screen on unselected Gaia light curves with measured blind spots (RUWE veto of blended events, colour test, short t_E) and a held-out GAPS null, Γ₉₅ < 1 × 10⁻⁵ /star/yr conditional (#228) | Goal rounds 1–10 (#152–#166, provisional); far-side image bound vs GRAVITY summary (#172, provisional); X1 S2: no astrometric far-side discriminant at 9 µas, a fragile flux-contrast one for a = 2.1 M (#177); X1 S3: far side Δχ² +23 to +37 on the 2018-07-22 flare (a ≤ 7 M, exploratory, #179); X1 S4: model-free bound gives only +6.8 (a ≤ 3 M) and ≤ 0 (a ≥ 5 M), so S3 rests on the traced centroid (#180); X2 S4: with the NED + phantom axial equations, GW150914 η_max = 1.52 at a − 2M = 3.6e-4, outside the test-field bracket (#194); X2 S5: the S4 background (SV, phantom + magnetic NED, a − 2M = 3.6e-4) is radially unstable, e-folding 9.9 M (#204); X2 S6a: 5D equal-J rotating Ellis wormholes remain unstable just past the mode merger (oscillating pair verified at J/A 0.121, 0.142; even parity; #208) | none (n0153 inconclusive, #61) | none of a candidate (SkyLLH D-080 is a benchmark reproduction); X1 tracer swept angles and image weights checked by second routes (#180, #181); X2 test-field scattering: Numerov vs DOP853, time domain vs frequency domain (#192); X2 S4 two-channel tool reproduces Yang–Lan's a = 2.20 M QNMs (#194); X2 S6a: re-derived 5D equations reproduce 1309.2448's static eigenvalue and merger point; dense spectrum vs Evans shooting agree to ≤ 1e-7 (#208) |

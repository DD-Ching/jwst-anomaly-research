# Pending conclusions index (needs-human PRs)

**Every entry is PROVISIONAL.** The results below sit on unmerged, owner-unreviewed PRs. Listing a PR here is
bookkeeping only: it does not approve, merge or endorse any scientific conclusion. Read the PR itself before
reusing a number.

**Purpose:** a session that orients on `main` must not redo work that these PRs already did or paused (D-083).

**Upkeep:** when a needs-human PR changes head, is merged or is closed, update its row in the same cycle. Head
SHAs are as of 2026-10-12 23:00 UTC (rows #199, #237, #279: 2026-10-18 15:20 UTC); base merges by other sessions move them without changing the conclusions.

| PR | Head | Scope of the provisional conclusion | Open questions | Reopen only if |
|---|---|---|---|---|
| #61 vet n0153 | 706caa0 | Sunrise n0153: a real nuclear brightening, inconclusive without a third epoch | Third-epoch imaging | A new JWST epoch of the field is public |
| #145 W3 gb13 | d14a48b | One chain survivor, explained by red noise; not a candidate | Owner decision on the red-noise explanation | — (field result; gb13 is not to be re-run) |
| #150 coordinator writer | 8ed4a29 | Proposal D-081: one coordinating writer for CHANGELOG/TASKS/exotic_limits | Owner decision | — |
| #152–#166 goal rounds 1–10 | see PRs | Communication-channel feasibility: published traversable constructions give no shortcut (r1); v-causal bounds (r2–3); beyond-Horndeski wormhole stability fails (r4); Krasnikov tube needs Planck-thin walls (r5–7); Einstein–æther 1-bit channel not practical in the computed case (r8); c13 ≠ 0 instantaneous response is model-internal (r9); no latency route found externally (r10) | Owner review of each round; NIST 2015 no-signalling check undetermined (r10) | A new mechanism or input named in the round's own "reopen" list |
| #167 S1 M2 re-screen | a59c255 | Frozen M2 method (stretch grid 1.05, free per-band scale); 747,709 pairs, 0 candidates; exploratory; conditional limit < 1.3 × 10⁻² deformed-twin pairs per eligible burst | Owner review | **Done; do not redo:** the s ≠ 1 chain and the full re-screen |
| #169 S1 blind spots | 95f0198 | Single-pulse bursts are not searchable by shape under M2. TTE timing helps about 3/30 bright short bursts; TTE expansion not run. Catalogue E_peak + M2 does not reach the pre-declared robust gain (ΔR +1.8/+1.4 pp optimistic/central, −0.1/−1.0 pp in the window/literature scenarios); **S1 energy upgrade paused** | Owner review | **Done; do not redo:** TTE for short/faint bursts, the E_peak joint constraint. Reopen only with a same-source E_peak systematic ≲ 0.10 shown on real data, or time-resolved spectra |
| #172 far-side light (X1 S0) | b7a5daa | SV far-side image bound b < b_c(a); against GRAVITY's single averaged loop (61 ± 9 µas), summary-level gaps only. a ≈ 7–9 M not distinguished without EHT; the observable is not yet sufficient | Owner review | X1 continued it and closed 2026-10-13 (docs/state/units/X1.md "Closure"; reopen inputs there) |
| #199 arcs flip scan + vet Abell 370 1565 | 0dbb6de | κ(z) scan of the 7 anti arcs in five fields: 1788 and MACS1149 1703 are κ ≈ 1 flips (explained); D-024 addendum 4 corrects addendum 3. 1565 vetted as inconclusive: z ≈ 2.9 photo-z vs NIRISS F150W; it sits in a z 0.386 MUSE [O II] nebula, with no redshift of its own; not a candidate | Owner review of the verdict; whether to run a grizli NIRISS extraction (best case ~4–5σ [O III]) | **Done; do not redo:** the κ(z) scan, five-model κ at 1565, the member subtraction, the MUSE-DEEP extraction, the NIRISS WFSS pipeline check and its line sensitivity, forced NIRCam/NIRISS photometry (colour undetermined), the row-artefact axis test. Reopen with a redshift for 1565 or imaging that fixes F150W/F200W to ≤ 0.15 |
| #237 two-epoch grids o052 + SMACS/VENUS | 3d5de60 | Sunrise o010/o052 grid: 0 candidates (−0.75 mag recovered 97 %). SMACS o001/VENUS: the pair on the system 5 arc is SN Atalanta (Fujimoto et al. arXiv:2609.30440), not a candidate; `g00791` most likely an unreported SN, `g02290` a nucleus brightening (AGN-like), neither a candidate | Owner review | **Done; do not redo:** both grids, the Atalanta identification and epoch check (`smacs_transient_epochs.py`). Reopen: a third epoch for `g02290`; a time delay for the ICLv2 third image only if VENUS publishes none |
| #279 clump-pair census, five other `arcs` fields | 882103e | Pooled two-knot rate f̂ 0.280 (95 % upper 0.370; 118 gate-passed rows); S1063 P(≥ 4 anti) 0.032 (0.056 at the upper); frozen outcome "excess remains" (exploratory). Re-run with #273's fixed `fit_models` (3 rows changed, outcome unchanged). On `main` since: the anti excess does not replicate in MACS0417, MACS1423 (clump-aware), WHL0137 | Owner call on the S1063 result | **Done; do not redo:** the five-field census with the fixed code. Reopen: 1316/1557 spectroscopy (116.28VK.001) or a new independent field |

Merged since the table was written: #273 (S1063 arcs vetting, 2026-10-16; record on `main` in
`results/arcs/abells1063-cats/vetting/`; owner decision: wait for 116.28VK.001, no stacking or whole-cube download).

# Pending conclusions index (needs-human PRs)

**Every entry is PROVISIONAL.** The results below sit on unmerged, owner-unreviewed PRs. Listing a PR here is
bookkeeping only: it does not approve, merge or endorse any scientific conclusion. Read the PR itself before
reusing a number.

**Purpose:** a session that orients on `main` must not redo work that these PRs already did or paused (D-083).

**Upkeep:** when a needs-human PR changes head, is merged or is closed, update its row in the same cycle. Head
SHAs are as of 2026-10-12 23:00 UTC; base merges by other sessions move them without changing the conclusions.

| PR | Head | Scope of the provisional conclusion | Open questions | Reopen only if |
|---|---|---|---|---|
| #61 vet n0153 | 706caa0 | Sunrise n0153: a real nuclear brightening, inconclusive without a third epoch | Third-epoch imaging | A new JWST epoch of the field is public |
| #145 W3 gb13 | d14a48b | One chain survivor, explained by red noise; not a candidate | Owner decision on the red-noise explanation | — (field result; gb13 is not to be re-run) |
| #150 coordinator writer | 8ed4a29 | Proposal D-081: one coordinating writer for CHANGELOG/TASKS/exotic_limits | Owner decision | — |
| #152–#166 goal rounds 1–10 | see PRs | Communication-channel feasibility: published traversable constructions give no shortcut (r1); v-causal bounds (r2–3); beyond-Horndeski wormhole stability fails (r4); Krasnikov tube needs Planck-thin walls (r5–7); Einstein–æther 1-bit channel not practical in the computed case (r8); c13 ≠ 0 instantaneous response is model-internal (r9); no latency route found externally (r10) | Owner review of each round; NIST 2015 no-signalling check undetermined (r10) | A new mechanism or input named in the round's own "reopen" list |
| #167 S1 M2 re-screen | a59c255 | Frozen M2 method (stretch grid 1.05, free per-band scale); 747,709 pairs, 0 candidates; exploratory; conditional limit < 1.3 × 10⁻² deformed-twin pairs per eligible burst | Owner review | **Done; do not redo:** the s ≠ 1 chain and the full re-screen |
| #169 S1 blind spots | 95f0198 | Single-pulse bursts are not searchable by shape under M2. TTE timing helps about 3/30 bright short bursts; TTE expansion not run. Catalogue E_peak + M2 does not reach the pre-declared robust gain (ΔR +1.8/+1.4 pp optimistic/central, −0.1/−1.0 pp in the window/literature scenarios); **S1 energy upgrade paused** | Owner review | **Done; do not redo:** TTE for short/faint bursts, the E_peak joint constraint. Reopen only with a same-source E_peak systematic ≲ 0.10 shown on real data, or time-resolved spectra |
| #172 far-side light (X1 S0) | b7a5daa | SV far-side image bound b < b_c(a); against GRAVITY's single averaged loop (61 ± 9 µas), summary-level gaps only. a ≈ 7–9 M not distinguished without EHT; the observable is not yet sufficient | Owner review | X1 continued it and closed 2026-10-13 (docs/state/units/X1.md "Closure"; reopen inputs there) |
| #273 arcs S1063 | (this PR) | 4 background anti vs 0.52; clump-aware P 0.026 is exploratory (simplified ordinary model; clump fraction uncertain; 2004's κ ≈ 1 not modelled); 2184 compatible with the tested noise model; 1316, 1557 pending vetting. Public spectroscopy: no catalogued z; the MUSE SODA cutouts give no line (excludes only extreme emitters); MUSE cannot resolve the clumps | Owner: stack the public 116.28VK OB cutouts, or wait for that program | **Done; do not redo:** steps 1–5 in `results/arcs/abells1063-cats/vetting/`. Reopen with new spectra or a second field's census |

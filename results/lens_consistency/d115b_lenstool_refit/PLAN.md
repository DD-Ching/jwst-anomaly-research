# D115b Lenstool re-fit: stage V, Lenstool's χ² of the pinned model with nothing freed (plan, frozen 2026-10-26 10:20Z)

Why: the D115b vetting record (`docs/candidates/macs0416-canucsref-D115b.md`) leaves one step: a Lenstool re-fit with
member 8515's σ free (and other parameters free), since Lenstool runs in the cloud (D-109). Before any parameter is
freed, Lenstool's own χ² evaluation of the pinned model must be validated; otherwise a later Δχ² is uninterpretable.

**Hypothesis (stage V):** Lenstool 8.15.5 `runmode inverse 0` (χ² evaluation, no optimisation) on `best.par`, with the
image section of the published `input.par` and the published image list, returns the χ² written in the `best.par`
header by the CANUCS optimisation (`#Chi2pos: 344.298308`, dof 354; `model_prediction`/`derived` by the authors).
Ordinary failure modes: version differences (the header was written 2024-09-05, version not stated), a missing
optimisation-only setting (e.g. per-family z or σ_pos rules not in `input.par`), a different image list than
`gold_canucs_v2.fk5.dat`, or a different χ² definition.

**Inputs** (pinned, sha256 in `lensmodel.MACS0416_CANUCS`; fetched via the MAST download API):
`best.par`, `input.par`, `arcs.dat` (= the published multiple-image file, 303 rows).

**Run file `v0.par`** (built by `scripts/d115b_lenstool_refit.py`, committed before the run): `best.par` unchanged
except (1) its `runmode` section becomes `reference 3 <best.par reference>`, `inverse 0`, `end`; (2) the `image`
section of `input.par` (`multfile 1 arcs.dat`, `forme -1`, `sigposArcsec 0.49`) is inserted after `runmode`.

**Gate V (frozen):** pass if Lenstool's reported position χ² is within 1 % (ASSUMPTION) of 344.298308.
- Pass → stage F (restricted re-fit, 8515 σ free) gets its own frozen plan before any run.
- Fail → record the value and the cause if identified; no re-fit Δχ² is read from Lenstool until a validated
  reference exists. Post hoc diagnosis is labelled as such.

Reported, not gated: the per-family χ² (`chires.dat` if written), Lenstool's dof and run time, and the port's
fiducial χ²_pos (330.81, record line "Grid step") for comparison.

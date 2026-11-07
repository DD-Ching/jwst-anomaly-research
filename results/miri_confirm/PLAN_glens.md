# Plan: vet MIRI's known NED `G_Lens` object (TASKS Then 5, 2026-11-07)

The first end-to-end run (`20261007T020124Z-4bfabaaa`) recovered one NED `G_Lens` object in the MIRI top 20; it was
never vetted. That run predates jwst 3.0.0, so the object is re-identified on the current MIRI config
(`configs/smacs_jwst3_miri_confirm.yaml`): every top-20 row whose cross-match lists NED `G_Lens` within the
classify radius, excluding the five already vetted (`docs/candidates/smacs_miri_f770w_top20_singles.md`).

**Hypothesis (ordinary, to test first).** The object ranks high because it is a lensed (magnified) background
galaxy of SMACS 0723 with PAH/dust-bright MIRI colours, or because of a MIRI processing or catalog effect (edge,
low weight, deblending, bright-neighbour light, single-band merge radius).

**Lensing-consistency question.** If it is a catalogued multiple image, does its position agree with the
Mahler+2022 model images of its system (within the model's rms)?

**What would change the plan.** An object that survives A–C of `/vet-candidate` goes to a needs-human candidate PR.
If no top-20 row carries `G_Lens` in the current run, record that the Oct-7 object no longer ranks and close the item.

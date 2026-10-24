# MACS J1931 stage 3c: wider background box, MACS0417 gate first, chance-controlled frame check F (frozen before any run)

**Why.** Stage 3 and its deviation 3b stopped at frame check F (18 / 51 < 20; `../macs1931_stage3/RESULT.md`). The
contact sheet there showed arcs 8.2, 9.1, 10.1 and 12.2 visible in F200W but not segmented at 3σ / 25 px after the
2″ background filter (`m1931_deblend.BKG_BOX` = 64 px at 0.0312″/px). Stage 3c changes only that box. Written after
stage 3's failure, before any build with the new box on either field.

## Rules (each a function in `scripts/m1931_stage3c.py` or the reused stage 3 code)
- **3c-0 builder:** `m1931_deblend.build` with `bkg_box` = 256 px (≈ 8.0″; ASSUMPTION: longer than the unsegmented
  arcs). Every other builder parameter is unchanged (3σ, 25 px, kernel, deblending, Kron columns).
- **3c-1 gate (MACS0417):** the stage V2 chain (`m1931_deblend.stage_v2`) with the new box on the same MACS0417
  section, output `../macs0417_ownref_3c/`. Passes iff `gate_C2_shared.passes` (as V2). **If it fails, stage 3c stops**
  and nothing is built on MACS J1931.
- **3c-2 frame check F (MACS J1931):** stage 3 rules 1–5 with the 3b section (all 51 used images) and the new box.
  F passes iff the stage 3 rule passes (≥ 20 of 51 within 0.5″, |median dx|, |median dy| ≤ 0.1″) **and** a chance
  control passes: the 51 positions moved by 5″ along each of 8 position angles (0°, 45°, …, 315°) are matched the
  same way; `n_matched ≥ 2 × mean(shifted counts)` (ASSUMPTIONs: 5″, 8, 2). The control is new: a deeper catalogue
  raises the chance of a 0.5″ match. If F fails: no reading, stop.
- **3c-3 (only if F passes):** stage 3 rules 6–10 unchanged (V2-b chain with and without the containing-source rule,
  footprint rule, local-noise rule, gate C2 with shared pairs, flags = `absent` predicted rows when C2 passes).
  Flags are unvetted until `/vet-candidate`; no reading without it.

## Run
`python scripts/m1931_stage3c.py gate`, then `python scripts/m1931_stage3c.py run --maps DIR` (refuses unless the
gate passed). Post hoc only: `... posthoc`.

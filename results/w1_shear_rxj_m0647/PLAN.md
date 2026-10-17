# W1 shear screen on RXJ2129 and MACS0647 (plan, frozen before the run)

Cloud routine, 2026-10-17 (claim PR #306). TASKS "JWST focus" 1 (W1): the next fields named in the handoff of
`results/w1_shear_more/RESULT.md`. Both have a JWST F200W pipeline catalogue, DJA v7.0 photo-z and map models that
pass `validate` (`results/arcs/rxj2129/`, `results/arcs/macs0647-clash/`):

| Field | Catalogue (sha256) | Photo-z (sha256) | Model (z_l) |
|---|---|---|---|
| RXJ2129 | `jw02767-o002_t001` F200W (`2ed1a2c7…`) | DJA v7.0 rxj2129 zout (`651c0b24…`) | `rxj2129-lenstool-maps` (0.234) |
| MACS0647 | `jw01433-o010_t001` F200W, epoch 1 (`473e32d2…`) | DJA v7.0 macs0647 zout (`ceb199c6…`) | `macs0647-clash-nfw` (0.591) |

The primary models of the `arcs` runs. Both catalogues are shallow (about 2 ks), so the samples are small. The
`arcs` alignment controls failed on both fields (67 %, one row short), but that control is not part of W1: the
shear screen's own known-case check is R (below). Inputs streamed to the gitignored cache, deleted after the run.

**Hypothesis (H).** A negative-mass (W1) lens of |M| ≥ 2 × 10¹² M☉ in a field gives a radial (E-mode positive in
the D-053 sign convention) aperture-mass peak.
**Ordinary explanations of a peak.** Diffraction spikes and streaks (two earlier hits, S1063 and MACS1423, sat on
linear streaks), PSF anisotropy and isophotal-moment noise (B-mode floor), blends near bright members, residual
cluster shear from model error (the CLASH model is HST-era), edge apertures.
**Rule (D-053 defaults, unchanged; ASSUMPTIONs).** Schirmer Q_TANH, R_ap 10″, 1″ grid, sources 1–10″, rotation
null; a peak counts only if p_rot < 0.05 **and** S > the field's max |S_×|. Grid half-width 100″ (as the other
shear-only fields; ASSUMPTION). Spike veto: catalogue segments plus Gaia DR3 G < 17 within 4′ of the model
reference point.
**Injections.** `inject_shear.py --fields <field>` seed 50, 200 per mass at the five D-053 masses, then an
eleven-field combination (nine-field headline plus these two) with `--combine-only`.
**What changes the plan.** A peak passing the rule: look at the F200W sheet; if it sits on a streak, run the
frozen `scripts/w1_peak_vet.py` test of `results/w1_shear_more/PLAN.md` addendum A unchanged (same outcome rule,
with this field's B-mode max), otherwise `/vet-candidate`. R outside (0, 1.5] or < 20 calibrating rows: the field
is dropped (D-053). A field that passes the rule is left out of the headline combination, as MACS1423 was.
**Known-case check.** R should fall near the 0.41–0.60 of the ten earlier fields (the cluster shear is seen).

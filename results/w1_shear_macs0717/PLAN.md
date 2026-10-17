# W1 shear screen (D-053) on MACS0717 with ASTRODEEP-FF photo-z: plan (2026-10-17, frozen before any run)

**Why.** TASKS "JWST focus" 1 (W1): MACS0717 was the one screened HFF cluster without W1, because DJA v7 has no
MACS0717 photo-z. HFF-DeepSpace (Shipley+18) downloads return HTTP 401 (2026-10-17). ASTRODEEP-FF (Di Criscienzo+17,
VizieR J/A+A/607/A30, tables `m0717cla` + `m0717clz`) is public: 4,082 rows, 2,866 with `relFlag` = 1.

**Hypothesis (W1, D-053).** A negative-mass (or Ellis) lens near the cluster gives a coherent radial (sign-reversed)
shear pattern in background-galaxy shapes; the screen looks for an aperture-mass peak S.
**Ordinary explanations, tested in this order:** shape noise (rotation null), cluster-model residual (CATS v4.1
shear removed; a peak is checked against the model's critical-curve region), PSF/streak/spike systematics (Gaia
spike veto, streak-mask veto of results/w1_streak_veto), blends with members (sheet looked at).
**What changes the plan:** a D-053 pass goes to `/vet-candidate`; a null adds MACS0717 to the W1 limit.

**Inputs (labels).** F200W VENUS `jw06882-o029_t063` catalogue (`observed`, sha256 in data/manifests/macs0717.ecsv);
ASTRODEEP-FF zbest → zout (`scripts/astrodeep_zout.py`, `derived`): z160/z840 = zbest ∓ 1.483 SIQR, rows with
`relFlag` ≠ 1 NaN (ASSUMPTIONs); match radius 0.3″ (`lens_consistency.attach_photoz`, unchanged). Model `macs0717-cats`
(CATS v4.1 maps, `model_prediction`). Gaia DR3 G < 17 within 4′ of the VENUS pointing (109.3864, +37.7550).

**Run (frozen; same settings as results/w1_streak_veto).** `inject_shear.py --fields macs0717 --streak-masks <masks>`
(seed 50, 200 lenses per mass, D-053 screen defaults, max radius 100″); streak mask from `w1_streak_mask.py --fields
macs0717` (unchanged code). Then `--combine-only` with the twelve veto-run fields (thirteen-field limit).

**Rule (computable; unchanged D-053, `inject_shear` / `exotic_screens`):** pass iff p_rot < 0.05 **and** S_max >
max |S_×|. Limits use the injections' efficiency per mass.

**Known limitation stated in advance.** Photo-z exists only over the HST WFC3 footprint (≈ 2′ × 2′ inside the NIRCam
mosaic), so rows outside it have no photo-z; the screen treats them as it treats rows without photo-z in other fields
(unchanged code). Expect a small effective area.

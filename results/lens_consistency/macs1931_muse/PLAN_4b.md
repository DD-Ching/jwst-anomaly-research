# MACS J1931 stage 4b (deviation from PLAN.md): register the MUSE cube to Gaia, rules otherwise unchanged

**Why (post hoc, not blind).** Stage 4 found no ≥ 5σ Lyα reference (best S/N 1.83, `summary.json`). A post hoc
check (`posthoc_gaia_offsets.json`, centroids on the z 5.079 band's median image, MUSE − Gaia DR3 at epoch 2015.5)
shows the cube frame is offset by a median (+1.30″, −0.88″) in (RA cos δ, Dec) on 6 stars with proper motions,
larger than the 0.8″ aperture. The image list is in the JWST (Gaia-tied) frame, so the stage 4 apertures sat beside
the images. A post hoc line search within 2″ of the system 6 and 100 / 101 images was also looked at before this
plan was written; this stage is therefore not blind. Gate L stays the quality gate.

## The only change
**Rule 0 (registration):** for each band cutout, Gaia DR3 sources inside the cube footprint with a proper motion
are moved to epoch 2015.5; each is centroided (`photutils.centroids.centroid_2dg`, 15 × 15 spaxel stamp minus its
median) on the cutout's per-spaxel median image; stamps touching the edge or a NaN are skipped. The offset is the
per-axis median of MUSE − Gaia. It is applied iff ≥ 4 stars are measured **and** each axis's MAD ≤ 0.2″
(ASSUMPTIONs); else the stage stops. Every listed image and target position is shifted by the offset into the cube
frame before rules 2–9, which are unchanged. The Gaia query (footprint box, `gaiadr3.gaia_source`) is stored as
`gaia_dr3.csv`.

Outputs go to `stage4b/`.

# Orphan-pair (dark-deflector) screen on Abell 2744: rules frozen before the run (2026-10-16)

Owner request 2026-10-16: one observational lensing test, at most 2 h, in a region not yet tested. TASKS
("orphan pairs … re-run on Abell 2744 (UNCOVER)") lists it as open. No orphan-pair statistic has been computed for
Abell 2744. The field itself has been used by other tests (counter-images, flux ratios, transients), so this is the
first run of this statistic here, **not** a blind sample of an unseen field.

## Question and ordinary model
- **Question.** Behind Abell 2744, are there more "orphan" image pairs than chance and redshift clustering predict?
  An orphan pair is two background sources 0.3–3″ apart with matching SEDs and photo-z, no catalogued galaxy
  between them and no published multiple image nearby.
- **Ordinary model.** The orphan count equals the D-055 nulls (e) and (f): chance SED matches plus same-redshift
  neighbours. Each top orphan is one of three things:
  - knots of one galaxy;
  - group members;
  - a cluster-scale pair the lens model explains (high |μ| or a parity flip near a critical curve).
- **Exotic alternative.** An unseen compact deflector, such as a dark point mass, a wormhole-like n = 2 lens or a
  negative-mass lens (D-047/D-051 types). It would add pairs that the cluster model does not explain and that have
  no visible lens.

## Data and tools (reuse; nothing new beyond configuration)
- **Catalogue.** DJA grizli v7.2-fix `abell2744clu` phot plus its eazy zout, both already pinned by sha256
  (`configs/abell2744.yaml`, `scripts/a2744_followups.py`). The phot file is 233.5 MB; the stated reason is in the
  config. In the cloud it is downloaded to the scratchpad and deleted after use.
- **Published images.** CATS v4.1 `arcs.txt` (`abell2744-cats`, 109 images) and the Bergamini+23 list
  (`abell2744-bergamini23`).
- **Lens check.** CATS v4.1 magnification and parity at both members.
- **Code.** `orphan_pairs.py` with every threshold unchanged (D-048, D-051, D-055 constants). The only change is a
  `FIELDS["abell2744"]` entry with the `dja` format and a plain (non-tar) zout file.

## Checks before reading the result
1. **Alignment.** Take the median offset between DJA sources and the published image positions they match within
   1″. It must be ≤ 0.3″, against the 1″ `IMAGE_RADIUS`; otherwise stop and record.
2. **Positive control.** SED-matched close pairs classed `published` are known lensed pairs that the screen
   recovers. Their number is reported. If it is 0, the field has no positive control at 0.3–3″, and that is stated.
3. **Injection recovery.** Run `inject_pairs.py` (point, W2, W1) on this catalogue if it accepts a cluster field
   with no new code beyond the field name. Otherwise record the gap and use the deep-field efficiencies only as a
   reference, not as this field's efficiency.
4. **Artifact burden.** Count the top orphans that the contact sheet shows as artifacts (spikes, edges, ICL
   residuals).

## Decision rules (frozen)
- **Count test (D-055's revisit criterion).** There is an excess only if P < 0.05 under **both** null (e) and
  null (f). Otherwise the field is null.
- **Escalation to `/vet-candidate`.** A top orphan (`--top 15`, `--cutouts`) goes on only if all of these hold:
  - it does not look like knots of one galaxy, or an artifact, on the sheet;
  - the CATS model gives |μ| < 3 and no parity flip at both members;
  - there is no visible deflector.
  Escalation is not a claim: vetting tests instruments and ordinary astrophysics first.
- **No changes after the run.** No threshold is changed after seeing the result. Any post hoc change is labelled
  and does not count toward the decision.

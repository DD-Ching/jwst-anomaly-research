# X2 S6a: prior-work (duplication) check, 2026-10-14

**Question.** Before any wording about novelty: has anyone already reported, for the 5D equal-J rotating Ellis
wormholes of arXiv:1309.2448, either of the following?
- an instability past the merger of the two purely imaginary modes (complex or oscillating modes), or
- an error in that paper's printed perturbation equations (V₁₀)?

**Method.**
- Read-only web search by a Claude subagent on 2026-10-14. No author or third party was contacted.
- Downloaded texts were kept in the session scratchpad and are not committed.

## Search scope (log)
**INSPIRE-HEP REST API.**
- `refersto:arxiv:1309.2448` is not parsed (25 145 unrelated hits). The paper's record id (1253547) was resolved
  and `refersto:recid:1253547` was used instead.
  - It returned **47 citing records; all 47 were examined**.
  - Titles were read for all; abstracts for the 12 relevant records and for every record without an arXiv ID that
    has "rotating" in its title.
- Citers were also screened (titles plus selected abstracts) for:
  - 2510.11406 (3 records);
  - 2403.08387 (12);
  - 2401.02898 (28);
  - 2301.05243 (25).
- Author queries:
  - `a Kunz, Jutta and date>2023`;
  - `a Blazquez-Salcedo and date>2022`;
  - `a Khoo and a Kunz`;
  - `a Kleihaus and date>2023`;
  - Folomeev, Chew and Azad, each combined with "wormhole".
- Erratum query `t erratum and a Dzhunushaliev and a Kleihaus`: 0 hits.

**arXiv.**
- Export API query `abs:wormhole AND (abs:rotating OR abs:spinning) AND (abs:unstable OR abs:instability OR
  abs:stability)`, newest first, 60 results screened.
- Version metadata for 1309.2448 (v1 only), 2510.11406, 2603.12748 and 2401.02898.

**Crossref.** The record for 10.1103/PhysRevD.88.124028 has no `updated-by` or relation entries, and
`filter=updates:` returns 0.

**Web search.** One extended query; only the papers already listed came back.

**Not reached:**
- Google Scholar;
- the APS article page itself (paywalled);
- full texts of the two 2017 proceedings that have no arXiv ID;
- theses and talks;
- papers that neither cite 1309.2448 nor are indexed yet.

## Relevant works
"Read" says how much of each work was checked: full text, or abstract/metadata only.

| arXiv | What it establishes (relevant part) | Read |
|---|---|---|
| 1309.2448 (PRD 88, 124028) | 5D, full rotation: two purely imaginary modes merge at (\|J\|/A)_cr = 0.10063 (symmetric). Beyond it, "no further unstable modes within the type of modes investigated". arXiv v1 only, no erratum | source (S6) |
| 2301.05243 (PLB 848, 138349) | 4D, second order in slow rotation: a second mode appears. The paper *conjectures* that the modes "merge and disappear" | full text (searched) |
| 2403.08387 (PRD 109, 124051) | 4D, detailed second-order slow rotation; branches truncated at their crossing. Reads the 5D result as "merge … and disappear beyond this critical value" | full text (searched) |
| 2510.11406 v2 (PLB 875, 140335) | Charged wormholes: the merged modes continue as a complex pair. *Conjectures* the same for rotating chargeless wormholes, "can only be confirmed by a non-perturbative study in the future" | full text |
| 2401.02898 (PRD 109, 084013) | 4D rapid rotation, spectral method: no instability in the M_z = 2, 3 sectors. The radial (M_z = 0) sector is not treated | abstract |
| 2509.22118 (Universe 11, 325) | Slow-rotation ℓ = 2 equations; no post-merger radial result | abstract |
| 2603.12748 | 5D equal-J backgrounds with large asymmetry. Says investigating the stability of these wormholes "is essential", i.e. still open | full text (searched) |
| 2607.18399 | Static charged wormholes: a non-radial polar ℓ = 2 instability | abstract |
| no arXiv ID (INSPIRE 1600650; 2017 proceedings by Dzhunushaliev et al.) | Restate 1309.2448: "when the rotation is sufficiently fast, the radial instability disappears" | INSPIRE abstract |

## Answers
**Complex continuation past the merger: not found as a computed result.**
- In the searched scope (arXiv- and INSPIRE-indexed literature up to 2026-10-14), no paper computes complex
  unstable modes, or any instability past the merger, for these 5D wormholes. None re-analyses the perturbation
  equations of 1309.2448.
- The behaviour **has been conjectured in print**. 2510.11406 conjectures it by analogy with the charged case.
  There "rotating" refers mainly to 4D, but that paper cites the 5D merger as its motivation.
- Several papers read 1309.2448 as showing stabilisation past J_cr: 2301.05243, 2403.08387 and the 2017
  proceedings. 1309.2448 itself limited its statement to the purely imaginary modes it searched.

**Erratum or typo in V₁₀: none found.**
- arXiv has v1 only, and Crossref and INSPIRE list no erratum.
- None of the citing full texts that were read (2510.11406, 2301.05243, 2403.08387, 2603.12748) mentions a
  misprint.
- Our evidence for the typo is our own re-derivation (`X2_s6a_5d.md`, Step 1). No one else has checked it.

**Wording consequence.**
- "First" is **not** used. The searched scope is limited, and the result answers a published conjecture rather
  than an unanticipated question.
- Permitted wording: "We found no prior computation of … in the literature indexed by INSPIRE and arXiv up to
  14 Oct 2026."

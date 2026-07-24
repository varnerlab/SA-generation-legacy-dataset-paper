# Handoff: npj Revision, Exact-Sampler Development, and Reviewer Audit

**Date:** 2026-07-24  
**Branch:** `revision/npj-major-revision`  
**Manuscript:** *Validated Synthetic Patient Generation for Small Longitudinal Cohorts: Coagulation Dynamics Across Pregnancy*

## Current Situation

The npj paper has been accepted. After acceptance, reviews of the underlying Stochastic Attention paper exposed an important mathematical simplification: the latent weighted Hopfield target used in the npj paper is exactly a finite Gaussian mixture and can be sampled ancestrally without Langevin dynamics.

The agreed strategy is **not** to reposition or reopen the accepted clinical paper. Instead, treat the exact sampler as a constructive outcome of extending the sampling-correctness and privacy studies requested by the npj reviewers.

The accepted ULA-generated cohort and all clinical analyses remain unchanged.

## Exact-Sampler Identity

For

\[
E_r(\xi)
=
\frac{1}{2}\|\xi\|^2
-
\frac{1}{\beta}
\log\sum_{k=1}^{K}r_k\exp(\beta m_k^\top\xi),
\]

the latent target is

\[
\pi_r(\xi)
=
\sum_{k=1}^{K}q_k\,
\mathcal N(m_k,\beta^{-1}I),
\qquad
q_k
\propto
r_k\exp\!\left(\frac{\beta}{2}\|m_k\|^2\right).
\]

Because the manuscript uses unit-normalized PCA memories,

\[
q_k=\frac{r_k}{\sum_jr_j}.
\]

An exact latent draw is therefore:

1. Select a stored patient according to the multiplicity weights.
2. Add isotropic Gaussian noise with covariance \(\beta^{-1}I\).
3. Apply the existing direction normalization, empirical PCA-magnitude draw, PCA reconstruction, and de-standardization.

The decoded assay-space distribution is a transformed mixture, but remains exactly sampleable through this pushforward.

## Preliminary Non-Mutating Diagnostic

A temporary Julia diagnostic compared the accepted ULA cohort against 100 independently generated exact-sampler cohorts. No diagnostic code or output was committed.

| Metric | Accepted ULA cohort | Exact sampler median | Exact 5th--95th percentile |
|---|---:|---:|---:|
| Marginal MRE | 1.235% | 1.38% | 0.98--1.84% |
| Cross-visit Frobenius error | 0.641 | 0.588 | 0.541--0.631 |
| Median novelty | 0.514 | 0.511 | 0.494--0.532 |
| Fraction novelty > 0.2 | 1.00 | 1.00 | 0.98--1.00 |
| Median DCR | 13.784 | 13.994 | 13.598--14.343 |

A seed-42 exact cohort evaluated through the existing BZ2012 harness produced:

- pooled mechanistic overlap: `0.892`
- KS \(D\): `0.0596`
- KS \(p\): `0.2726`

Interpretation: the clinical behavior appears to survive, while Langevin dynamics is one implementation rather than a necessary sampler.

## Agreed Camera-Ready Strategy

Keep the update contained:

- Keep the title.
- Keep the accepted ULA cohort.
- Do not rerun the main clinical, subgroup, simulation, or downstream-utility studies.
- Add an exact ancestral reference to the existing sampling-correctness/privacy analysis.
- Add no new main-text figure.
- Correct claims about interpolation, target interpretation, beta, multiplicity weights, and computational cost only where required.
- Notify the handling editor with a short note before incorporating final changes.
- Do not mention the confidential review from the other venue.

The detailed implementation and manuscript plan is:

- [`docs/superpowers/plans/2026-07-24-exact-sampler-camera-ready-update.md`](../superpowers/plans/2026-07-24-exact-sampler-camera-ready-update.md)

## npj Reviewer Audit

Source reviews:

- `peer-review-feedback/Reviewer-1-npj.docx`
- `peer-review-feedback/Reviewer-2-npj.docx`

Response:

- `peer-review-feedback/response-to-reviewers.md`
- `peer-review-feedback/response-to-reviewers.docx`

All 20 reviewer items appear in the response letter.

- **18/20 are substantively addressed.**
- **R1.2 and R2.3 are only partially resolved.** They are the same shared concern about the BZ2012 model being calibrated on the same cohort.

### Remaining R1.2/R2.3 contradiction

The Discussion correctly says the ODE agreement is consistency with a cohort-calibrated model rather than fully independent validation:

- `paper/sections/discussion.tex`, around lines 191--195, violet `\rboth{...}`.

However, Results still calls it independent:

- `paper/sections/results.tex`, around line 276:
  `an independent mechanistic model that knows nothing about the SA generation process?`
- `paper/sections/results.tex`, around lines 328--330:
  `An independent mechanistic model ... confirmed ...`

Required fix:

- Remove the erroneous question mark.
- Replace both independence claims with accurate cohort-calibrated wording.
- Wrap the corrections in `\rboth{...}` because they jointly address R1.2 and R2.3.

## Color-Coding Audit

The source color system is correct and active:

- `\rone{...}`: blue, Reviewer 1
- `\rtwo{...}`: red, Reviewer 2
- `\rboth{...}`: violet, both
- `\reviewmodetrue` is enabled.

Definitions:

- `paper/main.tex`, around lines 33--44
- `paper/supplementary.tex`, around lines 28--39
- `arxiv/main.tex`, around lines 32--41

`paper/sections` and `arxiv/sections` are synchronized.

Representative pages of `paper/main.pdf` were rendered and visually inspected. Blue, red, and violet are clearly visible and legible in:

- Abstract
- PCA ablation
- simulation/generalizability section
- extrapolation and privacy sections
- downstream-utility caveats
- tail/generalizability/PCA limitations
- Supplementary PCA, hull, MIA, and sampler-correctness material

No layout defects were observed on the inspected revised pages.

## Critical Packaging Defect

`paper/main.pdf` is current and contains the revised main text plus Supplement.

`paper/supplementary.pdf` is stale. It was last regenerated in April and does not contain the new PCA ablation, convex-hull, membership-inference, or sampling-correctness sections.

Cause:

```make
$(SUPP).pdf: $(SUPP).tex
```

in `paper/Makefile` does not depend on `sections/supplementary.tex`.

Required fix:

```make
$(SUPP).pdf: $(SUPP).tex sections/supplementary.tex
```

Apply the equivalent dependency fix to `arxiv/Makefile`, then rebuild and visually verify the standalone Supplement.

## Response-Letter Cleanup

1. Reviewer 1's original document numbers comments `1--9`, then `11--13`; it skips 10. The response letter renumbers the last three as R1.10--R1.12. Preserve or explicitly explain the original numbering to avoid confusion.
2. The R1.1 `Changes` line promises “tables/figure,” but the delivered PCA ablation is a Results paragraph plus Supplementary table. Correct the pointer.
3. R1.2/R2.3 currently claim that independence wording was tempered in Results; complete the Results correction described above.
4. The response date remains `[fill on submission]`.

## Relationship to the New Exact-Sampler Plan

The exact-sampler update will supersede parts of the current:

- “exact Metropolis-adjusted sampler” wording
- strongly multimodal/sharply peaked basin explanation
- interpretation of interpolation
- \(O(TKd)\) computational limitation
- beta/multiplicity explanation

The exact reference should be presented as an extension of R1.5:

> Extending the reviewer-requested sampling-correctness analysis to an analytic reference revealed that the unit-memory Hopfield target is exactly a finite isotropic Gaussian mixture. Exact ancestral draws passed through the same magnitude and decoding pipeline reproduced the ULA cohort's evaluated behavior, including the membership-inference signal.

## Recommended Restart Order

1. Read this handoff and the exact-sampler update plan.
2. Implement and test the general exact sampler.
3. Run the fixed-seed full-harness exact reference, 100-cohort stability study, and repeated-split exact MIA.
4. Freeze the resulting CSVs and decide whether “reproduced” is supported across the complete harness.
5. Fix the two R1.2/R2.3 Results sentences.
6. Apply the contained exact-sampler manuscript changes to both `paper/` and `arxiv/`.
7. Update the response letter and create the short handling-editor note.
8. Fix the Supplement Makefile dependencies and rebuild all PDFs.
9. Render and visually inspect every changed page.
10. Fill the response date and prepare marked and clean final versions.

## Repository State at Handoff

Before adding this handoff, the branch was clean and synchronized with:

`origin/revision/npj-major-revision`

New uncommitted planning artifacts:

- `docs/superpowers/plans/2026-07-24-exact-sampler-camera-ready-update.md`
- `docs/handoffs/2026-07-24-npj-exact-sampler-review-handoff.md`

No manuscript, analysis code, response letter, or PDF was modified during the exact-sampler discussion or reviewer audit.

## Environment Note

The two reviewer DOCX files were fully extracted with Pandoc. LibreOffice/`soffice` is not installed, so the DOCX files themselves could not be rendered to page images. Reviewer 2's OOXML comments part is empty; there are no hidden Word margin comments beyond the visible review text.

The manuscript PDF was rendered with Poppler and visually inspected successfully.

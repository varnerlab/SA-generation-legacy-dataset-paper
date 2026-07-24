# Exact-Sampler Camera-Ready Update Plan

**Manuscript:** *Validated Synthetic Patient Generation for Small Longitudinal Cohorts: Coagulation Dynamics Across Pregnancy*

**Objective:** Incorporate the newly recognized closed-form representation of the latent Stochastic Attention (SA) target as a focused strengthening of the reviewer-requested sampling-correctness and privacy analysis. Preserve the accepted clinical paper, its published ULA cohort, title, main results, and validation framework.

## Positioning Decision

The update will make one precise distinction:

- The **accepted cohort** was generated with the reported 2,000-step ULA implementation.
- The **SA-induced latent target** is analytically sampleable as a finite Gaussian mixture.
- An **exact ancestral reference** will test whether the ULA cohort's fidelity, novelty, mechanistic behavior, and privacy signal reflect the intended target distribution.

The result will be presented as an outcome of extending Reviewer 1's requested sampling-correctness and privacy investigation, not as a new method or a replacement study.

## Scope Guard

Keep the update small enough for final-manuscript handling:

- Keep the current title.
- Do not replace or regenerate the accepted synthetic cohort.
- Do not rerun the main clinical, subgroup, simulation, or downstream-utility analyses.
- Add no main-text figure.
- Add one analytic reference to the existing sampler-correctness evidence.
- Correct only claims directly affected by the identity.
- Do not import the unrelated transformer, DDPM, or theoretical issues raised in the underlying SA paper's review.
- Notify the handling editor before incorporating the final changes.

## Core Mathematical Statement

For the weighted Hopfield energy already given in the Methods,

\[
E_r(\xi)
=
\frac{1}{2}\|\xi\|^2
-
\frac{1}{\beta}
\log\sum_{k=1}^{K} r_k\exp(\beta m_k^\top\xi),
\]

the latent target is

\[
\pi_r(\xi)
\propto
\exp[-\beta E_r(\xi)]
=
\sum_{k=1}^{K}
r_k\exp\!\left(\frac{\beta}{2}\|m_k\|^2\right)
\exp\!\left[-\frac{\beta}{2}\|\xi-m_k\|^2\right].
\]

Therefore,

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

For the manuscript's unit-normalized PCA memories, \(q_k=r_k/\sum_jr_j\). Exact latent sampling is consequently:

1. Draw \(k\sim\mathrm{Categorical}(q)\).
2. Draw \(\xi=m_k+\beta^{-1/2}\varepsilon\), with \(\varepsilon\sim\mathcal N(0,I)\).
3. Apply the existing published pushforward unchanged: normalize \(\xi\) to a direction, draw an empirical PCA magnitude, reconstruct through PCA, and de-standardize.

The final decoded patient distribution is a transformed mixture rather than a Gaussian mixture in assay space, but it remains exactly sampleable by applying the same pushforward to exact latent draws.

---

## Phase 1: Add and Verify the Exact Latent Sampler

### Task 1.1: Implement a general exact sampler

**Modify:** `code/src/Compute.jl`

Add:

```julia
exact_sample(X, n; β, weights=ones(size(X, 2)), rng=Random.default_rng())
```

Requirements:

- Support non-unit memory columns with component probabilities
  `weights[k] * exp(β * dot(X[:,k], X[:,k]) / 2)`.
- Compute probabilities in log space for numerical stability.
- Return both the sampled latent matrix and selected component indices.
- Validate positive `β`, positive weights, finite inputs, and compatible dimensions.
- Do not change `sample`, `mala_sample`, `weighted_sample`, or any frozen generation result.

### Task 1.2: Add mathematical and statistical tests

**Create:** `code/experiments/test_exact_sampler.jl`

Tests:

1. **Density identity:** At random probe points, compare the Hopfield log target with the Gaussian-mixture log density up to one constant.
2. **Unit-norm reduction:** Confirm that equal-norm memories yield component probabilities proportional to the supplied multiplicity weights.
3. **General-norm weights:** Confirm the \(\exp(\beta\|m_k\|^2/2)\) correction.
4. **Moment check:** With a large synthetic draw, verify the empirical mixture mean and covariance against their analytic values.
5. **Seed reproducibility:** Identical RNG seeds must produce identical components and latent samples.

**Gate:** No manuscript editing begins until all five checks pass.

---

## Phase 2: Add the Exact Reference to the Existing Diagnostic Package

### Task 2.1: Generate one protocol-matched exact cohort

**Create:** `code/experiments/run_exact_sampler_reference.jl`

**Read:**

- `code/data/full_longitudinal_memory.jld2`
- `code/data/mcmc_ladder_results.csv`
- existing canonical real and ULA cohorts

**Write:**

- `code/data/exact_sampler_reference.csv`
- `code/data/exact_sampler_replicates.csv`

Canonical comparison:

- \(N=100\)
- \(\beta=\beta^*\)
- seed 42
- same unit memory, PCA model, empirical-magnitude distribution, reconstruction, and standardization as the published pipeline
- separate seeded random streams for latent sampling and magnitude draws
- where possible, reuse the magnitude sequence employed by the sampler ladder so that the direction sampler is the only varying element

Evaluate through the same harness as V0--V3:

- pooled median marginal MRE
- cross-visit correlation Frobenius error
- TF-only mechanistic cloud overlap and KS statistic
- median cosine novelty and fraction above 0.2
- median synthetic-to-real DCR
- real-to-real DCR reference

The output row should use the label `Exact` and sampler description `analytic ancestral`.

### Task 2.2: Quantify Monte Carlo stability

Generate 100 exact cohorts of \(N=100\) under independent seeds.

Report median and 5th--95th percentile for:

- MRE
- cross-visit Frobenius error
- novelty
- fraction above the novelty threshold
- DCR

The fixed-seed cohort remains the protocol-matched table row; the replicate distribution establishes that any agreement is not a favorable-seed artifact. The mechanistic ODE harness need only be run for the fixed-seed cohort unless runtime permits a predeclared smaller replicate set.

### Task 2.3: Extend the repeated-split membership-inference check

**Modify:** `code/experiments/run_mcmc_ladder_mia.jl`

Add scheme `:Exact`:

- Refit PCA and unit memories on each 15-patient training split exactly as V0/V3 do.
- Draw 100 exact latent samples at the split-specific \(\beta^*\).
- Apply the identical magnitude and decoding pipeline.
- Use the identical 40-split partition sequence and MIA scoring.

**Write:** append `Exact` records to `code/data/mcmc_ladder_mia.csv`, or write a dedicated `code/data/exact_sampler_mia.csv` if preserving the existing CSV byte-for-byte is preferable.

Report mean AUC, standard deviation, range, and fraction of perfect-separation splits.

### Task 2.4: Interpret outcomes without a pass/fail filter

Report all results. Use “reproduced” or “closely matched” only if the exact fixed-seed result and replicate distribution support it across the complete harness.

Predeclare the following as materially close for drafting purposes:

- MRE: within 0.75 percentage points
- cross-visit Frobenius error: within 0.10
- mechanistic overlap: within 0.05
- median novelty: within 0.05
- DCR: within 0.75 standardized units
- repeated-split MIA AUC: within 0.05

These thresholds control the language, not which results are disclosed. If any metric falls outside them, describe the discrepancy as finite-time or implementation dependence and do not claim sampler equivalence.

### Task 2.5: Optional conditioned-sampling sanity check

This task is not required for the camera-ready update.

If included, verify only that empirical exact-sampler component frequencies follow the prescribed multiplicity probabilities for the PCOS and Developed-PE weights. Do not reopen the full conditioned clinical analysis unless this check contradicts the analytical result.

---

## Phase 3: Make Targeted Manuscript Edits

Apply every prose edit to both:

- `paper/sections/*.tex`
- `arxiv/sections/*.tex`

The two trees are currently synchronized and must remain so.

### Task 3.1: Abstract - one contained clarification

**Modify:** `sections/abstract.tex`, current lines 4--10.

Keep SA, the Hopfield energy, and the fact that the accepted cohort was generated with Langevin dynamics. Replace the unsupported implication that Langevin dynamics is necessary or that every draw literally interpolates several patterns.

Add one compact clause stating that an exact analytic reference reproduced the ULA cohort's evaluated behavior. Do not add equations or turn the abstract into a methods correction.

Target meaning:

> SA defines a memory-centered latent distribution; the reported cohort was generated with ULA, and an exact ancestral reference reproduced its evaluated behavior.

### Task 3.2: Introduction - correct the method contrast without repositioning the paper

**Modify:** `sections/introduction.tex`, current lines 59--76 and 86--94.

Required corrections:

- Replace “rather than estimating a parametric distribution” with the narrower and accurate claim that SA avoids fitting a full covariance model or neural generator.
- Describe stored patients as mixture centers in the reduced latent space.
- Replace “Langevin dynamics generates samples that interpolate between patterns” with language about sampling a memory-centered distribution whose component overlap yields new directions.
- Describe multiplicity weighting as exact reweighting of stored-patient components at inference time.
- Preserve the small-\(n\), longitudinal, rank-deficiency, and subgroup-motivation narrative.

### Task 3.3: Methods - add the exact representation immediately after the energy

**Modify:** `sections/methods.tex`, immediately after Eq. `hopfield_energy` and before the ULA equation.

Add:

1. The completed-square derivation.
2. The general component probabilities \(q_k\).
3. The unit-norm simplification \(q_k\propto r_k\).
4. A statement that the accepted cohort was generated using the reported ULA implementation.
5. A statement that the exact sampler was added afterward as an analytic reference using the identical normalization, magnitude, and decoder pipeline.

Then revise the beta paragraph:

- \(\beta\) controls component variance/overlap and the angular novelty-fidelity trade-off.
- \(r_k\) controls component probability exactly for unit memories.
- \(\beta^*\) is an entropy-based geometric operating point, not evidence that Langevin dynamics is required.

Retain the ULA equation because it documents how the accepted cohort was actually generated.

### Task 3.4: Results - strengthen the privacy attribution

**Modify:** `sections/results.tex`, Re-identification Risk, current lines 186--215.

Replace “an exact Metropolis-adjusted sampler” with “Metropolis-adjusted sampling and an analytic exact reference.” MALA is asymptotically target-correct but is not itself a finite-draw exact sampler.

Add the exact-reference MIA, DCR, and full-harness results in one short paragraph. The conclusion should be:

> The exact reference reproduced the same membership signal, locating it in the target distribution rather than in anchor initialization or ULA discretization.

Avoid “K sharply peaked non-overlapping basins” unless the exact-mixture analysis directly supports that wording. Prefer the exact statement that every stored patient is a component center and that the decoded target retains measurable proximity to those centers.

### Task 3.5: Supplement - extend the existing correctness section

**Modify:** `sections/supplementary.tex`, current lines 580--658.

Changes:

- Rename “four-rung sampler ladder” to “sampler-correctness comparison” or retain V0--V3 as the ladder and call `Exact` an analytic reference outside the rung numbering.
- Add a compact completed-square derivation.
- Add the `Exact / analytic ancestral` row to Table `stab:sampling-ladder`.
- Report exact-reference MIA in the endpoint column.
- State the 100-cohort stability interval in the table note or adjacent paragraph.
- Keep the existing V0--V3 figure unless adding the exact point is visually trivial. A new figure is out of scope.
- Update the caption so it does not imply MALA is the strongest available correctness reference.

Also revise the earlier ULA-vs-MALA caption, current lines 184--209:

- MALA acceptance supports negligible ULA discretization bias.
- The analytic sampler, introduced later in the supplement, provides the direct target reference.

### Task 3.6: Discussion - turn the development into the computational win

**Modify:** `sections/discussion.tex`.

Required edits:

1. Replace statements that SA avoids a parametric distribution with the accurate distinction: it avoids estimating a full high-dimensional covariance or training a neural model, while retaining the observed patients as mixture centers.
2. Recast “interpolation” as transformed memory-centered sampling and component overlap.
3. Replace the current \(O(TKd)\) limitation:
   - the reported ULA implementation costs \(O(TKd)\);
   - the exact ancestral implementation removes the sequential \(T\) factor for future generation;
   - PCA fitting and decoding remain.
4. Strengthen the privacy interpretation using the explicit patient-centered mixture.
5. Preserve the accepted limitations concerning small \(K\), tails, PCA linearity, independence of the ODE validation, and non-DP status.

### Task 3.7: Conclusion - one sentence, no new positioning

**Modify:** `sections/conclusion.tex`.

Add at most one sentence:

> An analytic sampler for the same latent target reproduced the ULA implementation's evaluated behavior and provides a faster route for future use.

Do not change the title, clinical take-home message, or four-level validation summary.

---

## Phase 4: Update the Review Package and Prepare the Editor Note

### Task 4.1: Response letter

**Modify:** `peer-review-feedback/response-to-reviewers.md`

Update:

- General response, paragraph describing the privacy/sampling-correctness investigation.
- R1.5 response and its `Changes` line.
- R1.9 response only if it currently enumerates the relocated theoretical material.

Suggested framing:

> Extending the reviewer-requested sampling-correctness analysis to an analytic reference revealed that the unit-memory Hopfield target is exactly a finite isotropic Gaussian mixture. Exact ancestral draws passed through the same magnitude and decoding pipeline reproduced the ULA cohort's evaluated behavior, including the membership-inference signal. This both strengthens the attribution of the privacy behavior to the target law and provides a faster implementation for future applications; it does not alter the accepted cohort or the clinical conclusions.

Regenerate:

- `peer-review-feedback/response-to-reviewers.docx`

### Task 4.2: Handling-editor note

**Create:** `peer-review-feedback/editor-note-exact-sampler.md`

Keep it under 200 words. State:

- the finding arose while extending the reviewer-prompted correctness analysis;
- the accepted ULA cohort and all clinical conclusions are unchanged;
- the addition consists of a short identity, one analytic-reference comparison, and corrected computational language;
- ask permission to incorporate it in the final manuscript.

Do not mention another venue's confidential review.

---

## Phase 5: Verification and Submission Sequence

### Task 5.1: Numerical audit

Verify that every manuscript number is generated from the committed CSVs:

- exact fixed-seed full-harness row
- exact replicate intervals
- exact repeated-split MIA summary
- unchanged V0--V3 values

Search for hand-entered values and reconcile them with the canonical outputs.

### Task 5.2: Claim audit

Search both manuscript trees and the response letter for:

- `exact Metropolis`
- `rather than estimating a parametric distribution`
- `Langevin` near `required`, `necessary`, or exclusive mechanism language
- `interpolate between patterns`
- `phase transition`
- `sharply peaked basins`
- `O(TKd)` and `2,000`
- `four-rung`

Every remaining instance must be either historically descriptive of the accepted ULA cohort or mathematically compatible with the exact-mixture representation.

### Task 5.3: Cross-copy and build verification

1. Confirm `paper/sections` and `arxiv/sections` remain synchronized.
2. Build `paper/main.pdf` and `paper/supplementary.pdf`.
3. Build the arXiv version.
4. Check for undefined references, duplicate labels, overfull boxes, table overflow, and changed figure count.
5. Render the changed PDF pages and visually inspect the Methods derivation, Results paragraph, Discussion cost paragraph, and Supplement table.
6. Generate both a marked reviewer copy and a clean black-text final copy using the existing `reviewmode` toggle.

### Task 5.4: Submission order

1. Complete the exact analysis and freeze its outputs.
2. Draft the limited manuscript changes.
3. Send the short editor note with a concise change summary.
4. Incorporate any editor-requested boundary changes.
5. Build and visually verify the final package.
6. Fill the response-letter date and deliver the final manuscript, supplement, response DOCX, and change summary.

---

## Recommended Commit Sequence

1. `analysis(exact): add and verify analytic sampler for weighted Hopfield target`
2. `analysis(exact): evaluate exact reference on fidelity, privacy, and mechanistic harness`
3. `revision(exact): add analytic reference to Methods, Results, and Supplement`
4. `revision(exact): correct interpretation, cost, and privacy language`
5. `revision(response): document exact-sampler outcome and editor note`
6. `revision(final): synchronize mirrors, rebuild, and verify camera-ready package`

## Expected Manuscript Footprint

- Title: unchanged
- Abstract: one sentence revised, one short clause added
- Introduction: two localized paragraphs corrected
- Methods: one short derivation paragraph plus interpretation edits
- Results: one short exact-reference paragraph
- Discussion: sampling interpretation, cost, and privacy paragraphs corrected
- Conclusion: at most one sentence
- Supplement: one derivation, one table row, one stability sentence
- Main figures: unchanged
- Main clinical results: unchanged

This footprint keeps the development proportional to its role: an analytic clarification and computational strengthening of the already accepted reviewer-requested diagnostic work.

# npj Minor Revision (Round 2): Design

Date: 2026-11-13. Due: Monday morning, 2026-11-16. The author uploads the package to SNAPP.

## Why

The npj Systems Biology and Applications paper received a minor revision. Reviewer 1
recommends publication. Reviewer 3 has three comments
(`peer-review-feedback/npj-minor-comments-10-28-2026.png`):

1. Gaussian mixture and copula models match stochastic attention (SA) on several metrics.
   Clarify the advantage of SA, given that SA is exactly a finite Gaussian mixture. Discuss
   the relationship between SA and Gaussian kernel-density estimation (KDE).
2. The membership-inference AUC is about 0.97. Clarify what this means for sharing synthetic
   data. Separate exploratory computational use from privacy-preserving data release.
3. Only 54% of profiles pass the biological criteria, and 24% pass both the biological and
   strict mechanistic criteria. Emphasize this limitation. Discuss whether
   positivity-preserving transforms or constrained reconstruction are needed.

The goal is to close all three comments in one round, with no new requests.

## The argument for each comment

**Comment 1.** With equal multiplicity weights and unit-normalized memories, the SA sampling
distribution in Eq. (exact_mixture) is an isotropic Gaussian KDE on the 23 unit-normalized
PCA profiles, with bandwidth h = beta^(-1/2). At beta* = 2.94, h = 0.58. The E1 KDE baseline
was a different KDE: raw PCA coordinates, full empirical covariance, and Scott's factor
K^(-1/(d+4)) = 0.87. SA differs from a textbook KDE in three ways:

- The bandwidth comes from the attention-entropy transition, not from a large-sample rule.
- Multiplicity weights steer the mixture toward subgroups, and the participation ratio gives a
  diversity check.
- SA is defined by an energy. The energy can take a constraint term: E(x) + lambda C(x).

The author's framing: lambda = 0 gives the exact mixture. That case has a known answer, so it
checks the Langevin sampler, much as the ODE xdot = -kx checks an Euler solver. With a
non-quadratic C, the target is no longer a Gaussian mixture, and a sampler is needed.

Two limits on this claim must appear in the text:

- A quadratic C keeps the target a Gaussian mixture. Only a non-quadratic C breaks the
  mixture form.
- Rejection sampling from the mixture can sample a constrained target exactly. For a hard
  constraint, keep draws that pass. For a soft constraint C >= 0, keep each draw with
  probability exp(-beta lambda C). Rejection becomes inefficient when few draws pass. The
  energy form with Langevin does not need many draws to pass. Do not claim that a Gaussian
  mixture "cannot" handle constraints.

**Comment 2.** The source data are de-identified. The membership test can tell only whether
an already-known profile was stored, and it cannot identify a person. The membership signal
persisted across samplers and across beta (median distance 13.8 to 14.3 at beta/beta* down to
0.05, versus 15.9 between real patients; `code/data/beta_privacy_sweep.csv`). The new text
states that the synthetic cohorts support exploratory modeling. They are not a privacy
method. A cohort generated from data that cannot be shared should be shared under the same
terms as the source data. Releasing data from identifiable or sensitive cohorts would need a
formal privacy method, such as differential privacy, which was not tested. Do not suggest
lowering beta as a privacy fix.

**Comment 3.** All 46 failures are negative estradiol or progesterone values. Neither hormone
is a BZ2012 model input. The strict mechanistic check requires all 15 ratios (5 TGA features
by 3 visits) to fall inside each feature's real 5th-95th percentile band. Real patients may
also often fail this check. Analysis W2 measures how often. Positivity can be enforced by a
log transform, by rejection, or by a constraint term. The constrained-sampling demonstration
(W3) tests the constraint term.

## Approaches considered

- A. Prose only. This is fast and low risk, but it makes the constraint claim without a
  demonstration.
- B. Demonstration only. This shows that the method works, but it leaves comments 1 and 3
  without the discussion the reviewer asked for.
- **A + B (chosen by the author).** The prose answers all three comments. A small
  constrained-sampling demonstration goes in the Supplement and is discussed in the new
  Discussion paragraph.

## Work items

### W1. KDE bandwidth rules on the beta axis (comment 1)

Use the 23 unit-normalized PCA profiles in `code/data/full_longitudinal_memory.jld2`
(key `X̂`, d = 18). Let sigma be the square root of the mean per-coordinate variance of these
23 vectors. Compute these bandwidths:

- Scott: h = sigma K^(-1/(d+4)).
- Silverman: h = sigma (4/(d+2))^(1/(d+4)) K^(-1/(d+4)).
- Leave-one-out likelihood: the h that maximizes
  sum_i log[(1/(K-1)) sum_{j != i} N(m_i; m_j, h^2 I)].

Report each h, its equivalent beta = 1/h^2, and beta/beta*. Generate 100 profiles with the
exact sampler at each equivalent beta. Use the full canonical pipeline: normalize the
direction, draw a magnitude from the 23 empirical norms, and decode. Use seed 42. Report the
same metrics as the beta sweep: median MRE, cross-visit Frobenius error, pooled TF-only
mechanistic overlap, median novelty, and median distance to the nearest real profile. Include
beta* as a reference row.

Script: `code/experiments/run_kde_bandwidth_check.jl`. Output:
`code/data/kde_bandwidth_check.csv`. New supplement table: `stab:kde-bandwidth`.

### W2. Real-patient pass rate on the strict mechanistic check (comment 3)

In `code/experiments/run_convex_hull_analysis.jl`, `real_ratios` and the per-feature bands
already exist. Count the real patients whose 15 ratios all fall inside the bands. Use the same
in-sample bands that were applied to the synthetic profiles. Also run the biological check on
the 23 real patients; 23/23 should pass. Report both counts next to the synthetic counts (54,
41, 24) in the existing hull and plausibility supplement table. If leave-one-out bands are
cheap, add them as a secondary row. Do not change the canonical synthetic numbers. That script
asserts 54/100.

### W3. Constrained positivity demonstration (comments 1 and 3, Supplement only)

Purpose: show that a non-quadratic constraint term removes negative hormone values, and check
the constrained Langevin sampler against an exact reference.

- **Magnitude first.** The canonical pipeline draws each profile's magnitude r independently,
  uniformly from the 23 empirical PCA norms (`Patient.generate_patients_rescaled`). Draw r
  first for each profile. With r fixed, each decoded hormone entry is linear in the unit
  direction: x_j = a_j + r c_j^T (xi/||xi||). There are 6 entries: estradiol and
  progesterone at 3 visits. Take a_j and c_j from the PCA mean, the loadings, and the
  destandardization parameters. Verify this linear form against `decode_sample` numerically.
- **Constraint.** C(xi; r) = sum_j [max(0, -x_j / s_j)]^2, where s_j is the real-cohort SD of
  entry j. This C is non-quadratic. The builder may use a smooth softplus version instead, but
  must state the choice.
- **Target.** pi_lambda(xi | r) is proportional to exp(-beta*[E(xi) + lambda C(xi; r)]),
  with uniform multiplicities and beta* = 2.94.
- **Langevin.** Use the ULA update in Eq. (ula) with the gradient lambda grad C added. Use
  T = 2000, and reduce the step size if large lambda needs it. Report the step size used.
- **Exact reference at each lambda.** Draw xi with `exact_sample`, then keep the draw with
  probability exp(-beta* lambda C(xi; r)). Hard rejection (keep only draws with all 6 entries
  at or above 0) is the limit as lambda goes to infinity. Record acceptance rates.
- **Arms.** lambda in {0, 1, 10, 100}, plus hard rejection. Use 100 profiles per arm and
  seed 42. For the Langevin-versus-reference agreement, use 1000 draws and compare the decoded
  hormone distributions (KS statistic) and the metrics below.
- **Metrics per arm.** Report the biological pass rate (same definition as the 54/100),
  median MRE overall, hormone-only MRE, cross-visit Frobenius error, pooled TF-only
  mechanistic overlap, strict mechanistic pass rate, joint pass rate, median novelty, and the
  acceptance rate for the reference arms.
- **Expected outcome (not a requirement).** Pass rates approach 100% as lambda grows, Langevin
  agrees with the reference at each lambda, and the fidelity metrics change little. Report
  what happens, including any loss in fidelity.

Script: `code/experiments/run_constrained_positivity_demo.jl`. Output:
`code/data/constrained_positivity_demo.csv`. New supplement table
`stab:constrained-positivity`. The supplement holds results only. The method goes in the main
Methods (see W4). No main-text figure.

**Seeds.** Reuse the canonical seeds (latent 42, magnitude 424242) and the canonical
per-profile magnitudes. The lambda = 0 arm must then reproduce 54/100, and the beta* row of W1
must reproduce the existing exact-sampler row. Treat either mismatch as a bug.

### W3b. Log-transform check (comment 3, Supplement only; added 2026-11-13 at the author's request)

The constraint demonstration (W3) removed negative hormone values. It also raised hormone
error from about 5% to about 30%, because removing the negative tail pushes the hormone means
up. The author asked to also test a positive representation, so the answer is complete.

All 72 assays are strictly positive at every visit. The hormones are right-skewed at visits 1
and 2 (skewness 2-4). Test two representations against the untransformed one:

- **Log hormones.** Take the natural log of the 6 estradiol and progesterone columns.
- **Log all assays.** Take the natural log of all 216 concatenated columns.

Each representation follows the canonical pipeline (`run_full_longitudinal.jl`). The log
transform is applied to the selected columns of the 23 x 216 real matrix. Then: z-score, PCA
at 95% variance, unit-normalized memories, its own beta* from the entropy-transition rule,
100 profiles from the exact sampler (latent seed 42), and magnitudes from its own PCA norms
(seed 424242). To decode, invert the PCA, invert the z-score, and exponentiate the logged
columns. Evaluate in original units against the real cohort with the same metrics as W3.
Compute novelty in each representation's own PCA space. The untransformed row must reproduce
`exact_sampler_reference.csv`.

Script: `code/experiments/run_log_transform_check.jl`. Output:
`code/data/log_transform_check.csv`. Combine with the W3 results in one supplement table of
ways to keep hormone values positive. The method goes in the main Methods paragraph.

### W4. Manuscript text

All new text follows the author's prose rule: short, plain, direct sentences that say exactly
what happened. No vague words, no "rather than" contrasts, no strung-together clauses. See
`HOUSE_STYLE.md` and the "Checks" section below.

- **Methods, after Eq. (exact_mixture).** Add one or two sentences: with equal weights this
  distribution is an isotropic Gaussian KDE on the unit-normalized profiles with bandwidth
  beta^(-1/2).
- **Methods, new `\rthree` paragraph in Baseline and Validation Framework (about 150-200
  words).** npj does not permit Supplementary Methods (technical check, 19 August 2026, item
  11, `peer-review-feedback/npj-changes.pdf`), so every method lives in the main Methods. This
  paragraph describes the three bandwidth rules (W1), the real-patient strict-check count
  (W2), the constraint term, the constrained Langevin update, and the rejection reference
  (W3). The supplement holds only the result tables and a short results text, framed like the
  round-1 "Supplementary Sampling and Privacy Results" section (my call).
- **Results, E1 paragraph.** Add one sentence stating that the KDE baseline differs from SA in
  coordinates, covariance, and bandwidth (`stab:kde-bandwidth`).
- **Results, plausibility paragraph.** Add the real-patient pass rate from W2. Replace "could
  prevent negative values" with a pointer to the W3 supplement result.
- **No new Discussion paragraphs (author, 2026-11-13: "we are already long").** Insert a few
  sentences into existing paragraphs. Target about 180 added words across the whole Discussion.
  Do not trim legacy or round-1 text to make room, because deletions are invisible in the
  marked manuscript.
- **Comment 1 insert (2-3 sentences, about 50 words).** Put it in the existing PCA and beta
  paragraph (`discussion.tex` lines 29-51), right after "Its sampler still required a choice
  of inverse temperature." Cover these points: with equal weights, SA is a Gaussian KDE on the
  unit-normalized profiles; beta sets the kernel width; the W1 comparison of beta* with the
  standard bandwidth rules; and the energy form also accepts constraint terms. No
  self-references to other sections.
- **Comment 3 insert (4-5 sentences, about 80 words).** Put it in the existing representation
  limits paragraph (lines 115-138), after the PCA linearity passage ends ("...clinical
  trajectories are generally linear.") and before "Cohort size also affects computational
  cost." The insert says the following. Linear decoding does not keep values positive.
  Forty-six profiles had a negative estradiol or progesterone value, and 24 of 100 passed both
  checks. Give the W2 real-patient count for the strict check. A constraint term removed the
  negative values, and the Langevin sampler matched exact rejection (W3). Rejection was
  efficient here because about half the draws passed. A constraint term matters when few
  draws pass. With such a term, the target is no longer a Gaussian mixture, and the lambda = 0
  case checked the sampler.
- **Discussion, final limitations paragraph (comment 2).** Expand the membership sentences
  (current lines 158-161) with the sharing statement from "Comment 2" above.
- **Data Availability, in both `paper/main.tex` and `arxiv/main.tex`.** New sentence: "The
  de-identified real patient data and the generated synthetic datasets are available in the
  project repository at [URL]." This matches what the public repository holds.
- **Abstract.** No change (my call; the limitation goes in Discussion and Results).

### W5. Colors and trees

- Add a `\rthree` macro in a new color (dark green, rgb 0.05,0.50,0.25) to `paper/main.tex`
  and `paper/supplementary.tex`. All round-2 changes use `\rthree`, including the Data
  Availability sentence.
- In the round-2 marked manuscript, only round-2 changes are colored. `\rone`, `\rtwo`, and
  `\rboth` render black (my call). The clean manuscript renders everything black.
- `paper/sections/` is the source. `arxiv/sections/` holds the same prose without color
  wrappers, plus a few intentional differences. Apply each edit to `arxiv/` by hand. Do not
  copy files between the trees.

### W6. Response letter

Create a new round-2 letter in `peer-review-feedback/` that follows the format and tooling of
`response-to-reviewers.md`/`.tex`. Thank Reviewer 1 briefly. For R3.1-R3.3, quote the
comment, answer it directly, summarize the evidence, and give the manuscript location. Use
the same terms as the paper.

### W7. Build and package

Build `paper/` and `arxiv/` with `make all`. Build a `paper-revision-2.2/` package that
mirrors `paper-revision-2.1/` (read its `README.txt`, `upload-checklist.md`, and
`scripts/build_latex_source_zip.py`). The package holds clean and marked manuscripts, the
supplementary PDF, the LaTeX source zip, and the rebuttal letter, with root-level v2.2 copies
named like the v2.1 files.

## Constraints and anti-goals

- Do not change the canonical cohort or any existing reported number. New numbers come only
  from W1-W3 and live in the new tables or CSVs.
- Do not claim that a Gaussian mixture cannot handle constraints. Do not suggest lowering
  beta as a privacy fix. Do not cite unpublished constrained-SA work (the MEng PAI-1 project
  or the compositional-SA pilot).
- Do not mention any companion SA methods manuscript under review at another venue. Do not
  add new citations to `Alswaidan2026SA`; the existing citation stays.
- The standalone supplement must contain no Supplementary Methods section (npj rule; also in
  `paper-revision-2.1/upload-checklist.md`).
- Do not present future work as a promise. State the limitation and the type of method that
  addresses it.
- Do not commit, push, or upload without the author's go-ahead.

## Checks before handoff

- Both trees build with 0 undefined references, in both the main document and the
  supplement. `pdftotext ... | grep '??'` finds nothing. `git diff --check` is clean.
- Every new inline number has a home in a table or CSV.
- New prose: no em dashes, no "rather than", no sentence starting with an acronym, median
  sentence about 16 words, nothing over about 35 words.
- Each response-letter number matches its CSV.
- List anything not verified.

## Schedule

Friday: W1-W3 computations. Saturday: W4 and the supplement. Sunday: W6 and W7, then the
author's read-through. Monday morning: the author submits.

## Calls made without the author (check these)

- Abstract unchanged.
- Round-2 marked manuscript colors only round-2 changes, in a new green `\rthree`.
- Data Availability sentence colored `\rthree`.
- Methods for W1-W3 go in a short main-text Methods paragraph, because npj forbids
  Supplementary Methods. The demo results stay in the supplement.
- Comments 1 and 3 go into existing Discussion paragraphs as short inserts. The author decided
  against new paragraphs. Comment 2 extends the existing final limitations paragraph.

## Left to the builder

The hinge or softplus form of C, the step size for large lambda, whether to add leave-one-out
bands in W2, table layouts, the response-letter file names, and the exact wording within the
rules above.

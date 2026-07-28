# Reviewer-Response Closure Audit

Audit date: 2026-07-28  
Closure verification: 2026-07-28  
Status: **Complete**

## Bottom line

All 20 labeled reviewer comments are now acknowledged and substantively
addressed.

- Fully or substantively closed: **20 of 20**
- Partially closed: **0 of 20**
- Omitted: **0 of 20**
- Remaining reviewer-response blockers: **none**

This document supersedes the earlier pre-fix audit, which identified six
partially closed comments. The affected analyses were rerun, all resulting
numbers were synchronized across the code artifacts, manuscript, supplement,
cover letter, and response letter, and the final PDFs were rebuilt and
inspected.

## Closure of the previously partial comments

| Comment | Final resolution | Status |
|---|---|---|
| Reviewer 1, major comment 2 | Mechanistic validation is consistently described as a same-cohort, generator-blind consistency check rather than independent validation. The residual dependence is disclosed in the abstract, introduction, results, discussion, methods, conclusion, cover letter, and response. | Closed |
| Reviewer 1, major comment 3 | The subgroup analysis is explicitly limited to reproducing patterns implied by the three PCOS and five Developed PE patients. Synthetic cohorts do not increase the inferential sample size and are restricted to exploratory simulation, hypothesis generation, power-analysis workflows, and workflow testing. | Closed |
| Reviewer 1, major comment 9 | The inverse-temperature explanation now separates the retrieval landscape from the exact sampling law. For unit-normalized memories, exact mixture probabilities are the multiplicity weights and are independent of beta; beta controls the within-component covariance. | Closed |
| Reviewer 1, minor comment 1 | The Mann–Whitney procedure is now a descriptive, power-limited non-detection diagnostic and is not presented as an equivalence test. Counts agree on 20 of 24 passing and 4 of 24 failing pairs. | Closed |
| Reviewer 2, major comment 3 | The same-cohort calibration and the limits of the mechanistic evidence are disclosed consistently. Claims of independent mechanistic fidelity were removed. | Closed |
| Reviewer 2, major comment 4 | Progesterone was added to the non-negativity constraints, the analysis was rerun, and the biological-plausibility rate was corrected from 81/100 to 54/100. Violation magnitudes, the decoding mechanism, the no-clipping decision, an MVN comparison, the strict mechanistic rate, the joint rate, and the matched hull diagnostic are now reported and saved. | Closed |

## Final verified results

### Subgroup diagnostic

- The bootstrap Mann–Whitney analysis is described as a **descriptive
  non-detection diagnostic**, not an equivalence test.
- **20 of 24** feature–condition pairs met the reporting threshold.
- **4 of 24** pairs fell below the threshold.
- The median no-difference-detected fraction was **98.6%**.
- A `p > 0.05` result is explicitly interpreted only as no difference detected
  by that test at the available sample size.
- The conditioned cohorts remain based on **3 PCOS** and **5 Developed PE**
  patients and are not treated as additional independent clinical
  observations.

### Convex-hull and biological-plausibility analysis

- Synthetic patients inside the 23-patient hull: **0 of 100**
- Leave-one-out real patients inside the remaining 22-patient hull: **0 of 23**
- Severity-matched hull-distance medians:
  - Synthetic: **11.0**
  - Real leave-one-out: **11.9**
  - Mann–Whitney: **p = 0.0708**
- Biological constraints satisfied: **54 of 100**
- Strict all-features/all-visits mechanistic envelope satisfied: **41 of 100**
- Both biological and strict mechanistic criteria satisfied: **24 of 100**
- Per-feature mechanistic-envelope agreement remains **83–92%**.

The corrected biological failures comprise:

- Negative estradiol: **20 measurements across 19 patients**
- Negative progesterone: **37 measurements across 37 patients**
- Patients with both types of violation: **10**
- Total SA hormone violations: **57 measurements across 46 patients**
- Minimum estradiol: approximately **−2109.8**
- Minimum progesterone: approximately **−75.8**

The unconstrained concatenated MVN comparator produced **58 negative hormone
measurements across 43 patients**. The manuscript treats negative
concentrations as biologically invalid for both methods. It explains that the
values arise from Gaussian sampling in PCA space followed by unbounded linear
decoding. No post hoc clipping was applied because clipping would alter the
covariance structure being validated.

### Inverse-temperature interpretation

The manuscript and response now distinguish:

1. Low beta gives diffuse local attention and an average-like fixed point in
   the retrieval-landscape interpretation.
2. Exact mixture component probabilities equal the multiplicity weights for
   the unit-normalized memories used here and therefore do not vary with beta.
3. The exact sampler adds Gaussian noise with covariance
   beta-inverse times the identity, so low beta produces broader draws rather
   than samples collapsing to the cohort mean.

### Mechanistic-validation interpretation

The coagulation ODE model is consistently characterized as:

- separately specified;
- blind to the SA generation process;
- applied identically to real and synthetic profiles; and
- calibrated on the same real cohort.

The result is therefore a same-cohort mechanistic consistency check, not
validation against fully independent ground truth.

## Updated reproducibility artifacts

- [`code/data/bootstrap_mw_nondetection.csv`](code/data/bootstrap_mw_nondetection.csv)
  replaces the misleadingly named Mann–Whitney equivalence artifact.
- [`code/data/archive/tost_equivalence_results_legacy.csv`](code/data/archive/tost_equivalence_results_legacy.csv)
  archives and labels the stale TOST output as legacy.
- [`code/data/convex_hull_results.csv`](code/data/convex_hull_results.csv)
  contains the corrected 54/100 biological-plausibility classification.
- [`code/data/convex_hull_biological_violations.csv`](code/data/convex_hull_biological_violations.csv)
  records feature-level violation counts, magnitudes, and real-cohort scales.
- [`code/data/nonnegativity_baseline_comparison.csv`](code/data/nonnegativity_baseline_comparison.csv)
  records the SA and MVN non-negativity comparison.
- [`code/data/hull_severity_matched_summary.csv`](code/data/hull_severity_matched_summary.csv)
  saves the severity-matched hull-distance and radial-overshoot results.

## Updated submission documents

- [`paper/main.pdf`](paper/main.pdf)
- [`paper/supplementary.pdf`](paper/supplementary.pdf)
- [`paper/cover_letter.pdf`](paper/cover_letter.pdf)
- [`peer-review-feedback/response-to-reviewers.md`](peer-review-feedback/response-to-reviewers.md)
- [`peer-review-feedback/response-to-reviewers.docx`](peer-review-feedback/response-to-reviewers.docx)

The Markdown and DOCX response letters contain the same substantive content.
The stale “possible schematic” promise was removed, and the final
cross-reference heading now says “implemented change.”

## Final verification

The following checks passed:

- Exact-sampler correctness suite: **5 of 5**
- Canonical SA helper reproduced the released cohort.
- Canonical concatenated data construction reproduced the expected
  **23 × 216** matrix.
- Mechanistic regression reproduced the cached BZ2012 results:
  pooled overlap approximately **0.902**, KS **p ≈ 0.593**.
- PCA-space ablation reran successfully and reproduced the expected results.
- Convex-hull analysis and severity-matched radial diagnostic reran
  successfully.
- Biological-plausibility counts recomputed to **54/100**, **41/100**, and
  **24/100** for the biological, strict mechanistic, and joint criteria.
- Mann–Whitney artifact recomputed to **20/24 passing and 4/24 failing**.
- Main manuscript compiled to **48 pages**.
- Standalone supplement compiled to **17 pages**.
- Cover letter compiled to **2 pages**.
- LaTeX logs contain no overfull boxes, undefined references, or unresolved
  placeholders.
- Every manuscript and supplementary PDF page was visually inspected; the
  Supplementary Table S5 overflow and the smaller layout defects are resolved.
- The response DOCX passed ZIP/XML integrity and Markdown-to-DOCX plain-text
  round-trip checks. The only round-trip difference was the rendered width of
  the final Markdown table separators, not document content.
- `git diff --check` passed.

LibreOffice was not available in the execution environment, so the response
DOCX could not be rendered for Word-specific visual inspection. Its structure
and text were validated successfully. This is a tooling limitation, not a
remaining reviewer-response content issue.

## Scientific limitations retained in the paper

Closing the reviewer comments does not remove the study's disclosed
limitations:

- The clinical cohort contains only 23 patients.
- The rare-subgroup evidence remains limited to three and five real patients.
- The mechanistic model was calibrated on the same cohort.
- Unconstrained decoding can produce negative hormone concentrations.
- The target distribution exhibits a membership-inference/memorization signal
  at this cohort size.
- Independent-cohort and prospective clinical validation remain future work.

These limitations are now stated consistently and do not represent unresolved
response-letter promises.

## Final conclusion

From the reviewer-response audit standpoint, the revision is complete and the
previously identified issues are closed. The package is ready for the authors'
normal final read-through and submission checks.

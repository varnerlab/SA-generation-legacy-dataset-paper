# Simulated Peer Reviews — NPJ Digital Medicine

**Manuscript:** "Validated Synthetic Patient Generation for Small Longitudinal Cohorts: Coagulation Dynamics Across Pregnancy"

---

## Reviewer 1: Computational / Machine Learning Expert

**Overall Assessment:** The manuscript presents a novel application of modern Hopfield networks to synthetic patient generation from very small clinical cohorts (K=23). The four-level validation framework is a strength, and the mechanistic validation using an independent ODE model is creative and compelling. However, I have several concerns about the methodology, comparisons, and claims that need to be addressed before publication.

**Recommendation:** Major revisions

### Major Concerns

1. **Insufficient baselines.** The paper compares SA only against MVN, which is a straw-man baseline in 2025. The synthetic data generation literature has matured substantially, with methods like CTGAN, TVAE, and tabular diffusion models now widely available. While the authors argue that these methods require larger training sets, this claim is not empirically validated — the authors should at minimum attempt to run CTGAN on this dataset and show that it fails, rather than dismissing it with a citation. A reader familiar with the synthetic data literature will find the single-baseline comparison unconvincing.

2. **PCA truncation to 18 components discards information.** The authors retain 95% of variance, discarding 5%. But the discarded variance may contain clinically meaningful signal — rare features, outlier patterns, or condition-specific variation that lives in the tail components. The paper does not analyze what is lost. A sensitivity analysis varying the PCA threshold (e.g., 90%, 95%, 99%) would strengthen the claims. How do results change if 99% is retained (d_PCA ≈ 22)?

3. **Novelty metric is insufficient.** The reported novelty score of 0.44 (44% angular distance from nearest stored pattern) does not adequately address privacy/memorization concerns. In clinical synthetic data, privacy is a first-order concern. The authors should report: (a) the minimum distance to any real patient in feature space, (b) a formal membership inference attack or re-identification risk analysis, and (c) comparison of novelty distributions between SA and MVN. The current single-number summary is inadequate for a journal focused on digital medicine.

4. **No downstream utility validation.** The four validation levels assess statistical and mechanistic plausibility, but never test whether synthetic data actually improves a downstream task. Can a classifier trained on SA-augmented data outperform one trained on real data alone? Can a mechanistic model calibrated on synthetic data generalize to real patients? Without downstream utility evidence, the clinical value of the synthetic data remains hypothetical.

### Minor Concerns

5. The bootstrap Mann-Whitney test at Level 3 uses 1,000 replicates without correction for 24 simultaneous comparisons. The authors should apply Benjamini-Hochberg FDR correction or acknowledge the multiple testing issue.

6. The claim that "SA operates on the data manifold directly" is imprecise. SA operates in a PCA-reduced space, which is a linear approximation of the data manifold. The authors should temper this language.

7. Algorithm 1 shows the pipeline but does not specify how the initial state ξ_0 is sampled — is it drawn from the unit sphere uniformly? From a specific stored pattern? This affects reproducibility.

8. The paper states N=100 synthetic patients were generated but does not justify this choice. Is 100 sufficient? How do results change with N=50 or N=500?

9. Table 1 reports MRE for only 6 of 72 features. A complete table (all 72 features × 3 visits) should be provided in the supplement.

---

## Reviewer 2: Clinical / Coagulation Domain Expert

**Overall Assessment:** This manuscript addresses an important clinical problem — the scarcity of longitudinal coagulation data in pregnancy, particularly for rare complications like PCOS and preeclampsia. The introduction provides a thorough overview of pregnancy-related coagulopathy and correctly identifies the challenge of small cohorts. The mechanistic validation using the BZ2012 model is a particularly interesting approach. However, I have concerns about the clinical framing, the dataset description, and several claims that overstate the clinical implications.

**Recommendation:** Major revisions

### Major Concerns

1. **Clinical characterization of subgroups is insufficient.** The paper defines PCOS and PE in the introduction but provides no clinical characterization of the actual patients in the study. What were the gestational ages at each visit? Were PE patients diagnosed before or after enrollment? Were PCOS patients on any medications that might affect coagulation (e.g., metformin, oral contraceptives at baseline)? The reader needs to understand the clinical context of the real data before evaluating the synthetic data. A demographics table (age, BMI, gestational age, parity, medication use) is standard for clinical cohorts and is missing entirely.

2. **The "remaining patient with other complications" is unexplained.** The dataset has 23 patients: 14 Healthy, 3 PCOS, 5 PE, and 1 with "other complications." What is this patient's condition? Was this patient included in the unconditioned SA generation? This needs to be clarified, as a single unusual patient could disproportionately influence generation from K=23 stored patterns.

3. **Visit timing is vaguely described.** "Pre-pregnancy baseline, first trimester, and third trimester" spans a wide range. Pre-pregnancy could be months or years before conception. First trimester spans weeks 1-13. Third trimester spans weeks 28-40. Coagulation parameters change substantially within these windows. Were all patients sampled at comparable gestational ages within each visit? If not, the visit-level comparisons may be confounded by gestational age variation.

4. **The BZ2012 calibration raises clinical concerns.** The intrinsic tenase k_cat was reduced to 0.021× its literature value — a 50-fold reduction. This is a very large departure from the published biochemistry. The extrinsic tenase was increased 15.5×. These are not minor adjustments; they suggest that the model structure may be inappropriate for this patient population, or that the TGA experimental conditions differ from the model assumptions (e.g., TF concentration, phospholipid surface). The authors should discuss why such large deviations were necessary and whether this affects the interpretability of the mechanistic validation.

5. **PCOS is not primarily a coagulation disorder.** The introduction frames PCOS as associated with "elevated Factor VIII and impaired fibrinolysis," but PCOS is primarily a hormonal/metabolic disorder. The coagulation changes in PCOS are secondary and variable. With only 3 PCOS patients, it is premature to claim that SA "preserves condition-specific signatures" for PCOS — the statistical power is simply not there to characterize a PCOS coagulation phenotype from 3 patients, let alone validate that a synthetic cohort reproduces it.

### Minor Concerns

6. The manuscript uses "TF Initiator" and "TF + TM Initiator" TGA conditions without specifying the TF concentration used in the assay. This is essential for comparability with other studies.

7. The ROTEM parameters are mentioned in the dataset description but never analyzed in the results. Were they included in the 72 features? If so, why are they not discussed?

8. The claim that SA can "enable hypothesis generation and power analysis" for rare complications is speculative without demonstrating it. A concrete example — e.g., using the synthetic PCOS cohort to estimate the power needed for a prospective study — would strengthen this claim.

9. The acknowledgment cites NIH grants but does not mention IRB approval number. For a clinical dataset, the specific IRB protocol number should be provided.

---

## Reviewer 3: Statistical / Methodological Expert

**Overall Assessment:** The manuscript presents an interesting application of Hopfield-network-based sampling to a challenging small-sample problem. The mathematical framework is clearly presented, and the four-level validation is well-structured. The statistical analyses are generally sound, though I have concerns about several methodological choices and their interpretation. The paper would benefit from a more rigorous statistical treatment in several areas.

**Recommendation:** Minor revisions

### Major Concerns

1. **The KS test at Level 4 tests the wrong hypothesis.** The authors use a two-sample KS test to show that real and synthetic predicted/measured ratio distributions are "statistically indistinguishable" (all p > 0.20). However, with n_real = 69 and n_synth = 300, the KS test has limited power to detect differences — a non-significant result does not demonstrate equivalence. The authors correctly identified this issue at Level 3 and addressed it with bootstrap subsampling, but then abandoned this approach at Level 4 in favor of a standard KS test. Consistency would be improved by either: (a) using the same bootstrap subsampling approach at Level 4, or (b) providing a power analysis showing that the KS test at these sample sizes has adequate power to detect a clinically meaningful difference (e.g., a 0.1 shift in the ratio distribution).

2. **The direction-magnitude decomposition introduces a distributional assumption.** By drawing magnitudes from the empirical distribution of {r_k} (K=23 norms), the authors are implicitly assuming that the norm distribution is well-characterized by 23 samples. For conditioned generation, the PCOS norms are drawn from only 3 samples. This is essentially a 3-sample bootstrap of the magnitude — the resulting norm distribution will have very limited diversity. The authors should discuss this limitation and its potential impact on the dispersion of generated samples.

3. **The 95% PCA threshold creates a hard boundary that affects the eigenvalue spectrum comparison.** The SA spectrum drops at component 18 not because SA respects the data rank, but because PCA truncation mechanically zeroes all variance beyond component 18. This is a methodological artifact, not evidence of SA preserving data structure. The MVN spectrum extends beyond component 22 because Ledoit-Wolf regularization is designed to do exactly this — it is a feature, not a bug, as it provides a full-rank covariance estimate that avoids singularity. The authors should reframe this comparison more carefully: SA and MVN make different trade-offs (SA truncates, MVN regularizes), and the question is which trade-off is more appropriate for the downstream application, not which spectrum is "better."

### Minor Concerns

4. The bootstrap Mann-Whitney analysis at Level 3 subsamples the synthetic data to n_real but does not account for the correlation structure introduced by the SA generation process. If SA-generated patients are not independent (e.g., if the Langevin chain has not mixed), the effective sample size may be smaller than 100, and the subsampling approach would be anti-conservative. The authors should report mixing diagnostics (e.g., autocorrelation of the Langevin chain).

5. The f_target = 0.80 for conditional generation is not justified. Why 80%? How sensitive are the results to this choice? If f_target = 0.60 or 0.95, do the conditioned cohort properties change substantially?

6. The paper reports "median MRE ≈ 2%" but does not provide confidence intervals on this estimate. With K=23 and N=100, the MRE itself has sampling variability. Bootstrap confidence intervals on the MRE would strengthen the claim.

7. The Langevin dynamics uses T=2,000 iterations with α=0.05. Is this sufficient for convergence? The noise scale at β* is √(2×0.05/5.27) ≈ 0.138, which is relatively small. Have the authors verified that the chain has converged (e.g., by comparing results at T=1,000 and T=4,000)?

8. The cross-visit correlation analysis uses Pearson correlation, which assumes linear relationships. Given that the authors argue SA is superior precisely because it captures nonlinear structure, this metric may not capture SA's advantage. Rank-based (Spearman) or mutual-information-based metrics would be more appropriate.

9. The paper should report the computational time for generating 100 patients. If SA requires minutes while MVN requires milliseconds, this is relevant for practitioners considering adoption.

---

## Summary of Key Issues Across All Reviews

| Issue | R1 | R2 | R3 | Priority |
|-------|----|----|-----|----------|
| Insufficient baselines (CTGAN, TVAE) | X | | | High |
| No downstream utility validation | X | | | High |
| Missing demographics table | | X | | High |
| PCA sensitivity analysis | X | | X | High |
| Privacy/memorization analysis | X | | | Medium |
| Clinical characterization of subgroups | | X | | Medium |
| KS test power / equivalence testing | | | X | Medium |
| Direction-magnitude from 3 norms (PCOS) | | | X | Medium |
| BZ2012 calibration departures | | X | | Medium |
| f_target justification | | | X | Low |
| Langevin convergence diagnostics | | | X | Low |
| Computational cost | X | | X | Low |
| Multiple testing correction | X | | | Low |

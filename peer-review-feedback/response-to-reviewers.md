# Response to Reviewers

**Manuscript:** Validated Synthetic Patient Generation for Small Longitudinal Cohorts: Coagulation Dynamics Across Pregnancy
**Journal:** npj Systems Biology and Applications
**Decision:** Major revision
**Date:** 7 August 2026

---

## Note on this document

Each reviewer comment is reproduced below. The comments are paraphrased, but all key points are retained. Each comment is followed by our response and the location of the corresponding manuscript change. Added or revised manuscript text is color-coded by reviewer:

- **Reviewer 1 edits: blue** (`\rone{...}`)
- **Reviewer 2 edits: red** (`\rtwo{...}`)
- **Edits addressing both reviewers: violet** (`\rboth{...}`)

A `\reviewmode` setting in the preamble changes all colored text to black for the final version.

Reviewer 1's original report numbered the comments 1 through 9 and then 11 through 13, with no comment 10. We kept every comment in its original order and renumbered the final three as R1.10, R1.11, and R1.12. These correspond to comments 11, 12, and 13 in the original report. No comment was omitted.

---

## General response to the editor and reviewers

We thank both reviewers for their careful and constructive comments. Both recognized the value of combining statistical and mechanistic validation for very small longitudinal cohorts. They also identified four main concerns: generalizability beyond the single cohort of 23 patients, the separate contributions of principal-component analysis and stochastic attention, the limited independence of the mechanistic validation, and the generator's extrapolation and privacy behavior. We addressed each concern with new analyses and revised text.

We made four main additions. First, we compared stochastic attention with four alternative generators in the same 18-dimensional principal-component space (R1.1, R1.4). Second, we used a convex-hull analysis to ask whether generated profiles lay inside or beyond the region defined by the real cohort (R2.4). Third, we tested whether a synthetic cohort could reveal that an already-known de-identified profile had been stored in the memory matrix (R1.5). We then traced this membership signal to the sampling distribution and derived that distribution in closed form. Fourth, we tested robustness using populations with known covariance, repeated subsets of the clinical cohort, and a nonlinear simulated population (R1.8, R2.1). These analyses separated recovery of a population from reproduction of one small sample and showed the limits of the linear representation.

Two results changed how we interpreted the method. The baseline comparison showed that the shared principal-component space supplied much of the marginal and mechanistic fidelity. The mixture and copula baselines were competitive with stochastic attention on those measures. We therefore removed the claim that stochastic attention performed best overall. We instead evaluated covariance recovery in simulations where the population covariance was known.

The privacy analysis found a clear membership signal at the reported inverse-temperature value. Follow-up tests showed that the signal came from the sampling distribution at this cohort size, not from the numerical sampler or the selected inverse temperature. We also showed that the target distribution was exactly a finite Gaussian mixture centered on the stored profiles. It could therefore be sampled directly without Langevin dynamics. Direct samples reproduced the reported results, including the membership signal. The reported synthetic cohort and all clinical results remained unchanged. The new analysis explained their sampling distribution and showed that the method was not an anonymizer.

Several findings reflected the small source cohort and the choices used to represent and condition the data. These included the convex-hull result, the membership signal, the limits of conditional generation, and the modest downstream gain. We therefore framed the clinical study as a proof of concept. We removed claims of equivalence and limited the proposed uses to exploratory simulation and workflow testing. For statistical inference, the sample size remained the number of real patients, not the number of generated profiles.

We also clarified that the mechanistic validation was not independent of the source cohort (R1.2, R2.3). We added limitations concerning irregular visit timing, missing data, tail compression, and the linear representation. We simplified the Hopfield explanation and moved the detailed derivation to the Supplement (R1.9, R2.minor2). Finally, we reorganized the figures so that the main text still contained seven figures after the new analyses were added (R2.minor1).

---

# Reviewer 1

## Major comments

### R1.1: Disentangle PCA preprocessing from the SA energy landscape (ablation baselines in PCA space)
> The regularization may already arise from the PCA dimensionality-reduction step rather than SA itself. Baselines do not disentangle PCA from the SA energy landscape. Suggests baselines operating in the same PCA-reduced space: (i) PCA+GMM, (ii) PCA+KDE, (iii) PCA+nearest-neighbor interpolation, (iv) PCA+diffusion/interpolation, (v) PCA+copula.

*We agree that the original baselines did not separate the effects of principal-component analysis from those of stochastic attention.* We therefore compared stochastic attention with four generators that used the same 18-dimensional principal-component space: a two-component diagonal-covariance Gaussian mixture, a kernel-density estimator, $k$-nearest-neighbor interpolation, and a Gaussian copula. We decoded and evaluated all five methods in the same way.

The shared representation explained much of the marginal and mechanistic fidelity, but it did not guarantee good performance. The nearest-neighbor method produced profiles close to the stored records. Its median novelty was 0.05, where novelty was defined as $1-\max_k\cos(\hat{\xi},\hat{m}_k)$. Only 5% of its profiles exceeded the 0.2 cutoff, which required a cosine similarity below 0.8 to every stored profile. Marginal error was the median relative error between generated and real per-feature means. The kernel-density estimator had a marginal error of 1.9% and a mechanistic overlap of 0.82. An overlap of 0.82 meant that 82% of its pooled predicted-to-recorded thrombin-generation ratios fell within the 5th to 95th percentile range of the real ratios.

The mixture and copula were competitive with stochastic attention. Their marginal errors were 1.1% and 1.3%, compared with 1.2% for stochastic attention. Their mechanistic overlaps were 0.92 and 0.93, compared with 0.90. All profiles generated by stochastic attention and the copula exceeded the novelty cutoff, as did 94% of the mixture profiles. Stochastic attention had the largest cross-visit correlation error (0.64, compared with 0.42 for the mixture and copula). However, this metric measured how well each method reproduced the correlation matrix of the same 23 patients used to construct or estimate it. It did not measure population recovery. We tested population covariance recovery separately in simulations with a known population covariance (R1.8/R2.1). Stochastic attention had lower error than the sample-covariance Gaussian generator in 35 of 39 low-rank $n<p$ settings. We concluded that the shared principal-component representation supplied much of the marginal and mechanistic fidelity, and that no method performed best on every measure.

We did not include a neural diffusion model because estimating one from 23 observations would not provide a credible small-sample baseline. The nearest-neighbor method provided a non-neural interpolation comparison.

**Changes:** Results (Marginal Plausibility); Supplementary Methods (four baseline implementations); Supplementary Table S8 (five-generator comparison). Blue. `\rone{...}`

---

### R1.2: Mechanistic validation is not fully independent (shared cohort/variables/assumptions)
> The BZ2012 model operates on the same coagulation variables used for generation, is calibrated on the same cohort, and shares the same biological domain assumptions. Temper the "mechanistically indistinguishable" claim and clarify the residual coupling.

*We agree that the original independence claim was too strong.* The mechanistic model was calibrated on the same cohort used to construct the stochastic-attention memory matrix. The analysis therefore tested whether real and synthetic profiles were consistent with the same-cohort model, not whether the synthetic profiles agreed with independent biological evidence. The comparison remained informative because the model was blind to the generation process and applied the same fixed calculation to real and synthetic profiles. We revised the Abstract, Results, and Discussion to make this distinction clear. We did not add an independent calibration experiment because no suitable independent cohort was available, and we identified this as future work.

**Changes:** Abstract wording; Results (Mechanistic Consistency) clarification; Discussion limitations. Violet (shared with R2.3). `\rboth{...}`

---

### R1.3: Statistical reliability of very small subgroups (PCOS n=3, PE n=5)
> Subgroup distributions may reflect interpolation between a handful of trajectories rather than genuine subgroup variability. Add discussion of reliability, explicit acknowledgment of under-sampling, and caution on downstream clinical interpretation.

*We agree and made this limitation explicit.* The conditioned cohorts reproduced the patterns observed in the small subgroups, but generating 100 profiles did not create more independent subgroup evidence. Conditioning increased the sampling weights of the three profiles from patients with PCOS or the five profiles from patients who developed preeclampsia. Profile magnitudes were still drawn from all 23 patients. The participation ratios were 4.6 and 7.7. These values described how strongly repeated draws were concentrated around a few stored profiles; they were not independent patient counts. Treating the 100 generated profiles as independent observations would therefore understate uncertainty. We limited the proposed uses to exploratory hypothesis generation, workflow testing, and planning future studies. Conditional generation did not add evidence about either condition or replace recruitment of additional patients.

**Changes:** Results (Conditional Generation) caveat; Methods (participation-ratio explanation); Discussion. `\rone{...}`

---

### R1.4: Stronger and more diverse baselines
> CTGAN performs poorly at small n; TVAE was not configured for concatenated longitudinal generation; MVN is disadvantaged in rank-deficient settings. Recommend stronger baselines designed for low-sample tabular generation or manifold interpolation.

*We agree that stronger small-sample baselines were needed.* As described in R1.1, we added a Gaussian mixture, a kernel-density estimator, nearest-neighbor interpolation, and a Gaussian copula. All four used the same reduced space as stochastic attention. This avoided comparing stochastic attention only with methods that operated in the full-dimensional data or generated each visit separately.

The mixture and copula matched stochastic attention closely on marginal and mechanistic fidelity. Nearest-neighbor interpolation produced mostly near-copies, whereas the kernel-density estimator had weaker marginal and mechanistic fidelity. All profiles generated by stochastic attention and the copula exceeded the novelty cutoff, as did 94% of the mixture profiles. These results showed that the shared principal-component representation explained much of the fidelity and that no generator performed best on every measure.

**Changes:** Same as R1.1. `\rone{...}`

---

### R1.5: Privacy and memorization (membership inference, re-identification)
> Given K=23 and direct interpolation between stored memory patterns, evaluate membership inference risk and patient re-identification concerns.

*We agree that privacy and memorization required direct testing.* Each real patient's 18-dimensional principal-component profile was stored as one column of the memory matrix. We asked a specific question: if someone already had a candidate de-identified assay profile, could a synthetic cohort reveal whether that profile had been included in the memory matrix? This was not a test of personal identification. The source data contained no direct identifiers or links to named individuals, and neither the synthetic data nor our analysis identified a participant. When the first analysis showed a clear membership signal, we performed additional tests to determine its source.

The test showed that membership could usually be inferred for an already-known candidate profile. We repeated the analysis across 40 random splits. In each split, we constructed the memory matrix from 15 patients, held out 8, generated a synthetic cohort, and measured each candidate's distance to the nearest synthetic profile. A smaller distance provided stronger evidence that the candidate was stored in the memory matrix. Stored profiles were closer to the synthetic cohort than held-out profiles (median distances of 10.7 and 15.4). We summarized this separation using the area under the receiver operating characteristic curve. In this setting, the area was the probability that the test ranked a randomly chosen stored profile above a randomly chosen held-out profile. A value of 0.5 meant chance performance, and 1.0 meant perfect separation. The mean value was 0.97, and 9 of 40 splits reached 1.0. Thus, the test usually distinguished stored from held-out profiles, but it did not guarantee the correct answer for every record. A separate comparison supported this result: synthetic profiles were slightly closer to real records than real records were to one another (median distances of 13.8 and 15.9 in standardized space).

We next tested whether the membership signal came from the numerical sampling procedure. We compared the reported Langevin sampler with three alternatives that changed its starting points, discarded early samples and reduced dependence between retained samples, or added a Metropolis correction. The correction acceptance rate was approximately 0.999, which indicated negligible discretization error. These changes left the mean area under the curve nearly unchanged at 0.96 and only slightly increased the distance to the closest real record. Lowering the inverse temperature over a 40-fold range also failed to remove the signal and reduced cross-visit covariance fidelity.

We then sampled the target distribution directly, without a Markov chain. For the unit-normalized profiles used here, the weighted Hopfield energy defined a finite Gaussian mixture centered on the stored profiles. The probability of selecting each component was proportional to its multiplicity weight, and inverse temperature controlled the common within-component variance. Direct sampling therefore required only selecting a stored profile according to its weight and adding Gaussian noise. These samples reproduced the reported results: marginal error was 1.7% compared with 1.3%, cross-visit covariance error was 0.60 compared with 0.65, mechanistic overlap was 0.91 compared with 0.89, median novelty was 0.50 compared with 0.51, and median distance to the closest record was 13.9 compared with 14.0. Results from 100 independently seeded cohorts showed that this agreement did not depend on one favorable seed. Across the same 40 stored and held-out splits, the mean area under the curve was 0.98. Because direct sampling had no initialization, burn-in, or discretization error, this result showed that the signal came from the target distribution rather than the numerical sampler.

The membership signal therefore arose from the target distribution at this cohort size, where every stored profile was the center of a mixture component. It did not result from the numerical sampler or the selected inverse temperature. The analysis showed that inclusion could be inferred for an already-known de-identified record. It did not attach a name to that record, reconstruct an unknown patient's measurements, or identify a participant. We stated that stochastic attention was not differentially private and that synthetic generation provided no formal privacy guarantee beyond the de-identification of the source data. With 23 patients, neither changing the sampler nor lowering the temperature removed the membership signal while preserving joint fidelity. The intended uses, mechanistic calibration and hypothesis generation, did not require release of patient-level records. We added this interpretation to the manuscript and reported the full diagnostics in a new Supplement subsection.

**Changes:** Results and Discussion (revised privacy treatment); new Supplement subsection "Sampling correctness and privacy diagnostics" (E3, E3b, four sampling procedures, the inverse-temperature sweep, and direct sampling); Methods (closed-form sampling distribution); Introduction (finite-mixture explanation and the source of non-identical samples); corresponding revisions to the computational-cost and interpolation language. Blue. `\rone{...}`

---

### R1.6: Robustness to irregular visit timing, missingness, variable-length trajectories
> The concatenation strategy assumes visits are temporally and physiologically aligned across patients. Discuss robustness to irregular timing, incomplete trajectories, and variable-length records.

*We agree that the current data representation has this limitation.* Concatenating visits preserved relationships across visits, but it required every patient to have the same assays at the same aligned times. The current pipeline could not directly use irregular visit times, missing assays or visits, or variable-length records. Imputation, padding, or a sequence-based method could support these data, but each would add assumptions that would need validation in a small cohort. We added this limitation to the Discussion and presented it together with the limits of the linear principal-component representation.

**Changes:** Discussion (limitations). `\rone{...}`

---

### R1.7: Clinical importance of tail compression
> SA smooths distributional extremes, but rare/extreme states are often the phenomena of interest. Discuss effect on downstream mechanistic modeling, implications for rare-event hypothesis generation, and whether the framework systematically underestimates severe pathology.

*We agree that tail compression required a stronger warning.* Truncating the principal-component representation and drawing profile magnitudes from only 23 observed values may both have contributed to this compression. If a generated cohort under-represented severe states, it could not establish how often those states occurred or how severe they became. Although synthetic profiles extended beyond the real cohort's convex hull, that result did not show that they reproduced the true tails of the population. We therefore stated that these cohorts should not be used to estimate extreme risk without additional real data. We also clarified that the main mechanistic-overlap measure used the 5th to 95th percentile range of the real cohort and therefore did not validate tail behavior. The Kolmogorov–Smirnov tests compared the full distributions but were not designed specifically to test the tails.

**Changes:** Discussion (expand tail-compression paragraph); links to E2. `\rone{...}`

---

### R1.8 / R2.1: Generalizability beyond a single cohort
> Validation is on a single K=23 cohort from a single domain. Encourage external validation or explicit limitation; recommend proof-of-concept framing. (R2: a second dataset, or a synthetic benchmark / simulation study demonstrating robustness to varying data characteristics.)

*We agree that one clinical cohort could not establish generalizability.* We added three robustness tests and explained the question, comparison, and error measure for each in the Methods.

The first test asked whether each method could recover a known population covariance from a small, high-dimensional sample. We drew samples from low-rank Gaussian populations with known covariance. We varied the sample size ($n=8$, 15, 23, 40, or 80), dimension ($p=30$, 120, or 216), and latent rank ($r=3$, 5, or 10). For each setting, we constructed a stochastic-attention memory matrix and estimated a Gaussian generator from the same sample. Each method generated 100 profiles. We measured covariance error as the square root of the summed squared entry-by-entry differences between the generated and known covariance matrices, divided by the corresponding value for the known matrix. This relative Frobenius error allowed comparisons across settings and tested recovery of the population covariance rather than reproduction of the sampled profiles. Stochastic attention had lower error in 35 of the 39 rank-deficient $n<p$ settings. At $n=23$ and $p=216$, the Gaussian generator's error was approximately 1.5 times the stochastic-attention error after averaging across ranks.

The second test asked how stochastic attention changed when fewer real patients were stored in the memory matrix. We randomly selected 20, 15, 10, or 8 of the 23 patients and repeated the selection ten times at each size. We constructed a new memory matrix from each subset, generated 100 profiles, and compared every generated cohort with the same complete set of 23 real patients. This provided a common reference across subset sizes. As the subset size decreased from 20 to 8, mean cross-visit correlation error increased from 0.63 to 1.11 and pooled marginal error increased from about 2.0% to 4.9%. Mechanistic overlap remained near 0.90.

The third test asked how the linear representation performed when the population followed a curved pattern. We used samples of 23 points from a noisy S-curve embedded in 30 or 120 dimensions. We measured each method's distance from the known curve and its covariance error as curvature increased. The mean distance of stochastic-attention samples from the curve doubled relative to the straight-line case between curvatures 0.25 and 0.5. Stochastic-attention samples remained closer to the curve than Gaussian samples, but the Gaussian generator recovered covariance more accurately. Neither method performed better on both measures.

These tests showed good covariance recovery for stochastic attention in most small, linear, low-rank settings. They also showed worse performance with fewer patients and a clear limitation for the nonlinear population tested here. We also noted that the same stochastic-attention method had been used in separate studies of discrete protein-sequence generation and multiplicity-weighted conditional generation. We framed the present clinical study as a proof of concept and stated that validation in additional clinical cohorts remained necessary.

**Changes:** New Methods subsection "Robustness Tests"; new Results subsection "Robustness to Sample Size, Dimensionality, and Nonlinearity"; revised main-text figure (`fig:sim-recovery`) and caption; added cross-domain context in the Introduction and Discussion; explicit proof-of-concept framing in the Abstract, Introduction, Discussion, and Conclusion. The Discussion also stated that the subgroup test did not establish equivalence. Violet. `\rboth{...}`

---

### R1.9 / R2.minor2: Theory exposition for a biomedical audience
> The entropy-inflection (β) selection, participation-ratio interpretation, and geometric interpretation of the energy landscape are hard for a broad biomedical audience; add intuition or a schematic. (R2: streamline the Hopfield theory to focus on the applied contribution.)

*We agree that the theory needed a simpler explanation.* We revised the Methods to explain the two roles of β. In the Hopfield update, low β spread attention across many stored profiles and moved toward their weighted average. High β concentrated attention on individual profiles. In direct sampling, β did not change the probability of selecting each unit-normalized profile. Those probabilities were set by the multiplicity weights. Instead, β controlled how widely samples spread around the selected profile.

We chose β\* at the point where the attention-entropy curve bent downward most strongly. We calculated this as the most negative finite-difference second derivative of normalized entropy with respect to log β. We defined the participation ratio as the number of equally weighted mixture components that would give the same concentration of sampling probability. It was not an independent patient count. We moved the entropy definition and calculation details to a new Supplementary Methods subsection while keeping the explanation and resulting values in the main text. We did not add another figure because the main text already contained seven figures and Reviewer 2 requested figure consolidation.

**Changes:** Methods (simplified explanation); Supplement (detailed theory and calculation). Violet. `\rboth{...}`

---

### R1.10: Unvalidated feature categories (fibrinolytic, viscoelastic)
> Mechanistic validation covered only thrombin generation, not fibrinolytic or viscoelastic features, which are a substantial portion of the feature space. Discuss how unvalidated categories may influence the generated manifold and whether they contribute meaningfully to the energy landscape; expand future validation.

*We agree and clarified this limitation.* The rotational thromboelastometry and fibrinolytic features were part of the full profile used to construct the principal-component representation and memory matrix. They therefore affected the distribution from which synthetic profiles were drawn. Their median relative errors were comparable to those of the other feature groups, but we assessed them only statistically. The mechanistic model represented thrombin generation, not fibrinolysis or clot viscoelasticity. Mechanistic validation of these features would require models that represent those processes.

**Changes:** Results (Marginal Plausibility) + Discussion. `\rone{...}`

---

### R1.11: Downstream improvement may reflect smoothing/variance suppression
> The synthetic-calibrated model's improvement over the real-calibrated model may partly reflect smoothing/regularization that suppresses physiologically meaningful heterogeneity.

*We agree and made this interpretation more direct.* The modest downstream gain (0.94× ratio) did not show that synthetic data were superior. All 100 synthetic profiles were derived from the same 23 patients, so they did not increase the independent evidence about the broader patient population. Their value was computational. Averaging the calibration objective over a denser sample from the represented distribution may have made the optimization more stable. We also connected this interpretation to the tail-compression limitation and did not claim that we had established the cause of the numerical gain.

**Changes:** Results (Downstream Utility) + Discussion. `\rone{...}`

---

### R1.12: Computational scalability
> Discuss how the framework behaves as cohort size increases, whether the energy landscape becomes unstable at larger K, and how complexity scales with dimensionality.

*We added a direct description of computational cost.* In the reported Langevin implementation, every sampling step evaluated all stored profiles. The cost of generating one profile therefore scaled as $O(TKd)$, where $T$ was the number of steps, $K$ was the number of stored profiles, and $d$ was their dimension. The closed-form result removed the sequential factor $T$. The component probabilities could be calculated once. Each new profile then required selecting one component, drawing a $d$-dimensional Gaussian vector, and decoding the profile. We clarified that the inverse transformation and decoding were still required for every generated profile.

The simulation benchmark did not establish how the method would perform in larger clinical cohorts. Its advantage over the sample-covariance Gaussian generator was concentrated in the $n<p$ settings, and the clinical analyses included only 23 patients. We therefore removed the expectation that membership disclosure would necessarily decrease as $K$ increased. Stability, fidelity, and membership risk at larger cohort sizes remained open questions.

**Changes:** Discussion, including a brief connection to the simulation benchmark. `\rone{...}`

---

## Minor comments

### R1.M1: Temper strong wording
> Claims of "clinically useful synthetic cohorts" and "statistically indistinguishable" may be too strong given the sample size.

*We agree and revised these claims throughout the manuscript.* The subgroup analysis was a descriptive, low-power Mann–Whitney test. A result of `p > 0.05` meant only that the test did not detect a difference in that replicate; it did not establish equivalence. We also limited the utility claim to the exploratory modeling uses tested in this proof-of-concept study.

**Changes:** Abstract, Introduction, Discussion, Conclusion wording. `\rone{...}`

---

# Reviewer 2

## Major comments

### R2.1: Generalizability beyond a single dataset
> All validation is on a single cohort (K=23). Application to a second dataset, or a synthetic benchmark / simulation study demonstrating robustness to varying data characteristics, would significantly strengthen the manuscript.

*We agree that validation in one clinical cohort was a major limitation.* As described in **R1.8**, we added three robustness tests. They measured recovery of a known covariance across several sample sizes, dimensions, and latent ranks; changes in performance as fewer clinical profiles were stored in the memory matrix; and performance as a simulated population became more curved. We also added context from previous protein-sequence studies and framed the clinical study as a proof of concept. These simulations tested specific data properties, but they did not replace validation in another clinical cohort.

**Changes:** See R1.8. Violet. `\rboth{...}`

---

### R2.2: Role and limitations of PCA; nonlinear embeddings
> Discuss scenarios where PCA may be insufficient; comment on whether nonlinear embeddings could be incorporated and the expected trade-offs between interpretability, fidelity, and flexibility.

*We agree that the limits of the linear representation needed to be tested and discussed.* In the S-curve simulation, generated profiles moved farther from the known curve as curvature increased. Between curvatures 0.25 and 0.5, the mean stochastic-attention distance reached twice its value for a straight population. Stochastic-attention samples remained closer to the curve than Gaussian samples, but the Gaussian generator recovered covariance more accurately. Neither method performed best on both measures. We limited this conclusion to the tested S-curve and did not generalize it to all nonlinear populations.

We also stated that nonlinear representations could describe curved or branching relationships more accurately. However, they would require additional choices for estimation, interpretation, and decoding that would be difficult to support with 23 patients. The linear representation was a practical choice for this cohort, not evidence that clinical trajectories were generally linear. The baseline comparison also showed that this shared representation supplied much of the marginal and mechanistic fidelity.

**Changes:** Discussion (limits of the linear principal-component representation), with a connection to E1. Red. `\rtwo{...}`

---

### R2.3: Mechanistic validation: shared calibration dependence
> The ODE is calibrated on the same real dataset used to generate the synthetic data, creating partial dependence. Clarify that the validation shows consistency with a model fitted to the same cohort rather than fully independent mechanistic fidelity; if feasible, use independently calibrated parameters or an alternative model.

*We agree.* As described in **R1.2**, we stated that the mechanistic model was calibrated on the same cohort used to construct the memory matrix. The analysis therefore tested consistency with a same-cohort model, not agreement with independent biological evidence. We revised the independence language and identified calibration with independent data or an alternative model as future work.

**Changes:** See R1.2. Violet. `\rboth{...}`

---

### R2.4: Interpolation vs. extrapolation (convex hull, novel phenotypes)
> Does SA generate novel out-of-sample variation or primarily interpolate? To what extent do generated patients lie within the convex hull of observed data? Does it produce biologically plausible but previously unseen phenotypes?

*We addressed this question with geometric, biological, and mechanistic checks.* We first calculated the convex hull, the smallest convex region containing all 23 real profiles in the 18-dimensional principal-component space. All 100 synthetic profiles lay outside this hull. However, each real profile also lay outside the hull formed by the other 22 real profiles. Thus, every real profile was an outer corner of this sparse, high-dimensional region, and an inside-or-outside test alone could not show whether the synthetic profiles extrapolated unusually far.

We therefore compared distances from the hull. Our hypothesis was that excessive extrapolation would place synthetic profiles farther from the real cohort than real profiles treated as unseen. Both groups needed to be evaluated against hulls containing 22 real profiles. The leave-one-out calculation already did this for each real profile. For each synthetic profile, we removed its nearest real profile before constructing the comparison hull. The median distances were 11.0 for synthetic profiles and 11.9 for held-out real profiles. A descriptive Mann–Whitney test did not detect a difference ($p=0.07$). This result did not establish equivalence, but it provided no evidence that synthetic profiles lay farther from the observed cohort than held-out real profiles. It also showed why hull membership alone was not informative in this setting. The hull occupied little of the surrounding space, and the generator's normalization meant that a generated direction could lie inside it only by matching a stored direction exactly.

Distance from the hull did not show whether a profile was biologically valid, so we applied two additional checks. The first required values that should be positive to remain positive and required coagulation-factor and thrombin-generation values to remain within five real-cohort standard deviations of the mean. Fifty-four of 100 synthetic profiles passed. The other 46 contained a negative estradiol or progesterone value. The unconstrained multivariate-normal baseline had a similar result, with negative hormone values in 43 of 100 profiles. Neither method enforced positivity. The hormones were closer to zero relative to their observed variation than the coagulation factors, which explained why only hormones crossed zero in these samples. The positive factor values did not provide a general guarantee of positivity.

The second check used the BZ2012 model to test whether the generated coagulation-factor combinations produced thrombin-generation behavior consistent with the real cohort. Estradiol and progesterone were not model inputs, and all nine generated factor inputs were positive. We therefore included all 100 synthetic profiles. For each thrombin-generation feature, 83% to 92% of profiles had a model-predicted-to-recorded ratio within the 5th to 95th percentile range of the real cohort. Forty-one profiles met this criterion for every feature and visit, and 24 passed both the biological and mechanistic checks. These results showed that extrapolation and biological plausibility were separate questions. Positive-variable transformations or constrained decoding could enforce non-negativity in future implementations. The Supplement reported the hormone-specific counts, extreme values, and comparison with the real-cohort scale.

**Changes:** New Results paragraph + Supplement (E2). `\rtwo{...}`

---

### R2.5: Demonstration of real-world utility
> Does synthetic data improve predictive models beyond the original cohort? Enhance power for subgroup analyses? Enable concrete new hypotheses? Even a limited quantitative demonstration would strengthen translational significance.

*We agree that the paper needed a concrete account of its practical value.* The revised Results showed that synthetic profiles could support mechanistic-model calibration, exploratory comparisons of small subgroups, and factor-level hypotheses for testing in larger studies. However, confirming those hypotheses and establishing clinical value would require additional real patients or independent experimental data.

We also added the original motivation for developing the method to the Introduction and Discussion. Synthetic profiles could also be used, for example, to train deep learning models that predict TGA responses from coagulation factors or estimate factor profiles compatible with an observed TGA response. We have not yet evaluated either task. As noted in R1.11, the 0.94× downstream error ratio did not add independent evidence about the broader population. It may instead have reflected a more stable calibration objective.

**Changes:** Introduction (original motivating application); Results (Downstream Utility); Discussion (forward and inverse modeling use cases). `\rtwo{...}`

---

## Minor comments

### R2.minor1: Consolidate figures
> Several figures could move to the supplementary material to improve readability.

*We agree.* We moved the principal-component-by-visit comparison of stochastic attention and the multivariate-normal baseline to the Supplement. We placed the simulation benchmark (`fig:sim-recovery`) in its former position. The main text therefore retained seven figures despite the new analyses. We reported the baseline, convex-hull, membership, and sampling diagnostics in Supplementary tables and figures rather than adding more main-text panels.

**Changes:** Figure placement (floats.tex / supplementary.tex). `\rtwo{...}`

---

### R2.minor2: Streamline Hopfield theory exposition
> Some sections (particularly the Hopfield network theory) could be streamlined to focus on the applied contribution.

As described in **R1.9**, we moved the detailed theory to the Supplement and added a plain-language explanation to the Methods.

**Changes:** See R1.9. Violet. `\rboth{...}`

---

## Cross-reference of comments and implemented changes

| Comment | Type | Vehicle | Color |
|---|---|---|---|
| R1.1, R1.4 | New experiment | E1 PCA-space ablation suite | blue |
| R1.2, R2.3 | Prose (temper) | Mechanistic independence clarification | violet |
| R1.3 | Prose | Subgroup reliability caveats | blue |
| R1.5 | New experiment | E3 membership inference / privacy | blue |
| R1.6 | Prose | Irregular-timing / missingness limitations | blue |
| R1.7 | Prose | Tail-compression clinical importance | blue |
| R1.8, R2.1 | New experiment + framing | Simulation benchmark + cross-domain evidence | violet |
| R1.9, R2.minor2 | Writing | Theory streamlining + corrected intuition | violet |
| R1.10 | Prose | Fibrinolytic/viscoelastic features | blue |
| R1.11 | Prose | Downstream smoothing caveat | blue |
| R1.12 | Prose | Computational scalability | blue |
| R1.M1 | Writing | Temper wording | blue |
| R2.2 | Prose | PCA limits / nonlinear embeddings | red |
| R2.4 | New experiment | E2 convex-hull / extrapolation | red |
| R2.5 | Framing | Downstream utility translational framing | red |
| R2.minor1 | Layout | Figure consolidation | red |

# Response to Reviewers

**Manuscript:** Validated Synthetic Patient Generation for Small Longitudinal Cohorts: Coagulation Dynamics Across Pregnancy
**Journal:** npj Systems Biology and Applications
**Decision:** Major revision
**Date:** 7 August 2026

---

## Note on this document

Each reviewer comment is reproduced (paraphrased faithfully, key specifics retained), followed by our response and a pointer to the corresponding manuscript change. In the revised manuscript, added or changed text is color-coded by reviewer:

- **Reviewer 1 edits — blue** (`\rone{...}`)
- **Reviewer 2 edits — red** (`\rtwo{...}`)
- **Edits addressing both reviewers — violet** (`\rboth{...}`)

A `\reviewmode` toggle in the preamble reverts all colored text to black for the camera-ready version.

One note on numbering. Reviewer 1's comments are numbered 1 through 9 and then 11 through 13 in the original report, with no comment 10. We have kept every comment in its original order and renumbered the final three as R1.10, R1.11 and R1.12 so that the count is continuous. In terms of the original report, our R1.10 is comment 11 (unvalidated feature categories), R1.11 is comment 12 (downstream improvement and variance suppression), and R1.12 is comment 13 (computational scalability). No comment has been omitted.

---

## General response to the editor and reviewers

We thank both reviewers for their careful and constructive reading. We are grateful for the recognition that combining statistical with mechanistic validation is a distinguishing strength of the work, and that learning from very small longitudinal cohorts is an important and underserved problem. The reviewers converged on four priorities: generalizability beyond the single K=23 cohort, disentangling the contribution of PCA from the stochastic-attention energy landscape, the degree of independence of the mechanistic validation, and the interpolation-versus-extrapolation and privacy behavior of the generator. We have addressed each with new work rather than with argument alone.

Four additions carried the revision. A PCA-space ablation (R1.1, R1.4) fit four alternative generators in the identical d=18 subspace and evaluated them on the same harness as stochastic attention. A convex-hull analysis (R2.4) asked whether generated patients lay inside or beyond the real cohort. A membership-inference and distance-to-closest-record analysis (R1.5), followed by a sampling-correctness investigation that ended in a closed-form characterization of the target distribution, established what the generated cohort revealed about its training set and why. A three-part robustness analysis (R1.8, R2.1) tested recovery of a known low-rank covariance, degradation under repeated clinical-cohort subsampling, and sensitivity to a nonlinear population manifold. These tests allowed us to distinguish population recovery from reproduction of one small cohort and to identify limits of the linear representation.

Two of these did not come out the way we expected, and we have let the results stand rather than presenting them selectively. The ablation showed that the shared PCA subspace supplied much of the marginal and mechanistic fidelity, and that mixture and copula baselines fit in that subspace were competitive with stochastic attention on those axes. We therefore removed the claim that stochastic attention dominated and moved the covariance-generalization argument to the simulation benchmark, where reproducing the training sample could not substitute for recovering the population. The privacy analysis found a genuine memorization signal at the published operating point, and the follow-up investigation established that it was intrinsic to the target distribution at this cohort size rather than an artifact of the sampler or of the operating temperature. Pushed to its end, that investigation produced the sharpest result in the revision: the latent target our method defines is exactly a finite Gaussian mixture centered on the stored patients, so it can be drawn from analytically rather than by Langevin dynamics. Exact draws reproduced the reported cohort's behavior on every axis we evaluated, including the membership-inference signal, which settled the attribution and also removed the sequential cost factor from future generation. The reported cohort and every clinical result remained unchanged; the new analysis identified precisely which distribution they came from. We state plainly that the method is not an anonymizer.

A theme emerged that we have made explicit in the Discussion. Several of the behaviors the reviewers asked about — the convex-hull result, the membership-inference result, the effective pattern count under conditioning, and the modest downstream gain — were consequences of a single fact, that 23 patients in a 216-dimensional space formed a sparse sample, rather than four separate properties of the method. Recognizing this shaped the revision. We now frame the clinical study as a proof of concept, replace equivalence language with what the tests support at this sample size, and bound the translational claim to exploratory simulation and workflow testing whose inferential sample size remains the number of real patients.

We have also tempered the independence claim for the mechanistic validation (R1.2, R2.3), added the scope limits the reviewers identified around irregular visit timing, missingness, tail compression and the linear embedding, streamlined the Hopfield exposition with plain-language intuition while moving the denser derivation to the supplement (R1.9, R2.minor2), and consolidated the figure set so that the main text still carries seven figures despite the new experiments (R2.minor1).

---

# Reviewer 1

## Major comments

### R1.1 — Disentangle PCA preprocessing from the SA energy landscape (ablation baselines in PCA space)
> The regularization may already arise from the PCA dimensionality-reduction step rather than SA itself. Baselines do not disentangle PCA from the SA energy landscape. Suggests baselines operating in the same PCA-reduced space: (i) PCA+GMM, (ii) PCA+KDE, (iii) PCA+nearest-neighbor interpolation, (iv) PCA+diffusion/interpolation, (v) PCA+copula.

*Excellent and central question—this is the crux of attributing the result to stochastic attention rather than to dimensionality reduction, and our ablation answered it directly and candidly.* We ran a PCA-space ablation suite comprising a two-component diagonal-covariance Gaussian mixture, a kernel-density estimator, $k$-nearest-neighbor interpolation, and a Gaussian copula.
All four operated in the identical $d=18$ PCA subspace and were decoded through the same evaluation harness as stochastic attention. The result was nuanced. The shared representation appeared to account for much of the marginal and mechanistic fidelity but did not guarantee performance.
The nearest-neighbor samples remained close to stored profiles: their median novelty, defined as $1-\max_k\cos(\hat{\xi},\hat{m}_k)$, was 0.05, and only 5% exceeded the operational 0.2 cutoff, corresponding to cosine similarity below 0.8 to every stored pattern.
The kernel-density estimator had higher marginal error (1.9%) and lower pooled TF-only mechanistic cloud overlap (0.82): after pooling across patients, visits, and TGA features, 82% of the predicted-to-recorded TGA ratios for its synthetic patients fell within the real-data 5th–95th percentile range.
By contrast, the mixture and copula were competitive with stochastic attention on marginal error (1.1% and 1.3%, respectively, versus 1.2%) and mechanistic overlap (0.92 and 0.93 versus 0.90). All 100 stochastic-attention and copula samples exceeded the novelty cutoff, as did 94% of mixture samples.
Stochastic attention had the largest empirical cross-visit correlation error of the five methods (0.64, versus 0.42 for the mixture and copula), but this metric measured reproduction of the correlation matrix of the same 23 patients used to construct each generator rather than population generalization. We therefore assessed covariance generalization in the controlled simulation benchmark with known population covariance (R1.8/R2.1), where stochastic attention was more accurate than the sample-covariance Gaussian generator in 35 of 39 low-rank Gaussian $n<p$ settings.
The substantive answer was that the shared principal-component representation supplied much of the marginal and mechanistic fidelity. Within this comparison, no generator dominated across the evaluated metrics.
Regarding diffusion, we did not include a neural score model because fitting one to 23 observations would not provide a credible low-sample baseline; the nearest-neighbor method represented non-neural interpolation.

**Changes:** New combined Results paragraph (Marginal Plausibility) + Supplementary Table S8 reporting all five generators on the identical harness, violet where it also answers R2.2. `\rone{...}`

---

### R1.2 — Mechanistic validation is not fully independent (shared cohort/variables/assumptions)
> The BZ2012 model operates on the same coagulation variables used for generation, is calibrated on the same cohort, and shares the same biological domain assumptions. Temper the "mechanistically indistinguishable" claim and clarify the residual coupling.

*We appreciate this point and agree the independence claim was overstated.* We tempered the language throughout (abstract, results, discussion) and added an explicit statement of the residual dependence: the ODE parameters were calibrated on the same cohort (real → stochastic attention → synthetic → calibration), so the validation demonstrated *consistency with a model fitted to the same biological system*, not fully independent ground truth. We clarified that the strength of the test lay in the ODE being blind to the stochastic-attention generation process and processing real and synthetic patients through the identical fixed mapping, and we reported that no difference was detected by this model under the tested conditions. We considered independent or alternative calibration; we addressed the point in prose rather than with a new experiment and identified validation against independent data as future work.

**Changes:** Abstract wording; Results (Mechanistic Consistency) clarification; Discussion limitations. Violet (shared with R2.3). `\rboth{...}`

---

### R1.3 — Statistical reliability of very small subgroups (PCOS n=3, PE n=5)
> Subgroup distributions may reflect interpolation between a handful of trajectories rather than genuine subgroup variability. Add discussion of reliability, explicit acknowledgment of under-sampling, and caution on downstream clinical interpretation.

*We agree, and we have made this caveat explicit and quantitative.* The conditioned cohorts reproduced the observed subgroup signatures, but generating 100 profiles did not increase the amount of independent subgroup information. The PCOS and Developed PE signals still came from only three and five real patients; after accounting for the lower-weight background patterns, the corresponding participation ratios were 4.6 and 7.7. In the PCOS case, magnitudes were also drawn from only three empirical norms. The revised Results stated that treating the 100 synthetic profiles as independent observations would understate uncertainty. We therefore framed conditioning as providing a denser set of profiles around the distributions suggested by those few patients, which could support exploratory hypothesis generation, workflow testing, and simulations for planning future studies. It did not add evidence about the conditions or replace recruitment of additional patients. The Methods additionally explained that the participation ratio quantified the effective number of stored patient patterns contributing to generation, making the small underlying evidence base explicit.

**Changes:** Results (Conditional Generation) caveat; Methods (participation-ratio explanation); Discussion. `\rone{...}`

---

### R1.4 — Stronger and more diverse baselines
> CTGAN performs poorly at small n; TVAE was not configured for concatenated longitudinal generation; MVN is disadvantaged in rank-deficient settings. Recommend stronger baselines designed for low-sample tabular generation or manifold interpolation.

*A fair point.* The PCA-space ablation suite (R1.1) added exactly this class of fairer, low-sample-appropriate baselines—density estimation, mixture models, copulas, and manifold interpolation—all operating in the same reduced space as stochastic attention rather than being handicapped by the full-dimensional or per-visit formulations. We reported their results alongside the earlier multivariate-normal, CTGAN, and TVAE comparisons and presented them candidly.
The mixture and copula were competitive with stochastic attention on marginal and mechanistic fidelity. Nearest-neighbor interpolation produced mostly near-copies, whereas kernel-density estimation had weaker marginal and mechanistic fidelity. All 100 stochastic-attention and copula samples exceeded the operational novelty cutoff, and 94% of mixture samples did as well.
The substantive answer was that the shared principal-component representation supplied much of the marginal and mechanistic fidelity. Within this comparison, no generator dominated across the evaluated metrics.

**Changes:** Same as R1.1. `\rone{...}`

---

### R1.5 — Privacy and memorization (membership inference, re-identification)
> Given K=23 and direct interpolation between stored memory patterns, evaluate membership inference risk and patient re-identification concerns.

*An important omission, and pursuing it changed how we present the method. Thank you.* Each real patient's 18-dimensional principal-component profile formed one column of the memory matrix used during stochastic-attention sampling. We therefore considered a narrowly defined membership-inference question: if someone already possessed a candidate de-identified assay profile, could the synthetic cohort reveal whether that record had been included in the memory matrix? This was not personal identification: the source cohort was fully de-identified and contained no direct identifiers or link to named individuals, and neither the synthetic cohort nor our analysis identified a participant. Because the first result was a genuine membership signal rather than the reassuring one we had expected, we followed it with a full sampling-correctness investigation to attribute the cause and reported the complete picture.

The membership-inference test gave clear evidence that inclusion in the memory matrix could usually be inferred from an already-known candidate profile. Across 40 random splits, we constructed the memory matrix from 15 patients, held out the other 8, generated a synthetic cohort, and scored each candidate by its distance to the nearest synthetic profile; a smaller distance was treated as stronger evidence of inclusion. Memory-matrix members were closer than held-out records (median 10.7 versus 15.4). We summarized the separation using the area under the receiver operating characteristic curve (AUC), which here was the probability that the test ranked a randomly chosen member as more likely to be included than a randomly chosen held-out record. A value of 0.5 indicated no discrimination and 1.0 indicated perfect discrimination. The mean AUC was 0.97, and 9 of 40 splits reached 1.0. Thus, the test usually distinguished included from held-out profiles, although it did not guarantee a correct answer for every individual record. The initial distance-to-closest-record comparison was consistent with this result: synthetic profiles lay slightly closer to real records than real records lay to one another (median 13.8 versus 15.9 in the standardized space). We then established what drove the signal. First, it was not a numerical sampling artifact: a four-rung diagnostic ladder replaced the published anchor-initialized single-endpoint Langevin scheme with sphere-initialized, burn-in-controlled, thinned and pooled Metropolis-adjusted sampling. The acceptance rate was near 0.999, so discretization bias was negligible; the ladder left the membership-inference AUC essentially unchanged at 0.96 and raised the distance to the closest record only marginally, never reaching the real-to-real baseline. Second, the inverse temperature did not remove the signal: sweeping beta across a 40-fold range around the operating point raised this distance monotonically as beta fell but saturated well below the real-to-real baseline at every value tested, while cross-visit covariance fidelity degraded.

Third, and most directly, the signal was not a property of finite-time dynamics at all. Both checks above still compared one Markov chain against another, so we extended the correctness analysis to an analytic reference. Doing so showed that the latent target was not merely approximated by our sampler but was available in closed form: for the unit-normalized memories used in the manuscript, the weighted Hopfield energy induced exactly a finite isotropic Gaussian mixture whose components were centered on the stored patients, with covariance beta-inverse times the identity and mixture weights proportional to the multiplicities. That distribution could be sampled ancestrally, by choosing a stored patient according to its weight and adding isotropic noise, with no chain of any kind. Passing such exact draws through the identical magnitude, reconstruction, and de-standardization pipeline reproduced the reported cohort's evaluated behavior across the entire harness: marginal error 1.7% against 1.3%, cross-visit Frobenius 0.60 against 0.65, mechanistic cloud overlap 0.91 against 0.89, median novelty 0.50 against 0.51, and median distance to closest record 13.9 against 14.0, with each value falling inside the range spanned by the four ladder rungs. One hundred independently seeded exact cohorts placed these quantities in narrow intervals around those values, so the agreement did not depend on a favorable seed. Evaluated over the same 40 train/holdout splits, the membership-inference AUC under exact sampling was 0.98, if anything marginally higher than under the published scheme. Because ancestral sampling had no initialization, burn-in, or discretization error, this located the memorization in the target distribution itself rather than in any aspect of how that distribution was sampled.

The memorization was therefore intrinsic to the target at this cohort size, where every stored patient was the center of a mixture component, rather than a defect of the sampler or a poorly chosen operating point.

We therefore framed this correctly rather than as a favorable finding. The analysis demonstrated membership disclosure for an already-known de-identified record; it did not attach a name to that record, reconstruct an unknown patient's measurements, or re-identify a participant. We stated plainly that stochastic attention was not differentially private and that synthetic generation did not add a formal privacy guarantee beyond the source cohort's de-identification. At $K=23$, neither a better sampler nor a lower temperature yielded a configuration that eliminated the membership signal while preserving joint fidelity. The intended uses—mechanistic calibration and hypothesis generation—did not require releasing patient-level records. We added these results and this framing to the manuscript and detailed the diagnostics in a new Supplement subsection.

**Changes:** Rewritten privacy treatment in Results/Discussion + new Supplement subsection "Sampling correctness and privacy diagnostics" (E3, E3b, the MCMC ladder, the beta-privacy sweep, and the exact analytic reference) + the closed-form target in Methods + clarification in the Introduction that the induced distribution is a finite Gaussian mixture and that within-component variance produces non-identical samples + the corresponding correction to the computational-cost and interpolation language. Blue. `\rone{...}`

---

### R1.6 — Robustness to irregular visit timing, missingness, variable-length trajectories
> The concatenation strategy assumes visits are temporally and physiologically aligned across patients. Discuss robustness to irregular timing, incomplete trajectories, and variable-length records.

*A valuable point about real-world deployment.* We add a limitations paragraph to the Discussion. The concatenation that makes cross-visit structure available also constrains the data the method accepts: stacking three visits into one vector assumes every patient has every visit and that visits are aligned across patients, which holds in this cohort by design but not in most observational data. We state the three consequences separately. Irregular timing would place physiologically different states in the same coordinate, so the principal components would mix visit structure with between-patient variation. Missing assays or visits have no natural representation in a fixed-length vector and would require imputation before the embedding, which at this cohort size would itself be a substantial source of uncertainty. Variable-length records fall outside the formulation altogether and would need either padding to a common length or a sequence-aware extension in which attention operates over visits as well as over patients. We note that these are scope boundaries of the present formulation rather than properties of the energy landscape, which is indifferent to what the coordinates mean.

**Changes:** Discussion (limitations). `\rone{...}`

---

### R1.7 — Clinical importance of tail compression
> SA smooths distributional extremes, but rare/extreme states are often the phenomena of interest. Discuss effect on downstream mechanistic modeling, implications for rare-event hypothesis generation, and whether the framework systematically underestimates severe pathology.

*We agree this deserves more weight than we gave it.* We expanded the tail-compression discussion to cover its two sources (PCA truncation of higher moments and the finite empirical magnitude distribution of K=23 norms) and, more importantly, its clinical consequence: if the generated cohort under-represents the severe end of a distribution, then estimates of how often a severe phenotype occurs, or of how extreme it can become, are biased toward the center. We connected this to the extrapolation analysis but did not treat it as reassuring: stochastic attention placed patients outside the convex hull of the real cohort and therefore was not confined to the observed range, but producing out-of-range points was not the same as reproducing the frequency of extreme ones. We therefore recommended against using stochastic-attention cohorts to estimate extreme-quantile risk without real anchoring data. We also noted a scope limit of our evidence: the mechanistic checks compared distributions within a 5th-to-95th-percentile band and were therefore, by construction, a test of the central mapping rather than of the tails, so they were not interpreted as reassurance about tail behavior.

**Changes:** Discussion (expand tail-compression paragraph); links to E2. `\rone{...}`

---

### R1.8 / R2.1 — Generalizability beyond a single cohort
> Validation is on a single K=23 cohort from a single domain. Encourage external validation or explicit limitation; recommend proof-of-concept framing. (R2: a second dataset, or a synthetic benchmark / simulation study demonstrating robustness to varying data characteristics.)

*This is the most important shared concern, and we addressed it on two fronts.* (1) **Three robustness benchmarks:** we added a Methods subsection that defined each test, its reference target, and its error measure before the corresponding Results. First, we drew training cohorts from low-rank Gaussian populations with known covariance. For every combination of sample size (n=8, 15, 23, 40, or 80), dimensionality (p=30, 120, or 216), and latent rank (r=3, 5, or 10), stochastic attention and a Gaussian generator using the empirical sample covariance received the same training cohort and each generated 100 profiles. Relative Frobenius error against the known population covariance therefore measured population recovery rather than reproduction of the training draw. Stochastic attention had lower error in 35 of the 39 rank-deficient n<p settings; at n=23 and p=216, the Gaussian generator's error was about 1.5 times that of stochastic attention after averaging across ranks. Second, we drew ten random subsets of the clinical cohort at each of n=20, 15, 10, and 8, fit stochastic attention to each subset, and scored every generated cohort against the same full 23-patient reference. From n=20 to n=8, mean cross-visit correlation error rose from 0.63 to 1.11 and pooled marginal error rose from about 2.0% to 4.9%, whereas mechanistic overlap remained near 0.90. Third, we fit both generators to samples of 23 points from a noisy S-curve embedded in p=30 or 120 dimensions. We measured both distance from the known curve and covariance error as curvature increased. The off-manifold error of stochastic attention doubled relative to the straight-line case between curvatures 0.25 and 0.5. Its samples remained closer to the curve than the Gaussian samples, whereas the Gaussian generator recovered the population covariance more accurately; neither method dominated both criteria. (2) **Cross-domain evidence:** the identical stochastic-attention machinery, with β* predicted from PCA dimension alone, had been applied in separate studies to discrete protein-sequence generation from small family alignments and to multiplicity-weighted conditional steering—evidence that the geometric mechanism was not specific to clinical coagulation data. We also reframed the clinical study explicitly as a proof of concept and strengthened the limitations accordingly.

**Changes:** New Methods subsection "Robustness Benchmarks"; new Results subsection "Robustness to Sample Size, Dimensionality, and Nonlinearity"; revised main-text figure (`fig:sim-recovery`) and caption; Introduction and Discussion cross-domain framing elevated from aside to argument; proof-of-concept framing stated explicitly in the Abstract, Introduction, and Conclusion, with the Discussion carrying the matching bounded-claim language (exploratory modeling and workflow testing; the subgroup diagnostic is not an equivalence test). Violet. `\rboth{...}`

---

### R1.9 / R2.minor2 — Theory exposition for a biomedical audience
> The entropy-inflection (β) selection, participation-ratio interpretation, and geometric interpretation of the energy landscape are hard for a broad biomedical audience; add intuition or a schematic. (R2: streamline the Hopfield theory to focus on the applied contribution.)

*A helpful suggestion for accessibility.* We add plain-language intuition for both quantities the reviewers singled out, while separating the retrieval landscape from the exact sampling law. In the energy-landscape view, low β produces diffuse local attention and an average-like retrieval fixed point, whereas high β sharpens attention into distinct memory attractors. In the exact sampler, the component probabilities are the multiplicity weights and are independent of β for the unit-normalized memories used here; β controls the within-component covariance β⁻¹I, so low β produces broader draws rather than samples that collapse to the cohort mean. β\* is estimated from the transition in the retrieval landscape. For the participation ratio we pose the question it answers directly — if attention is spread unevenly across 23 patients, how many patients is the generator effectively drawing on? — and note that reporting it alongside a conditioned cohort states how much independent material that cohort rests on. We have moved the mechanical derivation (entropy definition, sweep parameters, finite-difference procedure) into a new Supplementary Methods subsection, keeping the concept, the domain-agnostic argument and the resulting values in the main text. We did not add a schematic figure: the main text carries seven figures and Reviewer 2 asks us to consolidate rather than expand, so we judged the prose intuition the better way to meet this request without displacing a results figure.

**Changes:** New/streamlined Methods exposition; denser theory moved to the Supplement. Violet. `\rboth{...}`

---

### R1.10 — Unvalidated feature categories (fibrinolytic, viscoelastic)
> Mechanistic validation covered only thrombin generation, not fibrinolytic or viscoelastic features, which are a substantial portion of the feature space. Discuss how unvalidated categories may influence the generated manifold and whether they contribute meaningfully to the energy landscape; expand future validation.

*Thank you—we clarified this.* We noted that ROTEM/fibrinolytic features were generated within the concatenated vector with comparable MREs but were not mechanistically validated. We expanded this statement to explain that these variables contributed through their loadings on the retained principal components to the memory patterns that defined the energy landscape. Their fidelity was therefore evaluated statistically, while mechanistic validation using models that explicitly represent fibrinolysis and clot viscoelasticity remains future work.

**Changes:** Results (Marginal Plausibility) + Discussion. `\rone{...}`

---

### R1.11 — Downstream improvement may reflect smoothing/variance suppression
> The synthetic-calibrated model's improvement over the real-calibrated model may partly reflect smoothing/regularization that suppresses physiologically meaningful heterogeneity.

*We agree and had already flagged the smoother loss landscape as the likely mechanism; we sharpened this.* We stated explicitly that the modest downstream gain (0.94× ratio) was consistent with variance reduction from N=100 versus K=23 training points and was not interpreted as stochastic attention adding biological information; it demonstrated that stochastic attention preserved *sufficient* structure for calibration, not that synthetic data was superior to real data. We connected this to the tail-compression caveat: suppressing heterogeneity was the same phenomenon that reduced optimization variance.

**Changes:** Results (Downstream Utility) + Discussion. `\rone{...}`

---

### R1.12 — Computational scalability
> Discuss how the framework behaves as cohort size increases, whether the energy landscape becomes unstable at larger K, and how complexity scales with dimensionality.

*A useful addition, and working through it corrected an expectation of our own.* On cost, the scaling follows directly from the update rule: each Langevin step evaluates the attention weights over the stored patterns and forms their weighted combination, so work per step is proportional to K×d and a sample costs that product times the number of steps T. The method is linear in both the cohort size and the retained dimension, and at very high feature counts the principal-component decomposition rather than the sampling becomes the practical bound.

On behavior we were more careful than we had planned to be. We had expected to argue that the simulation benchmark supported favorable scaling to larger cohorts. It did not show that. Instead, the advantage of stochastic attention over the sample-covariance Gaussian comparator was concentrated in the n<p settings, where the empirical covariance was rank deficient. The conclusion was that the method was positioned for the small-cohort, high-dimensional regime rather than as a general-purpose generator, and the Discussion now says so. We expect the memorization documented under R1.5 to ease as K grows, because basins that currently stand apart would begin to overlap and generation would interpolate across several patients rather than remaining near one, but we marked that as an expectation rather than a result: this cohort admitted only subsampling downward, and we did not test it.

**Changes:** Discussion (+ brief note tied to simulation benchmark). `\rone{...}`

---

## Minor comments

### R1.M1 — Temper strong wording
> Claims of "clinically useful synthetic cohorts" and "statistically indistinguishable" may be too strong given the sample size.

*Agreed.* We adopted more conservative wording throughout. The subgroup analysis was described as a power-limited, descriptive Mann–Whitney non-detection diagnostic: `p > 0.05` meant that this test did not detect a difference in that replicate and did not establish equivalence. We also qualified "clinically useful" as demonstrated only for the exploratory modeling use cases tested, within a proof-of-concept framing.

**Changes:** Abstract, Introduction, Discussion, Conclusion wording. `\rone{...}`

---

# Reviewer 2

## Major comments

### R2.1 — Generalizability beyond a single dataset
> All validation is on a single cohort (K=23). Application to a second dataset, or a synthetic benchmark / simulation study demonstrating robustness to varying data characteristics, would significantly strengthen the manuscript.

See **R1.8**—addressed jointly through three explicitly described robustness tests: recovery of a known covariance across combinations of sample size, dimensionality, and latent rank; degradation under repeated clinical-cohort subsampling against a common full-cohort reference; and departure from a known nonlinear manifold as curvature increased. These tests were accompanied by elevated cross-domain protein evidence and explicit proof-of-concept framing. *We appreciate the reviewer's recognition that additional clinical cohorts are hard to obtain; the simulation study was designed precisely to demonstrate robustness across data characteristics without one.*

**Changes:** See R1.8. Violet. `\rboth{...}`

---

### R2.2 — Role and limitations of PCA; nonlinear embeddings
> Discuss scenarios where PCA may be insufficient; comment on whether nonlinear embeddings could be incorporated and the expected trade-offs between interpretability, fidelity, and flexibility.

*A thoughtful point, and the simulation benchmark allowed us to answer it quantitatively rather than only in discussion.* On a nonlinear S-curve ground truth, the linear principal-component subspace used by stochastic attention began to place generated points off the true manifold as curvature grew, with the off-manifold residual crossing twice its linear-baseline value at a curvature of roughly 0.25 to 0.5. Stochastic attention stayed closer to the manifold than a multivariate Gaussian even past that point, whereas the Gaussian recovered the gross covariance better, so neither was uniformly preferable. This bounded where linearity started to matter. We added discussion of PCA's limits: it captures linear structure, so strongly nonlinear manifolds (branching trajectories, threshold effects) are imperfectly represented and can contribute to tail compression. We discussed nonlinear embeddings (kernel PCA, autoencoders, diffusion maps) as extensions, with the trade-off that they sacrifice the interpretability and the closed-form β* ≈ √d prediction that the linear PCA subspace provides, and can be unstable to fit at n=23, so linearity was a deliberate choice for the small-cohort regime, not an oversight. The E1 ablation results indicated which failure modes were PCA-driven versus sampler-driven.

**Changes:** Discussion (PCA limitations paragraph); ties to E1. Violet where shared with R1.1. `\rtwo{...}`

---

### R2.3 — Mechanistic validation: shared calibration dependence
> The ODE is calibrated on the same real dataset used to generate the synthetic data, creating partial dependence. Clarify that the validation shows consistency with a model fitted to the same cohort rather than fully independent mechanistic fidelity; if feasible, use independently calibrated parameters or an alternative model.

See **R1.2** — addressed jointly. *We really appreciate this precise articulation of the dependence structure* and have made it explicit in the Discussion, tempered the independence language, and flagged independent/alternative-model calibration as future work.

**Changes:** See R1.2. Violet. `\rboth{...}`

---

### R2.4 — Interpolation vs. extrapolation (convex hull, novel phenotypes)
> Does SA generate novel out-of-sample variation or primarily interpolate? To what extent do generated patients lie within the convex hull of observed data? Does it produce biologically plausible but previously unseen phenotypes?

*An excellent and precise question, and the analysis produced an instructive answer once we controlled for sample size.* We added a convex-hull/extrapolation analysis (E2) in the 18-dimensional principal-component space, followed by direct biological and mechanistic plausibility checks. The convex hull was the smallest convex region enclosing all 23 real patients. All 100 synthetic patients lay outside it. In a separate leave-one-out analysis, each real patient also lay outside the hull formed by the remaining 22 real patients, which indicated that every real patient formed an outer corner, or vertex, of this sparse, high-dimensional hull. This binary inside-or-outside result did not show whether synthetic patients lay just beyond the hull or much farther away. We therefore used distance to the hull to measure the degree of extrapolation. Our hypothesis was that unusually extrapolative generation would place synthetic patients farther from the hull than real patients treated as unseen. For a like-for-like comparison, the distances for both groups were calculated against hulls containing 22 real patients. The leave-one-out calculation already omitted the real patient being evaluated; for each synthetic patient, we omitted its nearest real patient. Median distance to the resulting hulls was 11.0 for synthetic patients and 11.9 for held-out real patients, and a descriptive Mann–Whitney test did not detect a difference (p = 0.07). Thus, the comparison provided no evidence that the synthetic profiles lay farther beyond the observed cohort than held-out real patients. This result did not establish equivalence, but it showed why binary hull membership alone was not an informative measure of extrapolation in this sparse, high-dimensional setting. The geometry further explained this limitation: the hull occupied little of the surrounding space, and the generator's normalization meant that a generated direction could remain inside it only by exactly matching a stored direction.

Because geometric proximity did not establish biological validity, we assessed plausibility in two ways. First, 54 of 100 synthetic patients satisfied all implemented biological constraints; the remaining 46 had at least one negative estradiol or progesterone concentration. The unconstrained concatenated multivariate-normal baseline had a comparable failure rate, with negative hormone concentrations in 43 of 100 patients. Thus, this limitation affected both unconstrained continuous generators rather than stochastic attention alone. Second, we asked whether the generated coagulation-factor combinations produced thrombin-generation behavior consistent with the real cohort under the BZ2012 model. Estradiol and progesterone were not inputs to this model, and all nine generated coagulation-factor inputs were positive. We therefore retained all 100 synthetic patients for the mechanistic check; the negative hormone values were not discarded or passed into the ODE system. For each thrombin-generation feature, 83–92% of synthetic patients had a model-predicted-to-recorded ratio within the 5th–95th percentile range observed in the real cohort. Forty-one patients met this range for every evaluated feature at every visit, and 24 met both this criterion and all biological constraints. The interpretation was therefore bounded: the distance comparison did not show that the synthetic profiles extrapolated farther than held-out real patients, but it did not establish equivalence; moreover, not every profile was biologically valid, and only 24% passed both the strict mechanistic and biological criteria. The supplementary table retained the hormone-specific counts, extreme values, and explanation of why negative decoded concentrations could occur. We also identified positive-variable transformations or constrained decoding as future extensions that could enforce non-negativity by construction.

**Changes:** New Results paragraph + Supplement (E2). `\rtwo{...}`

---

### R2.5 — Demonstration of real-world utility
> Does synthetic data improve predictive models beyond the original cohort? Enhance power for subgroup analyses? Enable concrete new hypotheses? Even a limited quantitative demonstration would strengthen translational significance.

*We appreciate the push toward translational concreteness.* We strengthened the framing of the downstream-utility experiment and connected it to the two use cases the reviewer named, while bounding the claim. The revised Results stated that the method supported mechanistic-model calibration without using the real cohort directly in that calibration, comparisons of subgroups containing only three or five patients, and factor-level hypotheses for testing in a larger study. We were equally explicit about what it did not do: it did not substitute for additional patients, and none of these results spoke to clinical deployment, which would require prospective validation against real clinical decisions. Per R1.11, we also stated that the 0.94× downstream ratio was a variance-reduction effect rather than evidence of new biological content.

**Changes:** Results (Downstream Utility) + Discussion framing. `\rtwo{...}`

---

## Minor comments

### R2.minor1 — Consolidate figures
> Several figures could move to the supplementary material to improve readability.

*Agreed, and this dovetails with the new material.* We moved the stochastic-attention-versus-multivariate-normal PCA-by-visit figure to the supplement and added the simulation-benchmark figure (`fig:sim-recovery`) in its place, so the main text carried seven figures, unchanged from the original submission despite the new experiments. The ablation, convex-hull, membership-inference, and sampling diagnostics were reported as supplementary tables and figures rather than as additional main-text panels.

**Changes:** Figure placement (floats.tex / supplementary.tex). `\rtwo{...}`

---

### R2.minor2 — Streamline Hopfield theory exposition
> Some sections (particularly the Hopfield network theory) could be streamlined to focus on the applied contribution.

See **R1.9** — addressed jointly by moving dense theory to the supplement and adding plain-language intuition.

**Changes:** See R1.9. Violet. `\rboth{...}`

---

## Cross-reference: comment → implemented change

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

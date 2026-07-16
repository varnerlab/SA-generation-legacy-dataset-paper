# Response to Reviewers

**Manuscript:** Validated Synthetic Patient Generation for Small Longitudinal Cohorts: Coagulation Dynamics Across Pregnancy
**Journal:** npj Systems Biology and Applications
**Decision:** Major revision
**Date:** [fill on submission]

---

## Note on this document

Each reviewer comment is reproduced (paraphrased faithfully, key specifics retained), followed by our response and a pointer to the corresponding manuscript change. In the revised manuscript, added or changed text is color-coded by reviewer:

- **Reviewer 1 edits — blue** (`\rone{...}`)
- **Reviewer 2 edits — red** (`\rtwo{...}`)
- **Edits addressing both reviewers — violet** (`\rboth{...}`)

A `\reviewmode` toggle in the preamble reverts all colored text to black for the camera-ready version.

**Status legend for our internal tracking (remove before submission):**
`[PLAN]` decided approach · `[PENDING]` awaiting experiment results · `[DRAFTED]` response written · `[DONE]` manuscript change made

---

## General response to the editor and reviewers

[DRAFTED LAST] We thank both reviewers for their careful and constructive reading. We are grateful for the recognition that the integration of statistical and independent mechanistic validation is a distinguishing strength of the work, and that the problem — learning from very small longitudinal cohorts — is important and underserved. The reviewers converged on a shared set of priorities: (1) generalizability beyond the single K=23 cohort, (2) disentangling the contribution of PCA from the stochastic-attention energy landscape, (3) the degree of independence of the mechanistic validation, and (4) interpolation-vs-extrapolation and privacy behavior of the generator. We have addressed each with new experiments and sharpened framing, summarized below and detailed in the point-by-point responses. [Summarize the three new experiments + simulation benchmark + cross-domain evidence once results are in.]

---

# Reviewer 1

## Major comments

### R1.1 — Disentangle PCA preprocessing from the SA energy landscape (ablation baselines in PCA space)
> The regularization may already arise from the PCA dimensionality-reduction step rather than SA itself. Baselines do not disentangle PCA from the SA energy landscape. Suggests baselines operating in the same PCA-reduced space: (i) PCA+GMM, (ii) PCA+KDE, (iii) PCA+nearest-neighbor interpolation, (iv) PCA+diffusion/interpolation, (v) PCA+copula.

**[DRAFTED — E1 complete]** *Excellent and central question — this is the crux of attributing the result to SA rather than to dimensionality reduction, and our ablation answers it directly (and candidly).* We ran a PCA-space ablation suite: PCA+GMM, PCA+KDE, PCA+k-NN interpolation, and PCA+Gaussian copula, all fit in the identical d=18 PCA subspace and decoded through the identical evaluation harness (marginal MRE, cross-visit covariance, mechanistic cloud overlap/KS, and novelty/memorization). The result is nuanced and we report it exactly as measured. The shared PCA subspace does carry much of the marginal and mechanistic fidelity: PCA+GMM and PCA+Gaussian copula are competitive with SA on marginal error (median 1.1–1.3% vs. SA 1.2%) and on the mechanistic envelope (overlap 0.92–0.93 vs. SA 0.90). Two baselines fail cleanly on a single axis — PCA+k-NN memorizes (median novelty 0.05; 5% of samples exceed the novelty threshold, vs. 100% for SA), and PCA+KDE degrades both marginals (1.9% MRE) and mechanistic plausibility (overlap 0.82, below the within-real band). The distinctive property of SA is that it is the **only method that combines full novelty (100% of samples, median 0.51 — roughly double GMM/copula and an order of magnitude above k-NN) with competitive marginal and mechanistic fidelity**: the memorization-leaning baselines obtain their fidelity in part by staying close to the 23 real patients, which SA does not do. Consistent with this, the empirical 216×216 cross-visit *correlation* distance is actually smaller for the memorization-leaning baselines than for SA — an expected artifact of that metric, which rewards reproducing the sample correlation matrix and therefore favors methods that generate less novel points. We therefore make the covariance-*generalization* argument where memorization cannot help — the controlled simulation benchmark with known population covariance (see R1.8/R2.1), in which SA recovers the true cross-visit structure in the n<p regime while parametric baselines do not. The honest reading for the reviewer's question is thus: the PCA subspace supplies much of the marginal/mechanistic fidelity, and SA's specific contribution is *principled novel generation at competitive fidelity* rather than dominance on every metric. Regarding (iv) diffusion: a score network cannot be trained at n=23; the k-NN interpolation baseline represents the interpolation family, and we note this explicitly.

**Changes:** New Results paragraph + Supplement section (E1 tables/figure), violet where it also answers R2.2. `\rone{...}`

---

### R1.2 — Mechanistic validation is not fully independent (shared cohort/variables/assumptions)
> The BZ2012 model operates on the same coagulation variables used for generation, is calibrated on the same cohort, and shares the same biological domain assumptions. Temper the "mechanistically indistinguishable" claim and clarify the residual coupling.

**[PLAN]** *We appreciate this point and agree the independence claim was overstated.* We temper the language throughout (abstract, results, discussion) and add an explicit statement of the residual dependence: the ODE parameters were calibrated on the same cohort (real → SA → synthetic → calibration), so the validation demonstrates *consistency with a model fitted to the same biological system*, not fully independent ground truth. We clarify that the strength of the test lies in the ODE being blind to the SA generation process and processing real and synthetic patients through the identical fixed mapping, and we soften "indistinguishable" to "not distinguishable by this model under these conditions." (We considered independent/alternative calibration; we address the point in prose rather than a new experiment, and flag independent mechanistic validation — including the fibrinolysis model under development — as future work.)

**Changes:** Abstract wording; Results (Mechanistic Consistency) clarification; Discussion limitations. Violet (shared with R2.3). `\rboth{...}`

---

### R1.3 — Statistical reliability of very small subgroups (PCOS n=3, PE n=5)
> Subgroup distributions may reflect interpolation between a handful of trajectories rather than genuine subgroup variability. Add discussion of reliability, explicit acknowledgment of under-sampling, and caution on downstream clinical interpretation.

**[PLAN]** *We agree and have made this caveat explicit.* We add discussion that at these sizes the conditioned manifolds are severely under-sampled: the effective pattern count under multiplicity weighting drops below 5 (K_eff ≈ 4.6 for PCOS), magnitudes are drawn from only 3 empirical norms, and the "calibration gap" between internal conditioning and decoded phenotype (characterized geometrically in our companion protein work) applies here. We frame subgroup amplification as hypothesis-generating and power-analysis-supporting, not as a substitute for recruiting real rare-subgroup patients, and add explicit cautions on downstream clinical interpretation.

**Changes:** Results (Conditional Generation) caveat; Discussion. `\rone{...}`

---

### R1.4 — Stronger and more diverse baselines
> CTGAN performs poorly at small n; TVAE was not configured for concatenated longitudinal generation; MVN is disadvantaged in rank-deficient settings. Recommend stronger baselines designed for low-sample tabular generation or manifold interpolation.

**[DRAFTED — E1 complete]** *A fair point.* The PCA-space ablation suite (R1.1) adds exactly this class of fairer, low-sample-appropriate baselines — density estimation, mixture models, copulas, and manifold interpolation — all operating in the same reduced space as SA rather than being handicapped by the full-dimensional or per-visit formulations. We report their results alongside MVN/CTGAN/TVAE and, per R1.1, present them candidly: given the shared PCA subspace the mixture-model and copula baselines are competitive with SA on marginal and mechanistic fidelity, while interpolation (k-NN) and kernel-density baselines each fail a specific axis (novelty and marginal/mechanistic fidelity, respectively). What distinguishes SA is generating fully novel patients (100% above the novelty threshold) while retaining that fidelity, rather than uniformly outperforming every baseline on every metric.

**Changes:** Same as R1.1. `\rone{...}`

---

### R1.5 — Privacy and memorization (membership inference, re-identification)
> Given K=23 and direct interpolation between stored memory patterns, evaluate membership inference risk and patient re-identification concerns.

**[PLAN → PENDING E3]** *Important omission — thank you.* We add a privacy analysis (E3): a full Distance-to-Closest-Record distribution (synthetic→real vs. real→real, with a holdout baseline) and a nearest-neighbor membership-inference attack reporting AUC. Results [pending] indicate [AUC ≈ 0.5, DCR comparable to real–real], i.e., low empirical re-identification risk consistent with the novelty scores already reported. We state honestly that SA is not differentially private and that K=23 permits no formal guarantee, and we discuss membership-inference and re-identification risk explicitly.

**Changes:** New Results/Discussion paragraph + Supplement (E3). `\rone{...}`

---

### R1.6 — Robustness to irregular visit timing, missingness, variable-length trajectories
> The concatenation strategy assumes visits are temporally and physiologically aligned across patients. Discuss robustness to irregular timing, incomplete trajectories, and variable-length records.

**[PLAN]** *A valuable point about real-world deployment.* We add a limitations/future-work discussion: the current concatenation assumes aligned, complete visit structure (which this cohort has by design); irregular timing and missingness would require [imputation prior to embedding / per-visit alignment / a time-aware embedding], and variable-length records would require padding or a sequence-aware extension. We frame this as a scope boundary of the present proof-of-concept and a concrete extension direction.

**Changes:** Discussion (limitations). `\rone{...}`

---

### R1.7 — Clinical importance of tail compression
> SA smooths distributional extremes, but rare/extreme states are often the phenomena of interest. Discuss effect on downstream mechanistic modeling, implications for rare-event hypothesis generation, and whether the framework systematically underestimates severe pathology.

**[PLAN]** *We agree this deserves more weight than we gave it.* We expand the tail-compression discussion: its two sources (PCA truncation of higher moments; finite empirical magnitude distribution of K=23 norms), its consequence that severe-pathology tails may be underrepresented, and the resulting caution for rare-event hypothesis generation. We connect this to the extrapolation analysis (E2): SA does produce some out-of-hull phenotypes, but tail coverage remains a genuine limitation, and we recommend against using SA to estimate extreme-quantile risk without real anchoring data.

**Changes:** Discussion (expand tail-compression paragraph); links to E2. `\rone{...}`

---

### R1.8 / R2.1 — Generalizability beyond a single cohort
> Validation is on a single K=23 cohort from a single domain. Encourage external validation or explicit limitation; recommend proof-of-concept framing. (R2: a second dataset, or a synthetic benchmark / simulation study demonstrating robustness to varying data characteristics.)

**[PLAN → PENDING SIM]** *This is the most important shared concern and we address it on two fronts.* (1) **Simulation benchmark:** we add a controlled study with known ground-truth structure (specified low-rank cross-visit covariance and marginal shapes, including heavy tails) sampled across varying (n, p) and intrinsic rank. SA recovers the known structure across the n<p regime where parametric baselines fail, and we map the regime boundaries where it degrades — demonstrating robustness to varying data characteristics without requiring a second clinical cohort. (2) **Cross-domain evidence:** the identical SA machinery, with β* predicted from PCA dimension alone, has been independently validated on discrete protein-sequence generation from small family alignments (median seed ≈ 22 sequences, matching our K=23 regime) and on multiplicity-weighted conditional steering — evidence that the geometric mechanism is domain-agnostic. We also reframe the clinical study explicitly as a proof-of-concept and strengthen the limitations accordingly.

**Changes:** New Results subsection (simulation benchmark, one main-text figure) + Supplement; Introduction/Discussion cross-domain framing elevated from aside to argument; proof-of-concept framing in Abstract/Discussion. Violet. `\rboth{...}`

---

### R1.9 / R2.minor2 — Theory exposition for a biomedical audience
> The entropy-inflection (β) selection, participation-ratio interpretation, and geometric interpretation of the energy landscape are hard for a broad biomedical audience; add intuition or a schematic. (R2: streamline the Hopfield theory to focus on the applied contribution.)

**[PLAN]** *A helpful suggestion for accessibility.* We [add a schematic figure illustrating the energy landscape, stored patterns, and Langevin sampling] and add plain-language intuition for β selection (the temperature at which stored patients stop blurring into one average and start behaving as distinct attractors) and the participation ratio (an effective count of contributing patterns). We move the denser theoretical derivations to the supplement and keep the main text focused on the applied contribution.

**Changes:** New/streamlined Methods exposition; possible schematic figure; denser theory → Supplement. Violet. `\rboth{...}`

---

### R1.10 — Unvalidated feature categories (fibrinolytic, viscoelastic)
> Mechanistic validation covered only thrombin generation, not fibrinolytic or viscoelastic features, which are a substantial portion of the feature space. Discuss how unvalidated categories may influence the generated manifold and whether they contribute meaningfully to the energy landscape; expand future validation.

**[PLAN]** *Thank you — we clarify this.* We already note that ROTEM/fibrinolytic features are generated within the concatenated vector with comparable MREs but are not mechanistically validated; we expand this to state that these features do contribute to the PCA-derived energy landscape (they load onto retained components) and therefore shape generation, so their statistical fidelity matters even without a mechanistic check. We point to the fibrinolysis mechanistic model under development as the route to closing this gap.

**Changes:** Results (Marginal Plausibility) + Discussion. `\rone{...}`

---

### R1.11 — Downstream improvement may reflect smoothing/variance suppression
> The synthetic-calibrated model's improvement over the real-calibrated model may partly reflect smoothing/regularization that suppresses physiologically meaningful heterogeneity.

**[PLAN]** *We agree and already flagged the smoother loss landscape as the likely mechanism; we sharpen this.* We state explicitly that the modest downstream gain (0.94× ratio) is consistent with variance reduction from N=100 vs. K=23 training points and should not be read as SA adding biological information; it demonstrates that SA preserves *sufficient* structure for calibration, not that synthetic data is superior to real data. We connect this to the tail-compression caveat (suppressing heterogeneity is the same phenomenon that reduces optimization variance).

**Changes:** Results (Downstream Utility) + Discussion. `\rone{...}`

---

### R1.12 — Computational scalability
> Discuss how the framework behaves as cohort size increases, whether the energy landscape becomes unstable at larger K, and how complexity scales with dimensionality.

**[PLAN]** *A useful addition.* We add a scalability discussion: per-sample cost scales with the K×d attention evaluation over T Langevin steps (linear in K and in d), so the method scales favorably to larger cohorts; the energy landscape becomes *more* stable, not less, as K grows relative to d (the small-K regime is the hard one). The simulation benchmark (R1.8) exercises larger n and provides empirical support. We note memory/PCA cost as the practical bound at very high p.

**Changes:** Discussion (+ brief note tied to simulation benchmark). `\rone{...}`

---

## Minor comments

### R1.M1 — Temper strong wording
> Claims of "clinically useful synthetic cohorts" and "statistically indistinguishable" may be too strong given the sample size.

**[PLAN]** *Agreed.* We adopt more conservative wording throughout: "statistically indistinguishable" → "not distinguishable by the tests applied at this sample size," and we qualify "clinically useful" as demonstrated for the modeling use-cases tested (mechanistic calibration, hypothesis generation), within a proof-of-concept framing.

**Changes:** Abstract, Introduction, Discussion, Conclusion wording. `\rone{...}`

---

# Reviewer 2

## Major comments

### R2.1 — Generalizability beyond a single dataset
> All validation is on a single cohort (K=23). Application to a second dataset, or a synthetic benchmark / simulation study demonstrating robustness to varying data characteristics, would significantly strengthen the manuscript.

**[PLAN → PENDING SIM]** See **R1.8** — addressed jointly via the simulation benchmark, elevated cross-domain protein evidence, and explicit proof-of-concept framing. *We appreciate the reviewer's recognition that additional clinical cohorts are hard to obtain; the simulation study is designed precisely to demonstrate robustness across data characteristics without one.*

**Changes:** See R1.8. Violet. `\rboth{...}`

---

### R2.2 — Role and limitations of PCA; nonlinear embeddings
> Discuss scenarios where PCA may be insufficient; comment on whether nonlinear embeddings could be incorporated and the expected trade-offs between interpretability, fidelity, and flexibility.

**[PLAN]** *A thoughtful point.* We add discussion of PCA's limits: it captures linear structure, so strongly nonlinear manifolds (e.g., branching trajectories, threshold effects) would be imperfectly represented and could contribute to tail compression. We discuss nonlinear embeddings (kernel PCA, autoencoders, diffusion maps) as extensions, with the trade-off that they sacrifice the interpretability and the closed-form β* ≈ √d prediction that the linear PCA subspace provides, and can be unstable to fit at n=23 — so linearity is a deliberate choice for the small-cohort regime, not an oversight. The E1 ablation results inform which failure modes are PCA-driven vs. sampler-driven.

**Changes:** Discussion (PCA limitations paragraph); ties to E1. Violet where shared with R1.1. `\rtwo{...}`

---

### R2.3 — Mechanistic validation: shared calibration dependence
> The ODE is calibrated on the same real dataset used to generate the synthetic data, creating partial dependence. Clarify that the validation shows consistency with a model fitted to the same cohort rather than fully independent mechanistic fidelity; if feasible, use independently calibrated parameters or an alternative model.

**[PLAN]** See **R1.2** — addressed jointly. *We really appreciate this precise articulation of the dependence structure* and have made it explicit in the Discussion, tempered the independence language, and flagged independent/alternative-model calibration as future work.

**Changes:** See R1.2. Violet. `\rboth{...}`

---

### R2.4 — Interpolation vs. extrapolation (convex hull, novel phenotypes)
> Does SA generate novel out-of-sample variation or primarily interpolate? To what extent do generated patients lie within the convex hull of observed data? Does it produce biologically plausible but previously unseen phenotypes?

**[PLAN → PENDING E2]** *An excellent and precise question.* We add a convex-hull / extrapolation analysis (E2): per-synthetic-patient convex-hull membership in the d=18 subspace (LP feasibility), reporting the inside/outside fractions, and — importantly — showing that the out-of-hull synthetic patients still satisfy biological constraints and the mechanistic envelope, i.e., genuine novel-but-plausible phenotypes rather than pure interpolation or implausible extrapolation. Results [pending]: SA [interpolates primarily, with controlled local extrapolation along principal directions that remains plausible].

**Changes:** New Results paragraph + Supplement (E2). `\rtwo{...}`

---

### R2.5 — Demonstration of real-world utility
> Does synthetic data improve predictive models beyond the original cohort? Enhance power for subgroup analyses? Enable concrete new hypotheses? Even a limited quantitative demonstration would strengthen translational significance.

**[PLAN]** *We appreciate the push toward translational concreteness.* We strengthen the framing of the existing downstream-utility experiment (synthetic-calibrated mechanistic model predicting held-out real outcomes as well as the real-calibrated one) and connect it to two concrete use-cases the reviewer names: subgroup power (conditional amplification of PCOS/PE enables comparisons impossible at n=3/n=5) and hypothesis generation (the preserved condition-specific signatures suggest testable factor-level differences). We are careful, per R1.11, not to overclaim: the utility is enabling analyses that small n forecloses, not surpassing real data.

**Changes:** Results (Downstream Utility) + Discussion framing. `\rtwo{...}`

---

## Minor comments

### R2.minor1 — Consolidate figures
> Several figures could move to the supplementary material to improve readability.

**[PLAN]** *Agreed, and this dovetails with the new material.* We consolidate the main-text figure set, moving [list once finalized] to the supplement, so that the new simulation-benchmark figure and the ablation summary can be added without increasing main-text figure count. Net main-text figures: [target ≈ current or fewer].

**Changes:** Figure placement (floats.tex / supplementary.tex). `\rtwo{...}`

---

### R2.minor2 — Streamline Hopfield theory exposition
> Some sections (particularly the Hopfield network theory) could be streamlined to focus on the applied contribution.

**[PLAN]** See **R1.9** — addressed jointly by moving dense theory to the supplement, adding plain-language intuition, and (optionally) a schematic.

**Changes:** See R1.9. Violet. `\rboth{...}`

---

## Cross-reference: comment → planned change

| Comment | Type | Vehicle | Color |
|---|---|---|---|
| R1.1, R1.4 | New experiment | E1 PCA-space ablation suite | blue |
| R1.2, R2.3 | Prose (temper) | Mechanistic independence clarification | violet |
| R1.3 | Prose | Subgroup reliability caveats | blue |
| R1.5 | New experiment | E3 membership inference / privacy | blue |
| R1.6 | Prose | Irregular-timing / missingness limitations | blue |
| R1.7 | Prose | Tail-compression clinical importance | blue |
| R1.8, R2.1 | New experiment + framing | Simulation benchmark + cross-domain evidence | violet |
| R1.9, R2.minor2 | Writing | Theory streamlining + schematic | violet |
| R1.10 | Prose | Fibrinolytic/viscoelastic features | blue |
| R1.11 | Prose | Downstream smoothing caveat | blue |
| R1.12 | Prose | Computational scalability | blue |
| R1.M1 | Writing | Temper wording | blue |
| R2.2 | Prose | PCA limits / nonlinear embeddings | red |
| R2.4 | New experiment | E2 convex-hull / extrapolation | red |
| R2.5 | Framing | Downstream utility translational framing | red |
| R2.minor1 | Layout | Figure consolidation | red |

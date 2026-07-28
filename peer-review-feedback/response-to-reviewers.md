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

---

## General response to the editor and reviewers

We thank both reviewers for their careful and constructive reading. We are grateful for the recognition that combining statistical with mechanistic validation is a distinguishing strength of the work, and that learning from very small longitudinal cohorts is an important and underserved problem. The reviewers converged on four priorities: generalizability beyond the single K=23 cohort, disentangling the contribution of PCA from the stochastic-attention energy landscape, the degree of independence of the mechanistic validation, and the interpolation-versus-extrapolation and privacy behavior of the generator. We have addressed each with new work rather than with argument alone.

Four additions carry the revision. A PCA-space ablation (R1.1, R1.4) fits four alternative generators in the identical d=18 subspace and evaluates them on the same harness as SA. A convex-hull analysis (R2.4) asks whether generated patients lie inside or beyond the real cohort. A membership-inference and distance-to-closest-record analysis (R1.5), followed by a sampling-correctness investigation that ended in a closed-form characterization of the target distribution, establishes what the generated cohort reveals about its training set and why. And a simulation benchmark with known ground truth (R1.8, R2.1) sweeps sample size, dimensionality, latent rank and manifold curvature, which lets us test generalization in a way no single real cohort can.

Two of these did not come out the way we expected, and we have let the results stand rather than presenting them selectively. The ablation shows that the shared PCA subspace supplies much of the marginal and mechanistic fidelity, and that mixture and copula baselines fit in that subspace are competitive with SA on those axes. We therefore no longer claim that SA dominates, and we have moved the covariance-generalization argument to the simulation benchmark, where reproducing the training sample cannot substitute for recovering the population. The privacy analysis found a genuine memorization signal at the published operating point, and the follow-up investigation established that it is intrinsic to the target distribution at this cohort size rather than an artifact of the sampler or of the operating temperature. Pushed to its end, that investigation produced the sharpest result in the revision: the latent target our method defines is exactly a finite Gaussian mixture centered on the stored patients, so it can be drawn from analytically rather than by Langevin dynamics. Exact draws reproduce the reported cohort's behavior on every axis we evaluate, including the membership-inference signal, which settles the attribution and also removes the sequential cost factor from future generation. The reported cohort and every clinical result are unchanged; what changes is that we can now say precisely what distribution they came from. We state plainly that the method is not an anonymizer.

A theme emerged that we have made explicit in the Discussion. Several of the behaviors the reviewers asked about — the convex-hull result, the membership-inference result, the effective pattern count under conditioning, and the modest downstream gain — are consequences of a single fact, that 23 patients in a 216-dimensional space are sparse, rather than four separate properties of the method. Recognizing this shaped the revision. We now frame the clinical study as a proof of concept, temper "indistinguishable" to what the tests support at this sample size, and bound the translational claim to enabling analyses that small n forecloses.

We have also tempered the independence claim for the mechanistic validation (R1.2, R2.3), added the scope limits the reviewers identified around irregular visit timing, missingness, tail compression and the linear embedding, streamlined the Hopfield exposition with plain-language intuition while moving the denser derivation to the supplement (R1.9, R2.minor2), and consolidated the figure set so that the main text still carries seven figures despite the new experiments (R2.minor1).

---

# Reviewer 1

## Major comments

### R1.1 — Disentangle PCA preprocessing from the SA energy landscape (ablation baselines in PCA space)
> The regularization may already arise from the PCA dimensionality-reduction step rather than SA itself. Baselines do not disentangle PCA from the SA energy landscape. Suggests baselines operating in the same PCA-reduced space: (i) PCA+GMM, (ii) PCA+KDE, (iii) PCA+nearest-neighbor interpolation, (iv) PCA+diffusion/interpolation, (v) PCA+copula.

*Excellent and central question — this is the crux of attributing the result to SA rather than to dimensionality reduction, and our ablation answers it directly (and candidly).* We ran a PCA-space ablation suite: PCA+GMM, PCA+KDE, PCA+k-NN interpolation, and PCA+Gaussian copula, all fit in the identical d=18 PCA subspace and decoded through the identical evaluation harness (marginal MRE, cross-visit covariance, mechanistic cloud overlap/KS, and novelty/memorization). The result is nuanced and we report it exactly as measured. The shared PCA subspace does carry much of the marginal and mechanistic fidelity: PCA+GMM and PCA+Gaussian copula are competitive with SA on marginal error (median 1.1–1.3% vs. SA 1.2%) and on the mechanistic envelope (overlap 0.92–0.93 vs. SA 0.90). Two baselines fail cleanly on a single axis — PCA+k-NN memorizes (median novelty 0.05; 5% of samples exceed the novelty threshold, vs. 100% for SA), and PCA+KDE degrades both marginals (1.9% MRE) and mechanistic plausibility (overlap 0.82, below the within-real band). The other two baselines do not fail, and we will not claim otherwise. PCA+Gaussian copula matches SA in producing 100% novel samples, and its median novelty (0.46), like PCA+GMM's (0.50), sits close to SA's (0.51) rather than far below it. No method dominates. The empirical 216×216 cross-visit *correlation* distance is in fact larger for SA (0.64) than for the mixture and copula baselines (0.42), but that metric measures reproduction of the correlation matrix of the same 23 patients every generator was fit to, which at this sample size is a noisy estimate of the population rather than a test of generalization. We therefore make the covariance-*generalization* argument where reproduction cannot help — the controlled simulation benchmark with known population covariance (see R1.8/R2.1), in which SA recovers the true cross-visit structure in the n<p regime while the multivariate Gaussian does not. The honest reading for the reviewer's question is thus: the shared PCA subspace supplies much of the marginal and mechanistic fidelity and is available to any generator that samples it sensibly, and what SA adds over the competitive subspace baselines is the energy-landscape formulation itself — an operating point fixed by the attention-entropy inflection rather than selected by cross-validation on 23 patients, together with the inference-time multiplicity steering used for conditional subgroup generation, which the subspace baselines do not provide directly. Regarding (iv) diffusion: a score network cannot be trained at n=23; the k-NN interpolation baseline represents the interpolation family, and we note this explicitly.

**Changes:** New Results paragraph + Supplement section (E1 tables/figure), violet where it also answers R2.2. `\rone{...}`

---

### R1.2 — Mechanistic validation is not fully independent (shared cohort/variables/assumptions)
> The BZ2012 model operates on the same coagulation variables used for generation, is calibrated on the same cohort, and shares the same biological domain assumptions. Temper the "mechanistically indistinguishable" claim and clarify the residual coupling.

*We appreciate this point and agree the independence claim was overstated.* We temper the language throughout (abstract, results, discussion) and add an explicit statement of the residual dependence: the ODE parameters were calibrated on the same cohort (real → SA → synthetic → calibration), so the validation demonstrates *consistency with a model fitted to the same biological system*, not fully independent ground truth. We clarify that the strength of the test lies in the ODE being blind to the SA generation process and processing real and synthetic patients through the identical fixed mapping, and we soften "indistinguishable" to "not distinguishable by this model under these conditions." (We considered independent/alternative calibration; we address the point in prose rather than a new experiment, and flag independent mechanistic validation — including the fibrinolysis model under development — as future work.)

**Changes:** Abstract wording; Results (Mechanistic Consistency) clarification; Discussion limitations. Violet (shared with R2.3). `\rboth{...}`

---

### R1.3 — Statistical reliability of very small subgroups (PCOS n=3, PE n=5)
> Subgroup distributions may reflect interpolation between a handful of trajectories rather than genuine subgroup variability. Add discussion of reliability, explicit acknowledgment of under-sampling, and caution on downstream clinical interpretation.

*We agree, and we have made this caveat explicit and quantitative.* At these sizes the conditioned landscapes are severely under-sampled. The participation ratio is 4.6 for the PCOS cohort and 7.7 for Developed PE, so each cohort of 100 generated patients is a recombination of three and five real individuals respectively, and in the PCOS case magnitudes are drawn from only three empirical norms. The Results now state the consequence directly: those cohorts carry no more information about the conditions than the underlying patients do, and any downstream analysis treating them as 100 independent observations will understate uncertainty. We frame subgroup amplification as supporting hypothesis generation and power calculation rather than as a substitute for recruiting real rare-subgroup patients. The Methods additionally explain that reporting the participation ratio alongside a conditioned cohort is precisely what states how much independent material that cohort rests on.

**Changes:** Results (Conditional Generation) caveat; Discussion. `\rone{...}`

---

### R1.4 — Stronger and more diverse baselines
> CTGAN performs poorly at small n; TVAE was not configured for concatenated longitudinal generation; MVN is disadvantaged in rank-deficient settings. Recommend stronger baselines designed for low-sample tabular generation or manifold interpolation.

*A fair point.* The PCA-space ablation suite (R1.1) adds exactly this class of fairer, low-sample-appropriate baselines — density estimation, mixture models, copulas, and manifold interpolation — all operating in the same reduced space as SA rather than being handicapped by the full-dimensional or per-visit formulations. We report their results alongside MVN/CTGAN/TVAE and, per R1.1, present them candidly: given the shared PCA subspace the mixture-model and copula baselines are competitive with SA on marginal and mechanistic fidelity, while interpolation (k-NN) and kernel-density baselines each fail a specific axis (novelty and marginal/mechanistic fidelity, respectively). What distinguishes SA is generating fully novel patients (100% above the novelty threshold) while retaining that fidelity, rather than uniformly outperforming every baseline on every metric.

**Changes:** Same as R1.1. `\rone{...}`

---

### R1.5 — Privacy and memorization (membership inference, re-identification)
> Given K=23 and direct interpolation between stored memory patterns, evaluate membership inference risk and patient re-identification concerns.

*An important omission, and pursuing it changed how we present the method. Thank you.* We added a privacy analysis, and because the first result was a genuine memorization signal rather than the reassuring one we had expected, we followed it with a full sampling-correctness investigation to attribute the cause. We report the complete picture honestly.

The membership-inference risk is real at the published operating point. A nearest-neighbor membership-inference attack evaluated over repeated train/holdout splits reaches a mean AUC of 0.97, and the distance-to-closest-record distribution places synthetic patients slightly closer to the real records than the real records are to one another (median DCR synthetic-to-real 13.8 against real-to-real 15.9 in the standardized space). Rather than explain this away, we established what drives it. First, it is not a sampling artifact: a four-rung diagnostic ladder that replaces the published anchor-initialized single-endpoint Langevin scheme with sphere-initialized, burn-in-controlled, thinned and pooled Metropolis-adjusted sampling (acceptance near 0.999, so discretization bias is negligible) leaves the attack AUC essentially unchanged at 0.96 and raises the DCR only marginally, never reaching the real-to-real baseline. Second, it is not removed by the inverse temperature: sweeping beta across a 40-fold range around the operating point raises the DCR monotonically as beta falls but saturates well below the real-to-real baseline at every value tested, while the cross-visit covariance fidelity degrades.

Third, and most directly, it is not a property of finite-time dynamics at all. Both of the checks above still compare one Markov chain against another, so we extended the correctness analysis to an analytic reference. Doing so showed that the latent target is not merely approximated by our sampler but is available in closed form: for the unit-normalized memories the manuscript uses, the weighted Hopfield energy induces exactly a finite isotropic Gaussian mixture whose components are centered on the stored patients, with covariance beta-inverse times the identity and mixture weights proportional to the multiplicities. That distribution can be sampled ancestrally, by choosing a stored patient according to its weight and adding isotropic noise, with no chain of any kind. Passing such exact draws through the identical magnitude, reconstruction and de-standardization pipeline reproduced the reported cohort's evaluated behavior across the entire harness: marginal error 1.7% against 1.3%, cross-visit Frobenius 0.60 against 0.65, mechanistic cloud overlap 0.91 against 0.89, median novelty 0.50 against 0.51, and median distance to closest record 13.9 against 14.0, with each value falling inside the range spanned by the four ladder rungs. One hundred independently seeded exact cohorts place these quantities in narrow intervals around those values, so the agreement is not a favorable seed. Evaluated over the same 40 train/holdout splits, the membership-inference AUC under exact sampling is 0.98, if anything marginally higher than under the published scheme. Because ancestral sampling has no initialization, no burn-in and no discretization error, this locates the memorization in the target distribution itself rather than in any aspect of how that distribution is sampled.

The memorization is therefore intrinsic to the target at this cohort size, where every stored patient is the center of a mixture component, not a defect of the sampler or a poorly chosen operating point.

We now frame this correctly rather than as a favorable finding. Stochastic attention at a usable operating point is not an anonymizer, and we say so plainly: it is not differentially private, K=23 admits no formal guarantee, and neither a better sampler nor a lower temperature yields a configuration that is simultaneously private and faithful. What bounds the practical risk is that the source cohort is fully de-identified, so no auxiliary dataset exists against which a 216-dimensional assay profile could be re-identified, and the intended uses (mechanistic calibration and hypothesis generation) do not require releasing patient-level records. We add these results and this framing to the manuscript and detail the diagnostics in a new Supplement subsection.

**Changes:** Rewritten privacy treatment in Results/Discussion + new Supplement subsection "Sampling correctness and privacy diagnostics" (E3, E3b, the MCMC ladder, the beta-privacy sweep, and the exact analytic reference) + the closed-form target in Methods and the corresponding correction to the computational-cost and interpolation language. `\rone{...}`

---

### R1.6 — Robustness to irregular visit timing, missingness, variable-length trajectories
> The concatenation strategy assumes visits are temporally and physiologically aligned across patients. Discuss robustness to irregular timing, incomplete trajectories, and variable-length records.

*A valuable point about real-world deployment.* We add a limitations paragraph to the Discussion. The concatenation that makes cross-visit structure available also constrains the data the method accepts: stacking three visits into one vector assumes every patient has every visit and that visits are aligned across patients, which holds in this cohort by design but not in most observational data. We state the three consequences separately. Irregular timing would place physiologically different states in the same coordinate, so the principal components would mix visit structure with between-patient variation. Missing assays or visits have no natural representation in a fixed-length vector and would require imputation before the embedding, which at this cohort size would itself be a substantial source of uncertainty. Variable-length records fall outside the formulation altogether and would need either padding to a common length or a sequence-aware extension in which attention operates over visits as well as over patients. We note that these are scope boundaries of the present formulation rather than properties of the energy landscape, which is indifferent to what the coordinates mean.

**Changes:** Discussion (limitations). `\rone{...}`

---

### R1.7 — Clinical importance of tail compression
> SA smooths distributional extremes, but rare/extreme states are often the phenomena of interest. Discuss effect on downstream mechanistic modeling, implications for rare-event hypothesis generation, and whether the framework systematically underestimates severe pathology.

*We agree this deserves more weight than we gave it.* We expand the tail-compression discussion to cover its two sources (PCA truncation of higher moments; the finite empirical magnitude distribution of K=23 norms) and, more importantly, its clinical consequence: if the generated cohort under-represents the severe end of a distribution, then estimates of how often a severe phenotype occurs, or of how extreme it can become, are biased toward the center. We connect this to the extrapolation analysis but are careful not to let it off the hook — SA does place patients outside the convex hull of the real cohort, so it is not confined to the observed range, but producing out-of-range points is not the same as reproducing the frequency of extreme ones, and we say so explicitly. We therefore recommend against using SA cohorts to estimate extreme-quantile risk without real anchoring data. We also note a scope limit of our own evidence: the mechanistic checks compare distributions within a 5th-to-95th-percentile band and are therefore, by construction, a test of the central mapping rather than of the tails, so they should not be read as reassurance about tail behavior.

**Changes:** Discussion (expand tail-compression paragraph); links to E2. `\rone{...}`

---

### R1.8 / R2.1 — Generalizability beyond a single cohort
> Validation is on a single K=23 cohort from a single domain. Encourage external validation or explicit limitation; recommend proof-of-concept framing. (R2: a second dataset, or a synthetic benchmark / simulation study demonstrating robustness to varying data characteristics.)

*This is the most important shared concern and we address it on two fronts.* (1) **Simulation benchmark (new main-text figure and Supplement):** we add a controlled study with known ground-truth structure sampled across varying sample size n, dimensionality p, and intrinsic rank. On a low-rank Gaussian ground truth SA recovers the known population covariance more accurately than a multivariate Gaussian throughout the rank-deficient n<p regime that contains the cohort's own shape (n=23, p=216, where the Gaussian's relative covariance error is about 1.5 times SA's, and SA is more accurate in 35 of the 39 Gaussian n<p cells); because the target covariance is fixed in advance, this measures generalization rather than reproduction of the particular training draw, a distinction our empirical cross-visit comparison could not make. Subsampling the real cohort to as few as eight patients degrades recovery of the full-cohort structure gracefully rather than catastrophically (pooled marginal error rising from roughly 2% at n=20 to 5% at n=8, mechanistic plausibility flat near 0.90). A nonlinear S-curve ground truth locates the curvature (roughly 0.25 to 0.5) at which the linear-PCA assumption begins to place points off the manifold, an honest applicability boundary; we report candidly that on that nonlinear target neither method dominates on every metric (SA stays closer to the manifold, the Gaussian recovers the gross covariance better). (2) **Cross-domain evidence:** the identical SA machinery, with β* predicted from PCA dimension alone, has been applied in separate studies to discrete protein-sequence generation from small family alignments and to multiplicity-weighted conditional steering — evidence that the geometric mechanism is not specific to clinical coagulation data. We also reframe the clinical study explicitly as a proof-of-concept and strengthen the limitations accordingly.

**Changes:** New Results subsection "Robustness Across Data Regimes" + main-text figure (`fig:sim-recovery`) + Supplement diagnostics; Introduction and Discussion cross-domain framing elevated from aside to argument; proof-of-concept framing in Abstract and Discussion. Violet. `\rboth{...}`

---

### R1.9 / R2.minor2 — Theory exposition for a biomedical audience
> The entropy-inflection (β) selection, participation-ratio interpretation, and geometric interpretation of the energy landscape are hard for a broad biomedical audience; add intuition or a schematic. (R2: streamline the Hopfield theory to focus on the applied contribution.)

*A helpful suggestion for accessibility.* We add plain-language intuition for both quantities the reviewers singled out. For β, the Methods now explain that it sets how sharply the sampler commits to one stored patient: at low β the attention weights are nearly uniform, every patient contributes about equally, and samples collapse toward the cohort average; at high β a single patient dominates and samples reproduce that individual; β\* is our estimate of the transition between the two, where stored patients stop blurring into one average and begin acting as distinct attractors. For the participation ratio we pose the question it answers directly — if attention is spread unevenly across 23 patients, how many patients is the generator effectively drawing on? — and note that reporting it alongside a conditioned cohort states how much independent material that cohort rests on. We have moved the mechanical derivation (entropy definition, sweep parameters, finite-difference procedure) into a new Supplementary Methods subsection, keeping the concept, the domain-agnostic argument and the resulting values in the main text. We did not add a schematic figure: the main text carries seven figures and Reviewer 2 asks us to consolidate rather than expand, so we judged the prose intuition the better way to meet this request without displacing a results figure.

**Changes:** New/streamlined Methods exposition; possible schematic figure; denser theory → Supplement. Violet. `\rboth{...}`

---

### R1.10 — Unvalidated feature categories (fibrinolytic, viscoelastic)
> Mechanistic validation covered only thrombin generation, not fibrinolytic or viscoelastic features, which are a substantial portion of the feature space. Discuss how unvalidated categories may influence the generated manifold and whether they contribute meaningfully to the energy landscape; expand future validation.

*Thank you — we clarify this.* We already note that ROTEM/fibrinolytic features are generated within the concatenated vector with comparable MREs but are not mechanistically validated; we expand this to state that these features do contribute to the PCA-derived energy landscape (they load onto retained components) and therefore shape generation, so their statistical fidelity matters even without a mechanistic check. We point to the fibrinolysis mechanistic model under development as the route to closing this gap.

**Changes:** Results (Marginal Plausibility) + Discussion. `\rone{...}`

---

### R1.11 — Downstream improvement may reflect smoothing/variance suppression
> The synthetic-calibrated model's improvement over the real-calibrated model may partly reflect smoothing/regularization that suppresses physiologically meaningful heterogeneity.

*We agree and already flagged the smoother loss landscape as the likely mechanism; we sharpen this.* We state explicitly that the modest downstream gain (0.94× ratio) is consistent with variance reduction from N=100 vs. K=23 training points and should not be read as SA adding biological information; it demonstrates that SA preserves *sufficient* structure for calibration, not that synthetic data is superior to real data. We connect this to the tail-compression caveat (suppressing heterogeneity is the same phenomenon that reduces optimization variance).

**Changes:** Results (Downstream Utility) + Discussion. `\rone{...}`

---

### R1.12 — Computational scalability
> Discuss how the framework behaves as cohort size increases, whether the energy landscape becomes unstable at larger K, and how complexity scales with dimensionality.

*A useful addition, and working through it corrected an expectation of our own.* On cost, the scaling follows directly from the update rule: each Langevin step evaluates the attention weights over the stored patterns and forms their weighted combination, so work per step is proportional to K×d and a sample costs that product times the number of steps T. The method is linear in both the cohort size and the retained dimension, and at very high feature counts the principal-component decomposition rather than the sampling becomes the practical bound.

On behavior we are more careful than we had planned to be. We had expected to argue that the simulation benchmark supports favorable scaling to larger cohorts. It does not show that. What it shows is that SA's advantage over a fitted multivariate Gaussian is largest when the sample is smallest and *narrows* as n approaches p, where the Gaussian becomes well specified. The honest conclusion is that the method is positioned for the small-cohort regime rather than as a general-purpose generator, and the Discussion now says so. We do expect the memorization documented under R1.5 to ease as K grows, because basins that currently stand apart would begin to overlap and generation would interpolate across several patients rather than remaining near one, but we mark that as an expectation rather than a result: this cohort admits only subsampling downward, and we did not test it.

**Changes:** Discussion (+ brief note tied to simulation benchmark). `\rone{...}`

---

## Minor comments

### R1.M1 — Temper strong wording
> Claims of "clinically useful synthetic cohorts" and "statistically indistinguishable" may be too strong given the sample size.

*Agreed.* We adopt more conservative wording throughout: "statistically indistinguishable" → "not distinguishable by the tests applied at this sample size," and we qualify "clinically useful" as demonstrated for the modeling use-cases tested (mechanistic calibration, hypothesis generation), within a proof-of-concept framing.

**Changes:** Abstract, Introduction, Discussion, Conclusion wording. `\rone{...}`

---

# Reviewer 2

## Major comments

### R2.1 — Generalizability beyond a single dataset
> All validation is on a single cohort (K=23). Application to a second dataset, or a synthetic benchmark / simulation study demonstrating robustness to varying data characteristics, would significantly strengthen the manuscript.

See **R1.8** — addressed jointly via the simulation benchmark (now complete: SA recovers a known population covariance across the n<p regime, degrades gracefully under real-cohort subsampling, and has a quantified nonlinear applicability boundary), elevated cross-domain protein evidence, and explicit proof-of-concept framing. *We appreciate the reviewer's recognition that additional clinical cohorts are hard to obtain; the simulation study is designed precisely to demonstrate robustness across data characteristics without one.*

**Changes:** See R1.8. Violet. `\rboth{...}`

---

### R2.2 — Role and limitations of PCA; nonlinear embeddings
> Discuss scenarios where PCA may be insufficient; comment on whether nonlinear embeddings could be incorporated and the expected trade-offs between interpretability, fidelity, and flexibility.

*A thoughtful point, and the simulation benchmark now lets us answer it quantitatively rather than only in discussion.* On a nonlinear S-curve ground truth, SA's linear principal-component subspace begins to place generated points off the true manifold as curvature grows, with the off-manifold residual crossing twice its linear-baseline value at a curvature of roughly 0.25 to 0.5; notably SA stays closer to the manifold than a multivariate Gaussian even past that point, while the Gaussian recovers the gross covariance better, so neither is uniformly preferable. This bounds where linearity starts to matter. We add discussion of PCA's limits: it captures linear structure, so strongly nonlinear manifolds (branching trajectories, threshold effects) are imperfectly represented and can contribute to tail compression. We discuss nonlinear embeddings (kernel PCA, autoencoders, diffusion maps) as extensions, with the trade-off that they sacrifice the interpretability and the closed-form β* ≈ √d prediction that the linear PCA subspace provides, and can be unstable to fit at n=23, so linearity is a deliberate choice for the small-cohort regime, not an oversight. The E1 ablation results inform which failure modes are PCA-driven versus sampler-driven.

**Changes:** Discussion (PCA limitations paragraph); ties to E1. Violet where shared with R1.1. `\rtwo{...}`

---

### R2.3 — Mechanistic validation: shared calibration dependence
> The ODE is calibrated on the same real dataset used to generate the synthetic data, creating partial dependence. Clarify that the validation shows consistency with a model fitted to the same cohort rather than fully independent mechanistic fidelity; if feasible, use independently calibrated parameters or an alternative model.

See **R1.2** — addressed jointly. *We really appreciate this precise articulation of the dependence structure* and have made it explicit in the Discussion, tempered the independence language, and flagged independent/alternative-model calibration as future work.

**Changes:** See R1.2. Violet. `\rboth{...}`

---

### R2.4 — Interpolation vs. extrapolation (convex hull, novel phenotypes)
> Does SA generate novel out-of-sample variation or primarily interpolate? To what extent do generated patients lie within the convex hull of observed data? Does it produce biologically plausible but previously unseen phenotypes?

*An excellent and precise question, and the analysis produced an instructive answer once we controlled for sample size.* We added a convex-hull / extrapolation analysis (E2): per-synthetic-patient convex-hull membership in the d=18 PCA subspace (solved as a feasibility/least-distance program), with out-of-hull points additionally checked for biological-constraint satisfaction and mechanistic-envelope consistency. Two findings. First, the criterion does not discriminate at this sample size. All 100 synthetic patients fall outside the hull of the 23 real patients, but so does every real patient with respect to the other 22. When the comparison is matched for severity, removing each point's nearest vertex on both sides so that a real and a synthetic patient face equally depleted hulls, the distance to the hull is not distinguishable between them (median 11.0 for synthetic patients against 11.9 for held-out real patients, Mann–Whitney p = 0.07). We flag that an unmatched comparison, which scores synthetic patients against all 23 vertices while each held-out real patient necessarily loses its own, overstates the difference in SA's favor; we do not rely on it. A geometric diagnostic explains why nothing of either kind falls inside: with 23 points in 18 dimensions the hull is nearly a simplex, reaching the patients' typical radius (13.8) only along the 23 vertex directions while extending about 2.0 in a generic direction, and the sampler additionally returns every sample to the unit sphere, where a point inside the hull of unit-norm patterns would have to be one of those patterns (all 100 samples have a radial exit factor below 1, maximum 0.23). In-hull membership is therefore structurally unavailable to this generator rather than merely uncommon, and we report the 0/100 as a limit on what the hull criterion resolves at K=23 rather than as evidence of runaway extrapolation. Second, for the extrapolating points we asked whether they remain plausible: 81% satisfy all biological constraints (the failures are one specific and interpretable mode — a single hormone assay, estradiol, going slightly negative under extrapolation, which we note as an honest limitation and a candidate for a non-negativity constraint in future work), and at the level of individual mechanistic features the calibrated ODE ratios remain within the real cohort's band for the large majority of checks (per-feature agreement 83–92%; the strict all-features-and-all-visits joint criterion is met by 41%, reflecting how demanding that conjunction is rather than gross mechanistic implausibility). The reading is that SA does not merely interpolate within the observed data — it generates genuinely novel phenotypes — but its extrapolation is bounded to the same neighborhood the real data themselves occupy and remains biologically and mechanistically plausible for the large majority of generated patients.

**Changes:** New Results paragraph + Supplement (E2). `\rtwo{...}`

---

### R2.5 — Demonstration of real-world utility
> Does synthetic data improve predictive models beyond the original cohort? Enhance power for subgroup analyses? Enable concrete new hypotheses? Even a limited quantitative demonstration would strengthen translational significance.

*We appreciate the push toward translational concreteness.* We strengthen the framing of the downstream-utility experiment and connect it to the two use-cases the reviewer names, while bounding the claim. The Results now state that the method enables analyses the original sample size forecloses — calibrating a mechanistic model without spending the real cohort on training, and comparing subgroups that hold three and five patients — and that the preserved condition-specific signatures point to factor-level differences a larger study could test. We are equally explicit about what it does not do: it does not substitute for additional patients, and none of these results speak to clinical deployment, which would require prospective validation against real clinical decisions. Per R1.11 we also state that the 0.94× downstream ratio is a variance-reduction effect rather than evidence of new biological content.

**Changes:** Results (Downstream Utility) + Discussion framing. `\rtwo{...}`

---

## Minor comments

### R2.minor1 — Consolidate figures
> Several figures could move to the supplementary material to improve readability.

*Agreed, and this dovetails with the new material.* We moved the SA-versus-MVN PCA-by-visit figure to the supplement and added the simulation-benchmark figure (`fig:sim-recovery`) in its place, so the main text carries seven figures, unchanged from the original submission despite the new experiments. The ablation, convex-hull, membership-inference and sampling diagnostics are reported as supplementary tables and figures rather than as additional main-text panels.

**Changes:** Figure placement (floats.tex / supplementary.tex). `\rtwo{...}`

---

### R2.minor2 — Streamline Hopfield theory exposition
> Some sections (particularly the Hopfield network theory) could be streamlined to focus on the applied contribution.

See **R1.9** — addressed jointly by moving dense theory to the supplement, adding plain-language intuition, and (optionally) a schematic.

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

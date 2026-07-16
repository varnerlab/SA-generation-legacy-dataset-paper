# Design: SA Paper Major Revision (npj Systems Biology and Applications)

**Date:** 2026-07-16
**Manuscript:** Validated Synthetic Patient Generation for Small Longitudinal Cohorts: Coagulation Dynamics Across Pregnancy
**Decision:** Major revision (two reviewers, both positive)
**Companion doc:** `peer-review-feedback/response-to-reviewers.md` (reviewer-facing, comment-by-comment)

This spec is the *execution blueprint*: the technical specs of the new experiments, the manuscript change-map, the figure plan, and the execution order. The response doc maps every comment to a vehicle; this doc says how each vehicle gets built.

---

## 1. Goal and posture

Both reviewers recommend major revision but praise the core (statistical + independent mechanistic validation). No reviewer attacks the method. The revision must satisfy four shared priorities and a set of prose/framing asks, without acquiring a second real clinical dataset. Overall posture: reframe as a rigorous proof-of-concept, add four experiments that are largely reuse of existing infrastructure, and temper overclaims.

## 2. Scope decisions

**In scope — four new experiments:**
- **E1** PCA-space ablation suite (R1.1, R1.4)
- **E2** Convex-hull / extrapolation analysis (R2.4)
- **E3** Membership-inference / privacy analysis (R1.5)
- **SIM** Simulation benchmark, three arms (R1.8/R2.1, with R2.2 + R1.12 + R1.7 coverage)

**In scope — prose/writing/layout:** temper mechanistic-independence and strength wording; subgroup-reliability caveats; irregular-timing/missingness limitations; expand tail-compression; PCA-limits/nonlinear-embeddings discussion; fibrinolytic/viscoelastic clarification; downstream-smoothing caveat; scalability; theory streamlining + schematic; figure consolidation; proof-of-concept framing.

**Explicitly OUT of scope (scope guard):**
- No second real clinical cohort (SIM benchmark substitutes).
- No alternative/independent mechanistic recalibration experiment — R1.2/R2.3 handled in prose only.
- No new nonlinear SA generator — nonlinearity appears only as a *ground-truth probe* in SIM arm 3.
- No differential-privacy mechanism — E3 is empirical only.

## 3. New experiment specs

All experiments reuse the canonical pipeline (α=0.01, β\*=2.94, T=2000, PCA 95%→d=18, N=100, seed=42) and the existing evaluation harness. Results feed the `[PENDING]` slots in the response doc.

### E1 — PCA-space ablation suite
**Purpose:** attribute joint-structure and mechanistic fidelity to the SA energy landscape, not PCA preprocessing.
**Setup:** reuse the canonical d=18 PCA embedding of the 23×216 concatenated real data. Each baseline fits *in this same subspace*, samples N=100, and decodes via inverse PCA to 216-D. Fair-by-construction: every method shares SA's preprocessing.
**Baselines:**
- PCA+GMM (1–3 components, regularized/shrunk covariance; expect near-degeneracy at n=23 — documented as illustrative of n<p).
- PCA+KDE (Gaussian kernel, Scott/Silverman bandwidth; expect bandwidth instability).
- PCA+k-NN interpolation (SMOTE-style: real point + random k-NN neighbor, interpolate at λ∈[0,1]; expect low novelty / memorization-adjacent).
- PCA+Gaussian copula (per-component marginal CDF + Gaussian correlation; note PCA components are near-decorrelated, itself instructive).
**Evaluation (reuse):** pooled marginal MRE (median + bootstrap CI); 216×216 cross-visit covariance relative Frobenius error; mechanistic cloud overlap + KS on BZ2012 ratio distributions (decode → run calibrated BZ2012); novelty score + NN distance (memorization).
**Output:** methods×metrics table (supplement) + compact main-text summary extending the current CTGAN/TVAE baseline paragraph; optional bar-chart figure (supplement) following style conventions.
**Expected narrative:** each baseline fails a *different* axis; SA is the only method passing marginals + covariance + mechanistic envelope + novelty simultaneously.

### E2 — Convex-hull / extrapolation
**Purpose:** quantify interpolation vs. extrapolation; identify novel-but-plausible phenotypes.
**Method:** in d=18 PCA space, per synthetic point test convex-hull membership w.r.t. the 23 reals via LP feasibility (nonneg weights summing to 1 reproducing the point). Report inside/outside fractions and signed distance-to-hull for outside points; low-D hull plots (PC1–PC2, PC1–PC3) for visualization.
**Plausibility of out-of-hull points:** run out-of-hull synthetics through (i) biological-constraint checks (existing `validate_biological_constraints`), (ii) mechanistic envelope (within-real ratio range / cloud overlap). Report fraction remaining plausible.
**Output:** supplement figure (hull projection + plausibility summary) + short Results paragraph.
**Expected:** majority inside (interpolation), meaningful minority just-outside along principal directions, out-of-hull points remain plausible → controlled novel extrapolation, not memorization or implausible extrapolation.

### E3 — Membership inference / privacy
**Purpose:** quantify re-identification / membership risk at K=23.
**Metrics:**
- DCR distribution: per-synthetic distance to closest real (synthetic→real) vs. real→real NN distances; optional holdout split (generate from a subset, compare DCR to held-out vs training reals — small at K=23, documented).
- Nearest-neighbor MIA: attacker labels a record "member" if distance to nearest synthetic < threshold; sweep threshold → ROC AUC over members (training reals) vs non-members (holdout/reference). AUC≈0.5 → no leakage.
**Output:** supplement figure (DCR histograms + ROC) + short Results/Discussion paragraph.
**Caveats (stated):** SA is not differentially private; K=23 admits no formal guarantee; this is empirical evidence.
**Expected:** AUC≈0.5, DCR ≥ real–real → low empirical re-identification risk.

### SIM — Simulation benchmark (three arms)
**Purpose:** demonstrate robustness across data characteristics without a second real cohort; empirically bound PCA's linear assumption; support scalability and tail claims.
- **Arm 1 — Linear low-rank Gaussian GT:** r-dim latent factor model → p-dim observations with specified cross-visit block covariance and marginal transforms (subset heavy-tailed, e.g. lognormal/t). For a grid over (n, p, r) and tail settings, draw a cohort of n, run SA and MVN, and measure recovery of the *known* population covariance (relative Frobenius) and marginals (MRE/KS) against a large fresh GT test draw. Produce a recovery phase-diagram over n/p and rank.
- **Arm 2 — Real subsampling:** subsample the real 23 → n∈{20,15,10,8}, regenerate with SA, measure recovery of the full-23 structure (covariance, marginals, mechanistic overlap). Degradation curve.
- **Arm 3 — Nonlinear manifold GT:** known nonlinear ground truth (e.g. S-curve / branching trajectory) embedded in p-dim; run SA; quantify residual off-manifold error vs. the linear case, locating where linear-PCA-based SA fidelity breaks.
**Output:** ONE main-text figure (phase diagram / recovery curves) + supplement (nonlinear-arm detail, subsampling curves). New Results subsection "Robustness across data regimes."
**Expected:** SA robust across linear regimes incl. n<p (MVN fails); graceful degradation with n; nonlinear arm shows the expected linear-subspace limitation (honest applicability boundary → evidences R2.2).

## 4. Manuscript change-map

Per-comment detail lives in the response doc. Section-level summary:

- **Abstract:** temper "indistinguishable"/"clinically useful"; add proof-of-concept + generalizability-evidence phrasing. (R1.M1, R1.2, R1.8)
- **Introduction:** elevate cross-domain SA evidence from aside to explicit generalizability argument; proof-of-concept framing. (R1.8/R2.1)
- **Results:** extend baseline paragraph with E1; new "Robustness across data regimes" subsection (SIM); short E2 and E3 paragraphs; subgroup-reliability caveat; fibrinolytic/viscoelastic clarification; downstream-smoothing caveat. (R1.1/1.3/1.4/1.5/1.10/1.11, R2.4/2.5, SIM)
- **Discussion:** mechanistic-independence clarification; PCA-limits/nonlinear-embeddings; expanded tail-compression; irregular-timing/missingness limitations; scalability; strengthened limitations + proof-of-concept. (R1.2/1.6/1.7/1.12, R2.2/2.3/2.5)
- **Methods:** streamline Hopfield/β exposition, add plain-language intuition; move dense theory to supplement; add schematic figure. (R1.9/R2.minor2)
- **Conclusion:** align tempered wording. (R1.M1)
- **Preamble:** add `\rone/\rtwo/\rboth` color macros + `\reviewmode` toggle.

## 5. Figure plan (reconcile R1 "add" vs. R2 "consolidate")

Current main-text figures (7): bio-correlations, pregnancy-progression, cross-visit-corr, pca-by-visit, conditioned-features, mechanistic-combined, downstream-utility.
- **Add to main text:** SIM phase-diagram (+1); theory schematic (+1, optional).
- **Demote to supplement:** `pca-by-visit` (dispersion point already carried by cross-visit-corr + eigenvalue spectrum); one further demotion if the schematic is added (candidate: fold pregnancy-progression into supplement or merge panels).
- **Supplement (new):** E1 ablation table/figure, E2 hull figure, E3 DCR/ROC figure, SIM nonlinear + subsampling detail.
- **Target:** main-text figure count ≤ 7 (satisfies R2.minor1).
- **Style (per author conventions):** consistent colors across figures, no plot titles, legends on all panels, gray scatter backgrounds, visible black error bars, statistical tests annotated on comparison figures.

## 6. Response-to-reviewers workflow

- Living doc: `peer-review-feedback/response-to-reviewers.md`, appreciative tone, comment→response→changes, color-coded.
- Fill `[PENDING]` slots as E1/E2/E3/SIM produce numbers; write general response last.
- Render final to `.docx` via `pandoc response-to-reviewers.md -o response-to-reviewers.docx`.
- Manuscript: colored tracked changes via `\rone` (blue), `\rtwo` (red), `\rboth` (violet); flip `\reviewmodefalse` for camera-ready.

## 7. Execution sequence and dependencies

1. **Prose-only edits** (no result dependency): mechanistic-independence tempering, subgroup caveats, irregular-timing, tail-compression, PCA-limits, fibrinolytic/viscoelastic, downstream-smoothing, scalability, wording, theory streamlining. Can start immediately, in parallel with experiments.
2. **E1, E2, E3** — reuse existing SA cohort + harness; E2/E3 are analysis-only, E1 adds four baseline samplers. Run first (cheap, parallelizable).
3. **SIM** — the only substantial new code (three arms); independent of E1–E3.
4. **Figures** — build new figures after SIM/E1 results; then consolidate/demote.
5. **Response doc** — fill pending numbers, write general response, render to docx.
6. **Final build** — add color macros; integrate all edits; rebuild `main.pdf` + `supplementary.pdf`; flip `\reviewmode` off for a clean camera-ready pass.

## 8. Success criteria

- Every R1 (13) and R2 (7) comment has a concrete response + a corresponding manuscript change.
- E1: SA uniquely passes all four fidelity axes; each PCA baseline fails ≥1.
- E2: quantified inside/outside-hull fractions; out-of-hull points demonstrably plausible.
- E3: MIA AUC ≈ 0.5; DCR ≥ real–real.
- SIM: robustness demonstrated across linear regimes incl. n<p; graceful degradation; honest nonlinear boundary.
- Main-text figure count not increased (R2.minor1 satisfied).
- Overclaims tempered; theory streamlined; proof-of-concept framing consistent.
- `response-to-reviewers.docx` generated; colored-diff manuscript compiles; `\reviewmode` toggle verified.

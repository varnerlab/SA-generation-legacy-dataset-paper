# SA Paper Major Revision — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Execute the npj Systems Biology & Applications major revision — four new experiments (E1 PCA-space ablation, E2 convex-hull/extrapolation, E3 membership-inference/privacy, SIM three-arm simulation benchmark), the figure reshuffle, the prose/framing edits, and the response-to-reviewers document — so every R1 (×13) and R2 (×7) comment has a concrete response and a corresponding manuscript change.

**Architecture:** Experiments reuse the canonical SA pipeline and evaluation harness already in `code/` (no method changes); each new experiment is a standalone script under `code/experiments/` that loads the canonical serialized memory, computes results, and writes a CSV (+ optional figure) into `code/data/` / `code/figs/`. Manuscript edits are color-coded LaTeX applied to both the `paper/` build and the `arxiv/` mirror. The response doc's `[PENDING]` slots are filled from the experiment CSVs, then rendered to `.docx`.

**Tech Stack:** Julia (Plots/GR, MultivariateStats, Distributions, HypothesisTests, JLD2; adds JuMP + HiGHS for E2 LP), LaTeX (pdflatex/bibtex via `paper/Makefile`), pandoc (response doc → docx). Verification uses the repo's own writing skills: `/style-check`, `/review-section`, `/audit-magic-numbers`.

## Global Constraints

- **Canonical SA hyperparameters (copy verbatim into every experiment):** `α = 0.01`, `β* = 2.94` (recomputed per cohort via `find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1,1000.0))`), `T = 2000`, PCA `pratio = 0.95` → `d = 18`, `N = 100`, `Random.seed!(42)`.
- **Run all Julia from `code/`.** Every experiment script begins with `include(joinpath(@__DIR__, "..", "Include.jl"))`. Scripts that need the mechanistic model additionally `using HockinMannModel, HypothesisTests, Optim`.
- **Reuse, do not reinvent, the harness.** Canonical artifact: `code/data/full_longitudinal_memory.jld2` with keys `X̂` (unit-norm PCA memory, `d×23`), `pca_concat` (fitted `PCA` model), `std_params_concat::StandardizationParams`, `pca_norms` (empirical magnitudes, len 23), `complete_subjects`, `subject_conditions`, `kept_cols`, `n_assays`, `df_clean`. Raw PCA coordinates of the 23 reals: `Y = X̂ .* pca_norms'` (`d×23`). Key reuse functions: `decode_sample` (`Patient.jl:364`), `feature_summary_comparison` (`Patient.jl:600`, gives `Mean_Rel_Error`), `sample_novelty` (`Utilities.jl:52`), `find_entropy_inflection` (`Utilities.jl:254`). Mechanistic envelope reuse: `CALIBRATED_OVERRIDES` + `run_patient` + `to_clinical` from `validate_mechanistic_plausibility.jl:79-126`; overlap/KS pattern from `regen_mechanistic_figures.jl:104-109` (`quantile(real_ratios,[0.05,0.95])` band → `frac_in_range`; `ApproximateTwoSampleKSTest`). Cross-visit Frobenius pattern from `validate_cross_visit_covariance.jl:108-109` with `safe_cor`. Pooled-median + bootstrap-CI pattern from `paper_summary_table.jl:100-112` (seed 42, 10 000 resamples, percentile indices 250/9750).
- **Every prose edit is applied to BOTH copies:** `paper/sections/<file>.tex` AND `arxiv/sections/<file>.tex` (they are independent copies, not symlinks). Line numbers in tasks below refer to `paper/`; find the matching text in `arxiv/` by content.
- **Tracked-changes convention:** `\rone{...}` = blue (R1), `\rtwo{...}` = red (R2), `\rboth{...}` = violet (both). A `\reviewmode` boolean toggles all three to black for camera-ready. Macros defined in `paper/main.tex`, `paper/supplementary.tex`, AND `arxiv/main.tex` (each has an independent preamble).
- **Figure style (author convention):** consistent per-visit palette across figures; NO plot titles; legends on all panels; gray scatter background `RGB(0.97,0.97,0.98)`; visible black error bars; statistical tests annotated on comparison figures. SA color `RGB(0.20,0.50,0.72)`, MVN color `RGB(0.85,0.45,0.25)`. Plots.jl (GR backend); Makie is not installed.
- **Main-text figure count must stay ≤ 7** (R2.minor1). New main-text figures are offset by demoting `fig:pca-by-visit` to the supplement.
- **Commit after every task.** Never commit on a task that failed its verification step.

---

## File map (what gets created / modified)

**New experiment scripts (`code/experiments/`):**
- `run_pca_ablation.jl` — E1: four PCA-space baselines + full-harness evaluation → `code/data/pca_ablation_results.csv`
- `run_convex_hull_analysis.jl` — E2: hull membership + out-of-hull plausibility → `code/data/convex_hull_results.csv`
- `run_privacy_analysis.jl` — E3: DCR + NN-MIA AUC → `code/data/privacy_results.csv`
- `sim_common.jl` — SIM shared helpers (ground-truth generators, recovery metrics)
- `run_sim_lowrank.jl` — SIM arm 1 (linear low-rank Gaussian GT sweep) → `code/data/sim_lowrank_results.csv`
- `run_sim_subsample.jl` — SIM arm 2 (real-23 subsampling) → `code/data/sim_subsample_results.csv`
- `run_sim_nonlinear.jl` — SIM arm 3 (nonlinear-manifold GT) → `code/data/sim_nonlinear_results.csv`

**New figure scripts (`code/experiments/`):**
- `fig_sim_phase_diagram.jl` → `code/figs/sim_phase_diagram.pdf` (main text)
- `fig_pca_ablation.jl` → `code/figs/pca_ablation_summary.pdf` (supplement)
- `fig_convex_hull.jl` → `code/figs/convex_hull.pdf` (supplement)
- `fig_privacy.jl` → `code/figs/privacy_dcr_roc.pdf` (supplement)
- `fig_theory_schematic.jl` (OPTIONAL) → `code/figs/theory_schematic.pdf` (main text)

**Manuscript (edit in both `paper/` and `arxiv/`):**
- `main.tex` (+ `paper/supplementary.tex`) preamble — color macros
- `sections/abstract.tex`, `introduction.tex`, `results.tex`, `discussion.tex`, `methods.tex`, `conclusion.tex` — prose
- `sections/floats.tex` — demote `fig:pca-by-visit`, add SIM (+ optional schematic) figures
- `sections/supplementary.tex` — new E1/E2/E3/SIM tables and figures

**Response doc:**
- `peer-review-feedback/response-to-reviewers.md` — fill `[PENDING]`, write general response
- `peer-review-feedback/response-to-reviewers.docx` — pandoc render (new)

---

## Phase 0 — Infrastructure

### Task 1: Tracked-change color macros + `\reviewmode` toggle

**Files:**
- Modify: `paper/main.tex` (after line 31, the `\norm` macro block)
- Modify: `paper/supplementary.tex` (preamble, lines 9–26 region)
- Modify: `arxiv/main.tex` (equivalent macro block)

**Interfaces:**
- Produces: `\rone{}`, `\rtwo{}`, `\rboth{}` text macros and a `\reviewmode` boolean, consumed by every later prose/figure/supplement task.

- [ ] **Step 1: Add the macro block.** Insert into all three preambles (xcolor is already loaded in `paper/main.tex:10`; add `\usepackage{xcolor}` first in any preamble that lacks it — check `paper/supplementary.tex` and `arxiv/main.tex`):

```latex
% --- Tracked-changes for revision (R1 blue, R2 red, both violet) ---
\newif\ifreviewmode
\reviewmodetrue   % set \reviewmodefalse for camera-ready (all edits render black)
\ifreviewmode
  \newcommand{\rone}[1]{\textcolor[rgb]{0.10,0.30,0.75}{#1}}   % Reviewer 1 — blue
  \newcommand{\rtwo}[1]{\textcolor[rgb]{0.75,0.15,0.15}{#1}}   % Reviewer 2 — red
  \newcommand{\rboth}[1]{\textcolor[rgb]{0.50,0.15,0.65}{#1}}  % both — violet
\else
  \newcommand{\rone}[1]{#1}
  \newcommand{\rtwo}[1]{#1}
  \newcommand{\rboth}[1]{#1}
\fi
```

- [ ] **Step 2: Smoke-test the macros compile.** Temporarily add `\rone{R1 test} \rtwo{R2 test} \rboth{both test}` to the top of `paper/sections/abstract.tex`, then build:

Run: `cd paper && make clean && make`
Expected: `main.pdf` and `supplementary.pdf` produced with no errors; the three test phrases render blue/red/violet. Then remove the test phrases.

- [ ] **Step 3: Verify the toggle.** Set `\reviewmodefalse` in `paper/main.tex`, rebuild, confirm test phrases render black; set back to `\reviewmodetrue`.

Run: `cd paper && make`
Expected: clean build both ways.

- [ ] **Step 4: Commit.**

```bash
git add paper/main.tex paper/supplementary.tex arxiv/main.tex
git commit -m "revision: add R1/R2 tracked-change color macros and reviewmode toggle"
```

---

## Phase 1 — Experiments E1 / E2 / E3 (reuse canonical harness; parallelizable)

### Task 2: E1 — PCA-space ablation suite

**Files:**
- Create: `code/experiments/run_pca_ablation.jl`
- Reads: `code/data/full_longitudinal_memory.jld2`, `code/data/synthetic_full_longitudinal.csv` (SA row for the comparison table)
- Writes: `code/data/pca_ablation_results.csv`

**Interfaces:**
- Consumes: canonical memory keys (Global Constraints); `decode_sample`, `feature_summary_comparison`, `sample_novelty`; mechanistic reuse block; cross-visit Frobenius pattern; pooled-median+bootstrap pattern.
- Produces: a methods×metrics CSV with columns `Method, Median_MRE, MRE_CI_lo, MRE_CI_hi, CrossVisit_Frob, Mech_Overlap_TFonly, Mech_KS_D, Mech_KS_p, Median_Novelty, Frac_Novel_gt_0p2`. `Method ∈ {SA, PCA+GMM, PCA+KDE, PCA+kNN, PCA+Copula}`.

All four baselines fit **in the raw PCA coordinate space** `Y = X̂ .* pca_norms'` (`d×23`), draw `N=100` samples as `d×100`, then decode each column with `Patient.decode_sample(y, pca_concat, std_params_concat)` to a 216-vector, reshaped to the 3-visit long format for evaluation.

- [ ] **Step 1: Write the four baseline samplers (fully specified, no packages beyond stdlib+Distributions).** In `run_pca_ablation.jl` define:

```julia
# Y :: d×23 raw PCA coords; returns d×N samples. rng passed for reproducibility.
using Random, LinearAlgebra, Statistics, Distributions

# (a) PCA + diagonal-covariance GMM (1–3 comps, variance floor); robust at n<p.
function sample_gmm(Y, N, rng; ncomp=2, floor_var=1e-3, iters=100)
    d, K = size(Y); pts = collect(eachcol(Y))
    # k-means++ init
    μ = [Y[:, rand(rng, 1:K)] for _ in 1:ncomp]
    σ2 = [fill(var(Y) + floor_var, d) for _ in 1:ncomp]; w = fill(1/ncomp, ncomp)
    for _ in 1:iters
        # E-step
        logr = zeros(K, ncomp)
        for k in 1:K, c in 1:ncomp
            logr[k,c] = log(w[c]) - 0.5*sum(log.(2π.*σ2[c])) - 0.5*sum((pts[k].-μ[c]).^2 ./ σ2[c])
        end
        r = exp.(logr .- maximum(logr, dims=2)); r ./= sum(r, dims=2)
        # M-step
        Nc = vec(sum(r, dims=1)); w = Nc ./ K
        for c in 1:ncomp
            Nc[c] < 1e-6 && continue
            μ[c] = sum(r[k,c].*pts[k] for k in 1:K) ./ Nc[c]
            σ2[c] = sum(r[k,c].*(pts[k].-μ[c]).^2 for k in 1:K) ./ Nc[c] .+ floor_var
        end
    end
    out = zeros(d, N)
    for n in 1:N
        c = rand(rng, Categorical(w)); out[:,n] = μ[c] .+ sqrt.(σ2[c]).*randn(rng, d)
    end
    out
end

# (b) PCA + Gaussian KDE sampler: pick a real point, add isotropic Gaussian noise (Scott bandwidth).
function sample_kde(Y, N, rng)
    d, K = size(Y)
    h = K^(-1/(d+4))                      # Scott's rule factor
    S = cov(Y'; corrected=true)           # d×d empirical cov of the 23 points
    L = cholesky(Symmetric(S) + 1e-8I).L
    out = zeros(d, N)
    for n in 1:N
        j = rand(rng, 1:K); out[:,n] = Y[:,j] .+ h .* (L*randn(rng, d))
    end
    out
end

# (c) PCA + k-NN interpolation (SMOTE-style): real point + λ·(neighbor − point), λ∼U(0,1).
function sample_knn(Y, N, rng; k=3)
    d, K = size(Y); out = zeros(d, N)
    D = [norm(Y[:,i]-Y[:,j]) for i in 1:K, j in 1:K]
    for n in 1:N
        i = rand(rng, 1:K)
        nbrs = sortperm(D[i,:])[2:k+1]     # exclude self
        j = rand(rng, nbrs); λ = rand(rng)
        out[:,n] = Y[:,i] .+ λ.*(Y[:,j].-Y[:,i])
    end
    out
end

# (d) PCA + Gaussian copula: empirical marginal CDF per component + Gaussian correlation.
function sample_copula(Y, N, rng)
    d, K = size(Y)
    U = similar(Y)                        # ranks → (0,1)
    for c in 1:d
        U[c,:] = (invperm(sortperm(Y[c,:])) .- 0.5) ./ K
    end
    Z = quantile.(Normal(), clamp.(U, 1/(2K), 1-1/(2K)))
    R = cor(Z'); R = Symmetric(R) + 1e-6I
    L = cholesky(R).L
    out = zeros(d, N)
    for n in 1:N
        z = L*randn(rng, d); u = cdf.(Normal(), z)
        for c in 1:d                       # inverse empirical marginal (nearest-rank quantile)
            out[c,n] = quantile(Y[c,:], clamp(u[c], 0.0, 1.0))
        end
    end
    out
end
```

- [ ] **Step 2: Add a self-check (the "failing test") for sampler shape/finiteness.** At the top of the run block:

```julia
rng = MersenneTwister(42)
@assert size(sample_gmm(Y, 100, rng)) == (size(Y,1), 100)
@assert all(isfinite, sample_kde(Y, 100, rng))
@assert all(isfinite, sample_knn(Y, 100, rng))
@assert all(isfinite, sample_copula(Y, 100, rng))
```

Run: `cd code && julia experiments/run_pca_ablation.jl`
Expected before samplers are correct: `AssertionError`. After: proceeds past the asserts.

- [ ] **Step 3: Wire the evaluation harness for each method.** For each method's `d×100` PCA samples: decode all columns to a `100×216` matrix, pivot into the long-format synthetic DataFrame (columns = `kept_cols`, plus `SyntheticID`, `Visit`), then compute:
  - **Median MRE + bootstrap CI** — build the pooled per-(visit,feature) `Mean_Rel_Error` via `feature_summary_comparison(df_clean_visit, df_synth_visit, kept_cols)` looped over the 3 visits, pool, then replicate `paper_summary_table.jl:100-112` (seed 42, 10 000 resamples, indices 250/9750).
  - **Cross-visit Frobenius** — `norm(safe_cor(synth_concat) - safe_cor(real_concat)) / norm(safe_cor(real_concat))`, `safe_cor` per `validate_cross_visit_covariance.jl:93`.
  - **Mechanistic overlap + KS (TF-only)** — run `run_patient` (TM=0) with `CALIBRATED_OVERRIDES` on all 100 decoded synth + 23 reals, form `ratios = Predicted./Measured` clamped `[0,5]`, compute `frac_in_range` (band `quantile(real_ratios,[0.05,0.95])`) and `ApproximateTwoSampleKSTest(real_ratios, synth_ratios)` → `.δ`, `pvalue`.
  - **Novelty** — for each sample's unit-normalized PCA direction, `sample_novelty(ξ̂, X̂)`; report median and fraction > 0.2.

- [ ] **Step 4: Run and verify the table is produced with sane invariants.**

Run: `cd code && julia experiments/run_pca_ablation.jl`
Expected: writes `data/pca_ablation_results.csv` with 5 rows; every `Median_MRE > 0`, `0 ≤ Mech_Overlap ≤ 1`, `0 ≤ Median_Novelty ≤ 1`; SA row matches the canonical numbers already in the paper (pooled median MRE ≈ 1.2%) within rounding — this cross-checks the harness wiring. Print the CSV to stdout for the reviewer.

- [ ] **Step 5: Commit.**

```bash
git add code/experiments/run_pca_ablation.jl code/data/pca_ablation_results.csv
git commit -m "revision(E1): PCA-space ablation suite (GMM/KDE/kNN/copula) vs SA"
```

---

### Task 3: E2 — Convex-hull / extrapolation analysis

**Files:**
- Create: `code/experiments/run_convex_hull_analysis.jl`
- Reads: `code/data/full_longitudinal_memory.jld2`, `code/data/synthetic_full_longitudinal.csv`
- Writes: `code/data/convex_hull_results.csv`
- Modify: `code/Project.toml` (add `JuMP`, `HiGHS`)

**Interfaces:**
- Consumes: `Y` (raw PCA reals `d×23`); synthetic PCA coords `Ysyn` (`d×100`) obtained by re-standardizing `synthetic_full_longitudinal.csv` with `std_params_concat` and `transform(pca_concat, Zsyn')`; biological-constraint checks from `validate_biological_constraints.jl`; mechanistic envelope reuse block.
- Produces: CSV with per-synthetic rows `SyntheticID, InHull::Bool, DistToHull, BioPlausible::Bool, MechInEnvelope::Bool` plus a summary header row (fractions).

- [ ] **Step 1: Add the LP/QP solver dependency.**

Run: `cd code && julia --project=. -e 'using Pkg; Pkg.add(["JuMP","HiGHS"]); using JuMP, HiGHS'`
Expected: both resolve and precompile; `Project.toml` gains `JuMP` and `HiGHS`.

- [ ] **Step 2: Write hull membership (LP feasibility) + distance (QP) and a failing self-check.**

```julia
using JuMP, HiGHS
# Is p in conv(cols of Y)? feasibility of Yw=p, sum(w)=1, w≥0. Returns (inhull, dist).
function hull_membership(Y, p; tol=1e-7)
    d, K = size(Y)
    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, w[1:K] >= 0); @variable(m, s[1:d])   # s = Yw − p (slack)
    @constraint(m, sum(w) == 1)
    @constraint(m, [i=1:d], sum(Y[i,k]*w[k] for k in 1:K) - p[i] == s[i])
    @objective(m, Min, sum(s.^2))                     # min squared distance to hull
    optimize!(m)
    dist = sqrt(max(objective_value(m), 0.0))
    (dist <= tol, dist)
end
# Self-check: a real vertex is in-hull (dist≈0); a far point is out.
let (inb, db) = hull_membership(Y, Y[:,1]); @assert inb && db < 1e-5 end
let (ino, _)  = hull_membership(Y, Y[:,1] .+ 100); @assert !ino end
```

Run: `cd code && julia experiments/run_convex_hull_analysis.jl`
Expected before correct: `AssertionError`; after: passes.

- [ ] **Step 3: Classify all 100 synthetic points and assess out-of-hull plausibility.** For each column of `Ysyn`: record `InHull`, `DistToHull`. For the out-of-hull subset, decode to measurement space and check (i) `BioPlausible` via the biological-constraint predicates in `validate_biological_constraints.jl` (non-negativity, <5σ outliers, expected factor↔TGA sign correlations preserved at the cohort level), (ii) `MechInEnvelope` via the mechanistic overlap band (ratio inside `quantile(real_ratios,[0.05,0.95])` for TF-only).

- [ ] **Step 4: Run and verify summary invariants.**

Run: `cd code && julia experiments/run_convex_hull_analysis.jl`
Expected: writes `data/convex_hull_results.csv`; `frac_in_hull + frac_out_hull == 1`; among out-of-hull points a reported `frac_plausible` in [0,1]; prints the three headline numbers (inside %, outside %, outside-but-plausible %).

- [ ] **Step 5: Commit.**

```bash
git add code/experiments/run_convex_hull_analysis.jl code/data/convex_hull_results.csv code/Project.toml code/Manifest.toml
git commit -m "revision(E2): convex-hull membership + out-of-hull plausibility (JuMP/HiGHS)"
```

---

### Task 4: E3 — Membership-inference / privacy analysis

**Files:**
- Create: `code/experiments/run_privacy_analysis.jl`
- Reads: `code/data/full_longitudinal_memory.jld2` (for `df_clean`, `std_params_concat`, `kept_cols`, `n_assays`), `code/data/synthetic_full_longitudinal.csv`
- Writes: `code/data/privacy_results.csv`

**Interfaces:**
- Consumes: standardized real matrix `Zreal` (`23×216`, via `standardize(df_clean_concat, kept_cols)` or `std_params_concat`), standardized synthetic matrix `Zsyn` (`100×216`); the inline concatenation+PCA+generation steps from `run_full_longitudinal.jl:53-135` (for the holdout retrain).
- Produces: CSV with a DCR block (`Kind ∈ {synth_to_real, real_to_real, holdout_to_train}`, per-record min distance) and a scalar `MIA_AUC`.

- [ ] **Step 1: DCR distributions (distance to closest record), Euclidean in standardized space.**

```julia
using LinearAlgebra, Statistics
# each column of A → Euclidean distance to its closest column in B
dcr(A, B) = [minimum(norm(a - b) for b in eachcol(B)) for a in eachcol(A)]
# synth→real, and leave-one-out real→real as the reference (Z_* are n×p; transpose so columns = records):
dcr_sr = dcr(permutedims(Zsyn), permutedims(Zreal))
rr = [minimum(norm(Zreal[i,:] - Zreal[j,:]) for j in 1:23 if j != i) for i in 1:23]
```

- [ ] **Step 2: NN-MIA with a holdout split + failing self-check.** Split the 23 reals into a train subset (15) and holdout (8); regenerate `N=100` SA synthetics **from the 15-subset only** (reuse the `run_full_longitudinal.jl:53-135` inline pipeline on the subset, `seed=42`). The attacker scores each real record by `−(distance to nearest synthetic)`; sweep the threshold to build an ROC of members (15 train) vs non-members (8 holdout); AUC via the rank statistic.

```julia
# score = negative NN-distance to synthetic set; AUC = P(score_member > score_nonmember)
auc(members, nonmembers) = mean([m > n for m in members, n in nonmembers])
@assert 0.0 <= mia_auc <= 1.0
```

Run: `cd code && julia experiments/run_privacy_analysis.jl`
Expected: proceeds past the assert; `mia_auc` computed.

- [ ] **Step 3: Run and verify.**

Run: `cd code && julia experiments/run_privacy_analysis.jl`
Expected: writes `data/privacy_results.csv`; prints `MIA AUC`, and median DCR for synth→real vs real→real. Sanity: `median(dcr_sr) ≳ median(rr)` (synthetics no closer to reals than reals are to each other) and `MIA_AUC ≈ 0.5` are the target findings — record whatever the data show honestly.

- [ ] **Step 4: Commit.**

```bash
git add code/experiments/run_privacy_analysis.jl code/data/privacy_results.csv
git commit -m "revision(E3): DCR distribution + nearest-neighbor membership-inference AUC"
```

---

## Phase 2 — SIM simulation benchmark (the only substantial new code)

### Task 5: SIM shared helpers + arm 1 (linear low-rank Gaussian ground truth)

**Files:**
- Create: `code/experiments/sim_common.jl`
- Create: `code/experiments/run_sim_lowrank.jl`
- Writes: `code/data/sim_lowrank_results.csv`

**Interfaces:**
- Produces (`sim_common.jl`): `make_lowrank_gt(p, r, rng; tail=:gaussian)` → `(sampler, Σpop)` where `sampler(n)` draws `n×p` observations from an `r`-dim latent factor model with a specified cross-visit block covariance and optional heavy-tailed (lognormal/t) marginal transforms; `sa_recover(Xtrain)` → generates `N=100` SA synthetics from `Xtrain` via the canonical inline pipeline (standardize → PCA 0.95 → unit-norm → `find_entropy_inflection` → rescaled ULA `T=2000, α=0.01, seed=42`); `mvn_recover(Xtrain)` → `MvNormal` fit + draw; `cov_frob(Ĉ, Σpop)` and `marginal_mre(Xsyn, Xtest)` recovery metrics against a large fresh GT test draw.
- Consumes: `find_entropy_inflection`, `MultivariateStats.fit(PCA,...)`, `Distributions.MvNormal`.

- [ ] **Step 1: Write `sim_common.jl` with a failing recovery self-check.** On a well-specified GT with large `n`, SA covariance-recovery error must be small:

```julia
rng = MersenneTwister(42)
sampler, Σpop = make_lowrank_gt(30, 5, rng)
Xtest = sampler(5000)
@assert cov_frob(cov(sampler(2000)'), Σpop) < 0.15   # sanity: sample cov ≈ population
```

Run: `cd code && julia experiments/run_sim_lowrank.jl`
Expected before helpers correct: `AssertionError`; after: passes.

- [ ] **Step 2: Sweep the (n, p, r, tail) grid.** For `n ∈ {8,15,23,40,80}`, `p ∈ {30,120,216}`, `r ∈ {3,5,10}`, `tail ∈ {:gaussian,:heavy}`: draw a cohort of `n`, run SA and MVN, measure `cov_frob` vs `Σpop` and `marginal_mre`/KS vs a fresh 5 000-row GT test set. Record one row per grid cell per method.

- [ ] **Step 3: Run and verify.**

Run: `cd code && julia experiments/run_sim_lowrank.jl`
Expected: writes `data/sim_lowrank_results.csv` with columns `n, p, r, tail, Method, Cov_Frob, Marginal_MRE, Marginal_KS`; in the `n<p` cells SA's `Cov_Frob` is finite while MVN's degrades/blows up (rank-deficient) — the intended contrast. Print a pivot of `Cov_Frob` by `(n,p)` for SA vs MVN.

- [ ] **Step 4: Commit.**

```bash
git add code/experiments/sim_common.jl code/experiments/run_sim_lowrank.jl code/data/sim_lowrank_results.csv
git commit -m "revision(SIM arm1): low-rank Gaussian ground-truth recovery sweep (SA vs MVN)"
```

---

### Task 6: SIM arm 2 (real-23 subsampling degradation)

**Files:**
- Create: `code/experiments/run_sim_subsample.jl`
- Reads: `code/data/full_longitudinal_memory.jld2`
- Writes: `code/data/sim_subsample_results.csv`

**Interfaces:**
- Consumes: `sim_common.jl` helpers; the full-23 real concatenated matrix and its `safe_cor` reference; mechanistic overlap reuse.
- Produces: CSV `n, Rep, Cov_Frob_vs_full23, Marginal_MRE, Mech_Overlap`.

- [ ] **Step 1: Subsample + regenerate.** For `n ∈ {20,15,10,8}` and several random subsets (`Rep ∈ 1:10`, seeds derived from index, not `rand`), subsample the 23 complete subjects to `n`, regenerate `N=100` SA synthetics via the canonical inline pipeline, and measure recovery of the **full-23** structure: cross-visit `Cov_Frob` vs `safe_cor` of all 23, pooled `Marginal_MRE` vs all 23, and mechanistic `frac_in_range`.

- [ ] **Step 2: Verify degradation curve is monotone-ish and finite.**

Run: `cd code && julia experiments/run_sim_subsample.jl`
Expected: writes `data/sim_subsample_results.csv`; all metrics finite; mean `Cov_Frob` increases as `n` drops (graceful degradation). Print the mean±sd per `n`.

- [ ] **Step 3: Commit.**

```bash
git add code/experiments/run_sim_subsample.jl code/data/sim_subsample_results.csv
git commit -m "revision(SIM arm2): real-23 subsampling degradation curve"
```

---

### Task 7: SIM arm 3 (nonlinear-manifold ground truth — applicability boundary)

**Files:**
- Create: `code/experiments/run_sim_nonlinear.jl`
- Writes: `code/data/sim_nonlinear_results.csv`

**Interfaces:**
- Consumes: `sim_common.jl` (`sa_recover`); a new `make_nonlinear_gt(p, rng; kind=:scurve)` producing a known nonlinear manifold (S-curve / branching trajectory) embedded in `p`-dim with an off-manifold residual metric `offmanifold_error(Xsyn)`.
- Produces: CSV `kind, curvature, Method, OffManifold_Err, Cov_Frob` comparing SA fidelity on nonlinear vs matched-variance linear GT.

- [ ] **Step 1: Nonlinear GT + off-manifold metric + failing self-check.** The metric must read ~0 for on-manifold points and grow with curvature:

```julia
@assert offmanifold_error(true_manifold_points) < 1e-6
@assert offmanifold_error(sa_recover(train)) >= offmanifold_error(linear_matched)  # nonlinear ≥ linear residual
```

- [ ] **Step 2: Sweep curvature.** For increasing curvature/branching, run SA, measure off-manifold residual vs the linear-GT baseline of matched variance; locate where residual becomes material (the honest linear-PCA boundary).

Run: `cd code && julia experiments/run_sim_nonlinear.jl`
Expected: writes `data/sim_nonlinear_results.csv`; residual rises with curvature; prints the curvature at which off-manifold error exceeds the linear-case residual by a stated factor.

- [ ] **Step 3: Commit.**

```bash
git add code/experiments/run_sim_nonlinear.jl code/data/sim_nonlinear_results.csv
git commit -m "revision(SIM arm3): nonlinear-manifold ground truth, linear-PCA applicability boundary"
```

---

## Phase 3 — Figures

### Task 8: SIM main-text figure (phase diagram / recovery curves)

**Files:**
- Create: `code/experiments/fig_sim_phase_diagram.jl`
- Reads: `sim_lowrank_results.csv`, `sim_subsample_results.csv`
- Writes: `code/figs/sim_phase_diagram.pdf` (+ `.png`), then copy to `paper/sections/figures/sim_phase_diagram.pdf` and `arxiv/sections/figures/`

**Interfaces:** consumes SIM CSVs; follows the figure-style Global Constraint.

- [ ] **Step 1: Build the figure** — a 2-panel main-text figure: (left) a recovery phase-diagram heatmap of SA `Cov_Frob` over the `(n, p)` grid with the MVN failure region annotated; (right) the arm-2 subsampling degradation curve (mean ± sd vs `n`) with SA vs MVN. No title; legends on both panels; gray background; black error bars; the `n<p` boundary annotated.
- [ ] **Step 2: Verify the PDF renders and matches style.**

Run: `cd code && julia experiments/fig_sim_phase_diagram.jl && ls -la figs/sim_phase_diagram.pdf`
Expected: PDF exists, non-empty; visually inspect palette/legend/no-title compliance.

- [ ] **Step 3: Copy into both build trees and commit.**

```bash
cp code/figs/sim_phase_diagram.pdf paper/sections/figures/ && cp code/figs/sim_phase_diagram.pdf arxiv/sections/figures/
git add code/experiments/fig_sim_phase_diagram.jl code/figs/sim_phase_diagram.* paper/sections/figures/sim_phase_diagram.pdf arxiv/sections/figures/sim_phase_diagram.pdf
git commit -m "revision(fig): SIM recovery phase-diagram main-text figure"
```

---

### Task 9: Supplement figures for E1 / E2 / E3

**Files:**
- Create: `code/experiments/fig_pca_ablation.jl` → `code/figs/pca_ablation_summary.pdf`
- Create: `code/experiments/fig_convex_hull.jl` → `code/figs/convex_hull.pdf`
- Create: `code/experiments/fig_privacy.jl` → `code/figs/privacy_dcr_roc.pdf`
- Copy each into `paper/sections/figures/` and `arxiv/sections/figures/`

**Interfaces:** consume `pca_ablation_results.csv`, `convex_hull_results.csv`, `privacy_results.csv`.

- [ ] **Step 1: E1 figure** — grouped bar chart of the five methods across the metric axes (MRE, cross-visit Frobenius, mechanistic overlap, novelty), each axis normalized so "SA passes all, each baseline fails one" is visible; consistent method colors.
- [ ] **Step 2: E2 figure** — PC1–PC2 and PC1–PC3 scatter of the 23 reals with hull outline, synthetic points colored by in/out-hull, out-but-plausible points marked; gray background.
- [ ] **Step 3: E3 figure** — DCR histograms (synth→real vs real→real overlaid) + the MIA ROC curve with AUC annotated.
- [ ] **Step 4: Verify all three PDFs render.**

Run: `cd code && julia experiments/fig_pca_ablation.jl && julia experiments/fig_convex_hull.jl && julia experiments/fig_privacy.jl && ls -la figs/pca_ablation_summary.pdf figs/convex_hull.pdf figs/privacy_dcr_roc.pdf`
Expected: three non-empty PDFs.

- [ ] **Step 5: Copy into build trees and commit.**

```bash
for f in pca_ablation_summary convex_hull privacy_dcr_roc; do cp code/figs/$f.pdf paper/sections/figures/; cp code/figs/$f.pdf arxiv/sections/figures/; done
git add code/experiments/fig_pca_ablation.jl code/experiments/fig_convex_hull.jl code/experiments/fig_privacy.jl code/figs/pca_ablation_summary.* code/figs/convex_hull.* code/figs/privacy_dcr_roc.* paper/sections/figures/pca_ablation_summary.pdf paper/sections/figures/convex_hull.pdf paper/sections/figures/privacy_dcr_roc.pdf arxiv/sections/figures/
git commit -m "revision(fig): E1/E2/E3 supplement figures"
```

---

### Task 10 (OPTIONAL): Theory schematic

**Files:** Create `code/experiments/fig_theory_schematic.jl` → `code/figs/theory_schematic.pdf`.

Only build this if, after Task 11, the main-text figure budget allows it (≤7). A schematic of the energy landscape (stored patterns as wells, Langevin trajectory, β controlling well separation) illustrating the β* intuition for R1.9/R2.minor2. If built, it demotes a second existing figure (candidate: fold `fig:pregnancy-progression` into supplement) to stay ≤7. Mark done or explicitly skipped in the response doc.

- [ ] **Step 1:** Decide build-or-skip based on figure budget; if skipped, note "schematic deferred; theory intuition handled in prose" and move on.
- [ ] **Step 2 (if built):** render, copy to both trees, commit.

---

### Task 11: Demote `fig:pca-by-visit` to supplement; reconcile figure count

**Files:**
- Modify: `paper/sections/floats.tex` (remove figure block lines 136–157), `arxiv/sections/floats.tex` (matching block)
- Modify: `paper/sections/supplementary.tex` (append figure after line 339 block), `arxiv/sections/supplementary.tex`
- Modify: `paper/sections/results.tex:95` (update `\ref{fig:pca-by-visit}`), `arxiv/sections/results.tex`

**Interfaces:** `fig:pca-by-visit` label is referenced only at `results.tex:95`.

- [ ] **Step 1:** Cut the `fig:pca-by-visit` figure block (graphics `sa_vs_mvn_pca_by_visit_v2.pdf`) from `floats.tex` and paste into the Supplementary Figures block of `supplementary.tex` with a new `sfig:` label; add the SIM figure block (graphics `sim_phase_diagram.pdf`, new label `fig:sim-recovery`) into `floats.tex` where `pca-by-visit` was, keeping 7 main-text figures. Do this in both `paper/` and `arxiv/`.
- [ ] **Step 2:** Update the prose reference at `results.tex:95` from `Fig.~\ref{fig:pca-by-visit}` to the new supplement label; add a `Fig.~\ref{fig:sim-recovery}` reference in the new SIM Results subsection (Task 13).
- [ ] **Step 3: Verify the build and figure count.**

Run: `cd paper && make clean && make`
Expected: builds clean; grep confirms exactly 7 `\begin{figure}` in `floats.tex`; no `??` undefined refs in `main.log`.

- [ ] **Step 4: Commit.**

```bash
git add paper/sections/floats.tex paper/sections/supplementary.tex paper/sections/results.tex arxiv/sections/floats.tex arxiv/sections/supplementary.tex arxiv/sections/results.tex
git commit -m "revision(fig): demote pca-by-visit to supplement, add SIM recovery figure (≤7 main)"
```

---

## Phase 4 — Prose edits (result-independent parts can run in parallel with Phase 1–3)

> Each task edits both `paper/sections/` and `arxiv/sections/`, wraps new/changed text in the correct color macro, and ends by rebuilding and running the repo's writing linters. Colors: R1→`\rone` (blue), R2→`\rtwo` (red), both→`\rboth` (violet), per the cross-reference table in the response doc.

### Task 12: Temper wording + mechanistic-independence + proof-of-concept framing (R1.M1, R1.2/R2.3, proof-of-concept)

**Files:** `sections/abstract.tex`, `introduction.tex`, `discussion.tex`, `conclusion.tex` (both trees).

- [ ] **Step 1: Temper the two flagged abstract phrases.** `abstract.tex:15` "mechanistically indistinguishable" → `\rone{not distinguishable by the mechanistic model under the conditions tested}`; `abstract.tex:20` and `conclusion.tex:14` "clinically useful" → `\rone{useful for the modeling use-cases demonstrated (mechanistic calibration, hypothesis generation), within a proof-of-concept}`. Add a proof-of-concept clause to the abstract closing (`\rboth`).
- [ ] **Step 2: Mechanistic-independence clarification** (violet, `\rboth`). In `discussion.tex` limitations (near lines 143–147, the TF+TM ODE-bias sentence): add the explicit residual-dependence statement (ODE calibrated on the same cohort → "consistency with a model fitted to the same biological system," not fully independent ground truth; strength = ODE blind to SA, identical fixed mapping for real and synthetic). Soften any "indistinguishable" in Results Mechanistic Consistency (`results.tex:143-203`) to the tempered phrasing.
- [ ] **Step 3: Proof-of-concept framing** in Introduction (elevate) and Discussion closing (`\rboth`).
- [ ] **Step 4: Verify build + style.**

Run: `cd paper && make` then invoke `/style-check` and `/review-abstract` on the edited files.
Expected: clean build; style-check passes (no em-dashes/subsection headings introduced; line numbers intact); abstract review OK.

- [ ] **Step 5: Commit.** `git commit -m "revision(prose): temper strength/independence wording, proof-of-concept framing (R1.M1,R1.2,R2.3)"`

---

### Task 13: Results edits — E1/E2/E3/SIM paragraphs + subgroup/fibrinolytic/downstream caveats

**Files:** `sections/results.tex` (both trees).

- [ ] **Step 1: Extend the baseline paragraph with E1** (`results.tex:57-65`, after CTGAN/TVAE) — one paragraph citing the new supplement E1 table, stating each PCA baseline fails a different axis while SA passes all four (numbers from `pca_ablation_results.csv`). `\rone`.
- [ ] **Step 2: New "Robustness across data regimes" subsection (SIM)** after the Cross-Visit subsection — summarize arm-1 phase diagram (SA robust across `n<p`, MVN fails), arm-2 degradation, arm-3 nonlinear boundary; reference `\ref{fig:sim-recovery}`; numbers from the three SIM CSVs. `\rboth`.
- [ ] **Step 3: E2 paragraph** (interpolation vs extrapolation; inside/outside-hull fractions; out-of-hull-but-plausible) from `convex_hull_results.csv`. `\rtwo`.
- [ ] **Step 4: E3 paragraph** (DCR, MIA AUC; not differentially private; K=23 no formal guarantee) from `privacy_results.csv`. `\rone`.
- [ ] **Step 5: Subgroup-reliability caveat** in Conditional Generation (`results.tex:100-141`): K_eff<5 for PCOS, magnitudes from 3 empirical norms, hypothesis-generating not substitute. `\rone`.
- [ ] **Step 6: Fibrinolytic/viscoelastic clarification** (`results.tex:51-55`): ROTEM/fibrinolytic features load onto retained PCs → shape generation even without mechanistic check. `\rone`.
- [ ] **Step 7: Downstream-smoothing caveat** (`results.tex:205-230`): 0.94× gain consistent with variance reduction (N=100 vs K=23), not SA adding biological information. `\rone`.
- [ ] **Step 8: Downstream-utility translational framing (R2.5)** (`results.tex:205-230` + a Discussion sentence): strengthen the framing of the existing downstream-utility experiment and connect it to the two use-cases R2 names — subgroup power (conditional amplification of PCOS/PE enables comparisons impossible at n=3/n=5) and hypothesis generation (preserved condition-specific signatures suggest testable factor-level differences) — while staying careful (per Step 7 / R1.11) not to claim synthetic surpasses real. `\rtwo`.
- [ ] **Step 9: Verify.**

Run: `cd paper && make` then `/style-check`, `/review-section results`, `/audit-magic-numbers results` (every new number must trace to a table/figure).
Expected: clean build; audit finds no untraceable magic numbers.

- [ ] **Step 10: Commit.** `git commit -m "revision(prose): Results — E1/E2/E3/SIM paragraphs + subgroup/fibrinolytic/downstream/translational (R2.5) framing"`

---

### Task 14: Discussion edits — PCA limits, tail compression, irregular timing, scalability

**Files:** `sections/discussion.tex` (both trees).

- [ ] **Step 1: PCA limits / nonlinear embeddings** (`\rtwo`, near the PCA discussion lines 37–56): linear-only structure, nonlinear extensions (kernel PCA/autoencoders/diffusion maps) trade interpretability + closed-form β*≈√d, unstable at n=23 → linearity is a deliberate small-cohort choice; tie to SIM arm-3 boundary and E1 (PCA-driven vs sampler-driven failures).
- [ ] **Step 2: Expand tail compression** (`\rone`, lines 131–136): two sources (PCA truncation; finite K=23 empirical magnitudes), severe-pathology tails underrepresented, caution for rare-event/extreme-quantile use; connect to E2 (some out-of-hull phenotypes but tail coverage still limited).
- [ ] **Step 3: Irregular timing / missingness / variable-length** (`\rone`, add to limitations): concatenation assumes aligned complete visits; deployment would need imputation/time-aware embedding/sequence extension; framed as scope boundary.
- [ ] **Step 4: Computational scalability** (`\rone`): per-sample cost linear in K and d over T steps; landscape more stable as K grows (small-K is the hard regime); SIM arm-1 large-n cells as evidence; PCA/memory as practical bound at high p.
- [ ] **Step 5: Elevate cross-domain evidence** (`\rboth`): promote `Varner2026SAProtein` (median seed ≈22 ≈ K=23) and `Varner2026Multiplicity` from aside to explicit generalizability argument (they are already cited at `discussion.tex:47,68,69`).
- [ ] **Step 6: Verify.** `cd paper && make` then `/style-check`, `/review-section discussion`.
- [ ] **Step 7: Commit.** `git commit -m "revision(prose): Discussion — PCA limits, tail compression, timing, scalability, cross-domain evidence"`

---

### Task 15: Methods theory streamlining + intuition (R1.9 / R2.minor2)

**Files:** `sections/methods.tex` (both trees); `sections/supplementary.tex` (receive demoted dense derivations).

- [ ] **Step 1:** Add plain-language intuition for β* (temperature at which stored patients stop blurring into one average and become distinct attractors) and the participation ratio K_eff (effective count of contributing patterns) near `methods.tex:95-117`. `\rboth`.
- [ ] **Step 2:** Move the denser entropy-inflection derivation to a new supplement subsection; keep the main-text Methods focused on the applied pipeline (`alg:sa-pipeline` stays in main text).
- [ ] **Step 3:** If the schematic (Task 10) was built, reference it here.
- [ ] **Step 4: Verify.** `cd paper && make` then `/style-check`, `/review-section methods`. Confirm no broken equation refs (`grep '??' paper/main.log`).
- [ ] **Step 5: Commit.** `git commit -m "revision(prose): Methods theory streamlined + plain-language β*/K_eff intuition"`

---

## Phase 5 — Response document + final build

### Task 16: Fill response-to-reviewers, write general response, render to docx

**Files:** Modify `peer-review-feedback/response-to-reviewers.md`; create `peer-review-feedback/response-to-reviewers.docx`.

- [ ] **Step 1:** Replace every `[PENDING]`/`[pending]` bracket with the real numbers from the experiment CSVs (E1 axes, E2 inside/outside/plausible fractions, E3 AUC + DCR, SIM regime boundaries). Flip `[PLAN]`/`[PLAN → PENDING …]` status tags to `[DONE]` as each corresponding manuscript change lands.
- [ ] **Step 2:** Write the "General response to the editor and reviewers" (marked `[DRAFTED LAST]`) summarizing the three experiments + simulation benchmark + cross-domain evidence.
- [ ] **Step 3:** Remove the internal status-legend block before rendering. Render:

Run: `cd peer-review-feedback && pandoc response-to-reviewers.md -o response-to-reviewers.docx`
Expected: `.docx` produced; open-check that tables render.

- [ ] **Step 4: Commit.** `git add peer-review-feedback/response-to-reviewers.md peer-review-feedback/response-to-reviewers.docx && git commit -m "revision: fill response-to-reviewers with results + general response, render docx"`

---

### Task 17: Final integrated build + camera-ready toggle check

**Files:** none (verification task); possibly touch `sections/supplementary.tex` intros for the new E1/E2/E3/SIM tables.

- [ ] **Step 1: Add supplement tables** for E1 (methods×metrics), E2 (hull summary), E3 (DCR/AUC), SIM (regime grid) into `supplementary.tex` (before the `\FloatBarrier` at line 336), both trees, with `\rone/\rtwo/\rboth` coloring matching the response doc.
- [ ] **Step 2: Full rebuild, review mode ON.**

Run: `cd paper && make clean && make`
Expected: `main.pdf` + `supplementary.pdf` build clean; no undefined refs/citations (`! grep -c '??' main.log`); all new figures/tables appear; colored edits visible.

- [ ] **Step 3: Camera-ready dry run, review mode OFF.** Set `\reviewmodefalse`, rebuild, confirm all edits render black and layout is intact; then restore `\reviewmodetrue` for the tracked-changes submission copy.

Run: `cd paper && make`
Expected: clean black build.

- [ ] **Step 4: arxiv mirror build.**

Run: `cd arxiv && make`
Expected: arxiv `main.pdf` builds clean with the mirrored edits.

- [ ] **Step 5: Final commit.** `git add -A paper arxiv && git commit -m "revision: supplement tables + final integrated build (review-mode + camera-ready verified)"`

---

## Success criteria (from the design spec §8)

- Every R1 (13) and R2 (7) comment has a concrete response + a corresponding manuscript change (verify against the response doc's cross-reference table).
- E1: SA uniquely passes all four fidelity axes; each PCA baseline fails ≥1 (visible in `pca_ablation_results.csv` + supplement figure).
- E2: quantified inside/outside-hull fractions; out-of-hull points demonstrably plausible.
- E3: MIA AUC ≈ 0.5; DCR(synth→real) ≥ DCR(real→real).
- SIM: robustness across linear regimes incl. `n<p`; graceful degradation with `n`; honest nonlinear boundary.
- Main-text figure count = 7 (not increased); `pca-by-visit` demoted; SIM figure added.
- Overclaims tempered; theory streamlined; proof-of-concept framing consistent (verified via `/style-check`, `/review-section`, `/audit-magic-numbers`).
- `response-to-reviewers.docx` generated; colored-diff `main.pdf` + `supplementary.pdf` compile; `\reviewmode` toggle verified in both directions; arxiv mirror builds.

---

## Notes / decisions carried from the design spec

- **Scope guards (do NOT violate):** no second real clinical cohort (SIM substitutes); no alternative/independent mechanistic recalibration experiment (R1.2/R2.3 = prose only); no new nonlinear SA generator (nonlinearity is only a GT probe in SIM arm 3); E3 is empirical only (no differential-privacy mechanism).
- **The `arxiv/` mirror is a separate copy** — every section edit must be duplicated there. This is the most likely source of drift; the final `cd arxiv && make` (Task 17.4) is the guard.
- **No single reusable "build concatenated longitudinal memory" function exists** — the concatenation+PCA+unit-norm+pca_norms logic is currently inlined and duplicated across `run_full_longitudinal.jl`, `paper_sa_vs_mvn.jl`, `paper_pca_by_visit.jl`, `paper_correlation_heatmaps.jl`. E1/E3/SIM re-use it by loading `full_longitudinal_memory.jld2` where possible; the holdout-retrain (E3) and SIM arms that must rebuild from a subset should copy the `run_full_longitudinal.jl:53-135` block. Optional refactor (out of scope for the revision, note for later): extract this into a `Patient.build_longitudinal_memory(df_clean, kept_cols, n_assays; pratio=0.95)` helper.

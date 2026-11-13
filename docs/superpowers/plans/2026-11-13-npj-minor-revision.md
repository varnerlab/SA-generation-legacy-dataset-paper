# npj Minor Revision (Round 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close Reviewer 3's three comments with three small computations, short manuscript
inserts, a round-2 response letter, and a v2.2 upload package, by Monday 2026-11-16.

**Architecture:** One helper file (`code/experiments/revision2_helpers.jl`) holds every new
function and is covered by one test file. Three thin scripts call the helpers and write CSVs.
Prose goes into `paper/sections/` in `\rthree`, mirrored by hand into `arxiv/sections/`
without color. The package copies the `paper-revision-2.1/` workflow.

**Tech Stack:** Julia (project env in `code/`, `include("Include.jl")`), HockinMannModel
(BZ2012 ODE), pdflatex/bibtex, xelatex/latexmk for the letter.

**Spec:** `docs/superpowers/specs/2026-11-13-npj-minor-revision-design.md`

## Global Constraints

- Do not change the canonical cohort or any existing reported number. New numbers come only
  from Tasks 2-4 and live in new tables or CSVs.
- Canonical values: beta* = 2.94 (from `find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))`),
  alpha = 0.01, T = 2000, N = 100, latent seed 42, magnitude seed 424242.
- Do not claim that a Gaussian mixture cannot handle constraints. Do not suggest lowering beta
  as a privacy fix. Do not cite unpublished constrained-SA work. Do not mention any companion
  SA methods manuscript. Do not add new citations to `Alswaidan2026SA`.
- npj forbids Supplementary Methods. Every method goes in `paper/sections/methods.tex`. The
  supplement holds results only.
- No new Discussion paragraphs. The Discussion grows by about 180 words in total. Do not trim
  legacy or round-1 text.
- All round-2 changes use `\rthree` (dark green, rgb 0.05,0.50,0.25). In the marked build,
  `\rone`, `\rtwo`, and `\rboth` render black.
- Prose: short, plain, direct sentences. No em dashes, no "rather than", no sentence that
  starts with an acronym, no self-references to sections. Median sentence about 16 words,
  nothing over about 35. Follow `HOUSE_STYLE.md`.
- Commit locally after each task. Never push, tag, or upload.

## Review Focus

1. Large lambda makes ULA unstable or biased. Expect Langevin to agree with the exact
   rejection reference at every lambda. The test lives in Task 4 (`agreement` check).
2. The copied metric code drifts from the published exact-sampler numbers. Expect an exact
   match with `code/data/exact_sampler_reference.csv`. The test lives in Task 1.
3. The re-implemented biological and strict mechanistic checks drift from the published
   counts. Expect 54, 41, and 24 out of 100 on the published cohort. The tests live in
   Tasks 1 and 3.
4. The leave-one-out bandwidth optimizer stops at a search boundary. Expect an interior
   optimum. The test lives in Task 1.
5. Manuscript drift: a number in the text that is missing from its CSV, a color macro leaking
   into `arxiv/`, or an undefined reference in the standalone supplement. Expect all three
   checks in Task 7 to pass.

---

### Task 1: Helper functions and tests

**Files:**
- Create: `code/experiments/revision2_helpers.jl`
- Create: `code/experiments/test_revision2_helpers.jl`

**Interfaces:**
- Consumes: `exact_sample(X, n; β, weights, rng)` (`code/src/Compute.jl:218`),
  `decode_sample(y, pca_model, std_params)` (`code/src/Patient.jl:364`),
  `bz2012_ratios`, `overlap_ks`, `TF_TGA_COLS` (`code/experiments/mechanistic_eval.jl`),
  `feature_summary_comparison`, `build_concat_matrix`, `sample_novelty`,
  `find_entropy_inflection`.
- Produces:
  - `load_reference() -> NamedTuple`. Fields: `X̂, pca_concat, std_params_concat, pca_norms, complete_subjects, kept_cols, n_assays, df_clean, β_star, C_real, Zreal, median_rr, df_real_complete, real_ratios, shared_magnitudes, bio_mu, bio_sigma`.
    `bio_mu` and `bio_sigma` are the `Dict`s of means and SDs over `outlier_features`,
    computed from `cleaned_full_data.csv` as at `run_convex_hull_analysis.jl:219-221`.
    It mirrors Step 0 of `run_exact_sampler_reference.jl:44-83` exactly, including
    `Random.seed!(424242)` before drawing `shared_magnitudes`.
  - `decode_cohort(ref, dirs::Matrix{Float64}, mags::Vector{Float64}) -> DataFrame`. A copy
    of `run_exact_sampler_reference.jl:98-116` that takes `ref` in place of globals.
  - `evaluate_cohort(ref, dirs, df; with_mech::Bool=true) -> NamedTuple`. A copy of
    `run_exact_sampler_reference.jl:119-150`, with the same field names.
  - `exact_directions(ref, β, rng, n) -> (dirs, k)`.
  - `bio_plausible_patient(rows::AbstractDataFrame, ref) -> Bool`. The same predicate as
    `run_convex_hull_analysis.jl:205-238`: the same `positive_cols`, the same
    `outlier_features`, and mu/sigma from `cleaned_full_data.csv`.
  - `strict_bands(ref) -> Dict`. Keys are the values of `ref.real_ratios.Feature`, and each
    value is the `(5th, 95th)` percentile tuple for that feature, as at
    `run_convex_hull_analysis.jl:265-267`.
  - `strict_mech_pass(ratios::DataFrame, bands, n_visits::Int) -> Bool`. True if
    `nrow(ratios) == n_visits * length(TF_TGA_COLS)` and every ratio lies inside its band.
  - `hormone_affine(ref) -> (idx::Vector{Int}, a::Vector{Float64}, c::Matrix{Float64}, s::Vector{Float64})`.
    It covers the 6 concatenated columns for `:Estradiol` and `:Progesterone` at visits 1-3.
    Then `decode_sample(r .* u, …)[idx] == a .+ r .* (c * u)` for any unit vector `u`, and
    `s` is `std_params_concat.σ[idx]`.
  - `positivity_penalty(ξ, r, A) -> Float64` and `positivity_grad(ξ, r, A) -> Vector{Float64}`,
    where `A` is the `hormone_affine` tuple and
    `C = Σ_j max(0, -(a_j + r c_j·ξ̂)/s_j)^2` with `ξ̂ = ξ/‖ξ‖`.
  - `constrained_ula(X̂, ξ₀, r, A; β, λ, α, T, rng) -> Vector{Float64}`. Returns the final
    state of
    `ξ ← ξ - α[(ξ - X̂ softmax(β X̂ᵀ ξ)) + λ ∇C(ξ; r)] + sqrt(2α/β) ε`.
    At λ = 0 this is Eq. (ula) with uniform weights.
  - `rejection_direction(X̂, r, A; β, λ, rng, max_tries=1_000_000) -> (ξ::Vector{Float64}, tries::Int)`.
    It draws `ξ` with `exact_sample(X̂, 1; β, rng)` and accepts with probability
    `exp(-β λ C(ξ; r))`. With `λ = Inf`, it accepts only when `C == 0`. It throws if it
    reaches `max_tries`.
  - `kde_bandwidths(X̂) -> NamedTuple{(:sigma, :scott, :silverman, :loo)}`.
    `sigma = sqrt(mean(var(X̂; dims=2)))`, `scott = sigma*K^(-1/(d+4))`,
    `silverman = sigma*(4/(d+2))^(1/(d+4))*K^(-1/(d+4))`. `loo` maximizes
    `Σ_i log[(1/(K-1)) Σ_{j≠i} N(m_i; m_j, h²I)]` over `h ∈ [0.01, 3]`, using
    log-sum-exp and a golden-section or grid-plus-refine search.

- [ ] **Step 1: Write `test_revision2_helpers.jl` with these testsets (it includes `Include.jl`, `mechanistic_eval.jl`, and the helpers)**

```julia
const ref = load_reference()   # once, at top level; every testset below uses it

@testset "metric copy reproduces exact_sampler_reference.csv" begin
    dirs, _ = exact_directions(ref, ref.β_star, MersenneTwister(42), 100)
    m = evaluate_cohort(ref, dirs, decode_cohort(ref, dirs, ref.shared_magnitudes); with_mech=true)
    row = CSV.read(joinpath(_PATH_TO_DATA, "exact_sampler_reference.csv"), DataFrame)[1, :]
    for f in (:Median_MRE, :CrossVisit_Frob, :Mech_Overlap_TFonly, :Median_Novelty, :DCR_synth_to_real_median)
        @test isapprox(getproperty(m, f), row[f]; rtol=1e-10)
    end
end
@testset "bio check reproduces 54/100 on the published cohort" begin
    # Build the published cohort's long rows exactly as run_convex_hull_analysis.jl does
    # (read synthetic_full_longitudinal.csv, project to PCA, decode).
    @test count_bio_pass_published(ref) == 54
end
@testset "hormone affine map matches decode_sample" begin
    A = hormone_affine(ref)
    for _ in 1:20
        u = normalize(randn(18)); r = rand(ref.pca_norms)
        @test isapprox(decode_sample(r .* u, ref.pca_concat, ref.std_params_concat)[A[1]], A[2] .+ r .* (A[3] * u); rtol=1e-10, atol=1e-10)
    end
end
@testset "penalty and gradient" begin
    # zero penalty when all six hormone values ≥ 0; central-difference gradient check rtol 1e-5
end
@testset "rejection with λ = Inf gives nonnegative hormones" begin
    # 200 draws: every decoded hormone value ≥ 0
end
@testset "loo bandwidth is interior" begin
    bw = kde_bandwidths(ref.X̂)
    @test 0.01 < bw.loo < 3.0
    @test loo_loglik(ref.X̂, bw.loo) ≥ max(loo_loglik(ref.X̂, 0.9bw.loo), loo_loglik(ref.X̂, 1.1bw.loo))
    @test isapprox(bw.scott, bw.sigma * 23^(-1/22); rtol=1e-12)
end
@testset "λ = 0 constrained ULA matches the exact mixture" begin
    # 2000 chains (T=2000, α=0.01) vs 2000 exact draws at β*: two-sample KS on the first two
    # coordinates of ξ̂, expect p > 0.01
end
```

`count_bio_pass_published(ref)` and `loo_loglik(X̂, h)` are also defined in the helper file.

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `cd code && julia experiments/test_revision2_helpers.jl`
Expected: FAIL with `UndefVarError: load_reference not defined`.

- [ ] **Step 3: Implement every function listed under Produces in `revision2_helpers.jl`**

Derive `hormone_affine` from `decode_sample = destandardize ∘ reconstruct`:
`x = μ + σ ⊙ (mean_pca + P y)`. So `a = μ[idx] + σ[idx] .* mean_pca[idx]` and
`c = σ[idx] .* P[idx, :]`. Take `P` and `mean_pca` from `MultivariateStats.projection` and
`mean` of `pca_concat`.

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `cd code && julia experiments/test_revision2_helpers.jl`
Expected: all testsets pass. The first testset runs the BZ2012 model (a few minutes).

- [ ] **Step 5: Commit**

```bash
git add code/experiments/revision2_helpers.jl code/experiments/test_revision2_helpers.jl
git commit -m "Add round-2 revision helpers and tests"
```

### Task 2: KDE bandwidth check (W1)

**Files:**
- Create: `code/experiments/run_kde_bandwidth_check.jl`
- Output: `code/data/kde_bandwidth_check.csv`

**Interfaces:**
- Consumes: `load_reference`, `kde_bandwidths`, `exact_directions`, `decode_cohort`, `evaluate_cohort`.
- Produces: a CSV with columns `Rule, h, Beta, Beta_over_BetaStar, Median_MRE, CrossVisit_Frob, Mech_Overlap_TFonly, Median_Novelty, DCR_synth_to_real_median`.
  It has 4 rows: `EntropyTransition` (β*), `Scott`, `Silverman`, `LeaveOneOut`.

- [ ] **Step 1: Write the script.** For each rule, set `β = 1/h^2`. Draw directions with `MersenneTwister(42)`, decode with `ref.shared_magnitudes`, and evaluate with `with_mech=true`. Assert that the `EntropyTransition` row equals the `Exact` row of `exact_sampler_reference.csv` to `rtol=1e-10`.
- [ ] **Step 2: Run it.** Run: `cd code && julia experiments/run_kde_bandwidth_check.jl`. Expected: the assertion passes and the CSV is written with 4 rows. Record the h and beta/beta* of each rule.
- [ ] **Step 3: Commit** the script and the CSV.

### Task 3: Real-patient pass rates (W2)

**Files:**
- Create: `code/experiments/run_real_strict_check.jl`
- Output: `code/data/real_strict_check.csv`

**Interfaces:**
- Consumes: `load_reference`, `strict_bands`, `strict_mech_pass`, `bio_plausible_patient`.
- Produces: a CSV with columns `Cohort, N, Bio_Pass, Strict_Mech_Pass, Both_Pass`.
  Rows: `Real (in-sample bands)`, `Published synthetic`, and optionally
  `Real (leave-one-out bands)`.

- [ ] **Step 1: Write the script.** For the 23 complete-case real patients, group `ref.real_ratios` by subject and apply `strict_mech_pass` with `strict_bands(ref)`. Apply `bio_plausible_patient` to each real patient's rows. Recompute the published synthetic row with the same functions. Assert `Published synthetic` = (54, 41, 24). Leave-one-out bands are optional (rebuild the bands without the patient being tested). Add them only if they take under 10 lines.
- [ ] **Step 2: Run it.** Run: `cd code && julia experiments/run_real_strict_check.jl`. Expected: the assertion passes. Record the real counts. Real Bio_Pass should be 23. If it is not 23, stop and report.
- [ ] **Step 3: Commit** the script and the CSV.

### Task 4: Constrained positivity demonstration (W3)

**Files:**
- Create: `code/experiments/run_constrained_positivity_demo.jl`
- Output: `code/data/constrained_positivity_demo.csv` and `code/data/constrained_positivity_agreement.csv`

**Interfaces:**
- Consumes: `load_reference`, `hormone_affine`, `constrained_ula`, `rejection_direction`, `decode_cohort`, `evaluate_cohort`, `bio_plausible_patient`, `strict_bands`, `strict_mech_pass`.
- Produces:
  - `constrained_positivity_demo.csv` with columns `Sampler, Lambda, Alpha, T, Bio_Pass, Hormone_Neg_Profiles, Median_MRE, Hormone_MRE, CrossVisit_Frob, Mech_Overlap_TFonly, Strict_Mech_Pass, Both_Pass, Median_Novelty, Acceptance_Rate`.
    It has 9 rows: Langevin and Reference at λ ∈ {0, 1, 10, 100}, plus `Reference, Inf`.
  - `constrained_positivity_agreement.csv` with columns `Lambda, Entry, KS_D, KS_p`. This is
    Langevin vs Reference on the 6 hormone entries, using 1000 draws per arm.

- [ ] **Step 1: Write the script.**
  - Magnitudes for the 100-profile arms come from `ref.shared_magnitudes`. The 1000-draw agreement arms use `Random.seed!(424243)`.
  - Langevin chains start at `normalize(X̂[:,k0] + 0.01 randn)` with `k0` uniform, which matches the canonical initialization.
  - Step size per arm: `α = min(0.01, 0.2/(λ*L))`, where `L = 2 * maximum(ref.pca_norms)^2 * maximum(sum(abs2, A[3] ./ A[4]; dims=2))`. Use `T = ceil(Int, 20/α)`.
  - Latent RNG `MersenneTwister(42)` for every arm.
  - `Hormone_MRE` is the median relative error over the 6 hormone columns only.
  - `Acceptance_Rate` = accepted draws / total tries for Reference rows. Leave it missing for Langevin rows.
- [ ] **Step 2: Add the agreement assertion.** For each finite λ and each entry, `KS_p > 0.01`. If an entry fails, halve `α` for that arm and rerun once. If it still fails, stop and report. Do not loosen the threshold.
- [ ] **Step 3: Run it.** Run: `cd code && julia experiments/run_constrained_positivity_demo.jl`. Expected: the agreement assertion passes and both CSVs are written. Runtime is about 30-45 minutes, mostly the BZ2012 runs. Record pass rates, acceptance rates, and any loss in fidelity.
- [ ] **Step 4: Commit** the script and both CSVs.

**Checkpoint after Task 4:** Show the author the three CSVs as a short table before any prose
is written. The prose depends on these numbers.

### Task 5: Supplement tables and results text

**Files:**
- Modify: `paper/sections/supplementary.tex`. Add tables near `stab:beta-sweep` (line 30) and `stab:hull-geometry` (line 391), and a results section after the sampling-ladder section (line 623 or later).
- Modify: `arxiv/sections/supplementary.tex` (the same edits without `\rthree`).

- [ ] **Step 1:** Add `stab:kde-bandwidth` (from `kde_bandwidth_check.csv`), `stab:constrained-positivity` (from `constrained_positivity_demo.csv`, with the KS agreement summarized in the caption), and the Task 3 rows added to `stab:hull-geometry`. Add a short results-only section headed "Supplementary Bandwidth and Constraint Results". The section and captions must not contain a method description beyond what a caption needs.
- [ ] **Step 2:** Run `cd paper && make all`, then `grep -n 'undefined' supplementary.log main.log`, then `pdftotext supplementary.pdf - | grep -c '??'`. Expected: no undefined references and a count of 0.
- [ ] **Step 3: Commit.**

### Task 6: Manuscript text and colors

**Files:**
- Modify: `paper/main.tex:53-61` and `paper/supplementary.tex:52-60`. Add `\rthree`; make `\rone`, `\rtwo`, and `\rboth` print `#1` uncolored in both modes. In clean mode, `\rthree` also prints `#1`.
- Modify: `paper/main.tex:102-107` (Data Availability) and `arxiv/main.tex:82-87` (no color).
- Modify: `paper/sections/methods.tex`. Add 1-2 sentences after line 104 ("...$q_k = r_k/\sum_j r_j$ exactly."). Add a new paragraph of 150-200 words at the end of "Baseline and Validation Framework" (before line 330).
- Modify: `paper/sections/results.tex`. Add one sentence to the ablation paragraph (lines 65-89) and one sentence to the plausibility paragraph (lines 176-192), and replace "Positive-variable transformations or constrained decoding could prevent negative values."
- Modify: `paper/sections/discussion.tex`. Insert after line 37 ("...choice of inverse temperature.}"), insert after line 128 ("...generally linear.}"), and extend lines 158-161.
- Mirror every section edit into `arxiv/sections/` without color wrappers.

- [ ] **Step 1:** Write each insert with the exact content listed in the spec under "W4. Manuscript text" and "The argument for each comment", using the numbers from Tasks 2-4. Use the Data Availability sentence exactly: "The de-identified real patient data and the generated synthetic datasets are available in the project repository at \url{https://github.com/varnerlab/SA-generation-legacy-dataset-paper}."
- [ ] **Step 2: Measure the new prose.** Extract the `\rthree` spans and report the word count added to the Discussion (target about 180), the median and maximum sentence length, and the hits for em dashes, "rather than", and acronym-initial sentences. Expected: Discussion additions ≤ 200 words, maximum sentence ≤ 35 words, and 0 hits.
- [ ] **Step 3: Check the trees.** Run `cd paper && make all` and `cd arxiv && make all`. Run `grep -c 'rthree' arxiv/sections/*.tex` (expect 0). Diff `paper/sections` against `arxiv/sections` after stripping color macros, and confirm that only the pre-existing differences remain.
- [ ] **Step 4: Commit.**

### Task 7: Response letter, build, and package

**Files:**
- Create: `peer-review-feedback/response-to-reviewers-round2.md` and `peer-review-feedback/response-to-reviewers-round2.tex`. Reuse the preamble from `response-to-reviewers.tex:1-73`.
- Modify: `peer-review-feedback/Makefile`. Add a target for the round-2 PDF.
- Create: `paper-revision-2.2/`, copied from `paper-revision-2.1/` with `sections/` replaced from `paper/sections/`. Update every `v2.1` string to `v2.2` (in `Makefile`, `README.txt`, `upload-checklist.md`, and `scripts/build_latex_source_zip.py`), and point the response target at the round-2 letter.
- Root outputs: `Varner-Manuscript-Clean-v2.2.pdf`, `Varner-Manuscript-Marked-v2.2.pdf`, `Varner-Supplementary-Information-v2.2.pdf`, `Varner-LaTeX-Source-Clean-v2.2.zip`.

- [ ] **Step 1: Write the letter.** Thank Reviewer 1 in two sentences. Then write sections R3.1-R3.3. Each section quotes the comment, answers it directly in the first sentence, summarizes the evidence with numbers from the CSVs, and names the manuscript location. Add a short note that changes appear in green and that the Data Availability statement was corrected.
- [ ] **Step 2: Build.** Run `cd peer-review-feedback && make` and `cd paper-revision-2.2 && make all && make package latex-source` (use the target names in its Makefile). Expected: every PDF is built, `grep -n undefined *.log` finds nothing, and `pdftotext` of each PDF contains no `??`.
- [ ] **Step 3: Verify the package.**
  - The marked PDF contains the `\rthree` color operator and none of the round-1 color operators (decompress the streams and count the RGB ops).
  - The clean PDF contains no reviewer colors.
  - The standalone supplement has no "Supplementary Methods" heading.
  - Every number in the letter matches its CSV.
  - `git diff --check` is clean.
  - List what was not verified.
- [ ] **Step 4: Commit locally.** Hand off to the author for the read-through. Do not push or upload.

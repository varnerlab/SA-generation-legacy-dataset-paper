# SA Sampling-Correctness Investigation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Attribute the E3 membership-inference leak (MIA ≈ 0.97, DCR synth→real < real→real) to a specific cause — chain initialization, single-endpoint extraction, or ULA discretization bias — or show it is intrinsic to the target distribution π(ξ) ∝ exp(−βE) at β\*, via a mixing diagnostic and a four-rung MCMC ladder (V0→V3).

**Architecture:** Two new standalone experiment scripts reuse the existing samplers and evaluation harness. The mixing diagnostic derives a principled burn-in `B` and thinning `τ`; the ladder generates four N=100 cohorts that differ only in the direction sampler and runs each through one identical fidelity + privacy harness. The published cohort and all prior results stay frozen — these are diagnostic probes, not a method replacement.

**Tech Stack:** Julia (existing `Compute.sample`/`mala_sample`, `find_entropy_inflection`, `decode_sample`, `mechanistic_eval.jl`, the E3 DCR + E3b MIA machinery, `run_sampling_diagnostics.jl` ESS/ACF primitives), Plots/GR.

**Design doc:** `docs/superpowers/specs/2026-07-17-sa-sampling-investigation-design.md`.

## Global Constraints

- **Held constant across all rungs (only the direction sampler varies):** β = β\* recomputed via `find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1,1000.0))` on the 23-memory X̂ (≈ 2.94); α = 0.01; the 18-dim unit-normalized memory `X̂` from `code/data/full_longitudinal_memory.jld2`; the decode = normalize each sampled state to a unit direction, multiply by a magnitude **drawn from the 23 empirical `pca_norms`**, then `decode_sample(ξ_rescaled, pca_concat, std_params_concat)`; **N = 100** synthetic patients per rung; reproducible seeds.
- **Run all Julia from `code/`**, scripts begin `include(joinpath(@__DIR__, "..", "Include.jl"))`; mechanistic/MIA steps additionally `using HockinMannModel, HypothesisTests` and `include("mechanistic_eval.jl")`.
- **Frozen references — do NOT modify:** `run_full_longitudinal.jl`, `src/Patient.jl:sa_generate_from_matrix`, `mechanistic_eval.jl`, `run_privacy_analysis.jl`, `run_privacy_split_sensitivity.jl`, `run_sampling_diagnostics.jl`. Task S1 ADDS a helper to `src/Patient.jl` without changing existing behavior; new ladder code uses it (the E3/E3b copies stay as-is — their duplication is a noted Minor, not re-opened here).
- **Uniform-sphere init** (V1+): `ξ₀ = randn(d); ξ₀ ./= norm(ξ₀)`.
- **β fixed** — never sweep β in this investigation (that is a separate axis, out of scope).
- **Privacy metric coverage (approved default):** DCR on all four rungs; repeated-split MIA (E3b machinery) on **V0 and the best-DCR rung only** (best = largest median DCR synth→real).
- **V2/V3 pooling (approved default):** 20 chains × 5 thinned post-burn-in draws = 100.
- **Prose edits (Task S6) apply to BOTH `paper/sections/` and `arxiv/sections/`** and the response doc; color macros `\rone`/`\rboth` (privacy is R1.5 → `\rone`).
- **Commit after every task; never commit on a failed verification.**

## File map

**New scripts (`code/experiments/`):**
- `run_mcmc_mixing_diagnostic.jl` → `code/data/mcmc_mixing_diagnostic.csv`, `code/figs/mcmc_mixing_diagnostic.pdf`
- `run_mcmc_ladder.jl` → `code/data/mcmc_ladder_results.csv`
- `run_mcmc_ladder_mia.jl` → `code/data/mcmc_ladder_mia.csv`
- `run_beta_privacy_sweep.jl` → `code/data/beta_privacy_sweep.csv` (Task Sβ)
- `fig_mcmc_ladder.jl` → `code/figs/mcmc_ladder.pdf`

**Modified:**
- `code/src/Patient.jl` — add `build_concat_matrix` helper (Task S1).
- `paper/sections/{results,discussion,supplementary}.tex` + `arxiv/` mirror + `peer-review-feedback/response-to-reviewers.md` (Task S6).

---

## Task S1: Extract `build_concat_matrix` helper

**Files:**
- Modify: `code/src/Patient.jl` (add near `sa_generate_from_matrix`)
- Create: `code/experiments/test_build_concat.jl`

**Interfaces:**
- Produces: `build_concat_matrix(df, id_col::Symbol, id_values::Vector, kept_cols::Vector{Symbol}, n_assays::Int; visit_col::Symbol=:Visit) -> Matrix{Float64}` returning a `length(id_values) × (n_assays·3)` matrix, row `i` = the `[V1|V2|V3]` concatenation of `id_values[i]`'s three visits in `kept_cols` order, `offset=(v-1)·n_assays` (the layout used across the repo).
- Consumed by: Tasks S2, S3, S4.

- [ ] **Step 1: Write the failing test.** In `test_build_concat.jl`: load `full_longitudinal_memory.jld2` (`df_clean`, `complete_subjects`, `kept_cols`, `n_assays`), build the real 23×216 matrix with the new helper, and assert it equals the inline construction used by `run_full_longitudinal.jl:57-65` (rebuild that inline block locally and compare):

```julia
X_ref = <inline [V1|V2|V3] build over complete_subjects, as in run_full_longitudinal.jl:57-65>
X_new = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)
@assert size(X_new) == (length(complete_subjects), n_assays*3)
@assert isapprox(X_new, X_ref; atol=1e-12) "build_concat_matrix diverges from the canonical inline layout"
println("build_concat_matrix reproduces canonical concat ✓")
```

Run: `cd code && julia experiments/test_build_concat.jl`
Expected before helper: `UndefVarError: build_concat_matrix`.

- [ ] **Step 2: Implement `build_concat_matrix`** in `Patient.jl` (general over `id_col`/`visit_col`, so it works for both real `SubjectID` and synthetic `SyntheticID` cohorts).

- [ ] **Step 3: Run the test → passes** (`build_concat_matrix reproduces canonical concat ✓`).

- [ ] **Step 4: Commit.**
```bash
git add code/src/Patient.jl code/experiments/test_build_concat.jl
git commit -m "revision(sampling): extract build_concat_matrix helper (reproduces canonical layout)"
```

---

## Task S2: Mixing diagnostic

**Files:**
- Create: `code/experiments/run_mcmc_mixing_diagnostic.jl`
- Reads: `code/data/full_longitudinal_memory.jld2`
- Writes: `code/data/mcmc_mixing_diagnostic.csv`, `code/figs/mcmc_mixing_diagnostic.pdf`

**Interfaces:**
- Consumes: `X̂` (via memory JLD2), `find_entropy_inflection`, `Compute.sample`, `Compute.mala_sample`, and `integrated_acf` (from `run_sampling_diagnostics.jl:126` — copy the small function or include the definition; do not run that whole script).
- Produces: `B` (burn-in) and `τ` (thinning) written to the CSV summary for Task S3 to read.

- [ ] **Step 1: Run M sphere-init chains and compute R̂ + ESS.** `M=20` chains, `ξ₀` uniform on the sphere, ULA `sample(X̂, ξ₀, T; β=β_star, α=0.01)` for `T=5000`, seeded per chain. Record the energy trace `E(ξ_t)` per chain via `hopfield_energy(ξ_t, X̂, β_star)`.

```julia
# split-R̂ on the energy scalar over a growing window → find B where R̂ < 1.05
function split_rhat(traces::Matrix{Float64})   # traces: M × L (post-warmup window)
    M, L = size(traces); Lh = L ÷ 2
    chains = vcat(traces[:,1:Lh], traces[:,Lh+1:2Lh])   # 2M half-chains of length Lh
    m = size(chains,1)
    means = vec(mean(chains, dims=2)); vars = vec(var(chains, dims=2))
    B = Lh * var(means); W = mean(vars)
    V = ((Lh-1)/Lh)*W + B/Lh
    sqrt(V / W)
end
```

Compute R̂ over windows `[t0:T]` for increasing `t0`; set `B` = smallest `t0` with `R̂ < 1.05` (fallback `B = 1000` if never, and log it). `τ = ceil(mean integrated_acf of post-B energy traces))`. ESS = `(T-B)/τ`.

- [ ] **Step 2: Anchor-distance trace + init-dependence.** For one representative chain, record `min_k ‖ξ_t − X̂[:,k]‖` and `‖ξ_t − X̂[:,k0]‖` vs t. Then run `M` **anchor-init** chains and `M` **sphere-init** chains at the same `T`, and compare their endpoint distributions with a KS test on (a) endpoint energies and (b) endpoint nearest-memory distances. Report the KS D + p (via `ApproximateTwoSampleKSTest`) — small D / large p ⇒ init-independent (mixed); large D ⇒ init still matters at this T.

- [ ] **Step 3: MALA acceptance** (context): run `M` MALA chains, report mean acceptance (reuse `mala_sample.accept_rate`).

- [ ] **Step 4: Self-check + run.** Assert `B < T`, `τ ≥ 1`, `0 ≤ MALA acceptance ≤ 1`, KS D ∈ [0,1].

Run: `cd code && julia experiments/run_mcmc_mixing_diagnostic.jl`
Expected: writes the CSV (R̂-vs-window, B, τ, ESS, anchor-distance summary, init-dependence KS D/p, MALA acceptance) and the figure (energy traces + burn-in marker, R̂ curve, anchor-distance curve). Prints `B`, `τ`, and the init-dependence verdict (mixed / not-mixed).

- [ ] **Step 5: Commit.**
```bash
git add code/experiments/run_mcmc_mixing_diagnostic.jl code/data/mcmc_mixing_diagnostic.csv code/figs/mcmc_mixing_diagnostic.*
git commit -m "revision(sampling): MCMC mixing diagnostic (R-hat, ESS, anchor-distance, init-dependence)"
```

---

## Task S3: The V0–V3 ladder — generation + fidelity + DCR

**Files:**
- Create: `code/experiments/run_mcmc_ladder.jl`
- Reads: `full_longitudinal_memory.jld2`, `mcmc_mixing_diagnostic.csv` (for the starting τ estimate only), `cleaned_full_data.csv`/`df_clean` (real reference), `synthetic_full_longitudinal.csv` (V0 published-cohort cross-check)
- Writes: `code/data/mcmc_ladder_results.csv`, `code/data/mcmc_ladder_burnin.csv`

**Interfaces:**
- Consumes: `build_concat_matrix` (S1); `Compute.sample`, `Compute.mala_sample`; `hopfield_energy`; `decode_sample`; `sample_novelty`; `feature_summary_comparison`; `safe_cor` (local, per `validate_cross_visit_covariance.jl:93`); `mechanistic_eval.jl` (`bz2012_ratios`,`overlap_ks`); the E3 DCR computation.
- Produces: `mcmc_ladder_results.csv` with one row per rung `{Rung, Median_MRE, CrossVisit_Frob, Mech_Overlap_TFonly, Mech_KS_D, Median_Novelty, DCR_synth_to_real_median, DCR_real_to_real_median, best_rung::Bool}`; plus `mcmc_ladder_burnin.csv` with the derived `B`, settled-segment `τ`, per-chain settling times, and Geweke z-scores.

- [ ] **Step 1: Derive and VALIDATE the burn-in B and thinning τ (do NOT hardcode 1000).** The S2 mechanical B=0 is a known artifact (rising-R̂, multimodal basin-locking). Instead derive B from settling times: run the 20 sphere-init ULA chains long (`T=6000`), record each chain's energy `hopfield_energy(ξ_t, X̂, β_star)` trace, and for each chain find its **settling step** = first `t` after which the trailing-window (width 500) mean energy stays within a small ε of the chain's final-basin mean energy (chain has entered and stabilized in its basin). Set **B = ceil(1.5 × the 90th-percentile settling step)** (starting expectation ≈1000; let the data move it). Then **recompute τ** as `ceil(mean integrated_acf of the POST-B energy segments)` (the settled-segment autocorrelation — may differ from S2's transient-inflated 107). **Validate** with a Geweke-style stationarity check on each chain's post-B energy: compare mean of the first 10% vs the last 50%, report the z-scores; assert the fraction of chains with `|z|<2` is high (e.g. ≥0.8) — if not, increase B and re-derive, logging it. Write `B`, τ, the settling-time distribution, and the Geweke z-scores to `mcmc_ladder_burnin.csv`. Print them.

- [ ] **Step 2: Four direction samplers (shared decode).** Implement a `draw_directions(rung; ...) -> d×100 matrix of unit directions`:
  - **V0:** `k0=rand(1:K); ξ₀=X̂[:,k0]+0.01·randn; normalize; sample(...,T=2000); normalize endpoint` — 100 chains×1 (reproduces the published scheme).
  - **V1:** identical but `ξ₀` uniform-sphere.
  - **V2:** 20 sphere-init chains, ULA, run `B+5τ` steps, take states at `B+τ,…,B+5τ`, normalize each → 100 directions (B, τ from Step 1).
  - **V3:** same pooling as V2 but `mala_sample`.
  Then decode uniformly for every rung: for each direction, draw a magnitude from `pca_norms` (seeded), `decode_sample(dir .* mag, pca_concat, std_params_concat)` → 216-vector; assemble the 100×3-visit long DataFrame (`kept_cols`+`SyntheticID`+`Visit`).

- [ ] **Step 3: Self-check (the failing test).** Assert each rung yields a finite `100×216` decoded cohort, and that **V0's pooled median MRE ≈ the published ~1.2%** (V0 must reproduce the published fidelity — it is the same scheme; if not, the ladder harness is mis-wired).

```julia
@assert isapprox(v0_median_mre, 0.012; atol=0.004) "V0 does not reproduce published fidelity — harness mis-wired"
```

Run: `cd code && julia experiments/run_mcmc_ladder.jl` → fails here until the harness is correct.

- [ ] **Step 4: Fidelity + DCR per rung.** For each rung compute: pooled median MRE (per-visit `feature_summary_comparison` → pool); cross-visit Frobenius (`safe_cor`); mechanistic overlap + KS (TF-only, via `bz2012_ratios`/`overlap_ks`); median novelty (`sample_novelty` on the unit directions); DCR synth→real + real→real medians (E3 computation, canonical standardized space). Mark `best_rung=true` for the rung with the largest `DCR_synth_to_real_median`.

- [ ] **Step 5: Run + verify invariants.**

Run: `cd code && julia experiments/run_mcmc_ladder.jl`
Expected: writes `mcmc_ladder_burnin.csv` + `mcmc_ladder_results.csv` (4 rows); every `Median_MRE>0`, `0≤overlap≤1`, `0≤novelty≤1`; V0 MRE ≈ 1.2%; exactly one `best_rung`; the Geweke check passes (≥80% of chains stationary post-B). Print the burn-in derivation (B, τ, settling times, Geweke) and the rung × metric table. Report honestly whether/where DCR rises above the real→real baseline (privacy improving) and at what fidelity cost.

- [ ] **Step 6: Commit.**
```bash
git add code/experiments/run_mcmc_ladder.jl code/data/mcmc_ladder_results.csv code/data/mcmc_ladder_burnin.csv
git commit -m "revision(sampling): V0-V3 MCMC ladder — derived/validated burn-in, fidelity + DCR per rung"
```

---

## Task S4: Repeated-split MIA on V0 and the best-DCR rung

**Files:**
- Create: `code/experiments/run_mcmc_ladder_mia.jl`
- Reads: `mcmc_ladder_results.csv` (to identify `best_rung`), `full_longitudinal_memory.jld2`
- Writes: `code/data/mcmc_ladder_mia.csv`

**Interfaces:**
- Consumes: the E3b repeated-split MIA machinery (`run_privacy_split_sensitivity.jl` — its `split_mia_auc` logic) generalized to accept a **direction-sampler for the holdout retrain** (so it can generate the 15-subset synthetics with V0's anchor-init scheme vs the best rung's scheme), `build_concat_matrix` (S1).
- Produces: `mcmc_ladder_mia.csv` with the AUC distribution (mean, SD, p05/p50/p95, R, frac==1) for `V0` and `best_rung`.

- [ ] **Step 1: Generalize the split-MIA to a sampler argument.** Read `run_privacy_split_sensitivity.jl`; factor its per-split retrain so the direction sampler for the 15-subset holdout retrain is a parameter (V0 anchor-init endpoint vs the best rung's init/pooling/sampler). Keep gen_seed=42, the adaptive stop (SD/√R<0.01, floor 40, cap 200), and the same scoring.

- [ ] **Step 2: Self-check.** The V0 path must reproduce E3b's mean AUC ≈ 0.975 (same scheme):

```julia
@assert isapprox(v0_mean_auc, 0.975; atol=0.03) "V0 MIA does not reproduce E3b — sampler wiring wrong"
```

Run: `cd code && julia experiments/run_mcmc_ladder_mia.jl` → fails until the V0 path matches E3b.

- [ ] **Step 3: Run for V0 and best rung.**

Run: `cd code && julia experiments/run_mcmc_ladder_mia.jl`
Expected: writes `mcmc_ladder_mia.csv`; V0 mean AUC ≈ 0.975 (cross-check); best-rung AUC reported honestly (≈0.5 ⇒ leak fixed by that sampler change; still high ⇒ intrinsic to π at β\*). Prints both distributions.

- [ ] **Step 4: Commit.**
```bash
git add code/experiments/run_mcmc_ladder_mia.jl code/data/mcmc_ladder_mia.csv
git commit -m "revision(sampling): repeated-split MIA on V0 and best-DCR rung"
```

---

## Task Sβ: β–privacy diagnostic curve (DCR vs β)

**Decision (user, 2026-07-18):** the ladder verdict is "leak intrinsic to π at β\*" — fixing the sampler (V1→V3) barely moved DCR. The design's R1.5 branch (design §6, line 73) says the honest answer is then "privacy requires lowering β (a real trade-off)." That trade-off must be **measured, not asserted**: `run_beta_sweep.jl` already gives novelty/MRE vs β but never touches DCR. This task adds the privacy side. **β\* = 2.94 stays the canonical operating point — this is a diagnostic sweep *around* it, NOT a change to the paper's operating β** (moving the operating point was explicitly rejected as cascading). Stays within the "diagnostic probes, not method replacement" scope guard, exactly like the ladder.

**Files:**
- Create: `code/experiments/run_beta_privacy_sweep.jl`
- Reads: `code/data/full_longitudinal_memory.jld2` (`X̂`, `pca_concat`, `std_params_concat`, `pca_norms`, `df_clean`, `kept_cols`, `n_assays`, `complete_subjects`), `code/data/cleaned_full_data.csv` (real reference for DCR + fidelity)
- Writes: `code/data/beta_privacy_sweep.csv`

**Interfaces:**
- Consumes: the **published V0 decode scheme** (anchor-init single-endpoint ULA, `T=2000`, magnitude drawn from `pca_norms`) — the SAME scheme as the ladder V0 rung / `run_full_longitudinal.jl`, but with β swept; `find_entropy_inflection` (for the β\* reference value); the **E3 DCR computation** (Euclidean in the canonical standardized space, `dcr(A,B)` from `run_privacy_analysis.jl`); `sample_novelty`; per-visit `feature_summary_comparison` → pooled median MRE. Reuse `build_concat_matrix` (S1) for the real reference matrix.
- Produces: `beta_privacy_sweep.csv` with one row per β: `beta, beta_frac, DCR_synth_to_real_median, DCR_real_to_real_median, Median_Novelty, Median_MRE, Mech_Overlap_TFonly` (mechanistic overlap optional if cheap).

- [ ] **Step 1: Mirror the existing sweep grid.** Read the β grid semantics from `run_beta_sweep.jl` (β_frac × β\*) so this privacy curve aligns on the same x-axis as the paper's existing novelty/fidelity curve. Grid: `β_frac ∈ {0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0}`, `β = β_frac · β_star`. `β_frac = 1.0` is the canonical operating point.

- [ ] **Step 2: Generate + measure per β.** At each β: `Random.seed!(42)`, generate `N=100` synthetics with the V0 scheme at that β, decode, compute DCR synth→real median (+ the constant real→real baseline), median novelty, pooled median MRE. Optional: mechanistic overlap (TF-only) via `mechanistic_eval.jl` if runtime allows.

- [ ] **Step 3: Self-check (cross-wiring anchor).** At `β_frac = 1.0` the row MUST reproduce the ladder V0 rung within tolerance (the same scheme at the same β):
```julia
@assert isapprox(row_betastar.DCR_synth_to_real_median, 13.97; atol=0.5) "β* row does not reproduce ladder V0 DCR — harness mis-wired"
@assert isapprox(row_betastar.Median_MRE, 0.013; atol=0.004) "β* row does not reproduce V0 fidelity"
```
Run: `cd code && julia experiments/run_beta_privacy_sweep.jl` → fails until wired correctly.

- [ ] **Step 4: Run + verify the trade-off is visible.**
Run: `cd code && julia experiments/run_beta_privacy_sweep.jl`
Expected: writes `beta_privacy_sweep.csv`; as β **decreases**, DCR synth→real **rises** (toward/past the real→real baseline = better privacy) while MRE and novelty rise (worse fidelity / more dispersion) — the frontier. β\* row matches V0. Print the β × {DCR, novelty, MRE} table and state where DCR first crosses the real→real baseline.

- [ ] **Step 5: (OPTIONAL) MIA vs β at the extremes only.** Repeated-split MIA is expensive; if included, run it only at the two bracketing β (lowest and β\*) to bracket the AUC, not the full grid. Default: skip — DCR + novelty + MRE already draw the frontier. Log the decision either way.

- [ ] **Step 6: Commit.**
```bash
git add code/experiments/run_beta_privacy_sweep.jl code/data/beta_privacy_sweep.csv
git commit -m "revision(sampling): β–privacy diagnostic curve (DCR vs β; frontier around frozen β*)"
```

> **S5/S6 now incorporate Sβ:** S5's ladder figure gains a β-privacy panel (DCR + fidelity vs β, β\* marked); S6's R1.5 response quantifies "lowering β trades fidelity for privacy" with the crossing point from `beta_privacy_sweep.csv` instead of asserting it.

---

## Task S5: Ladder + diagnostic supplement figures

**Files:**
- Create: `code/experiments/fig_mcmc_ladder.jl` → `code/figs/mcmc_ladder.pdf`
- Copy `mcmc_ladder.pdf` and `mcmc_mixing_diagnostic.pdf` into `paper/sections/figures/` and `arxiv/sections/figures/`

**Interfaces:** consumes `mcmc_ladder_results.csv`, `mcmc_ladder_mia.csv`, `mcmc_mixing_diagnostic.csv`; follows the figure-style convention (consistent colors, no title, legends, gray scatter bg, black error bars, tests annotated).

- [ ] **Step 1: Build the ladder figure** — a rung (V0→V3) × metric panel: privacy (DCR synth→real vs the real→real baseline line; MIA points for V0/best) on one axis and fidelity (MRE, mechanistic overlap) on another, so the privacy-improvement-and-its-fidelity-cost is visible across rungs.
- [ ] **Step 2: Verify the PDF renders.** Run: `cd code && julia experiments/fig_mcmc_ladder.jl && ls -la figs/mcmc_ladder.pdf` → non-empty PDF.
- [ ] **Step 3: Copy to both trees and commit.**
```bash
cp code/figs/mcmc_ladder.pdf paper/sections/figures/ && cp code/figs/mcmc_ladder.pdf arxiv/sections/figures/
cp code/figs/mcmc_mixing_diagnostic.pdf paper/sections/figures/ && cp code/figs/mcmc_mixing_diagnostic.pdf arxiv/sections/figures/
git add code/experiments/fig_mcmc_ladder.jl code/figs/mcmc_ladder.* paper/sections/figures/mcmc_ladder.pdf paper/sections/figures/mcmc_mixing_diagnostic.pdf arxiv/sections/figures/mcmc_ladder.pdf arxiv/sections/figures/mcmc_mixing_diagnostic.pdf
git commit -m "revision(sampling): ladder + mixing-diagnostic supplement figures"
```

---

## Task S6: Write-up — R1.5 response + supplement section + Methods/Discussion note

**Files:**
- Modify: `peer-review-feedback/response-to-reviewers.md` (R1.5), `paper/sections/supplementary.tex` (+ arxiv), `paper/sections/discussion.tex` (+ arxiv), `paper/sections/methods.tex` (+ arxiv, brief sampling note)

**Interfaces:** consumes the results of S2–S4. Coordinates with the main-plan Task 13 (Results privacy paragraph) and Task 16 (response doc) — this task owns the *sampling-investigation* narrative; those tasks reference it.

- [ ] **Step 1: Determine the verdict from the data** — read `mcmc_mixing_diagnostic.csv` (init-dependence) + `mcmc_ladder_results.csv`/`mcmc_ladder_mia.csv` and state which rung (if any) closes the leak and at what fidelity cost, per the design's decision logic.
- [ ] **Step 2: R1.5 response** (`response-to-reviewers.md`): update the privacy response with (a) the leak (E3/E3b: MIA ≈ 0.97, DCR), (b) the **de-identified framing** (method property, not a PII breach; no realistic linkage for 216-dim assays), (c) the **attribution** from the ladder (sampler-driven vs intrinsic to β\*), and (d) the honest scope: SA at the published operating point is not a privacy mechanism; if a rung fixes it, note it as a sampling recommendation. Numbers from the CSVs.
- [ ] **Step 3: Supplement section** (`supplementary.tex`, both trees): new subsection "Sampling diagnostics and the memorization–privacy relationship" — the mixing diagnostic (R̂, ESS, init-dependence) + the ladder table + the two figures (S5). `\rone`.
- [ ] **Step 4: Methods/Discussion note** (`methods.tex` + `discussion.tex`, both trees): a short, honest note that the published cohort uses anchor-initialized single-endpoint ULA, that a proper-MCMC sampling of π was examined (supplement), and the resulting privacy characterization. `\rone`.
- [ ] **Step 5: Verify.** `cd paper && make` (clean build); run `/style-check` and `/audit-magic-numbers` on the edited sections; confirm every new number traces to a CSV/figure.
- [ ] **Step 6: Commit.**
```bash
git add peer-review-feedback/response-to-reviewers.md paper/sections/supplementary.tex paper/sections/discussion.tex paper/sections/methods.tex arxiv/sections/supplementary.tex arxiv/sections/discussion.tex arxiv/sections/methods.tex
git commit -m "revision(sampling): R1.5 response + supplement section + methods/discussion note"
```

---

## Success criteria

- Mixing diagnostic reports R̂, τ_int/ESS, MALA acceptance, and a clear anchor-vs-sphere init-dependence verdict; B and τ derived and consumed by the ladder.
- Ladder yields a complete rung × {MRE, cross-visit Frobenius, mechanistic overlap/KS, novelty, DCR} table; V0 reproduces published fidelity (~1.2% MRE) and E3b MIA (~0.975) as harness cross-checks.
- MIA on V0 + best-DCR rung reported; the leak is attributed to init / pooling / discretization, or declared intrinsic to π at β\*, with the fidelity cost of any privacy gain quantified.
- Published cohort untouched; results framed as diagnostic probes; R1.5 response + supplement section + methods/discussion note drafted, build clean, magic-number audit clean.

## Self-review notes

- Spec coverage: diagnostic (S2), ladder V0–V3 (S3), MIA endpoints (S4), figures (S5), write-up + R1.5 + de-identified framing (S6), DRY helper (S1) — all spec sections mapped.
- Frozen references honored: S1 only adds to `Patient.jl`; ladder code is new; E3/E3b copies untouched.
- Cross-check anchors prevent silent mis-wiring: V0 must reproduce published MRE (S3) and E3b MIA (S4).

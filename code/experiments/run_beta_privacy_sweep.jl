# ──────────────────────────────────────────────────────────────────────────────
# run_beta_privacy_sweep.jl
#
# Task Sβ — β–privacy diagnostic curve. The V0→V3 MCMC ladder (Task S3) and the
# S4 MIA both showed the DCR / membership-inference leak is intrinsic to
# π(ξ) ∝ exp(−βE) at the FROZEN operating point β*≈2.94 — not a sampler
# artifact. `run_beta_sweep.jl` already sweeps β and measures novelty/MRE but
# NEVER touches DCR. This script adds the privacy side of that trade-off: DCR
# vs β, on the SAME β_frac grid, so the paper can show the fidelity/privacy
# frontier with numbers instead of asserting it.
#
# Task Sβ-extend — the original marginal-MRE metric barely moves across β
# (1.24%→1.45%), which a reviewer could read as "lower β buys privacy almost
# for free." To make the claim airtight, this extends the grid downward
# (adding β_frac=0.05, 0.1) and adds two JOINT-fidelity columns — CrossVisit_Frob
# and Mech_Overlap_TFonly — mirroring the exact patterns `run_mcmc_ladder.jl`
# already uses per rung (safe_cor Frobenius; BZ2012 TF-only mechanistic cloud
# overlap via mechanistic_eval.jl), so the sweep can show privacy rising WHILE
# joint structure degrades, not just a flat marginal-fidelity cost.
#
# Scope guard: β*=2.94 stays the paper's canonical operating point. This is a
# DIAGNOSTIC sweep AROUND it — the published cohort, `run_full_longitudinal.jl`,
# and every frozen reference are read-only here. Only β varies.
#
# Generator at every grid point: `sa_generate_from_matrix` (Task S1 helper,
# src/Patient.jl) — the exact published V0 decode scheme (anchor init, single
# T=2000 ULA chain per sample, magnitude drawn from pca_norms, seed=42) — called
# once per β. DCR uses the E3 `dcr(A,B)` computation from run_privacy_analysis.jl
# verbatim, in the canonical standardized space (`std_params_concat`).
# CrossVisit_Frob mirrors `run_mcmc_ladder.jl`'s per-rung `safe_cor` Frobenius
# block (itself from `validate_cross_visit_covariance.jl`); Mech_Overlap_TFonly
# mirrors `run_mcmc_ladder.jl`'s per-rung `bz2012_ratios`/`overlap_ks` block
# (`mechanistic_eval.jl`, TF-only, TM=0.0). Both frozen scripts are read-only —
# only their patterns are copied here.
#
# Reads:  data/full_longitudinal_memory.jld2 (df_clean, complete_subjects,
#         kept_cols, n_assays, pca_concat, std_params_concat, pca_norms, X̂)
#         data/cleaned_full_data.csv       (not re-read; df_clean from the JLD2
#                                            is the same cleaned frame)
#         data/mcmc_ladder_results.csv     (V0 row — self-check + soft cross-check)
# Writes: data/beta_privacy_sweep.csv
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))

using Random, LinearAlgebra, Statistics

# ══════════════════════════════════════════════════════════════════════════════
# Step 0: Load canonical pipeline memory + real reference (held fixed across β)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 0: Loading canonical pipeline memory …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
JLD2.@load mem_path X̂ pca_concat std_params_concat pca_norms complete_subjects kept_cols n_assays df_clean

d, K = size(X̂)
n_visits = 3
@info "  X̂: $(size(X̂))  (d=$d, K=$K real patients),  n_assays=$n_assays"

# Real reference matrix — via the Task-S1 helper (not a re-copied inline loop).
X_real = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)
@info "  X_real (build_concat_matrix): $(size(X_real))"

Zreal = (X_real .- std_params_concat.μ') ./ std_params_concat.σ'

# E3 DCR computation (run_privacy_analysis.jl:115, verbatim): records as COLUMNS.
dcr(A, B) = [minimum(norm(a - b) for b in eachcol(B)) for a in eachcol(A)]

# DCR real→real: leave-one-out median (CONSTANT across β — computed, not hardcoded).
rr = [minimum(norm(Zreal[i, :] .- Zreal[j, :]) for j in 1:K if j != i) for i in 1:K]
median_rr = median(rr)
@info "  DCR real→real median (reference, constant across β) = $(round(median_rr, digits=4))"

df_real_complete = df_clean[in.(df_clean.SubjectID, Ref(complete_subjects)), :]

# ── Cross-visit Frobenius reference (held fixed across β): safe_cor pattern
# copied verbatim from run_mcmc_ladder.jl:84-89 / validate_cross_visit_covariance.jl:93-97 ──
function safe_cor(X)
    C = cor(X)
    C[isnan.(C)] .= 0.0
    return C
end
C_real = safe_cor(X_real)

# ── Mechanistic (TF-only) reference (held fixed across β): bz2012_ratios on the
# K real complete-case patients, once — copied pattern from run_mcmc_ladder.jl:99-101 ──
@info "  Running BZ2012 (TF-only, calibrated) on the $K real complete-case patients …"
real_ratios = bz2012_ratios(df_real_complete, TF_TGA_COLS; TM=0.0)
@info "  real_ratios: $(nrow(real_ratios)) rows"

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Mirror the existing sweep grid — β_frac × β*.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 1: Deriving β* via find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1,1000.0)) …"
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  β* = $(round(β_star, digits=6))"

# Task Sβ-extend: add β_frac=0.05, 0.1 at the low end (grid kept sorted ascending).
const β_FRACS = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
const N_SAMPLES = 100
const T_STEPS = 2000

@info "  Grid: β_frac ∈ $β_FRACS  →  β ∈ $(round.(β_FRACS .* β_star, digits=3))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Generate + measure per β
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 2: Generating + measuring per β …"

results = DataFrame(
    beta = Float64[],
    beta_frac = Float64[],
    DCR_synth_to_real_median = Float64[],
    DCR_real_to_real_median = Float64[],
    Median_Novelty = Float64[],
    Median_MRE = Float64[],
    CrossVisit_Frob = Float64[],
    Mech_Overlap_TFonly = Float64[],
)

for frac in β_FRACS
    β = frac * β_star
    @info "  β_frac=$frac  β=$(round(β, digits=4)) …"

    Random.seed!(42)
    synth_df, mem = sa_generate_from_matrix(X_real, kept_cols, n_assays, N_SAMPLES, T_STEPS;
                                             β=β, seed=42)

    # ── DCR synth→real: E3 computation, canonical standardized space ──
    X_synth = build_concat_matrix(synth_df, :SyntheticID, collect(1:N_SAMPLES), kept_cols, n_assays)
    Zsynth = (X_synth .- std_params_concat.μ') ./ std_params_concat.σ'
    dcr_sr = dcr(permutedims(Zsynth), permutedims(Zreal))
    med_dcr_sr = median(dcr_sr)

    # ── Median novelty: recover each sample's pre-decode unit direction by
    # inverting decode_sample (standardize with mem.std_params, then
    # MultivariateStats.transform via mem.pca_model — the exact adjoint of the
    # reconstruct() used inside sa_generate_from_matrix's decode step, so this
    # recovers ξ_dir up to floating-point round-off without touching the
    # frozen helper's internals), then `sample_novelty(dir, mem.X̂)`. ──
    Z_synth_std = (X_synth .- mem.std_params.μ') ./ mem.std_params.σ'
    ξ_pca = MultivariateStats.transform(mem.pca_model, permutedims(Z_synth_std))  # d_pca × N_SAMPLES
    novelties = [sample_novelty(ξ_pca[:, i] ./ (norm(ξ_pca[:, i]) + 1e-12), mem.X̂) for i in 1:N_SAMPLES]
    med_nov = median(novelties)

    # ── Median MRE: pooled per-(visit,feature) Mean_Rel_Error, looped over
    # the 3 visits (same pattern as the ladder / E1). ──
    all_mre = Float64[]
    for v in 1:n_visits
        real_v = df_real_complete[df_real_complete.Visit .== v, :]
        synth_v = synth_df[synth_df.Visit .== v, :]
        comp = feature_summary_comparison(real_v, synth_v, kept_cols)
        append!(all_mre, comp.Mean_Rel_Error)
    end
    med_mre = median(all_mre)

    # ── Cross-visit Frobenius: safe_cor pattern, copied from run_mcmc_ladder.jl:498-500 ──
    C_synth = safe_cor(X_synth)
    frob = norm(C_synth - C_real) / norm(C_real)

    # ── Mechanistic (TF-only) overlap: bz2012_ratios/overlap_ks pattern, copied
    # from run_mcmc_ladder.jl:502-504 ──
    @info "    BZ2012 mechanistic eval (TF-only) on $N_SAMPLES synthetic patients …"
    synth_ratios = bz2012_ratios(synth_df, TF_TGA_COLS; TM=0.0)
    mech_overlap, mech_ksD, mech_ksp = overlap_ks(real_ratios.Ratio, synth_ratios.Ratio)

    push!(results, (β, frac, med_dcr_sr, median_rr, med_nov, med_mre, frob, mech_overlap))
    @info "    DCR(s→r)=$(round(med_dcr_sr,digits=3))  [DCR(r→r)=$(round(median_rr,digits=3))]  " *
          "Novelty=$(round(med_nov,digits=3))  MRE=$(round(100*med_mre,digits=2))%  " *
          "Frob=$(round(frob,digits=4))  MechOverlap=$(round(mech_overlap,digits=3))  KS D=$(round(mech_ksD,digits=3))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Self-check — β_frac=1.0 MUST reproduce the ladder V0 rung (same
# scheme at the same β; this is the wiring guard for the whole harness).
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 3: Self-check against ladder V0 …"

ladder_path = joinpath(_PATH_TO_DATA, "mcmc_ladder_results.csv")
df_ladder = CSV.read(ladder_path, DataFrame)
v0_row = df_ladder[findfirst(==("V0"), df_ladder.Rung), :]
@info "  Ladder V0 target: DCR(s→r)=$(round(v0_row.DCR_synth_to_real_median,digits=4))  MRE=$(round(v0_row.Median_MRE,digits=5))"

row_betastar = results[findfirst(isapprox(1.0), results.beta_frac), :]
@info "  Sweep β_frac=1.0 row:  DCR(s→r)=$(round(row_betastar.DCR_synth_to_real_median,digits=4))  " *
      "MRE=$(round(row_betastar.Median_MRE,digits=5))"

@assert isapprox(row_betastar.DCR_synth_to_real_median, 13.975; atol=0.6) "β* row does not reproduce ladder V0 DCR — harness mis-wired"
@assert isapprox(row_betastar.Median_MRE, 0.013; atol=0.004) "β* row does not reproduce V0 fidelity"
@info "  Self-check PASSED: β_frac=1.0 row reproduces ladder V0 within tolerance."

# Range-sanity for the new joint-fidelity columns, across ALL rows (Task Sβ-extend).
@assert all(0 .<= results.Mech_Overlap_TFonly .<= 1) "Mech_Overlap_TFonly must be in [0,1] for every row"
@assert all(isfinite, results.CrossVisit_Frob) && all(results.CrossVisit_Frob .>= 0) "CrossVisit_Frob must be finite and non-negative for every row"
@info "  Range-sanity PASSED: CrossVisit_Frob finite & ≥0, Mech_Overlap_TFonly ∈ [0,1] for all $(nrow(results)) rows."

# SOFT cross-check (report, do NOT hard-assert): the sweep's β_frac=1.0 uses the
# canonical sa_generate_from_matrix cohort, whereas ladder V0 uses a separate
# magnitude-RNG stream, so expect CLOSE (same neighborhood) not identical.
@info "\n  Soft cross-check — β_frac=1.0 joint-fidelity columns vs ladder V0 (different magnitude-RNG stream, expect close not identical):"
@info "    Sweep β_frac=1.0:  CrossVisit_Frob=$(round(row_betastar.CrossVisit_Frob,digits=4))  Mech_Overlap_TFonly=$(round(row_betastar.Mech_Overlap_TFonly,digits=4))"
@info "    Ladder V0:         CrossVisit_Frob=$(round(v0_row.CrossVisit_Frob,digits=4))  Mech_Overlap_TFonly=$(round(v0_row.Mech_Overlap_TFonly,digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Write CSV + report the frontier honestly.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Writing data/beta_privacy_sweep.csv …"
CSV.write(joinpath(_PATH_TO_DATA, "beta_privacy_sweep.csv"), results)

println("\n" * "="^100)
println("β–privacy diagnostic sweep (Task Sβ-extend) — β × {DCR, novelty, MRE, CrossVisit_Frob, Mech_Overlap_TFonly}")
println("="^100)
pretty_table(results)

println("\n" * "="^100)
println("Honest read-out: does DCR(synth→real) cross the real→real baseline ($(round(median_rr,digits=3))) " *
        "inside the grid?")
println("="^100)
crossing_idx = findfirst(r -> r.DCR_synth_to_real_median >= median_rr, eachrow(results))
if crossing_idx === nothing
    println("  NO crossing within β_frac ∈ [$(minimum(β_FRACS)), $(maximum(β_FRACS))]×β*: " *
            "DCR(s→r) stays BELOW the real→real baseline at every grid point tested.")
    println("  Highest DCR(s→r) observed: $(round(maximum(results.DCR_synth_to_real_median),digits=3)) " *
            "at β_frac=$(results.beta_frac[argmax(results.DCR_synth_to_real_median)]) " *
            "(gap to baseline: $(round(median_rr - maximum(results.DCR_synth_to_real_median),digits=3)))")
else
    crow = results[crossing_idx, :]
    println("  Crosses at β_frac=$(crow.beta_frac) (β=$(round(crow.beta,digits=3))): " *
            "DCR(s→r)=$(round(crow.DCR_synth_to_real_median,digits=3)) ≥ baseline=$(round(median_rr,digits=3))")
end
for row in eachrow(results)
    moved = row.DCR_synth_to_real_median >= median_rr ? "ABOVE real→real — privacy improves" :
            "still below real→real"
    println("  β_frac=$(row.beta_frac)  β=$(round(row.beta,digits=3))  DCR(s→r)=$(round(row.DCR_synth_to_real_median,digits=3))  " *
            "Novelty=$(round(row.Median_Novelty,digits=3))  MRE=$(round(100*row.Median_MRE,digits=2))%  " *
            "Frob=$(round(row.CrossVisit_Frob,digits=4))  MechOverlap=$(round(row.Mech_Overlap_TFonly,digits=3))  — $moved")
end
println("="^100)

println("\n" * "="^100)
println("Honest read-out: does JOINT fidelity (CrossVisit_Frob, Mech_Overlap_TFonly) degrade toward low β? " *
        "(this is the whole point of Task Sβ-extend — privacy up, joint fidelity down, no sweet spot)")
println("="^100)
# Sorted by beta ascending (β_FRACS is sorted ascending and results are pushed in that order).
frob_lowest_beta, frob_highest_beta = results.CrossVisit_Frob[1], results.CrossVisit_Frob[end]
overlap_lowest_beta, overlap_highest_beta = results.Mech_Overlap_TFonly[1], results.Mech_Overlap_TFonly[end]
println("  CrossVisit_Frob:      β_frac=$(results.beta_frac[1]) → $(round(frob_lowest_beta,digits=4))   " *
        "β_frac=$(results.beta_frac[end]) → $(round(frob_highest_beta,digits=4))   " *
        "($(frob_lowest_beta > frob_highest_beta ? "RISES toward low β (worse cross-visit fidelity) ✓ cliff" : "does NOT rise toward low β"))")
println("  Mech_Overlap_TFonly:  β_frac=$(results.beta_frac[1]) → $(round(overlap_lowest_beta,digits=4))   " *
        "β_frac=$(results.beta_frac[end]) → $(round(overlap_highest_beta,digits=4))   " *
        "($(overlap_lowest_beta < overlap_highest_beta ? "DROPS toward low β (worse mechanistic fidelity) ✓ cliff" : "does NOT drop toward low β"))")
println("="^100)

# ── Invariant checks ──
@assert nrow(results) == length(β_FRACS) "expected $(length(β_FRACS)) rows, got $(nrow(results))"
@assert all(results.DCR_synth_to_real_median .> 0) "every DCR_synth_to_real_median must be > 0"
@assert all(results.DCR_real_to_real_median .== median_rr) "DCR_real_to_real_median must be constant across β"
@assert all(0 .<= results.Median_Novelty .<= 1) "Median_Novelty must be in [0,1]"
@assert all(results.Median_MRE .> 0) "every Median_MRE must be > 0"
@assert issorted(results.beta) "β grid should be increasing (β_FRACS is sorted ascending)"
@assert all(0 .<= results.Mech_Overlap_TFonly .<= 1) "Mech_Overlap_TFonly must be in [0,1]"
@assert all(isfinite, results.CrossVisit_Frob) && all(results.CrossVisit_Frob .>= 0) "CrossVisit_Frob must be finite and non-negative"
@info "\nAll self-checks and invariants passed. Wrote data/beta_privacy_sweep.csv ($(nrow(results)) rows)."

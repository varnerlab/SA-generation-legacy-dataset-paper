# ──────────────────────────────────────────────────────────────────────────────
# run_sim_subsample.jl
#
# SIM benchmark, arm 2: real-23 subsampling degradation curve.
#
# The full 23-patient cohort is already small. This arm asks: if only n<23 of
# those patients had been available, how gracefully does SA's recovery of the
# FULL-23 cohort's structure degrade? Every subset (n ∈ {20,15,10,8}, 10 reps
# each) is scored against the SAME full-23 reference — safe_cor(real_full23)
# for the cross-visit metric, the full-23 clinical cohort for pooled MRE, and
# the full-23 TF-only mechanistic ratio band for the mechanistic overlap — so
# the numbers measure recovery of the full cohort's structure from a smaller
# draw, not each subset recovering itself. This complements arm 1 (synthetic
# low-rank ground truth) with the real data, and evidences SA is usable at
# (and somewhat below) the operating cohort size. Report the truth regardless
# of whether degradation turns out to be monotone.
#
# Reads:  data/full_longitudinal_memory.jld2  (df_clean, complete_subjects,
#                                               kept_cols, n_assays)
# Writes: data/sim_subsample_results.csv      (n, Rep, Cov_Frob_vs_full23,
#                                               Marginal_MRE, Mech_Overlap)
#
# Real/synthetic concat-matrix reconstruction uses `build_concat_matrix`
# (code/src/Patient.jl, read-only). `safe_cor` is copied locally from
# validate_cross_visit_covariance.jl:93 (also read-only). The mechanistic
# TF-only axis uses `mechanistic_eval.jl`'s `bz2012_ratios`/`overlap_ks`,
# mirroring run_mcmc_ladder.jl's mechanistic block (read-only) exactly.
#
# Runtime: 40 cells × (100-sample SA regeneration with per-subset β* search +
# 100×3-row BZ2012 ODE mechanistic re-evaluation). The mechanistic ODE solves
# dominate (40×100 solves plus one 23-patient reference solve) — expect
# ~15-30 minutes total. Do NOT shrink N, T, or the grid to speed this up.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "sim_common.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 0: Load the real cohort + build the FULL-23 reference (held constant
# across every subset/rep — the point is recovery of THIS structure).
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 0: Loading canonical pipeline memory …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
JLD2.@load mem_path complete_subjects kept_cols n_assays df_clean

K_full = length(complete_subjects)
n_visits = 3
@info "  complete_subjects: $K_full real patients, n_assays=$n_assays"

function safe_cor(X)
    # Copied per validate_cross_visit_covariance.jl:93 (local, as directed by the brief).
    C = cor(X)
    C[isnan.(C)] .= 0.0
    return C
end

X_real_full23 = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)
C_real_full23 = safe_cor(X_real_full23)
@info "  X_real_full23 (build_concat_matrix): $(size(X_real_full23))"

df_real_full23 = df_clean[in.(df_clean.SubjectID, Ref(complete_subjects)), :]
@info "  Running BZ2012 (TF-only, calibrated) on the $K_full real full-23 patients (reference, once) …"
real_ratios_full23 = bz2012_ratios(df_real_full23, TF_TGA_COLS; TM=0.0)
@info "  real_ratios_full23: $(nrow(real_ratios_full23)) rows"

# ══════════════════════════════════════════════════════════════════════════════
# Deterministic per-cell seeding (no `rand()`), anchored to the canonical seed
# 42, mirroring sim_common.jl's `cell_seed` pattern but sized for this task's
# 2-D (n-index, Rep) grid.
# ══════════════════════════════════════════════════════════════════════════════
"""
    subsample_seed(i_n, rep; base=42) -> Int

Deterministic seed from 1-based (n-index, Rep) grid-cell indices. Used both
to seed the WITHOUT-replacement subsampling RNG and as `sa_generate_from_matrix`'s
`seed=` (its own independent `Random.seed!` call), so the two seeded streams
never alias each other despite sharing an integer value.
"""
subsample_seed(i_n::Int, rep::Int; base::Int=42) = base + 1_000 * rep + 10 * i_n

n_grid = [20, 15, 10, 8]

"""
    run_cell(n, i_n, rep) -> NamedTuple

Subsample `n` of the 23 complete subjects (seeded, without replacement),
regenerate `N=100` SA synthetics from that subset, and score all three
metrics against the FULL-23 reference built in Step 0.
"""
function run_cell(n::Int, i_n::Int, rep::Int)
    seed = subsample_seed(i_n, rep)

    rng_sub = MersenneTwister(seed)
    subset_ids = Random.shuffle(rng_sub, complete_subjects)[1:n]

    X_subset = build_concat_matrix(df_clean, :SubjectID, subset_ids, kept_cols, n_assays)

    synth_df, _mem = sa_generate_from_matrix(X_subset, kept_cols, n_assays, 100, 2000;
                                              β=nothing, α=0.01, pratio=0.95, seed=seed)

    # ── Metric 1: cross-visit Frobenius vs safe_cor(real_full23) ──
    synth_concat = build_concat_matrix(synth_df, :SyntheticID, collect(1:100), kept_cols, n_assays)
    C_synth = safe_cor(synth_concat)
    cf = norm(C_synth - C_real_full23) / norm(C_real_full23)

    # ── Metric 2: pooled median per-(visit,feature) MRE vs the full-23 real cohort ──
    all_mre = Float64[]
    for v in 1:n_visits
        real_v = df_real_full23[df_real_full23.Visit .== v, :]
        synth_v = synth_df[synth_df.Visit .== v, :]
        comp = feature_summary_comparison(real_v, synth_v, kept_cols)
        append!(all_mre, comp.Mean_Rel_Error)
    end
    mre = median(all_mre)

    # ── Metric 3: TF-only mechanistic cloud overlap vs the full-23 real ratio band ──
    synth_ratios = bz2012_ratios(synth_df, TF_TGA_COLS; TM=0.0)
    frac, _ksD, _ksp = overlap_ks(real_ratios_full23.Ratio, synth_ratios.Ratio)

    return (n=n, Rep=rep, Cov_Frob_vs_full23=cf, Marginal_MRE=mre, Mech_Overlap=frac)
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: FAST self-check (TDD) — one subset, n=20 Rep=1, wiring sanity
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: self-check — one subset (n=20, Rep=1) …"
t0 = time()
row1 = run_cell(20, 1, 1)
@info "  Step 1 cell finished in $(round(time()-t0,digits=1))s: " *
      "cf=$(row1.Cov_Frob_vs_full23)  mre=$(row1.Marginal_MRE)  mech=$(row1.Mech_Overlap)"

cf, mre, mech = row1.Cov_Frob_vs_full23, row1.Marginal_MRE, row1.Mech_Overlap
@assert isfinite(cf) && cf >= 0
@assert isfinite(mre) && 0 < mre < 0.10        # clinical MRE, ~1-5%
@assert 0 <= mech <= 1
@info "  Step 1 self-check PASSED ✓"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: full sweep — n ∈ {20,15,10,8} × Rep 1:10 (row1 above is cell
# (n=20,Rep=1); reused here rather than recomputed).
# ══════════════════════════════════════════════════════════════════════════════
total_cells = length(n_grid) * 10
@info "Step 2: sweeping $total_cells cells …"

rows = NamedTuple[row1]
cell_i = 1
t_start = time()

for (i_n, n) in enumerate(n_grid), rep in 1:10
    (n == 20 && rep == 1) && continue   # already computed in Step 1
    global cell_i += 1
    elapsed = round(time() - t_start, digits=1)
    @info "[$cell_i/$total_cells] n=$n Rep=$rep seed=$(subsample_seed(i_n, rep))  (elapsed $(elapsed)s)"
    push!(rows, run_cell(n, i_n, rep))
end

@info "Step 2 complete: $(length(rows)) rows in $(round(time()-t_start,digits=1))s"

results = DataFrame(rows)
sort!(results, [:n, :Rep], rev=[true, false])

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: write CSV + verify graceful degradation
# ══════════════════════════════════════════════════════════════════════════════
out_path = joinpath(_PATH_TO_DATA, "sim_subsample_results.csv")
CSV.write(out_path, results)
@info "Step 3: wrote $out_path ($(nrow(results)) rows)"

@assert nrow(results) == total_cells "expected $total_cells rows, got $(nrow(results))"
@assert all(isfinite, results.Cov_Frob_vs_full23) && all(results.Cov_Frob_vs_full23 .>= 0)
@assert all(isfinite, results.Marginal_MRE) && all(results.Marginal_MRE .> 0)
@assert all(isfinite, results.Mech_Overlap) && all(0 .<= results.Mech_Overlap .<= 1)
@info "  Invariant checks passed: all metrics finite and in-range."

agg = combine(groupby(results, :n),
              :Cov_Frob_vs_full23 => mean => :Mean_Cov_Frob,
              :Cov_Frob_vs_full23 => std  => :SD_Cov_Frob,
              :Marginal_MRE => mean => :Mean_Marginal_MRE,
              :Marginal_MRE => std  => :SD_Marginal_MRE,
              :Mech_Overlap => mean => :Mean_Mech_Overlap,
              :Mech_Overlap => std  => :SD_Mech_Overlap)
sort!(agg, :n, rev=true)   # n=20 → n=8

println()
println("═"^90)
println("Mean ± SD per n (n=20 → n=8, descending; full-23 reference held fixed)")
println("═"^90)
pretty_table(agg)

cf_by_n = agg.Mean_Cov_Frob  # already sorted n descending: [n=20, n=15, n=10, n=8]
monotone = all(diff(cf_by_n) .> 0)
@info "  Mean Cov_Frob_vs_full23 by n (20→8): $(round.(cf_by_n, digits=4))"
if monotone
    @info "  Degradation is MONOTONE: mean Cov_Frob_vs_full23 strictly increases as n drops. ✓"
else
    @warn "  Degradation is NON-MONOTONE somewhere in the n=20→8 curve — reporting as-is (no massaging)."
end

println()
println("Reference (published, full n=23 cohort): pooled median MRE ≈1.2%, " *
        "TF-only mechanistic overlap ≈0.90, cross-visit Frobenius ≈0.64.")
println("n=20 (nearest to full cohort): Cov_Frob=$(round(agg.Mean_Cov_Frob[1],digits=4)), " *
        "MRE=$(round(100*agg.Mean_Marginal_MRE[1],digits=2))%, " *
        "Mech_Overlap=$(round(agg.Mean_Mech_Overlap[1],digits=3))")

@info "run_sim_subsample.jl complete."

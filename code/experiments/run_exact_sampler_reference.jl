# ──────────────────────────────────────────────────────────────────────────────
# run_exact_sampler_reference.jl
#
# Task 2.1/2.2 of the exact-sampler update. The V0-V3 ladder (run_mcmc_ladder.jl)
# asked whether the published cohort's behavior is an artifact of anchor
# initialization, missing burn-in, or ULA discretization, by varying the
# direction sampler while holding everything else fixed. All four rungs are
# still Markov chains. This script adds the analytic reference: the weighted
# Hopfield target is exactly the Gaussian mixture
#
#     pi(xi) = sum_k q_k N(m_k, beta^-1 I),   q_k = r_k / sum_j r_j   (unit memories)
#
# so it can be sampled ancestrally with no chain at all (`exact_sample`,
# verified by test_exact_sampler.jl). If the exact reference reproduces the
# ladder's behavior, the fidelity/novelty/privacy signal belongs to the TARGET
# rather than to any sampler's finite-time dynamics.
#
# Protocol match with the ladder: same X̂, same beta*, same N=100, same decode
# (direction -> empirical PCA magnitude -> reconstruct -> de-standardize), and
# the SAME seeded magnitude draw (MAGNITUDE_SEED = 424_242) so that the
# direction sampler is the only element that differs from V0-V3.
#
# Reads:  data/full_longitudinal_memory.jld2
#         data/mcmc_ladder_results.csv        (V0-V3 rows, for the comparison print)
# Writes: data/exact_sampler_reference.csv    (one protocol-matched row, label "Exact")
#         data/exact_sampler_replicates.csv   (100 independent cohorts, stability)
#
# The published ULA cohort is NOT regenerated or replaced.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))

using Random, LinearAlgebra, Statistics

const N_SYNTH       = 100
const MAGNITUDE_SEED = 424_242      # identical to run_mcmc_ladder.jl:107
const EXACT_SEED     = 42           # canonical latent-draw seed
const N_REPLICATES   = 100
const NOVELTY_THRESH = 0.2

unitize(v::Vector{Float64}) = v ./ (norm(v) + 1e-12)

# ══════════════════════════════════════════════════════════════════════════════
# Step 0: canonical memory + real reference (held fixed, mirrors the ladder)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 0: Loading canonical pipeline memory …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
JLD2.@load mem_path X̂ pca_concat std_params_concat pca_norms complete_subjects kept_cols n_assays df_clean

d, K = size(X̂)
n_visits = 3
@info "  X̂: $(size(X̂))  (d=$d, K=$K real patients),  n_assays=$n_assays"

# Unit-norm memories are what make q_k = r_k / sum(r); assert rather than assume.
mem_norms = [norm(X̂[:, k]) for k in 1:K]
@assert all(isapprox.(mem_norms, 1.0; atol=1e-8)) "memories are not unit-norm (max dev $(maximum(abs.(mem_norms .- 1))))"
@info "  memory columns unit-norm ✓ (max deviation $(maximum(abs.(mem_norms .- 1.0))))"

@info "Step 0b: β* via find_entropy_inflection …"
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  β* = $(round(β_star, digits=4))"

X_real = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)

safe_cor(X) = (C = cor(X); C[isnan.(C)] .= 0.0; C)
C_real = safe_cor(X_real)

Zreal = (X_real .- std_params_concat.μ') ./ std_params_concat.σ'
rr = [minimum(norm(Zreal[i, :] .- Zreal[j, :]) for j in 1:K if j != i) for i in 1:K]
median_rr = median(rr)
@info "  DCR real→real median (reference): $(round(median_rr, digits=4))"

df_real_complete = df_clean[in.(df_clean.SubjectID, Ref(complete_subjects)), :]
@info "  Running BZ2012 (TF-only, calibrated) on the $K real complete-case patients …"
real_ratios = bz2012_ratios(df_real_complete, TF_TGA_COLS; TM=0.0)

# Shared magnitude channel — byte-identical draw to the ladder's.
Random.seed!(MAGNITUDE_SEED)
shared_magnitudes = [rand(pca_norms) for _ in 1:N_SYNTH]
@info "  Shared magnitude draw (seed=$MAGNITUDE_SEED): range " *
      "[$(round(minimum(shared_magnitudes),digits=3)), $(round(maximum(shared_magnitudes),digits=3))]"

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: exact direction draws + shared decode
# ══════════════════════════════════════════════════════════════════════════════
"""
    exact_directions(β, rng) -> (dirs, components)

Ancestral draw from π, then the published normalization to a unit direction.
`dirs` is `d × N_SYNTH`; `components` records which stored patient seeded each draw.
"""
function exact_directions(β::Float64, rng::Random.AbstractRNG)
    out = exact_sample(X̂, N_SYNTH; β=β, rng=rng)
    dirs = Matrix{Float64}(undef, d, N_SYNTH)
    for i in 1:N_SYNTH
        dirs[:, i] = unitize(out.Ξ[i, :])
    end
    return dirs, out.k
end

function decode_cohort(dirs::Matrix{Float64}, magnitudes::Vector{Float64})
    _, n = size(dirs)
    concat = Matrix{Float64}(undef, n, n_assays * n_visits)
    for i in 1:n
        concat[i, :] = decode_sample(dirs[:, i] .* magnitudes[i], pca_concat, std_params_concat)
    end
    df_visits = DataFrame[]
    for v in 1:n_visits
        offset = (v - 1) * n_assays
        df_v = DataFrame(concat[:, (offset+1):(offset+n_assays)], kept_cols)
        df_v.SyntheticID = 1:n
        df_v.Visit = fill(v, n)
        push!(df_visits, df_v)
    end
    df = vcat(df_visits...)
    sort!(df, [:SyntheticID, :Visit])
    return df
end

"""Fidelity/novelty/privacy metrics. `with_mech=false` skips the BZ2012 ODE run."""
function evaluate_cohort(dirs::Matrix{Float64}, df_cohort::DataFrame; with_mech::Bool=true)
    n = size(dirs, 2)

    all_mre = Float64[]
    for v in 1:n_visits
        comp = feature_summary_comparison(
            df_real_complete[df_real_complete.Visit .== v, :],
            df_cohort[df_cohort.Visit .== v, :], kept_cols)
        append!(all_mre, comp.Mean_Rel_Error)
    end
    med_mre = median(all_mre)

    X_coh = build_concat_matrix(df_cohort, :SyntheticID, collect(1:n), kept_cols, n_assays)
    frob = norm(safe_cor(X_coh) - C_real) / norm(C_real)

    novs = [sample_novelty(dirs[:, i], X̂) for i in 1:n]
    med_nov = median(novs)
    frac_nov = count(>(NOVELTY_THRESH), novs) / n

    Zcoh = (X_coh .- std_params_concat.μ') ./ std_params_concat.σ'
    dcr_sr = [minimum(norm(Zcoh[i, :] .- Zreal[j, :]) for j in 1:K) for i in 1:n]
    med_dcr = median(dcr_sr)

    overlap = ksD = ksp = NaN
    if with_mech
        synth_ratios = bz2012_ratios(df_cohort, TF_TGA_COLS; TM=0.0)
        overlap, ksD, ksp = overlap_ks(real_ratios.Ratio, synth_ratios.Ratio)
    end

    return (Median_MRE=med_mre, CrossVisit_Frob=frob,
            Mech_Overlap_TFonly=overlap, Mech_KS_D=ksD, Mech_KS_p=ksp,
            Median_Novelty=med_nov, Frac_Novelty_gt_thresh=frac_nov,
            DCR_synth_to_real_median=med_dcr, DCR_real_to_real_median=median_rr)
end

@info "\nStep 1: protocol-matched exact cohort (latent seed $EXACT_SEED, shared magnitudes) …"
dirs_exact, comps_exact = exact_directions(β_star, MersenneTwister(EXACT_SEED))
df_exact = decode_cohort(dirs_exact, shared_magnitudes)

X_check = build_concat_matrix(df_exact, :SyntheticID, collect(1:N_SYNTH), kept_cols, n_assays)
@assert size(X_check) == (N_SYNTH, n_assays * n_visits) "decoded cohort wrong shape $(size(X_check))"
@assert all(isfinite, X_check) "decoded cohort has non-finite values"
@info "  decoded cohort finite $(N_SYNTH)×$(n_assays*n_visits) ✓; " *
      "$(length(unique(comps_exact))) of $K memories used as component centers"

@info "  Running BZ2012 on the exact cohort (this takes a few minutes) …"
ref = evaluate_cohort(dirs_exact, df_exact; with_mech=true)
@info "  [Exact] MRE=$(round(100*ref.Median_MRE,digits=2))%  " *
      "Frob=$(round(ref.CrossVisit_Frob,digits=4))  " *
      "overlap=$(round(ref.Mech_Overlap_TFonly,digits=4))  " *
      "nov=$(round(ref.Median_Novelty,digits=4))  " *
      "DCR(s→r)=$(round(ref.DCR_synth_to_real_median,digits=3))  " *
      "[DCR(r→r)=$(round(median_rr,digits=3))]"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Monte Carlo stability — 100 independent exact cohorts
# ══════════════════════════════════════════════════════════════════════════════
# Both the latent stream AND the magnitude stream vary per replicate, so the
# interval reflects the full generation process rather than the direction
# sampler alone. The fixed-seed row above is the protocol-matched comparison.
@info "\nStep 2: $N_REPLICATES independent exact cohorts (fidelity/novelty/DCR; no ODE) …"
rep_rows = NamedTuple[]
for i in 1:N_REPLICATES
    dirs_i, _ = exact_directions(β_star, MersenneTwister(1_000 + i))
    Random.seed!(2_000 + i)
    mags_i = [rand(pca_norms) for _ in 1:N_SYNTH]
    df_i = decode_cohort(dirs_i, mags_i)
    m = evaluate_cohort(dirs_i, df_i; with_mech=false)
    push!(rep_rows, (Replicate=i, LatentSeed=1_000 + i, MagnitudeSeed=2_000 + i,
                     Median_MRE=m.Median_MRE, CrossVisit_Frob=m.CrossVisit_Frob,
                     Median_Novelty=m.Median_Novelty,
                     Frac_Novelty_gt_thresh=m.Frac_Novelty_gt_thresh,
                     DCR_synth_to_real_median=m.DCR_synth_to_real_median))
    i % 20 == 0 && @info "  … $i/$N_REPLICATES"
end
df_reps = DataFrame(rep_rows)

pct(v, p) = quantile(v, p)
summarize(col) = (median(df_reps[!, col]), pct(df_reps[!, col], 0.05), pct(df_reps[!, col], 0.95))

@info "\n  Replicate distribution (median [5th–95th]):"
for (label, col, scale) in [("MRE (%)", :Median_MRE, 100.0),
                            ("cross-visit Frobenius", :CrossVisit_Frob, 1.0),
                            ("median novelty", :Median_Novelty, 1.0),
                            ("fraction novelty > $(NOVELTY_THRESH)", :Frac_Novelty_gt_thresh, 1.0),
                            ("median DCR (s→r)", :DCR_synth_to_real_median, 1.0)]
    m, lo, hi = summarize(col)
    @info "    $label: $(round(scale*m,digits=3)) [$(round(scale*lo,digits=3))–$(round(scale*hi,digits=3))]"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: write CSVs and print the ladder comparison
# ══════════════════════════════════════════════════════════════════════════════
df_ref = DataFrame(
    Rung=["Exact"], Sampler=["analytic ancestral"],
    Median_MRE=[ref.Median_MRE], CrossVisit_Frob=[ref.CrossVisit_Frob],
    Mech_Overlap_TFonly=[ref.Mech_Overlap_TFonly],
    Mech_KS_D=[ref.Mech_KS_D], Mech_KS_p=[ref.Mech_KS_p],
    Median_Novelty=[ref.Median_Novelty],
    Frac_Novelty_gt_thresh=[ref.Frac_Novelty_gt_thresh],
    DCR_synth_to_real_median=[ref.DCR_synth_to_real_median],
    DCR_real_to_real_median=[ref.DCR_real_to_real_median],
    LatentSeed=[EXACT_SEED], MagnitudeSeed=[MAGNITUDE_SEED], N=[N_SYNTH], Beta=[β_star],
)
CSV.write(joinpath(_PATH_TO_DATA, "exact_sampler_reference.csv"), df_ref)
CSV.write(joinpath(_PATH_TO_DATA, "exact_sampler_replicates.csv"), df_reps)
@info "\n  Wrote exact_sampler_reference.csv and exact_sampler_replicates.csv"

@assert ref.Median_MRE > 0 "Median_MRE must be positive"
@assert 0 <= ref.Median_Novelty <= 1 "novelty out of range"
@assert isfinite(ref.Mech_Overlap_TFonly) "mechanistic overlap not computed"

ladder_path = joinpath(_PATH_TO_DATA, "mcmc_ladder_results.csv")
println("\n", "="^92)
println("EXACT ANALYTIC REFERENCE vs THE V0–V3 SAMPLER LADDER  (N=$N_SYNTH, β*=$(round(β_star,digits=4)))")
println("="^92)
row(rung, sampler, mre, frob, ov, nov, dcr) =
    println(rpad(rung, 7), rpad(sampler, 32),
            lpad(mre,  9), lpad(frob, 9), lpad(ov, 9), lpad(nov, 9), lpad(dcr, 10))

row("Rung", "Sampler", "MRE(%)", "Frob", "overlap", "novelty", "DCR(s→r)")
println("-"^92)
if isfile(ladder_path)
    for r in eachrow(CSV.read(ladder_path, DataFrame))
        desc = r.Rung == "V3" ? "sphere+MALA+pool" :
               r.Rung == "V2" ? "sphere+ULA+pool" :
               r.Rung == "V1" ? "sphere+ULA+endpoint" :
                                "anchor+ULA+endpoint (published)"
        row(r.Rung, desc,
            round(100*r.Median_MRE, digits=2), round(r.CrossVisit_Frob, digits=4),
            round(r.Mech_Overlap_TFonly, digits=4), round(r.Median_Novelty, digits=4),
            round(r.DCR_synth_to_real_median, digits=3))
    end
end
row("Exact", "analytic ancestral",
    round(100*ref.Median_MRE, digits=2), round(ref.CrossVisit_Frob, digits=4),
    round(ref.Mech_Overlap_TFonly, digits=4), round(ref.Median_Novelty, digits=4),
    round(ref.DCR_synth_to_real_median, digits=3))
println("-"^92)
println("real→real DCR reference: ", round(median_rr, digits=3))
println("="^92)

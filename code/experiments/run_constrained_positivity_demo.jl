# ──────────────────────────────────────────────────────────────────────────────
# run_constrained_positivity_demo.jl   (npj round 2, Reviewer 3 comments 1 and 3; spec W3)
#
# Adds a positivity constraint term to the stochastic-attention energy:
#
#     pi_lambda(xi | r)  ∝  exp(-beta* [E(xi) + lambda C(xi; r)]),
#     C(xi; r) = sum_j max(0, -x_j / s_j)^2,
#
# where x_j are the six decoded estradiol and progesterone values (3 visits) of
# a profile with magnitude r and direction xi/|xi|, and s_j is the real-cohort
# SD of entry j. C is not quadratic, so for lambda > 0 the target is no longer a
# Gaussian mixture. Two samplers target the same distribution:
#
#   Langevin   ULA on E + lambda C (constrained_ula)
#   Reference  exact mixture draw, kept with probability exp(-beta* lambda C)
#              (rejection_direction); lambda = Inf keeps only C = 0
#
# Each profile's magnitude r is drawn first (canonical magnitudes, seed 424242),
# then its direction is sampled with r fixed. The lambda = 0 arms are new draws
# from the unconstrained target; they are not the published ULA cohort.
#
# Reads:  data/full_longitudinal_memory.jld2, data/cleaned_full_data.csv
# Writes: data/constrained_positivity_demo.csv
#         data/constrained_positivity_agreement.csv
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))
include(joinpath(@__DIR__, "revision2_helpers.jl"))

const LAMBDAS = [0.0, 1.0, 10.0, 100.0]
const N_AGREE = 1000
const FAMILY_ALPHA = 0.01

ref = load_reference()
A = hormone_affine(ref)
bands = strict_bands(ref)
d, K = size(ref.X̂)

# Step size: α = min(0.01, 0.2 / (λ L)), T = ceil(20 / α), so α T = 20 in every arm.
L = 2 * maximum(ref.pca_norms)^2 * maximum(vec(sum(abs2, A.c ./ A.s; dims=2)))
step_size(λ) = λ == 0 ? 0.01 : min(0.01, 0.2 / (λ * L))
n_steps(α) = ceil(Int, 20 / α)
@info "Penalty curvature bound L = $(round(L, digits=3))"

"""Langevin directions (d × n) for magnitudes `mags`, started near a random stored profile."""
function langevin_dirs(λ::Float64, α::Float64, T::Int, mags::Vector{Float64}, rng)
    n = length(mags)
    dirs = Matrix{Float64}(undef, d, n)
    for i in 1:n
        k0 = rand(rng, 1:K)
        ξ0 = r2_unitize(ref.X̂[:, k0] .+ 0.01 .* randn(rng, d))
        ξ = constrained_ula(ref.X̂, ξ0, mags[i], A; β=ref.β_star, λ=λ, α=α, T=T, rng=rng)
        dirs[:, i] = r2_unitize(ξ)
    end
    return dirs
end

"""Reference directions (d × n) by rejection from the exact mixture; also returns total tries."""
function reference_dirs(λ::Float64, mags::Vector{Float64}, rng)
    n = length(mags)
    dirs = Matrix{Float64}(undef, d, n)
    total_tries = 0
    for i in 1:n
        ξ, tries = rejection_direction(ref.X̂, mags[i], A; β=ref.β_star, λ=λ, rng=rng)
        dirs[:, i] = r2_unitize(ξ)
        total_tries += tries
    end
    return dirs, total_tries
end

"""Decoded hormone values (6 × n) for unit directions and magnitudes."""
function hormone_values(dirs::Matrix{Float64}, mags::Vector{Float64})
    n = size(dirs, 2)
    H = Matrix{Float64}(undef, length(A.idx), n)
    for i in 1:n
        H[:, i] = A.a .+ mags[i] .* (A.c * dirs[:, i])
    end
    return H
end

"""All reported metrics for one 100-profile arm."""
function arm_metrics(dirs::Matrix{Float64}, mags::Vector{Float64})
    df = decode_cohort(ref, dirs, mags)
    m = evaluate_cohort(ref, dirs, df; with_mech=true)
    n = size(dirs, 2)

    H = hormone_values(dirs, mags)
    neg_profiles = 0
    for i in 1:n
        if any(H[:, i] .< 0)
            neg_profiles += 1
        end
    end

    hormone_mre = Float64[]
    for v in 1:R2_N_VISITS
        comp = feature_summary_comparison(
            ref.df_real_complete[ref.df_real_complete.Visit .== v, :],
            df[df.Visit .== v, :], [:Estradiol, :Progesterone])
        append!(hormone_mre, comp.Mean_Rel_Error)
    end

    ratios = bz2012_ratios(df, TF_TGA_COLS; TM=0.0)
    bio_pass = 0
    strict_pass = 0
    both_pass = 0
    for id in 1:n
        bio = bio_plausible_patient(df[df.SyntheticID .== id, :], ref)
        strict = strict_mech_pass(ratios[ratios.PatientID .== string(id), :], bands, R2_N_VISITS)
        bio_pass += bio
        strict_pass += strict
        both_pass += (bio && strict)
    end

    return (Bio_Pass=bio_pass, Hormone_Neg_Profiles=neg_profiles,
            Median_MRE=m.Median_MRE, Hormone_MRE=median(hormone_mre),
            CrossVisit_Frob=m.CrossVisit_Frob, Mech_Overlap_TFonly=m.Mech_Overlap_TFonly,
            Strict_Mech_Pass=strict_pass, Both_Pass=both_pass, Median_Novelty=m.Median_Novelty)
end

# ── 100-profile arms (canonical magnitudes, latent seed 42) ────────────────────
mags = ref.shared_magnitudes
rows = NamedTuple[]
for λ in LAMBDAS
    α = step_size(λ)
    T = n_steps(α)
    @info "λ=$λ: Langevin (α=$(round(α, sigdigits=3)), T=$T) …"
    dirs_l = langevin_dirs(λ, α, T, mags, MersenneTwister(42))
    ml = arm_metrics(dirs_l, mags)
    push!(rows, merge((Sampler="Langevin", Lambda=λ, Alpha=α, T=T), ml, (Acceptance_Rate=missing,)))
    @info "  Langevin: bio=$(ml.Bio_Pass) neg=$(ml.Hormone_Neg_Profiles) strict=$(ml.Strict_Mech_Pass) both=$(ml.Both_Pass) " *
          "MRE=$(round(100 * ml.Median_MRE, digits=2))% hMRE=$(round(100 * ml.Hormone_MRE, digits=1))% " *
          "Frob=$(round(ml.CrossVisit_Frob, digits=3)) overlap=$(round(ml.Mech_Overlap_TFonly, digits=3)) nov=$(round(ml.Median_Novelty, digits=3))"

    dirs_r, tries = reference_dirs(λ, mags, MersenneTwister(42))
    mr = arm_metrics(dirs_r, mags)
    push!(rows, merge((Sampler="Reference", Lambda=λ, Alpha=missing, T=missing), mr,
                      (Acceptance_Rate=length(mags) / tries,)))
    @info "  Reference: bio=$(mr.Bio_Pass) neg=$(mr.Hormone_Neg_Profiles) strict=$(mr.Strict_Mech_Pass) both=$(mr.Both_Pass) " *
          "MRE=$(round(100 * mr.Median_MRE, digits=2))% hMRE=$(round(100 * mr.Hormone_MRE, digits=1))% " *
          "acc=$(round(length(mags) / tries, digits=3))"
end
dirs_h, tries_h = reference_dirs(Inf, mags, MersenneTwister(42))
mh = arm_metrics(dirs_h, mags)
push!(rows, merge((Sampler="Reference", Lambda=Inf, Alpha=missing, T=missing), mh,
                  (Acceptance_Rate=length(mags) / tries_h,)))
@info "Hard rejection: bio=$(mh.Bio_Pass) neg=$(mh.Hormone_Neg_Profiles) strict=$(mh.Strict_Mech_Pass) both=$(mh.Both_Pass) " *
      "MRE=$(round(100 * mh.Median_MRE, digits=2))% acc=$(round(length(mags) / tries_h, digits=3))"
df_demo = DataFrame(rows)
CSV.write(joinpath(_PATH_TO_DATA, "constrained_positivity_demo.csv"), df_demo)
@info "Wrote data/constrained_positivity_demo.csv"

# ── Agreement: Langevin vs exact reference, 1000 draws per arm ────────────────
# Family-wise threshold: Bonferroni over (number of λ) × (6 hormone entries).
entry_names = String[]
for v in 1:R2_N_VISITS
    push!(entry_names, "Estradiol_V$v")
    push!(entry_names, "Progesterone_V$v")
end
p_threshold = FAMILY_ALPHA / (length(LAMBDAS) * length(entry_names))

Random.seed!(424_243)
mags_agree = [rand(ref.pca_norms) for _ in 1:N_AGREE]

agree_rows = NamedTuple[]
for λ in LAMBDAS
    dirs_r, _ = reference_dirs(λ, mags_agree, MersenneTwister(2042))
    H_r = hormone_values(dirs_r, mags_agree)
    α = step_size(λ)
    passed = false
    attempt = 1
    local results
    while !passed && attempt <= 2
        T = n_steps(α)
        @info "Agreement λ=$λ attempt $attempt: Langevin α=$(round(α, sigdigits=3)), T=$T, n=$N_AGREE …"
        dirs_l = langevin_dirs(λ, α, T, mags_agree, MersenneTwister(1042))
        H_l = hormone_values(dirs_l, mags_agree)
        results = NamedTuple[]
        for j in eachindex(entry_names)
            t = ApproximateTwoSampleKSTest(H_l[j, :], H_r[j, :])
            push!(results, (Lambda=λ, Entry=entry_names[j], Alpha=α, KS_D=t.δ, KS_p=pvalue(t)))
        end
        min_p = minimum(r.KS_p for r in results)
        @info "  min KS p = $(round(min_p, sigdigits=3)) (threshold $(round(p_threshold, sigdigits=3)))"
        if min_p > p_threshold
            passed = true
        else
            α = α / 2
            attempt += 1
        end
    end
    if !passed
        error("Langevin and the exact reference disagree at λ=$λ after halving the step size")
    end
    append!(agree_rows, results)
end
df_agree = DataFrame(agree_rows)
CSV.write(joinpath(_PATH_TO_DATA, "constrained_positivity_agreement.csv"), df_agree)
@info "  ✓ Langevin agrees with the exact reference at every λ (family-wise α = $FAMILY_ALPHA)"
@info "Wrote data/constrained_positivity_agreement.csv"

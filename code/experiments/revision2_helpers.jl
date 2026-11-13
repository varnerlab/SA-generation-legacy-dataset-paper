# ──────────────────────────────────────────────────────────────────────────────
# revision2_helpers.jl
#
# Shared functions for the npj round-2 minor revision (Reviewer 3):
#   W1  KDE bandwidth rules placed on the beta axis   (run_kde_bandwidth_check.jl)
#   W2  real-patient pass rates on the strict check   (run_real_strict_check.jl)
#   W3  constrained positivity demonstration          (run_constrained_positivity_demo.jl)
#
# The reference loader, decode, and metric functions are copies of Step 0-1 of
# run_exact_sampler_reference.jl that take a `ref` NamedTuple in place of
# script globals. The biological and strict mechanistic checks are copies of
# the predicates in run_convex_hull_analysis.jl. test_revision2_helpers.jl pins
# both copies to the published numbers.
#
# Requires: Include.jl, HockinMannModel, HypothesisTests, mechanistic_eval.jl.
# ──────────────────────────────────────────────────────────────────────────────

using Random, LinearAlgebra, Statistics

const R2_N_SYNTH        = 100
const R2_MAGNITUDE_SEED = 424_242     # identical to run_exact_sampler_reference.jl
const R2_N_VISITS       = 3
const R2_NOVELTY_THRESH = 0.2

r2_unitize(v::AbstractVector{Float64}) = v ./ (norm(v) + 1e-12)

function r2_safe_cor(X::AbstractMatrix{Float64})
    C = cor(X)
    C[isnan.(C)] .= 0.0
    return C
end

"""
    load_reference() -> NamedTuple

Load the canonical memory, beta*, the real-cohort references, the shared
magnitude draw (seed 424242), and the biological-check parameters.
Mirrors Step 0 of run_exact_sampler_reference.jl.
"""
function load_reference()
    mem = JLD2.load(joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2"))
    X̂ = mem["X̂"]
    pca_concat = mem["pca_concat"]
    std_params_concat = mem["std_params_concat"]
    pca_norms = mem["pca_norms"]
    complete_subjects = mem["complete_subjects"]
    kept_cols = mem["kept_cols"]
    n_assays = mem["n_assays"]
    df_clean = mem["df_clean"]

    d, K = size(X̂)
    for k in 1:K
        @assert isapprox(norm(X̂[:, k]), 1.0; atol=1e-8) "memory column $k is not unit-norm"
    end

    phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
    β_star = phase.β_star

    X_real = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)
    C_real = r2_safe_cor(X_real)
    Zreal = (X_real .- std_params_concat.μ') ./ std_params_concat.σ'
    rr = Float64[]
    for i in 1:K
        push!(rr, minimum(norm(Zreal[i, :] .- Zreal[j, :]) for j in 1:K if j != i))
    end
    median_rr = median(rr)

    df_real_complete = df_clean[in.(df_clean.SubjectID, Ref(complete_subjects)), :]
    real_ratios = bz2012_ratios(df_real_complete, TF_TGA_COLS; TM=0.0)

    Random.seed!(R2_MAGNITUDE_SEED)
    shared_magnitudes = [rand(pca_norms) for _ in 1:R2_N_SYNTH]

    # Biological-check parameters (run_convex_hull_analysis.jl, Step 3b (i))
    positive_cols = [:II, :V, :VII, :VIII, :IX, :X, :AT, :PC, :Fbgn,
                     :Estradiol, :Progesterone,
                     Symbol("TF Initiator Peak (nM)"), Symbol("TF Initiator ETP (nM·min)"),
                     Symbol("TF Initiator Lagtime (min)")]
    outlier_features = [:II, :V, :VII, :VIII, :IX, :X, :AT, :PC, :Fbgn,
                        Symbol("TF Initiator Peak (nM)"), Symbol("TF Initiator ETP (nM·min)")]
    positive_cols = filter(c -> c in kept_cols, positive_cols)
    outlier_features = filter(c -> c in kept_cols, outlier_features)
    df_real_full = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
    bio_mu = Dict(col => mean(collect(skipmissing(df_real_full[:, col]))) for col in outlier_features)
    bio_sigma = Dict(col => std(collect(skipmissing(df_real_full[:, col]))) for col in outlier_features)

    return (X̂=X̂, pca_concat=pca_concat, std_params_concat=std_params_concat,
            pca_norms=pca_norms, complete_subjects=complete_subjects, kept_cols=kept_cols,
            n_assays=n_assays, df_clean=df_clean, β_star=β_star, C_real=C_real,
            Zreal=Zreal, median_rr=median_rr, df_real_complete=df_real_complete,
            real_ratios=real_ratios, shared_magnitudes=shared_magnitudes,
            positive_cols=positive_cols, outlier_features=outlier_features,
            bio_mu=bio_mu, bio_sigma=bio_sigma)
end

"""
    exact_directions(ref, β, rng, n) -> (dirs, k)

Draw `n` points from the exact mixture at `β` and normalize each to a unit
direction. `dirs` is `d × n`; `k` holds the selected component of each draw.
"""
function exact_directions(ref, β::Float64, rng::Random.AbstractRNG, n::Int)
    out = exact_sample(ref.X̂, n; β=β, rng=rng)
    d = size(ref.X̂, 1)
    dirs = Matrix{Float64}(undef, d, n)
    for i in 1:n
        dirs[:, i] = r2_unitize(out.Ξ[i, :])
    end
    return dirs, out.k
end

"""
    decode_cohort(ref, dirs, mags) -> DataFrame

Scale each unit direction by its magnitude, decode to measurement space, and
return the long format (one row per profile and visit, `SyntheticID`, `Visit`).
"""
function decode_cohort(ref, dirs::Matrix{Float64}, mags::Vector{Float64})
    n = size(dirs, 2)
    n_assays = ref.n_assays
    concat = Matrix{Float64}(undef, n, n_assays * R2_N_VISITS)
    for i in 1:n
        concat[i, :] = decode_sample(dirs[:, i] .* mags[i], ref.pca_concat, ref.std_params_concat)
    end
    df_visits = DataFrame[]
    for v in 1:R2_N_VISITS
        offset = (v - 1) * n_assays
        df_v = DataFrame(concat[:, (offset + 1):(offset + n_assays)], ref.kept_cols)
        df_v.SyntheticID = 1:n
        df_v.Visit = fill(v, n)
        push!(df_visits, df_v)
    end
    df = vcat(df_visits...)
    sort!(df, [:SyntheticID, :Visit])
    return df
end

"""
    evaluate_cohort(ref, dirs, df; with_mech=true) -> NamedTuple

Median feature MRE, cross-visit correlation Frobenius error, pooled TF-only
mechanistic overlap and KS, novelty, and distance to the closest real profile.
Copy of `evaluate_cohort` in run_exact_sampler_reference.jl.
"""
function evaluate_cohort(ref, dirs::Matrix{Float64}, df_cohort::DataFrame; with_mech::Bool=true)
    n = size(dirs, 2)
    K = size(ref.X̂, 2)

    all_mre = Float64[]
    for v in 1:R2_N_VISITS
        comp = feature_summary_comparison(
            ref.df_real_complete[ref.df_real_complete.Visit .== v, :],
            df_cohort[df_cohort.Visit .== v, :], ref.kept_cols)
        append!(all_mre, comp.Mean_Rel_Error)
    end
    med_mre = median(all_mre)

    X_coh = build_concat_matrix(df_cohort, :SyntheticID, collect(1:n), ref.kept_cols, ref.n_assays)
    frob = norm(r2_safe_cor(X_coh) - ref.C_real) / norm(ref.C_real)

    novs = [sample_novelty(dirs[:, i], ref.X̂) for i in 1:n]
    med_nov = median(novs)
    frac_nov = count(>(R2_NOVELTY_THRESH), novs) / n

    Zcoh = (X_coh .- ref.std_params_concat.μ') ./ ref.std_params_concat.σ'
    dcr_sr = [minimum(norm(Zcoh[i, :] .- ref.Zreal[j, :]) for j in 1:K) for i in 1:n]
    med_dcr = median(dcr_sr)

    overlap = NaN
    ksD = NaN
    ksp = NaN
    if with_mech
        synth_ratios = bz2012_ratios(df_cohort, TF_TGA_COLS; TM=0.0)
        overlap, ksD, ksp = overlap_ks(ref.real_ratios.Ratio, synth_ratios.Ratio)
    end

    return (Median_MRE=med_mre, CrossVisit_Frob=frob,
            Mech_Overlap_TFonly=overlap, Mech_KS_D=ksD, Mech_KS_p=ksp,
            Median_Novelty=med_nov, Frac_Novelty_gt_thresh=frac_nov,
            DCR_synth_to_real_median=med_dcr, DCR_real_to_real_median=ref.median_rr)
end

"""
    bio_plausible_patient(rows, ref) -> Bool

True when no `positive_cols` value is negative and no `outlier_features`
value lies more than 5 real-cohort SDs from the real-cohort mean, across all
of the patient's visit rows. Copy of run_convex_hull_analysis.jl Step 3b (i).
"""
function bio_plausible_patient(rows::AbstractDataFrame, ref)
    for col in ref.positive_cols
        for v in rows[:, col]
            if v < 0
                return false
            end
        end
    end
    for col in ref.outlier_features
        σ = ref.bio_sigma[col]
        if σ <= 1e-12
            continue
        end
        for v in rows[:, col]
            if abs(v - ref.bio_mu[col]) > 5 * σ
                return false
            end
        end
    end
    return true
end

"""
    published_cohort_long(ref) -> DataFrame

The published 100-profile cohort (synthetic_full_longitudinal.csv), projected
to PCA space and decoded, exactly as run_convex_hull_analysis.jl does before
its plausibility checks.
"""
function published_cohort_long(ref)
    df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
    ids = sort(unique(df_synth.SyntheticID))
    n = length(ids)
    n_assays = ref.n_assays
    concat = Matrix{Float64}(undef, n, n_assays * R2_N_VISITS)
    for (i, id) in enumerate(ids)
        for v in 1:R2_N_VISITS
            row = findfirst((df_synth.SyntheticID .== id) .& (df_synth.Visit .== v))
            if row === nothing
                error("published_cohort_long: missing record for id=$id, visit=$v")
            end
            offset = (v - 1) * n_assays
            for (j, col) in enumerate(ref.kept_cols)
                concat[i, offset + j] = Float64(df_synth[row, col])
            end
        end
    end
    Z = (concat .- ref.std_params_concat.μ') ./ ref.std_params_concat.σ'
    Y = MultivariateStats.transform(ref.pca_concat, Z')
    decoded = Matrix{Float64}(undef, n, n_assays * R2_N_VISITS)
    for i in 1:n
        decoded[i, :] = decode_sample(Y[:, i], ref.pca_concat, ref.std_params_concat)
    end
    df_visits = DataFrame[]
    for v in 1:R2_N_VISITS
        offset = (v - 1) * n_assays
        df_v = DataFrame(decoded[:, (offset + 1):(offset + n_assays)], ref.kept_cols)
        df_v.SyntheticID = ids
        df_v.Visit = fill(v, n)
        push!(df_visits, df_v)
    end
    df = vcat(df_visits...)
    sort!(df, [:SyntheticID, :Visit])
    return df
end

"""
    count_bio_pass_published(ref) -> Int

Number of published synthetic profiles that pass `bio_plausible_patient`.
"""
function count_bio_pass_published(ref)
    df = published_cohort_long(ref)
    n_pass = 0
    for id in unique(df.SyntheticID)
        if bio_plausible_patient(df[df.SyntheticID .== id, :], ref)
            n_pass += 1
        end
    end
    return n_pass
end

"""
    strict_bands(ref) -> Dict

Per-feature (5th, 95th) percentile band of the real-cohort predicted/measured
ratios, keyed by the `Feature` values of `ref.real_ratios`.
"""
function strict_bands(ref)
    bands = Dict{String,Tuple{Float64,Float64}}()
    for f in unique(ref.real_ratios.Feature)
        r = ref.real_ratios.Ratio[ref.real_ratios.Feature .== f]
        bands[f] = (quantile(r, 0.05), quantile(r, 0.95))
    end
    return bands
end

"""
    strict_mech_pass(ratios, bands, n_visits) -> Bool

True when the patient has all `n_visits × 5` ratios and every ratio lies
inside its feature's band. A missing ratio (ODE failure) counts as a failure.
"""
function strict_mech_pass(ratios::DataFrame, bands, n_visits::Int)
    if nrow(ratios) != n_visits * length(TF_TGA_COLS)
        return false
    end
    for row in eachrow(ratios)
        lo, hi = bands[row.Feature]
        if row.Ratio < lo || row.Ratio > hi
            return false
        end
    end
    return true
end

"""
    hormone_affine(ref) -> NamedTuple (idx, a, c, s)

For the 6 concatenated columns of estradiol and progesterone (visits 1-3),
the decoded value of a profile with magnitude `r` and unit direction `u` is
`a .+ r .* (c * u)`. `s` is the real-cohort SD of each column.
"""
function hormone_affine(ref)
    idx = Int[]
    for v in 1:R2_N_VISITS
        for h in (:Estradiol, :Progesterone)
            j = findfirst(==(h), ref.kept_cols)
            push!(idx, (v - 1) * ref.n_assays + j)
        end
    end
    P = MultivariateStats.projection(ref.pca_concat)
    μpca = mean(ref.pca_concat)
    μ = ref.std_params_concat.μ
    σ = ref.std_params_concat.σ
    a = μ[idx] .+ σ[idx] .* μpca[idx]
    c = σ[idx] .* P[idx, :]
    s = σ[idx]
    return (idx=idx, a=a, c=c, s=s)
end

"""
    positivity_penalty(ξ, r, A) -> Float64

C(ξ; r) = Σ_j max(0, -x_j / s_j)^2, where x = a + r c ξ̂ are the decoded
hormone values and ξ̂ = ξ / ‖ξ‖.
"""
function positivity_penalty(ξ::AbstractVector{Float64}, r::Float64, A)
    u = ξ ./ norm(ξ)
    x = A.a .+ r .* (A.c * u)
    C = 0.0
    for j in eachindex(x)
        z = -x[j] / A.s[j]
        if z > 0
            C += z^2
        end
    end
    return C
end

"""
    positivity_grad(ξ, r, A) -> Vector{Float64}

Gradient of `positivity_penalty` with respect to ξ, including the
normalization ξ̂ = ξ / ‖ξ‖.
"""
function positivity_grad(ξ::AbstractVector{Float64}, r::Float64, A)
    nξ = norm(ξ)
    u = ξ ./ nξ
    x = A.a .+ r .* (A.c * u)
    dCdx = zeros(length(x))
    for j in eachindex(x)
        z = -x[j] / A.s[j]
        if z > 0
            dCdx[j] = -2 * z / A.s[j]
        end
    end
    g_u = r .* (A.c' * dCdx)
    return (g_u .- u .* dot(u, g_u)) ./ nξ
end

"""
    constrained_ula(X̂, ξ₀, r, A; β, λ, α, T, rng) -> Vector{Float64}

Unadjusted Langevin on E(ξ) + λ C(ξ; r) with uniform multiplicities:
ξ ← ξ - α[(ξ - X̂ softmax(β X̂ᵀ ξ)) + λ ∇C] + sqrt(2α/β) ε.
At λ = 0 this is Eq. (ula) of the manuscript. Returns the final state.
"""
function constrained_ula(X̂::Matrix{Float64}, ξ₀::Vector{Float64}, r::Float64, A;
                         β::Float64, λ::Float64, α::Float64, T::Int,
                         rng::Random.AbstractRNG)
    d = size(X̂, 1)
    ξ = copy(ξ₀)
    noise_scale = sqrt(2 * α / β)
    for t in 1:T
        logits = β .* (X̂' * ξ)
        logits .-= maximum(logits)
        p = exp.(logits)
        p ./= sum(p)
        grad = ξ .- X̂ * p
        if λ > 0
            grad .+= λ .* positivity_grad(ξ, r, A)
        end
        ξ = ξ .- α .* grad .+ noise_scale .* randn(rng, d)
    end
    return ξ
end

"""
    rejection_direction(X̂, r, A; β, λ, rng, max_tries=1_000_000) -> (ξ, tries)

Exact sampler for π(ξ) exp(-β λ C(ξ; r)): draw ξ from the exact mixture and
keep it with probability exp(-β λ C). With λ = Inf, keep only draws with
C = 0 (all six hormone values nonnegative). Throws after `max_tries` draws.
"""
function rejection_direction(X̂::Matrix{Float64}, r::Float64, A;
                             β::Float64, λ::Float64, rng::Random.AbstractRNG,
                             max_tries::Int=1_000_000)
    for tries in 1:max_tries
        ξ = exact_sample(X̂, 1; β=β, rng=rng).Ξ[1, :]
        C = positivity_penalty(ξ, r, A)
        if isinf(λ)
            if C == 0.0
                return ξ, tries
            end
        else
            if rand(rng) < exp(-β * λ * C)
                return ξ, tries
            end
        end
    end
    error("rejection_direction: no draw accepted in $max_tries tries (r=$r, λ=$λ)")
end

"""
    loo_loglik(X̂, h) -> Float64

Leave-one-out log-likelihood of an isotropic Gaussian KDE with bandwidth h on
the columns of X̂: Σ_i log[(1/(K-1)) Σ_{j≠i} N(m_i; m_j, h² I)].
"""
function loo_loglik(X̂::Matrix{Float64}, h::Float64)
    d, K = size(X̂)
    total = 0.0
    for i in 1:K
        terms = Float64[]
        for j in 1:K
            if j != i
                push!(terms, -sum(abs2, X̂[:, i] .- X̂[:, j]) / (2 * h^2))
            end
        end
        m = maximum(terms)
        lse = m + log(sum(exp.(terms .- m)))
        total += lse - log(K - 1) - (d / 2) * log(2π * h^2)
    end
    return total
end

"""
    kde_bandwidths(X̂) -> NamedTuple (sigma, scott, silverman, loo)

Standard bandwidth rules for an isotropic Gaussian KDE on the unit-normalized
profiles. `sigma` is the root mean per-coordinate variance. `loo` maximizes
`loo_loglik` over h in [0.01, 3] (log grid, then golden-section refinement).
"""
function kde_bandwidths(X̂::Matrix{Float64})
    d, K = size(X̂)
    sigma = sqrt(mean(var(X̂; dims=2)))
    scott = sigma * K^(-1 / (d + 4))
    silverman = sigma * (4 / (d + 2))^(1 / (d + 4)) * K^(-1 / (d + 4))

    # coarse log grid
    hs = exp.(range(log(0.01), log(3.0), length=400))
    ll = [loo_loglik(X̂, h) for h in hs]
    ib = argmax(ll)
    lo = log(hs[max(ib - 1, 1)])
    hi = log(hs[min(ib + 1, length(hs))])

    # golden-section refinement in log h
    gr = (sqrt(5) - 1) / 2
    x1 = hi - gr * (hi - lo)
    x2 = lo + gr * (hi - lo)
    f1 = loo_loglik(X̂, exp(x1))
    f2 = loo_loglik(X̂, exp(x2))
    for _ in 1:80
        if f1 > f2
            hi = x2
            x2 = x1
            f2 = f1
            x1 = hi - gr * (hi - lo)
            f1 = loo_loglik(X̂, exp(x1))
        else
            lo = x1
            x1 = x2
            f1 = f2
            x2 = lo + gr * (hi - lo)
            f2 = loo_loglik(X̂, exp(x2))
        end
    end
    loo = exp((lo + hi) / 2)
    return (sigma=sigma, scott=scott, silverman=silverman, loo=loo)
end

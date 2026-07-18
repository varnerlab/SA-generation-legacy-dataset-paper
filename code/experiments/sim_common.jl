# ──────────────────────────────────────────────────────────────────────────────
# sim_common.jl
#
# Shared helpers for the SIM simulation benchmark (Tasks 5, 6, 7). Included
# (not run standalone) by driver scripts AFTER Include.jl, e.g.:
#   include(joinpath(@__DIR__, "..", "Include.jl"))
#   include(joinpath(@__DIR__, "sim_common.jl"))
#
# Orientation convention (binding everywhere in this file): every data matrix
# is n×p — rows are observations, columns are features/variables. `cov(X)`
# for such a matrix already treats rows as observations, so NONE of the
# functions below transpose before calling `cov`. `Σpop` and any covariance
# estimate `Ĉ` are p×p.
# ──────────────────────────────────────────────────────────────────────────────

using HypothesisTests: ApproximateTwoSampleKSTest

# ══════════════════════════════════════════════════════════════════════════════
# Ground truth: low-rank latent factor model
# ══════════════════════════════════════════════════════════════════════════════

"""
    _heavy_transform(X; ν=3.0) -> Matrix

Fixed monotone heavy-tailed marginal transform applied column-wise: each
column is mapped through the Gaussian CDF then the quantile function of a
Student-t(ν) distribution (a Gaussian-copula → Student-t marginal warp).
Monotone per column ⇒ preserves the copula (rank) structure of the latent
Gaussian factor draw while inflating marginal kurtosis and changing the
linear (Pearson) covariance — which is exactly why the heavy-tail Σpop must
be estimated empirically rather than taken analytically (see
`make_lowrank_gt`).
"""
function _heavy_transform(X::AbstractMatrix; ν::Float64=3.0)
    u = cdf.(Normal(), X)
    u = clamp.(u, 1e-10, 1 - 1e-10)   # guard quantile() blow-up at the tails
    return quantile.(TDist(ν), u)
end

"""
    make_lowrank_gt(p, r, rng; tail=:gaussian, ψ=0.1) -> (sampler, Σpop)

Latent factor ground truth in the n<p regime. Loadings `L` (p×r) and
idiosyncratic noise variance `ψ` define
    Z = randn(rng, n, r); E = sqrt(ψ) .* randn(rng, n, p); X = Z*L' + E
`sampler(n)` returns an n×p matrix (rows = observations), consistent with the
orientation convention above.

- `tail=:gaussian`: `sampler` returns the raw linear-Gaussian draw; the
  population covariance is analytic, `Σpop = L*L' + ψ*I`.
- `tail=:heavy`: `sampler` additionally applies `_heavy_transform` to every
  draw (so every call, including the caller's later `sampler(n)` calls for
  Xtrain/Xtest, emits heavy-tailed data). Because that nonlinear transform
  changes the covariance, `Σpop` is NOT `L*L'+ψ*I` — it is estimated as the
  sample covariance of a fresh 50,000-row draw of the SAME (transformed)
  process, computed once here before returning, so it is the true population
  covariance of what `sampler` emits. (No transpose: `cov(big_draw)` on the
  50,000×p draw already treats rows as observations, per the orientation
  convention at the top of this file.)

Note: for `tail=:heavy`, computing the 50,000-row reference draw consumes
`rng` state before returning. This is deterministic given `rng`'s incoming
state (same seed ⇒ same Σpop ⇒ same subsequent draws), so reproducibility is
preserved; it merely means `sampler`'s *first* call by the caller continues
the RNG stream from *after* that internal reference draw, not from `rng`'s
original state.
"""
function make_lowrank_gt(p::Int, r::Int, rng::AbstractRNG; tail::Symbol=:gaussian, ψ::Float64=0.1)
    L = randn(rng, p, r)

    raw_sampler(n::Int) = Matrix{Float64}(randn(rng, n, r) * L' .+ sqrt(ψ) .* randn(rng, n, p))

    if tail == :gaussian
        Σpop = Matrix{Float64}(L * L' + ψ * I)
        return raw_sampler, Σpop
    elseif tail == :heavy
        sampler(n::Int) = _heavy_transform(raw_sampler(n))
        big_draw = sampler(50_000)
        Σpop = Matrix{Float64}(cov(big_draw))
        return sampler, Σpop
    else
        throw(ArgumentError("tail must be :gaussian or :heavy, got $(repr(tail))"))
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Deterministic per-cell seeding (no `rand()`)
# ══════════════════════════════════════════════════════════════════════════════

"""
    cell_seed(i_n, i_p, i_r, i_tail; base=42) -> Int

Deterministic seed derived from 1-based grid-cell indices (not `rand()`), so
sweeps are reproducible and every cell gets an independent, collision-free
seed. `base` anchors the whole grid to the canonical seed 42.
"""
cell_seed(i_n::Int, i_p::Int, i_r::Int, i_tail::Int; base::Int=42) =
    base + 100_000 * i_tail + 10_000 * i_r + 1_000 * i_p + 10 * i_n

# ══════════════════════════════════════════════════════════════════════════════
# Recovery methods
# ══════════════════════════════════════════════════════════════════════════════

"""
    sa_recover(Xtrain; T=2000, seed=42) -> Xsyn::Matrix (100×p)

Generate `N=100` SA synthetics from the n×p training cohort by reusing the
frozen helpers `sa_generate_from_matrix` and `build_concat_matrix`
(`code/src/Patient.jl`, read-only). Every grid `p` is divisible by 3, so the
p columns are treated as `n_assays = p÷3` features across 3 formal "visits"
— purely a reshape convenience to fit the concatenated-visit SA pipeline;
`sa_generate_from_matrix`'s internal PCA/ULA operate on all p columns jointly
regardless of this reshape, so no genuine longitudinal structure is assumed
or required. Canonical hyperparameters: α=0.01, pratio=0.95, β=nothing
(recomputed per cohort).
"""
function sa_recover(Xtrain::AbstractMatrix; T::Int=2000, seed::Int=42)
    Xtrain_f = Matrix{Float64}(Xtrain)
    n, p = size(Xtrain_f)
    na = p ÷ 3
    na * 3 == p || throw(ArgumentError("p=$p must be divisible by 3"))
    kept = [Symbol("f", i) for i in 1:na]
    synth_df, _ = sa_generate_from_matrix(Xtrain_f, kept, na, 100, T;
                                           β=nothing, α=0.01, pratio=0.95, seed=seed)
    ids = unique(synth_df.SyntheticID)
    return build_concat_matrix(synth_df, :SyntheticID, ids, kept, na)
end

"""
    mvn_recover(Xtrain; seed=42, n_samples=100) -> Xsyn::Matrix (100×p)

Fit a Gaussian to the training cohort honestly: `μ = mean(Xtrain, dims=1)`,
`Σ = cov(Xtrain)`. This mirrors the fitting style of
`paper_sa_vs_mvn.jl:137-176` (symmetrize, jitter, `MvNormal`, draw + transpose
back to n×p) but — unlike that reference, which uses a fixed 1e-6 jitter
unconditionally — escalates jitter only as far as strictly needed for a
numerically valid Cholesky draw, starting from the brief's minimal 1e-8. No
shrinkage toward Σpop and no rescue of the rank deficiency: when n≤p, Σ is
rank-deficient and the recovered covariance is expected to be poor — that
degradation is the point of the benchmark, so it is never hidden. If even
escalated jitter fails to produce a valid Cholesky factor, falls back to an
eigen-based draw retaining all available (nonnegative-eigenvalue) directions
so a result is always recorded, never silently skipped.
"""
function mvn_recover(Xtrain::AbstractMatrix; seed::Int=42, n_samples::Int=100)
    Xtrain_f = Matrix{Float64}(Xtrain)
    n, p = size(Xtrain_f)
    μ = vec(mean(Xtrain_f, dims=1))
    Σ = cov(Xtrain_f)
    Σ = (Σ + Σ') / 2   # symmetrize: numerical hygiene only, not shrinkage

    rng = MersenneTwister(seed)

    jitter = 1e-8
    max_jitter = 1e-2
    Xsyn = nothing
    last_err = nothing
    while jitter <= max_jitter
        try
            dist = MvNormal(μ, Σ + jitter * I)
            draws = rand(rng, dist, n_samples)   # p × n_samples
            Xsyn = Matrix{Float64}(draws')       # n_samples × p
            break
        catch e
            last_err = e
            jitter *= 10
        end
    end

    if Xsyn === nothing
        @warn "mvn_recover: jittered Cholesky failed up to jitter=$max_jitter " *
              "(n=$n, p=$p, last error: $(typeof(last_err))) — falling back to " *
              "eigen-based draw retaining all available directions."
        F = eigen(Symmetric(Σ))
        λ = max.(F.values, 0.0)   # clip negative/round-off eigenvalues, keep the rest as-is
        V = F.vectors
        Z = randn(rng, p, n_samples)
        draws = V * Diagonal(sqrt.(λ)) * Z .+ μ   # p × n_samples
        Xsyn = Matrix{Float64}(draws')
        n_pos = count(>(1e-10), λ)
        @warn "mvn_recover: eigen fallback retained $n_pos / $p positive-eigenvalue directions."
    end

    return Xsyn
end

# ══════════════════════════════════════════════════════════════════════════════
# Metrics
# ══════════════════════════════════════════════════════════════════════════════

"""
    cov_frob(Ĉ, Σpop) -> Float64

Relative Frobenius covariance-recovery error, `norm(Ĉ - Σpop) / norm(Σpop)`.
`norm(::AbstractMatrix)` in Julia's LinearAlgebra is the entrywise (Frobenius)
norm by default (NOT the operator/spectral norm — that's `opnorm`), so no
extra `sqrt(sum(abs2, ·))` is needed.
"""
cov_frob(Ĉ::AbstractMatrix, Σpop::AbstractMatrix) = norm(Ĉ - Σpop) / norm(Σpop)

"""
    marginal_mre(Xsyn, Xtest) -> Float64

Median over columns j of the relative marginal mean/std recovery error
against a fresh large GT test draw `Xtest`:
    err_j = 0.5*(|mean(Xsyn[:,j])-mean(Xtest[:,j])| / (|mean(Xtest[:,j])|+1e-8)
              + |std(Xsyn[:,j])-std(Xtest[:,j])| / (std(Xtest[:,j])+1e-8))
"""
function marginal_mre(Xsyn::AbstractMatrix, Xtest::AbstractMatrix)
    p = size(Xsyn, 2)
    size(Xtest, 2) == p || throw(DimensionMismatch("Xsyn has $p cols, Xtest has $(size(Xtest,2))"))
    errs = Vector{Float64}(undef, p)
    for j in 1:p
        μs, μt = mean(view(Xsyn, :, j)), mean(view(Xtest, :, j))
        σs, σt = std(view(Xsyn, :, j)), std(view(Xtest, :, j))
        errs[j] = 0.5 * (abs(μs - μt) / (abs(μt) + 1e-8) + abs(σs - σt) / (σt + 1e-8))
    end
    return median(errs)
end

"""
    marginal_ks(Xsyn, Xtest) -> Float64

Mean over columns of the two-sample KS statistic
`ApproximateTwoSampleKSTest(Xsyn[:,j], Xtest[:,j]).δ`.
"""
function marginal_ks(Xsyn::AbstractMatrix, Xtest::AbstractMatrix)
    p = size(Xsyn, 2)
    size(Xtest, 2) == p || throw(DimensionMismatch("Xsyn has $p cols, Xtest has $(size(Xtest,2))"))
    δs = Vector{Float64}(undef, p)
    for j in 1:p
        δs[j] = ApproximateTwoSampleKSTest(view(Xsyn, :, j), view(Xtest, :, j)).δ
    end
    return mean(δs)
end

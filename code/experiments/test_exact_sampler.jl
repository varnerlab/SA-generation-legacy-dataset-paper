# ──────────────────────────────────────────────────────────────────────────────
# test_exact_sampler.jl — correctness checks for the analytic (exact) sampler
#
# Verifies that `exact_sample` draws from the same target that the ULA sampler
# in `sample`/`weighted_sample` targets, namely
#
#     π_r(ξ) ∝ exp[-β E_r(ξ)] = Σ_k q_k N(m_k, β⁻¹ I),
#     q_k ∝ r_k exp(β‖m_k‖²/2).
#
# Gate for the exact-sampler manuscript update: all five checks must pass.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

const TOL_DENSITY = 1e-8     # spread of the log-density offset across probes
const TOL_WEIGHT  = 1e-12    # component-probability identities
const TOL_MEAN    = 0.02     # Monte Carlo mean error (n = 400_000)
const TOL_COV     = 0.03     # Monte Carlo covariance error (n = 400_000)

results = Tuple{String,Bool,String}[]
record!(name, ok, detail) = push!(results, (name, ok, detail))

# ── helpers ───────────────────────────────────────────────────────────────────

"""Unnormalized log target implied by the weighted Hopfield energy: -βE_r(ξ)."""
function log_hopfield_target(X, ξ, β, r)
    logits = β .* (X' * ξ) .+ log.(r)
    m = maximum(logits)
    lse = m + log(sum(exp.(logits .- m)))
    return -0.5 * β * dot(ξ, ξ) + lse
end

"""Unnormalized log density of the Gaussian mixture Σ_k q_k N(m_k, β⁻¹I)."""
function log_mixture_density(X, ξ, β, q)
    K = size(X, 2)
    terms = Vector{Float64}(undef, K)
    for k in 1:K
        diff = ξ .- view(X, :, k)
        terms[k] = log(q[k]) - 0.5 * β * dot(diff, diff)
    end
    m = maximum(terms)
    return m + log(sum(exp.(terms .- m)))
end

# ── Test 1: density identity ──────────────────────────────────────────────────
# The two log densities must differ by a single constant, independent of ξ.

let
    rng = MersenneTwister(1)
    d, K, β = 6, 9, 2.94
    X = randn(rng, d, K)                       # deliberately NON-unit columns
    r = 0.5 .+ rand(rng, K)                    # non-uniform multiplicities

    q = exact_sample(X, 1; β=β, weights=r, rng=MersenneTwister(0)).q

    offsets = Float64[]
    for _ in 1:200
        ξ = 2.0 .* randn(rng, d)
        push!(offsets, log_hopfield_target(X, ξ, β, r) - log_mixture_density(X, ξ, β, q))
    end
    spread = maximum(offsets) - minimum(offsets)
    record!("1. density identity (general norms)", spread < TOL_DENSITY,
            "offset spread over 200 probes = $(spread)")
end

# ── Test 2: unit-norm reduction  q_k = r_k / Σ r_j ────────────────────────────

let
    rng = MersenneTwister(2)
    d, K, β = 8, 12, 2.94
    X = randn(rng, d, K)
    for k in 1:K
        X[:, k] ./= norm(X[:, k])              # unit-norm memories
    end
    r = 0.5 .+ rand(rng, K)

    q = exact_sample(X, 1; β=β, weights=r, rng=MersenneTwister(0)).q
    err = maximum(abs.(q .- r ./ sum(r)))
    record!("2. unit-norm reduction q ∝ r", err < TOL_WEIGHT,
            "max |q_k - r_k/Σr| = $(err)")
end

# ── Test 3: general-norm correction  q_k ∝ r_k exp(β‖m_k‖²/2) ─────────────────

let
    rng = MersenneTwister(3)
    d, K, β = 5, 7, 1.7
    X = randn(rng, d, K)
    for k in 1:K                                # spread the norms widely
        X[:, k] .*= (0.3 + 0.4 * k)
    end
    r = 0.5 .+ rand(rng, K)

    q = exact_sample(X, 1; β=β, weights=r, rng=MersenneTwister(0)).q
    expected = [r[k] * exp(0.5 * β * dot(X[:, k], X[:, k])) for k in 1:K]
    expected ./= sum(expected)
    err = maximum(abs.(q .- expected))
    record!("3. general-norm weight correction", err < TOL_WEIGHT,
            "max |q_k - expected| = $(err)")
end

# ── Test 4: moment check against analytic mixture mean and covariance ─────────
# mean = Σ q_k m_k;  cov = β⁻¹I + Σ q_k (m_k - μ)(m_k - μ)ᵀ

let
    rng = MersenneTwister(4)
    d, K, β = 4, 6, 2.94
    X = randn(rng, d, K)
    for k in 1:K
        X[:, k] ./= norm(X[:, k])
    end
    r = 0.5 .+ rand(rng, K)
    n = 400_000

    out = exact_sample(X, n; β=β, weights=r, rng=MersenneTwister(11))
    q = out.q

    μ_analytic = X * q
    Σ_analytic = Matrix(I * (1.0 / β), d, d)
    for k in 1:K
        δ = X[:, k] .- μ_analytic
        Σ_analytic .+= q[k] .* (δ * δ')
    end

    μ_emp = vec(mean(out.Ξ, dims=1))
    Σ_emp = cov(out.Ξ)

    mean_err = maximum(abs.(μ_emp .- μ_analytic))
    cov_err  = maximum(abs.(Σ_emp .- Σ_analytic))
    record!("4. Monte Carlo moments (n = $n)",
            mean_err < TOL_MEAN && cov_err < TOL_COV,
            "max mean err = $(round(mean_err, digits=5)), " *
            "max cov err = $(round(cov_err, digits=5))")
end

# ── Test 5: seed reproducibility ──────────────────────────────────────────────

let
    rng = MersenneTwister(5)
    d, K, β = 7, 10, 2.94
    X = randn(rng, d, K)
    for k in 1:K
        X[:, k] ./= norm(X[:, k])
    end
    r = 0.5 .+ rand(rng, K)

    a = exact_sample(X, 500; β=β, weights=r, rng=MersenneTwister(42))
    b = exact_sample(X, 500; β=β, weights=r, rng=MersenneTwister(42))
    c = exact_sample(X, 500; β=β, weights=r, rng=MersenneTwister(43))

    same = (a.Ξ == b.Ξ) && (a.k == b.k)
    differs = (a.Ξ != c.Ξ)
    record!("5. seed reproducibility", same && differs,
            "seed 42 twice identical = $(same); seed 43 differs = $(differs)")
end

# ── report ────────────────────────────────────────────────────────────────────

println("\n", "="^74)
println("EXACT SAMPLER — CORRECTNESS CHECKS")
println("="^74)
for (name, ok, detail) in results
    println(ok ? "  PASS  " : "  FAIL  ", name)
    println("        ", detail)
end
n_pass = count(r -> r[2], results)
println("-"^74)
println("$(n_pass)/$(length(results)) checks passed")
println("="^74)

if n_pass != length(results)
    error("exact sampler correctness gate FAILED — do not proceed to manuscript edits")
end

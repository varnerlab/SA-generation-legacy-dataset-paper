"""
    sample(X, ξ₀, T; β=1.0, α=0.1, seed=nothing) -> (t, Ξ)

Run the Stochastic Attention Sampler (Algorithm 1) on a memory matrix `X`.

The sampler implements the Unadjusted Langevin Algorithm (ULA) on the modern
Hopfield energy ``E(\\boldsymbol{\\xi}) = \\tfrac{1}{2}\\|\\boldsymbol{\\xi}\\|^2
- \\tfrac{1}{\\beta}\\log\\sum_k \\exp(\\beta\\,\\mathbf{m}_k^\\top \\boldsymbol{\\xi})``,
yielding the update

```math
\\boldsymbol{\\xi}_{t+1}
= (1-\\alpha)\\,\\boldsymbol{\\xi}_t
  + \\alpha\\,\\mathbf{X}\\,\\operatorname{softmax}(\\beta\\,\\mathbf{X}^\\top\\boldsymbol{\\xi}_t)
  + \\sqrt{2\\alpha/\\beta}\\;\\boldsymbol{\\epsilon}_t,
\\qquad \\boldsymbol{\\epsilon}_t \\sim \\mathcal{N}(\\mathbf{0}, \\mathbf{I}).
```

# Arguments (positional)
- `X::Matrix{Float64}`: Memory matrix of size `d × K`, where each column is a
  stored pattern of dimension `d`.
- `ξ₀::Vector{Float64}`: Initial state vector of length `d`.
- `T::Int`: Number of iterations to run.

# Keyword Arguments
- `β::Float64=1.0`: Inverse temperature.
- `α::Float64=0.1`: Step size (learning rate) in `(0, 1)`.
- `seed::Union{Int, Nothing}=nothing`: Optional random seed.

# Returns
A named tuple `(t, Ξ)` where:
- `t::Vector{Int}`: Time indices `[0, 1, …, T]` of length `T + 1`.
- `Ξ::Matrix{Float64}`: State trajectory of size `(T + 1) × d`.
"""
function sample(X::Matrix{Float64}, ξ₀::Vector{Float64}, T::Int;
    β::Float64 = 1.0, α::Float64 = 0.1,
    seed::Union{Int, Nothing} = nothing)

    # --- input validation ---
    d, K = size(X)
    length(ξ₀) == d || throw(DimensionMismatch(
        "Initial state has length $(length(ξ₀)) but X has $d rows"))
    T > 0   || throw(ArgumentError("T must be positive, got T = $T"))
    β > 0   || throw(ArgumentError("β must be positive, got β = $β"))
    0 < α < 1 || throw(ArgumentError("α must be in (0,1), got α = $α"))

    # --- seed the RNG if requested ---
    if seed !== nothing
        Random.seed!(seed)
    end

    # --- pre-allocate output ---
    Ξ = Matrix{Float64}(undef, T + 1, d)
    Ξ[1, :] .= ξ₀

    # --- pre-compute constants ---
    noise_scale = sqrt(2.0 * α / β)

    # --- run the Langevin chain ---
    ξ = copy(ξ₀)
    for t in 1:T
        # attention weights  a = softmax(β X⊤ ξ)
        logits = β .* (X' * ξ)
        w = softmax(logits)

        # Gaussian noise  ε ~ N(0, I_d)
        ε = randn(d)

        # Langevin update
        ξ .= (1.0 - α) .* ξ .+ α .* (X * w) .+ noise_scale .* ε

        Ξ[t + 1, :] .= ξ
    end

    t_vec = collect(0:T)
    return (t = t_vec, Ξ = Ξ)
end


"""
    mala_sample(X, ξ₀, T; β=1.0, α=0.1, seed=nothing) -> (t, Ξ, accept_rate)

Run the Metropolis-Adjusted Langevin Algorithm (MALA) on the modern Hopfield energy.

Uses the same Langevin proposal as `sample` but adds a Metropolis–Hastings
acceptance correction to eliminate discretization bias.

# Returns
A named tuple `(t, Ξ, accept_rate)` where `accept_rate` is the fraction of
accepted proposals.
"""
function mala_sample(X::Matrix{Float64}, ξ₀::Vector{Float64}, T::Int;
    β::Float64 = 1.0, α::Float64 = 0.1,
    seed::Union{Int, Nothing} = nothing)

    d, K = size(X)
    length(ξ₀) == d || throw(DimensionMismatch(
        "Initial state has length $(length(ξ₀)) but X has $d rows"))
    T > 0   || throw(ArgumentError("T must be positive, got T = $T"))
    β > 0   || throw(ArgumentError("β must be positive, got β = $β"))
    0 < α < 1 || throw(ArgumentError("α must be in (0,1), got α = $α"))

    if seed !== nothing
        Random.seed!(seed)
    end

    Ξ = Matrix{Float64}(undef, T + 1, d)
    Ξ[1, :] .= ξ₀

    noise_scale = sqrt(2.0 * α / β)
    noise_var   = 2.0 * α / β
    log_norm_const = -0.5 * d * log(2π * noise_var)
    inv_2var = 1.0 / (2.0 * noise_var)

    n_accept = 0

    function langevin_mean(ξ::Vector{Float64})
        logits = β .* (X' * ξ)
        w = softmax(logits)
        return (1.0 - α) .* ξ .+ α .* (X * w)
    end

    function neg_log_target(ξ::Vector{Float64})
        logits = β .* (X' * ξ)
        logits_max = maximum(logits)
        lse = logits_max + log(sum(exp.(logits .- logits_max)))
        return 0.5 * β * dot(ξ, ξ) - lse
    end

    function log_proposal(y::Vector{Float64}, μx::Vector{Float64})
        diff = y .- μx
        return log_norm_const - inv_2var * dot(diff, diff)
    end

    ξ = copy(ξ₀)
    μ_curr = langevin_mean(ξ)
    nlt_curr = neg_log_target(ξ)

    for t in 1:T
        ε = randn(d)
        ξ_prop = μ_curr .+ noise_scale .* ε

        nlt_prop = neg_log_target(ξ_prop)
        μ_prop   = langevin_mean(ξ_prop)

        log_accept = (-nlt_prop + log_proposal(ξ, μ_prop)) -
                     (-nlt_curr + log_proposal(ξ_prop, μ_curr))

        if log(rand()) < log_accept
            ξ .= ξ_prop
            μ_curr  = μ_prop
            nlt_curr = nlt_prop
            n_accept += 1
        end

        Ξ[t + 1, :] .= ξ
    end

    t_vec = collect(0:T)
    accept_rate = n_accept / T

    return (t = t_vec, Ξ = Ξ, accept_rate = accept_rate)
end


"""
    weighted_sample(X, ξ₀, T, ρ; β=1.0, α=0.1, seed=nothing) -> (t, Ξ)

Run the multiplicity-weighted SA sampler. Adds a log-multiplicity bias `ρ` to the
softmax logits to tilt generation toward a designated subset of stored patterns.

The weighted update replaces the attention logits with:
```math
\\ell_k = \\beta\\,\\mathbf{m}_k^\\top \\boldsymbol{\\xi}_t + \\log \\rho_k
```

# Arguments
- `X::Matrix{Float64}`: Memory matrix `d × K`.
- `ξ₀::Vector{Float64}`: Initial state.
- `T::Int`: Number of iterations.
- `ρ::Vector{Float64}`: Multiplicity weights (length `K`), all positive.

# Keyword Arguments
- `β`, `α`, `seed`: Same as `sample`.
"""
function weighted_sample(X::Matrix{Float64}, ξ₀::Vector{Float64}, T::Int,
    ρ::Vector{Float64};
    β::Float64 = 1.0, α::Float64 = 0.1,
    seed::Union{Int, Nothing} = nothing)

    d, K = size(X)
    length(ξ₀) == d || throw(DimensionMismatch(
        "Initial state has length $(length(ξ₀)) but X has $d rows"))
    length(ρ) == K || throw(DimensionMismatch(
        "ρ has length $(length(ρ)) but X has $K columns"))
    T > 0   || throw(ArgumentError("T must be positive"))
    β > 0   || throw(ArgumentError("β must be positive"))
    0 < α < 1 || throw(ArgumentError("α must be in (0,1)"))
    all(ρ .> 0) || throw(ArgumentError("All multiplicity weights must be positive"))

    if seed !== nothing
        Random.seed!(seed)
    end

    Ξ = Matrix{Float64}(undef, T + 1, d)
    Ξ[1, :] .= ξ₀

    noise_scale = sqrt(2.0 * α / β)
    log_ρ = log.(ρ)

    ξ = copy(ξ₀)
    for t in 1:T
        logits = β .* (X' * ξ) .+ log_ρ
        w = softmax(logits)
        ε = randn(d)
        ξ .= (1.0 - α) .* ξ .+ α .* (X * w) .+ noise_scale .* ε
        Ξ[t + 1, :] .= ξ
    end

    t_vec = collect(0:T)
    return (t = t_vec, Ξ = Ξ)
end

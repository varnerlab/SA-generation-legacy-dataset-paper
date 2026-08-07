"""
    nearest_cosine_similarity(ξ, X) -> Float64

Maximum cosine similarity between state `ξ` and columns of memory matrix `X`.
"""
function nearest_cosine_similarity(ξ::Vector{Float64}, X::Matrix{Float64})::Float64
    nξ = norm(ξ)
    nξ < 1e-12 && return 0.0
    max_sim = -Inf
    for k in 1:size(X, 2)
        mk = @view X[:, k]
        sim_k = dot(mk, ξ) / (norm(mk) * nξ)
        max_sim = max(max_sim, sim_k)
    end
    return max_sim
end

"""
    hopfield_energy(ξ, X, β) -> Float64

Modern Hopfield energy at state `ξ`:
E(ξ) = ½‖ξ‖² - (1/β) log Σ_k exp(β mₖᵀξ)
"""
function hopfield_energy(ξ::Vector{Float64}, X::Matrix{Float64}, β::Float64)::Float64
    logits = β .* (X' * ξ)
    logits_max = maximum(logits)
    lse = logits_max + log(sum(exp.(logits .- logits_max)))
    return 0.5 * dot(ξ, ξ) - lse / β
end

"""
    attention_entropy(ξ, X, β) -> Float64

Shannon entropy of the attention weights p = softmax(β Xᵀξ).
"""
function attention_entropy(ξ::Vector{Float64}, X::Matrix{Float64}, β::Float64)::Float64
    p = NNlib.softmax(β .* (X' * ξ))
    H = 0.0
    for pk in p
        if pk > 1e-30
            H -= pk * log(pk)
        end
    end
    return H
end

"""
    sample_novelty(ξ, X) -> Float64

Novelty = 1 - max_k cos(ξ, mₖ). Zero means exact copy of a stored pattern.
"""
function sample_novelty(ξ::Vector{Float64}, X::Matrix{Float64})::Float64
    return 1.0 - nearest_cosine_similarity(ξ, X)
end

"""
    sample_diversity(samples) -> Float64

Mean pairwise cosine distance across a collection of generated samples.
"""
function sample_diversity(samples::Vector{Vector{Float64}})::Float64
    S = length(samples)
    S < 2 && return 0.0
    total = 0.0
    count = 0
    for i in 1:(S-1)
        nᵢ = norm(samples[i])
        nᵢ < 1e-12 && continue
        for j in (i+1):S
            nⱼ = norm(samples[j])
            nⱼ < 1e-12 && continue
            cos_ij = dot(samples[i], samples[j]) / (nᵢ * nⱼ)
            total += (1.0 - cos_ij)
            count += 1
        end
    end
    return count > 0 ? total / count : 0.0
end

# ──────────────────────────────────────────────────────────────────────────────
# Multiplicity-weighted conditioning functions
# (ported from SA-Binding-Generation-Study)
# ──────────────────────────────────────────────────────────────────────────────

"""
    multiplicity_vector(K, subset_indices; ρ=10.0) -> Vector{Float64}

Construct a multiplicity weight vector for K patterns where designated
patterns get multiplicity ρ and all others get multiplicity 1.

The multiplicity ratio ρ = r_designated / r_background controls conditioning:
  - ρ = 1: standard SA (uniform weighting)
  - ρ → ∞: hard curation (only designated patterns contribute)
  - intermediate ρ: soft conditioning with tunable strength

For unit-normalized memories, the designated finite-mixture component
probability is:
  f_eff = K_d · ρ / (K_d · ρ + K_bg)
"""
function multiplicity_vector(K::Int, subset_indices::Vector{Int}; ρ::Float64=10.0)
    ρ > 0 || throw(ArgumentError("Multiplicity ratio ρ must be positive"))
    r = ones(K)
    for idx in subset_indices
        r[idx] = ρ
    end
    return r
end

"""
    effective_subset_fraction(r, subset_indices) -> Float64

Fraction of finite-mixture component probability assigned to designated
patterns for unit-normalized memories:
  f_eff = Σ_{k ∈ subset} r_k / Σ_k r_k
"""
function effective_subset_fraction(r::Vector{Float64}, subset_indices::Vector{Int})
    total = sum(r)
    subset_total = sum(r[idx] for idx in subset_indices)
    return subset_total / total
end

"""
    effective_num_patterns(r) -> Float64

Participation ratio of the multiplicity vector:
  K_eff(r) = (Σ_k r_k)² / Σ_k r_k²

This is the inverse concentration of the component probabilities q_k = r_k/Σ_j r_j
across repeated draws. When all r_k are equal, K_eff = K. When one r_k dominates,
K_eff → 1. It is not a count of independent patients.
"""
function effective_num_patterns(r::Vector{Float64})
    return sum(r)^2 / sum(r .^ 2)
end

"""
    ρ_for_target_fraction(K_sub, K_bg, f_target) -> Float64

Compute the multiplicity ratio ρ needed to assign a target finite-mixture
component-probability fraction f_target to K_sub designated patterns, with
K_bg background patterns, when memories are unit normalized.

From f_eff = K_sub · ρ / (K_sub · ρ + K_bg), solving for ρ:
  ρ = f_target · K_bg / (K_sub · (1 - f_target))
"""
function ρ_for_target_fraction(K_sub::Int, K_bg::Int, f_target::Float64)
    0 < f_target < 1 || throw(ArgumentError("f_target must be in (0, 1)"))
    return f_target * K_bg / (K_sub * (1 - f_target))
end

"""
    weighted_hopfield_energy(ξ, X, β, r) -> Float64

Multiplicity-weighted Hopfield energy:
  E_r(ξ) = ½‖ξ‖² - (1/β) log Σ_k r_k exp(β m_k^T ξ)
"""
function weighted_hopfield_energy(ξ::Vector{Float64}, X::Matrix{Float64},
                                   β::Float64, r::Vector{Float64})::Float64
    logits = β .* (X' * ξ) .+ log.(r)
    logits_max = maximum(logits)
    lse = logits_max + log(sum(exp.(logits .- logits_max)))
    return 0.5 * dot(ξ, ξ) - lse / β
end

"""
    weighted_attention_entropy(ξ, X, β, r) -> Float64

Shannon entropy of the multiplicity-weighted attention distribution:
  p_k = r_k exp(β m_k^T ξ) / Σ_j r_j exp(β m_j^T ξ)
"""
function weighted_attention_entropy(ξ::Vector{Float64}, X::Matrix{Float64},
                                     β::Float64, r::Vector{Float64})::Float64
    logits = β .* (X' * ξ) .+ log.(r)
    p = NNlib.softmax(logits)
    H = 0.0
    for pk in p
        pk > 1e-30 && (H -= pk * log(pk))
    end
    return H
end

"""
    find_weighted_entropy_inflection(X̂, r; α, n_betas, β_range) -> NamedTuple

Select the weighted attention-entropy operating point β*(r) at the greatest
downward curvature with respect to log β. The historical function name is
retained for compatibility. Uses weighted attention entropy.
"""
function find_weighted_entropy_inflection(X̂::Matrix{Float64}, r::Vector{Float64};
                                           α::Float64=0.01, n_betas::Int=50,
                                           β_range::Tuple{Float64,Float64}=(0.1, 500.0))
    d, K = size(X̂)
    βs = 10 .^ range(log10(β_range[1]), log10(β_range[2]), length=n_betas)

    n_probes = min(K, 20)
    Hs = zeros(n_betas)
    for (bi, β) in enumerate(βs)
        H_sum = 0.0
        for k in 1:n_probes
            H_sum += weighted_attention_entropy(X̂[:, k], X̂, β, r)
        end
        Hs[bi] = H_sum / n_probes
    end

    log_βs = log.(βs)
    dH = diff(Hs) ./ diff(log_βs)
    d2H = diff(dH) ./ diff(log_βs[1:end-1])

    transition_idx = 1
    min_d2H = Inf
    for i in eachindex(d2H)
        if d2H[i] < min_d2H
            min_d2H = d2H[i]
            transition_idx = i + 1
        end
    end

    β_star = βs[transition_idx]
    K_eff = effective_num_patterns(r)

    @info "  Weighted attention-entropy transition (d=$d, K=$K, K_eff=$(round(K_eff, digits=1))):"
    @info "    β*(ρ) = $(round(β_star, digits=2))"

    return (β_star=β_star, K_eff=K_eff, βs=βs, Hs=Hs)
end

"""
    fisher_separation_index(X̂, subset_indices) -> Float64

Compute the Fisher separation index between subset and background patterns
in PCA space. Measures how well the two groups are separated:
  S = ‖μ_sub - μ_bg‖² / (σ²_sub + σ²_bg)

High S → multiplicity weighting will transfer the phenotype.
Low S → PCA intermixes the groups; consider hard curation instead.
"""
function fisher_separation_index(X̂::Matrix{Float64}, subset_indices::Vector{Int})
    d, K = size(X̂)
    bg_indices = setdiff(1:K, subset_indices)

    μ_sub = vec(mean(X̂[:, subset_indices], dims=2))
    μ_bg = vec(mean(X̂[:, bg_indices], dims=2))

    # within-group variance (mean squared distance from group centroid)
    σ²_sub = mean(sum((X̂[:, k] .- μ_sub).^2) for k in subset_indices)
    σ²_bg = mean(sum((X̂[:, k] .- μ_bg).^2) for k in bg_indices)

    denom = σ²_sub + σ²_bg
    denom < 1e-12 && return 0.0
    return sum((μ_sub .- μ_bg).^2) / denom
end

"""
    find_entropy_inflection(X̂; α, n_betas, β_range) -> NamedTuple

Select the attention-entropy operating point β* for the memory matrix X̂ at
the greatest downward curvature with respect to log β. The historical function
name is retained for compatibility.
Returns β*, SNR*, theoretical predictions, and the full β/H curves.
"""
function find_entropy_inflection(X̂::Matrix{Float64};
                                  α::Float64=0.01,
                                  n_betas::Int=50,
                                  β_range::Tuple{Float64,Float64}=(0.1, 500.0))

    d, K = size(X̂)
    βs = 10 .^ range(log10(β_range[1]), log10(β_range[2]), length=n_betas)

    # compute mean entropy at each β using a few random probes
    n_probes = min(K, 20)
    Hs = zeros(n_betas)
    for (bi, β) in enumerate(βs)
        H_sum = 0.0
        for k in 1:n_probes
            H_sum += attention_entropy(X̂[:, k], X̂, β)
        end
        Hs[bi] = H_sum / n_probes
    end

    # normalize: H/log(K) so curves are in [0,1]
    H_max = log(K)
    Hs_norm = Hs ./ H_max

    # numerical second derivative in log-β space
    log_βs = log.(βs)
    dH = diff(Hs_norm) ./ diff(log_βs)
    d2H = diff(dH) ./ diff(log_βs[1:end-1])

    # operating point: where d2H is most negative (greatest downward curvature)
    transition_idx = 1
    min_d2H = Inf
    for i in 1:length(d2H)
        if d2H[i] < min_d2H
            min_d2H = d2H[i]
            transition_idx = i + 1
        end
    end

    β_star = βs[transition_idx]
    snr_star = sqrt(α * β_star / (2 * d))

    # theoretical prediction for random unit-norm patterns
    β_star_theory = sqrt(d)
    snr_star_theory = sqrt(α / (2 * sqrt(d)))

    @info "  Attention-entropy transition analysis (d=$d, K=$K):"
    @info "    Greatest downward curvature: β* = $(round(β_star, digits=2)),  SNR* = $(round(snr_star, digits=4))"
    @info "    Theoretical (√d):      β* = $(round(β_star_theory, digits=2)),  SNR* = $(round(snr_star_theory, digits=4))"

    return (β_star=β_star, snr_star=snr_star,
            β_star_theory=β_star_theory, snr_star_theory=snr_star_theory,
            βs=βs, Hs=Hs, Hs_norm=Hs_norm)
end

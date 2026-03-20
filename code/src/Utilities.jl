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

"""
    find_entropy_inflection(X̂; α, n_betas, β_range) -> NamedTuple

Compute the entropy inflection point β* for the memory matrix X̂.
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

    # numerical second derivative in log-β space to find inflection
    log_βs = log.(βs)
    dH = diff(Hs_norm) ./ diff(log_βs)
    d2H = diff(dH) ./ diff(log_βs[1:end-1])

    # inflection: where d2H is most negative
    inflection_idx = 1
    min_d2H = Inf
    for i in 1:length(d2H)
        if d2H[i] < min_d2H
            min_d2H = d2H[i]
            inflection_idx = i + 1
        end
    end

    β_star = βs[inflection_idx]
    snr_star = sqrt(α * β_star / (2 * d))

    # theoretical prediction for random unit-norm patterns
    β_star_theory = sqrt(d)
    snr_star_theory = sqrt(α / (2 * sqrt(d)))

    @info "  Phase transition analysis (d=$d, K=$K):"
    @info "    Empirical inflection:  β* = $(round(β_star, digits=2)),  SNR* = $(round(snr_star, digits=4))"
    @info "    Theoretical (√d):      β* = $(round(β_star_theory, digits=2)),  SNR* = $(round(snr_star_theory, digits=4))"

    return (β_star=β_star, snr_star=snr_star,
            β_star_theory=β_star_theory, snr_star_theory=snr_star_theory,
            βs=βs, Hs=Hs, Hs_norm=Hs_norm)
end

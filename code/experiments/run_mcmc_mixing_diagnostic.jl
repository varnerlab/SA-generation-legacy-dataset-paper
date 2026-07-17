# ──────────────────────────────────────────────────────────────────────────────
# run_mcmc_mixing_diagnostic.jl
#
# Task S2 — MCMC mixing diagnostic for the ULA-based Stochastic Attention (SA)
# sampler. The published generation pipeline (run_full_longitudinal.jl) draws
# each synthetic patient from a single anchor-initialized ULA chain run for
# T=2000 steps. This script checks two things about that design choice:
#
#   (1) Does the anchor-initialized ULA chain actually reach the Hopfield
#       energy stationary distribution by the time it stops, or is it still
#       burning in? We estimate a burn-in B and thinning τ via split-R̂ /
#       integrated autocorrelation time on M=20 sphere-initialized chains run
#       out to T=5000, well beyond the published T=2000.
#   (2) Is the sampler's output init-dependent at this horizon? We compare
#       endpoint energies and endpoint nearest-memory distances between
#       M sphere-init chains and M anchor-init chains (same T) via a
#       two-sample KS test.
#
# We additionally run M MALA chains (Metropolis-corrected) for context on
# discretization bias via acceptance rate.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HypothesisTests

# ══════════════════════════════════════════════════════════════════════════════
# Step 0: Load memory matrix + β*
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 0: Loading memory matrix …"
memory_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
JLD2.@load memory_path X̂
d, K = size(X̂)
@info "  Memory: d=$d × K=$K patients"

@info "Step 0b: Finding β* …"
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  β* = $(round(β_star, digits=3))"

# ══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════════════════════

"""
    integrated_acf(e; max_lag=500, threshold=0.05) -> (τ_int, acf)

Estimate the integrated autocorrelation time from an energy trace using a
windowed estimator (stops when ACF drops below threshold).

Copied verbatim from experiments/run_sampling_diagnostics.jl:126 (that script
is a standalone top-to-bottom run and is not `include`-able here).
"""
function integrated_acf(e::Vector{Float64}; max_lag::Int=500, threshold::Float64=0.05)
    n = length(e)
    μ = mean(e)
    var_e = var(e)
    var_e < 1e-30 && return (τ_int=1.0, acf=Float64[1.0])

    acf_vals = Float64[]
    τ = 0.5  # initial contribution from lag 0
    for lag in 1:min(max_lag, n-1)
        c = mean((e[1:n-lag] .- μ) .* (e[lag+1:n] .- μ)) / var_e
        push!(acf_vals, c)
        if c < threshold
            break
        end
        τ += c
    end
    return (τ_int=2.0 * τ, acf=acf_vals)
end

"""
    split_rhat(traces) -> R̂

Split-R̂ (Gelman-Rubin) diagnostic on a scalar quantity recorded across `M`
chains. `traces` is `M × L` (post-warmup window of length `L`); each chain is
split into two halves, giving `2M` half-chains used to compare within-chain
vs. between-chain variance. R̂ ≈ 1 indicates the chains agree on the target
distribution; R̂ ≫ 1 indicates they have not mixed. (Formula given verbatim
in the task-S2 brief.)
"""
function split_rhat(traces::Matrix{Float64})
    M, L = size(traces); Lh = L ÷ 2
    chains = vcat(traces[:,1:Lh], traces[:,Lh+1:2Lh])   # 2M half-chains of length Lh
    m = size(chains,1)
    means = vec(mean(chains, dims=2)); vars = vec(var(chains, dims=2))
    B = Lh * var(means); W = mean(vars)
    V = ((Lh-1)/Lh)*W + B/Lh
    sqrt(V / W)
end

# ══════════════════════════════════════════════════════════════════════════════
# Self-test (TDD): verify split_rhat / integrated_acf behave correctly on
# synthetic data BEFORE trusting them on the real chains below. A broken
# implementation of either helper should fail loudly here, not silently
# produce a plausible-looking but wrong diagnostic downstream.
# ══════════════════════════════════════════════════════════════════════════════
@info "Self-test: verifying split_rhat / integrated_acf on synthetic chains …"
let
    Random.seed!(999)

    # Well-mixed synthetic chains (same stationary N(5,1)) → R̂ ≈ 1
    mixed = 5.0 .+ randn(10, 2000)
    r_mixed = split_rhat(mixed)
    @assert r_mixed < 1.05 "Self-test FAILED: well-mixed synthetic chains should give R̂ < 1.05 (got $r_mixed)"

    # Poorly-mixed synthetic chains (different, non-overlapping means per chain) → R̂ ≫ 1
    poorly_mixed = Matrix{Float64}(undef, 10, 2000)
    for c in 1:10
        poorly_mixed[c, :] = fill(5.0 * c, 2000) .+ 0.01 .* randn(2000)
    end
    r_bad = split_rhat(poorly_mixed)
    @assert r_bad > 1.05 "Self-test FAILED: poorly-mixed synthetic chains should give R̂ > 1.05 (got $r_bad)"

    # White noise has no autocorrelation → τ_int ≈ 1
    white = randn(5000)
    τ_white = integrated_acf(white).τ_int
    @assert 0.5 < τ_white < 3.0 "Self-test FAILED: white-noise τ_int should be ≈ 1 (got $τ_white)"

    @info "  Self-test passed (R̂_mixed=$(round(r_mixed,digits=3)), R̂_poorly_mixed=$(round(r_bad,digits=1)), τ_white=$(round(τ_white,digits=2)))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: M sphere-init ULA chains → split-R̂, burn-in B, τ, ESS
# ══════════════════════════════════════════════════════════════════════════════
M = 20
T_chain = 5000
α_step = 0.01

@info "\nStep 1: Running M=$M sphere-init ULA chains (T=$T_chain, α=$α_step) …"
sphere_seed_base = 1000
sphere_energies = Matrix{Float64}(undef, M, T_chain + 1)
sphere_endpoints = Matrix{Float64}(undef, d, M)

for c in 1:M
    seed = sphere_seed_base + c
    Random.seed!(seed)
    ξ₀ = randn(d)
    ξ₀ ./= norm(ξ₀)
    res = sample(X̂, ξ₀, T_chain; β=β_star, α=α_step, seed=seed)
    for t in 0:T_chain
        sphere_energies[c, t+1] = hopfield_energy(res.Ξ[t+1, :], X̂, β_star)
    end
    sphere_endpoints[:, c] = res.Ξ[end, :]
end
@info "  Done. Energy range at t=0: [$(round(minimum(sphere_energies[:,1]),digits=2)), $(round(maximum(sphere_energies[:,1]),digits=2))]" *
      "  at t=T: [$(round(minimum(sphere_energies[:,end]),digits=2)), $(round(maximum(sphere_energies[:,end]),digits=2))]"

# split-R̂ over growing windows [t0:T_chain]
window_starts = collect(0:50:(T_chain - 1000))
rhat_vals = Float64[split_rhat(sphere_energies[:, (t0+1):end]) for t0 in window_starts]

B_idx = findfirst(r -> r < 1.05, rhat_vals)
B_fallback_used = false
local B
if B_idx === nothing
    B = 1000
    B_fallback_used = true
    @warn "  R̂ never dropped below 1.05 over tested windows (min R̂ = $(round(minimum(rhat_vals),digits=3))); falling back to B=$B"
else
    B = window_starts[B_idx]
end
@assert B < T_chain "Burn-in B=$B must be less than T_chain=$T_chain"

post_B_energies = sphere_energies[:, (B+1):end]
τ_per_chain = [integrated_acf(post_B_energies[c, :]).τ_int for c in 1:M]
τ = ceil(mean(τ_per_chain))
@assert τ >= 1 "τ=$τ must be ≥ 1"

ESS = (T_chain - B) / τ

@info "  B (burn-in)   = $B"
@info "  τ (thinning)  = $τ  (per-chain mean integrated ACF; fallback used = $B_fallback_used)"
@info "  ESS           = $(round(ESS, digits=1))  (post-burn-in samples = $(T_chain - B))"

# Sanity check on the shape of the R̂-vs-window curve: for a genuinely mixing
# chain we expect R̂ to fall (or stay flat) as t0 grows and more of the early
# transient is discarded. If R̂ instead *rises* with t0, the low R̂ found near
# t0=0 is likely an artifact of transient-inflated within-chain variance
# (a huge early swing masks real between-chain differences), not evidence of
# genuine convergence — flag this loudly rather than silently trusting B.
if rhat_vals[end] > rhat_vals[1]
    @warn "  R̂ INCREASES with window start t0 (from $(round(rhat_vals[1],digits=3)) at t0=0 " *
          "to $(round(rhat_vals[end],digits=3)) at t0=$(window_starts[end])). This is atypical: " *
          "it suggests persistent/growing between-chain heterogeneity (chains settling into " *
          "different memory basins) rather than genuine convergence, and that the B=$B found " *
          "by the threshold rule is likely an artifact of transient-inflated within-chain " *
          "variance at small t0, not true mixing. See run_mcmc_mixing_diagnostic report."
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: M anchor-init ULA chains (same T) + representative anchor-distance trace
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 2: Running M=$M anchor-init ULA chains (T=$T_chain, α=$α_step) …"
anchor_seed_base = 2000
anchor_energies = Matrix{Float64}(undef, M, T_chain + 1)
anchor_endpoints = Matrix{Float64}(undef, d, M)

representative_Ξ = Matrix{Float64}(undef, T_chain + 1, d)
representative_k0 = 0

for c in 1:M
    seed = anchor_seed_base + c
    Random.seed!(seed)
    k0 = mod1(c, K)
    ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d)
    ξ₀ ./= norm(ξ₀)
    res = sample(X̂, ξ₀, T_chain; β=β_star, α=α_step, seed=seed)
    for t in 0:T_chain
        anchor_energies[c, t+1] = hopfield_energy(res.Ξ[t+1, :], X̂, β_star)
    end
    anchor_endpoints[:, c] = res.Ξ[end, :]
    if c == 1
        representative_Ξ .= res.Ξ
        global representative_k0 = k0
    end
end

# anchor-distance trace for the representative (first) anchor-init chain:
#   dist_nearest[t]    = min_k ‖ξ_t − X̂[:,k]‖   (distance to nearest memory)
#   dist_own_anchor[t] = ‖ξ_t − X̂[:,k0]‖         (distance to its own initial anchor)
n_frames = T_chain + 1
dist_nearest = Vector{Float64}(undef, n_frames)
dist_own_anchor = Vector{Float64}(undef, n_frames)
for t in 1:n_frames
    ξt = @view representative_Ξ[t, :]
    dist_nearest[t] = minimum(norm(ξt .- X̂[:, k]) for k in 1:K)
    dist_own_anchor[t] = norm(ξt .- X̂[:, representative_k0])
end
@info "  Representative chain: k0=$representative_k0, dist_nearest(0)=$(round(dist_nearest[1],digits=4)), dist_nearest(T)=$(round(dist_nearest[end],digits=4))"
@info "                        dist_own_anchor(0)=$(round(dist_own_anchor[1],digits=4)), dist_own_anchor(T)=$(round(dist_own_anchor[end],digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Init-dependence — KS test on endpoint distributions
#         (sphere-init from Step 1 vs. anchor-init from Step 2, same T_chain)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 3: Init-dependence check (sphere-init vs anchor-init endpoints, T=$T_chain) …"

sphere_endpoint_energy = sphere_energies[:, end]
anchor_endpoint_energy = anchor_energies[:, end]

nearest_dist(ξ) = minimum(norm(ξ .- X̂[:, k]) for k in 1:K)
sphere_endpoint_dist = [nearest_dist(sphere_endpoints[:, c]) for c in 1:M]
anchor_endpoint_dist = [nearest_dist(anchor_endpoints[:, c]) for c in 1:M]

ks_energy = ApproximateTwoSampleKSTest(sphere_endpoint_energy, anchor_endpoint_energy)
ks_energy_D = ks_energy.δ
ks_energy_p = pvalue(ks_energy)

ks_dist = ApproximateTwoSampleKSTest(sphere_endpoint_dist, anchor_endpoint_dist)
ks_dist_D = ks_dist.δ
ks_dist_p = pvalue(ks_dist)

@assert 0.0 <= ks_energy_D <= 1.0 "KS D (energy) out of [0,1]: $ks_energy_D"
@assert 0.0 <= ks_dist_D <= 1.0 "KS D (nearest-memory distance) out of [0,1]: $ks_dist_D"

init_dependent = (ks_energy_p < 0.05) || (ks_dist_p < 0.05)
verdict = init_dependent ? "NOT MIXED (init-dependent at T=$T_chain)" : "MIXED (init-independent at T=$T_chain)"

@info "  KS (endpoint energy):            D=$(round(ks_energy_D,digits=4)), p=$(round(ks_energy_p,digits=4))"
@info "  KS (endpoint nearest-mem. dist.): D=$(round(ks_dist_D,digits=4)), p=$(round(ks_dist_p,digits=4))"
@info "  Verdict: $verdict"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: MALA acceptance rate (discretization-bias context)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Running M=$M MALA chains (T=$T_chain, α=$α_step) for acceptance-rate context …"
mala_seed_base = 3000
mala_accept_rates = Vector{Float64}(undef, M)

for c in 1:M
    seed = mala_seed_base + c
    Random.seed!(seed)
    k0 = mod1(c, K)
    ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d)
    ξ₀ ./= norm(ξ₀)
    res = mala_sample(X̂, ξ₀, T_chain; β=β_star, α=α_step, seed=seed)
    mala_accept_rates[c] = res.accept_rate
end

@assert all(0.0 .<= mala_accept_rates .<= 1.0) "MALA acceptance rates out of [0,1]"
mean_mala_accept = mean(mala_accept_rates)
std_mala_accept = std(mala_accept_rates)
@info "  MALA acceptance rate: $(round(mean_mala_accept,digits=4)) ± $(round(std_mala_accept,digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Final self-checks (re-asserted together for a single loud failure point)
# ══════════════════════════════════════════════════════════════════════════════
@assert B < T_chain
@assert τ >= 1
@assert 0.0 <= mean_mala_accept <= 1.0
@assert 0.0 <= ks_energy_D <= 1.0 && 0.0 <= ks_dist_D <= 1.0
@info "\nStep 5: All self-checks passed."

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Save results (CSV)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 6: Saving results …"

df_summary = DataFrame(
    row_type = ["summary"],
    M = [M], T = [T_chain], alpha = [α_step], beta_star = [β_star],
    B = [B], B_fallback_used = [B_fallback_used],
    tau = [τ], ESS = [ESS],
    ks_energy_D = [ks_energy_D], ks_energy_p = [ks_energy_p],
    ks_dist_D = [ks_dist_D], ks_dist_p = [ks_dist_p],
    init_dependent = [init_dependent],
    mala_accept_mean = [mean_mala_accept], mala_accept_std = [std_mala_accept],
    representative_k0 = [representative_k0],
    anchor_dist_nearest_t0 = [dist_nearest[1]], anchor_dist_nearest_T = [dist_nearest[end]],
    anchor_dist_own_t0 = [dist_own_anchor[1]], anchor_dist_own_T = [dist_own_anchor[end]],
)

df_rhat = DataFrame(
    row_type = fill("rhat_window", length(window_starts)),
    t0 = window_starts,
    rhat = rhat_vals,
)

df_anchor = DataFrame(
    row_type = fill("anchor_distance", n_frames),
    t = 0:T_chain,
    dist_nearest = dist_nearest,
    dist_own_anchor = dist_own_anchor,
)

df_out = vcat(df_summary, df_rhat, df_anchor; cols=:union)
CSV.write(joinpath(_PATH_TO_DATA, "mcmc_mixing_diagnostic.csv"), df_out)
@info "  Saved data/mcmc_mixing_diagnostic.csv ($(nrow(df_out)) rows: 1 summary + $(length(window_starts)) rhat_window + $n_frames anchor_distance)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Figure (energy traces + burn-in marker, R̂ curve, anchor-distance curve)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 7: Plotting …"

bg_color = RGB(0.97, 0.97, 0.98)
sphere_color = RGB(0.10, 0.45, 0.70)   # steelblue-ish
anchor_color = RGB(0.85, 0.35, 0.10)   # coral-ish
common_attrs = (background_color_inside=bg_color, grid=false, framestyle=:box,
                guidefontsize=11, tickfontsize=9, legendfontsize=8, margin=6Plots.mm)

# Panel 1: energy traces (subset of chains for legibility) + burn-in marker
p_energy = plot(; xlabel="Iteration t", ylabel="Hopfield energy E(ξ_t)",
                 legend=:topright, common_attrs...)
n_show = min(M, 8)
for c in 1:n_show
    plot!(p_energy, 0:T_chain, sphere_energies[c, :],
          color=sphere_color, alpha=0.5, lw=1,
          label = c == 1 ? "Sphere-init" : nothing)
    plot!(p_energy, 0:T_chain, anchor_energies[c, :],
          color=anchor_color, alpha=0.5, lw=1,
          label = c == 1 ? "Anchor-init" : nothing)
end
vline!(p_energy, [B], color=:gray30, ls=:dash, lw=1.5, label="Burn-in B=$B")

# Panel 2: split-R̂ vs window start
p_rhat = plot(window_starts, rhat_vals; xlabel="Window start t0", ylabel="Split R̂",
              color=sphere_color, lw=2, label="R̂(t0:T)", legend=:topright, common_attrs...)
hline!(p_rhat, [1.05], color=:gray30, ls=:dash, lw=1.5, label="R̂ = 1.05")
vline!(p_rhat, [B], color=anchor_color, ls=:dot, lw=1.5, label="B=$B")

# Panel 3: anchor-distance trace for representative chain
p_dist = plot(0:T_chain, dist_nearest; xlabel="Iteration t", ylabel="Distance",
              color=sphere_color, lw=1.5, label="min_k ‖ξ_t − X̂[:,k]‖",
              legend=:topright, common_attrs...)
plot!(p_dist, 0:T_chain, dist_own_anchor; color=anchor_color, lw=1.5,
      label="‖ξ_t − X̂[:,k₀]‖  (own anchor, k₀=$representative_k0)")

p_all = plot(p_energy, p_rhat, p_dist; layout=(1, 3), size=(1800, 480), left_margin=10Plots.mm)
savefig(p_all, joinpath(_PATH_TO_FIG, "mcmc_mixing_diagnostic.pdf"))
savefig(p_all, joinpath(_PATH_TO_FIG, "mcmc_mixing_diagnostic.png"))
@info "  Saved figs/mcmc_mixing_diagnostic.pdf"

# ══════════════════════════════════════════════════════════════════════════════
# Final report
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^78)
println("MCMC Mixing Diagnostic — Summary")
println("="^78)
println("  M = $M chains,  T = $T_chain,  α = $α_step,  β* = $(round(β_star,digits=3))")
println("-"^78)
println("  B (burn-in)  = $B" * (B_fallback_used ? "  [fallback: R̂ never < 1.05]" : ""))
println("  τ (thinning) = $τ")
println("  ESS          = $(round(ESS,digits=1))  (post-burn-in samples = $(T_chain - B))")
println("-"^78)
println("  Init-dependence KS test (sphere-init vs anchor-init, T=$T_chain):")
println("    Endpoint energy:            D=$(round(ks_energy_D,digits=4)), p=$(round(ks_energy_p,digits=4))")
println("    Endpoint nearest-mem. dist.: D=$(round(ks_dist_D,digits=4)), p=$(round(ks_dist_p,digits=4))")
println("    Verdict: $verdict")
println("-"^78)
println("  MALA acceptance rate: $(round(mean_mala_accept,digits=4)) ± $(round(std_mala_accept,digits=4))")
println("="^78)

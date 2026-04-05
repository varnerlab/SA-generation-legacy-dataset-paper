# ──────────────────────────────────────────────────────────────────────────────
# run_sampling_diagnostics.jl
#
# ULA vs MALA control experiment: compares the Unadjusted Langevin Algorithm
# with the Metropolis-Adjusted Langevin Algorithm on the pregnancy coagulation
# dataset. Computes energy traces, autocorrelation, effective sample size,
# and MALA acceptance rates to confirm that ULA discretization bias is
# negligible at α = 0.01.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load data and build memory matrix (same as run_full_longitudinal.jl)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Loading dataset and building memory matrix …"
xlsx_path = joinpath(_PATH_TO_ORIGINAL_DATA, "20200604_R61_Legacy_Biochemical_Measurements.xlsx")
df_raw = load_legacy_data(xlsx_path)
df_clean, kept_cols = clean_assay_data(df_raw; max_missing_frac=0.3)
n_assays = length(kept_cols)
n_visits = 3

subjects = unique(df_clean.SubjectID)
complete_subjects = Int[]
for s in subjects
    visits = sort(unique(df_clean[df_clean.SubjectID .== s, :Visit]))
    if visits == [1, 2, 3]
        push!(complete_subjects, s)
    end
end
K = length(complete_subjects)
d_concat = n_assays * n_visits

X_concat = Matrix{Float64}(undef, K, d_concat)
for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_clean.SubjectID .== s) .& (df_clean.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_concat[i, offset + j] = df_clean[row, col]
        end
    end
end

# Standardize + PCA
μ_concat = vec(mean(X_concat, dims=1))
σ_concat = vec(std(X_concat, dims=1))
for i in eachindex(σ_concat)
    σ_concat[i] < 1e-12 && (σ_concat[i] = 1.0)
end
Z_concat = (X_concat .- μ_concat') ./ σ_concat'
Z_t = Z_concat'
pca_concat = MultivariateStats.fit(PCA, Z_t; pratio=0.95)
d_pca = MultivariateStats.outdim(pca_concat)
X_pca = MultivariateStats.transform(pca_concat, Z_t)

X̂ = copy(X_pca)
for k in 1:K
    X̂[:, k] ./= (norm(X̂[:, k]) + 1e-12)
end
@info "  Memory matrix: $d_pca × $K (PCA dimensions × patients)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Find β*
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 2: Finding β* …"
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  β* = $β_star"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Run diagnostic chains (ULA and MALA)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 3: Running diagnostic chains …"

n_chains = 10
T_chain = 5000
α_step = 0.01
T_burn = 1000  # burn-in to discard

# storage
ula_energies = Matrix{Float64}(undef, n_chains, T_chain + 1)
mala_energies = Matrix{Float64}(undef, n_chains, T_chain + 1)
mala_accept_rates = Vector{Float64}(undef, n_chains)

# energy difference statistics
ula_final_states = Matrix{Float64}(undef, n_chains, d_pca)
mala_final_states = Matrix{Float64}(undef, n_chains, d_pca)

for c in 1:n_chains
    # same initialization for both samplers
    k0 = mod1(c, K)
    ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d_pca)
    ξ₀ ./= norm(ξ₀)

    # ULA
    res_ula = sample(X̂, copy(ξ₀), T_chain; β=β_star, α=α_step, seed=100 + c)
    for t in 0:T_chain
        ula_energies[c, t+1] = hopfield_energy(res_ula.Ξ[t+1, :], X̂, β_star)
    end
    ula_final_states[c, :] = res_ula.Ξ[end, :]

    # MALA
    res_mala = mala_sample(X̂, copy(ξ₀), T_chain; β=β_star, α=α_step, seed=100 + c)
    for t in 0:T_chain
        mala_energies[c, t+1] = hopfield_energy(res_mala.Ξ[t+1, :], X̂, β_star)
    end
    mala_final_states[c, :] = res_mala.Ξ[end, :]
    mala_accept_rates[c] = res_mala.accept_rate

    @info "  Chain $c: MALA acceptance rate = $(round(res_mala.accept_rate, digits=4))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Compute autocorrelation and effective sample size
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Computing autocorrelation diagnostics …"

"""
    integrated_acf(e; max_lag=500, threshold=0.05) -> (τ_int, acf_values)

Estimate the integrated autocorrelation time from an energy trace using
a windowed estimator (stops when ACF drops below threshold).
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

# post burn-in traces
ula_τ = Float64[]
mala_τ = Float64[]
ula_ess = Float64[]
mala_ess = Float64[]

for c in 1:n_chains
    e_ula = ula_energies[c, (T_burn+1):end]
    e_mala = mala_energies[c, (T_burn+1):end]
    n_post = length(e_ula)

    res_ula = integrated_acf(e_ula)
    res_mala = integrated_acf(e_mala)

    push!(ula_τ, res_ula.τ_int)
    push!(mala_τ, res_mala.τ_int)
    push!(ula_ess, n_post / res_ula.τ_int)
    push!(mala_ess, n_post / res_mala.τ_int)
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Summary statistics
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 5: Summary …"

# post burn-in mean energies
ula_mean_E = [mean(ula_energies[c, (T_burn+1):end]) for c in 1:n_chains]
mala_mean_E = [mean(mala_energies[c, (T_burn+1):end]) for c in 1:n_chains]

println("\n" * "="^70)
println("ULA vs MALA Sampling Diagnostics")
println("="^70)
println("  β* = $(round(β_star, digits=2)),  α = $α_step,  T = $T_chain,  burn-in = $T_burn")
println("  Chains: $n_chains,  Memory: $d_pca × $K")
println("-"^70)
println("MALA acceptance rate:  $(round(mean(mala_accept_rates), digits=4)) ± $(round(std(mala_accept_rates), digits=4))")
println("-"^70)
println("Post-burn-in mean energy:")
println("  ULA:  $(round(mean(ula_mean_E), digits=4)) ± $(round(std(ula_mean_E), digits=4))")
println("  MALA: $(round(mean(mala_mean_E), digits=4)) ± $(round(std(mala_mean_E), digits=4))")
println("  Δ(ULA − MALA): $(round(mean(ula_mean_E) - mean(mala_mean_E), digits=6))")
println("-"^70)
println("Integrated autocorrelation time (τ_int):")
println("  ULA:  $(round(mean(ula_τ), digits=2)) ± $(round(std(ula_τ), digits=2))")
println("  MALA: $(round(mean(mala_τ), digits=2)) ± $(round(std(mala_τ), digits=2))")
println("-"^70)
println("Effective sample size (post-burn-in):")
println("  ULA:  $(round(mean(ula_ess), digits=1)) ± $(round(std(ula_ess), digits=1))")
println("  MALA: $(round(mean(mala_ess), digits=1)) ± $(round(std(mala_ess), digits=1))")
println("="^70)

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Energy trace plot (ULA vs MALA, single representative chain)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 6: Plotting energy traces …"

p_energy = plot(0:T_chain, ula_energies[1, :],
    label="ULA", color=:steelblue, alpha=0.8, lw=1,
    xlabel="Iteration", ylabel="Hopfield Energy",
    title="ULA vs MALA Energy Traces (Chain 1)")
plot!(p_energy, 0:T_chain, mala_energies[1, :],
    label="MALA", color=:coral, alpha=0.8, lw=1)
vline!(p_energy, [T_burn], label="Burn-in", color=:gray, ls=:dash, lw=1)
savefig(p_energy, joinpath(_PATH_TO_FIG, "ula_vs_mala_energy_trace.pdf"))
savefig(p_energy, joinpath(_PATH_TO_FIG, "ula_vs_mala_energy_trace.png"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Save results to CSV
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 7: Saving results …"
df_diag = DataFrame(
    chain = 1:n_chains,
    mala_accept_rate = mala_accept_rates,
    ula_mean_energy = ula_mean_E,
    mala_mean_energy = mala_mean_E,
    ula_tau_int = ula_τ,
    mala_tau_int = mala_τ,
    ula_ess = ula_ess,
    mala_ess = mala_ess,
)
CSV.write(joinpath(_PATH_TO_DATA, "ula_vs_mala_diagnostics.csv"), df_diag)

@info "\n✓ Done! Results saved to data/ula_vs_mala_diagnostics.csv"
@info "  Figure saved to figs/ula_vs_mala_energy_trace.pdf"

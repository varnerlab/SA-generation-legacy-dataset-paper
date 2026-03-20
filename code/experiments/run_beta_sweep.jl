# ──────────────────────────────────────────────────────────────────────────────
# run_beta_sweep.jl
#
# Sweep β from well below β* to above β* and measure how dispersion,
# novelty, and feature fidelity change. Goal: find the operating point
# that produces representative synthetic patients with full dispersion.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load saved pipeline
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading saved pipeline …"
JLD2.@load joinpath(_PATH_TO_DATA, "patient_memory.jld2") X̂ pca_model std_params df_clean kept_cols

d_pca, K = size(X̂)
@info "  Memory: d=$d_pca, K=$K"

# compute real data PCA projections for reference
Z_real, _ = standardize(df_clean, kept_cols)
X_real_pca = MultivariateStats.transform(pca_model, Z_real')  # d_pca × K

# real data dispersion metrics (target to match)
real_pc1_std = std(X_real_pca[1, :])
real_pc2_std = std(X_real_pca[2, :])
real_pc1_range = maximum(X_real_pca[1, :]) - minimum(X_real_pca[1, :])
real_pc2_range = maximum(X_real_pca[2, :]) - minimum(X_real_pca[2, :])
@info "  Real data PC1: std=$(round(real_pc1_std, digits=3)), range=$(round(real_pc1_range, digits=3))"
@info "  Real data PC2: std=$(round(real_pc2_std, digits=3)), range=$(round(real_pc2_range, digits=3))"

# real data: compute dispersion as mean distance from centroid in PCA space
real_centroid = vec(mean(X_real_pca, dims=2))
real_dists = [norm(X_real_pca[:, k] - real_centroid) for k in 1:K]
real_mean_dist = mean(real_dists)
real_std_dist = std(real_dists)
@info "  Real mean distance from centroid: $(round(real_mean_dist, digits=3)) ± $(round(real_std_dist, digits=3))"

# ══════════════════════════════════════════════════════════════════════════════
# β sweep
# ══════════════════════════════════════════════════════════════════════════════
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  β* = $(round(β_star, digits=2))"

# sweep from 0.1×β* to 3×β*
β_fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0]
β_values = β_fractions .* β_star

n_synth = 200
T_steps = 2000
α_step = 0.05

results = DataFrame(
    β_frac = Float64[],
    β = Float64[],
    noise_scale = Float64[],
    pc1_std = Float64[],
    pc2_std = Float64[],
    pc1_std_ratio = Float64[],
    pc2_std_ratio = Float64[],
    mean_dist_from_centroid = Float64[],
    dist_ratio = Float64[],
    mean_novelty = Float64[],
    diversity = Float64[],
    mean_rel_error = Float64[],
    median_rel_error = Float64[],
    correlation_mae = Float64[]
)

synth_pca_collections = Dict{Float64, Matrix{Float64}}()

for (i, β) in enumerate(β_values)
    frac = β_fractions[i]
    noise = sqrt(2.0 * α_step / β)
    @info "  β = $(round(β, digits=2)) ($(frac)×β*), noise_scale = $(round(noise, digits=4))"

    # generate
    df_synth = generate_patients(X̂, pca_model, std_params, n_synth, T_steps;
                                  β=β, α=α_step, seed=42)

    # project to PCA
    Z_synth = Matrix{Float64}(df_synth[!, kept_cols])
    Z_synth_std = (Z_synth .- std_params.μ') ./ std_params.σ'
    X_synth_pca = MultivariateStats.transform(pca_model, Z_synth_std')

    synth_pca_collections[frac] = X_synth_pca

    # dispersion
    pc1_s = std(X_synth_pca[1, :])
    pc2_s = std(X_synth_pca[2, :])
    synth_centroid = vec(mean(X_synth_pca, dims=2))
    synth_dists = [norm(X_synth_pca[:, k] - synth_centroid) for k in 1:size(X_synth_pca, 2)]
    mean_dist = mean(synth_dists)

    # novelty & diversity
    samples_vec = [X_synth_pca[:, j] for j in 1:size(X_synth_pca, 2)]
    novelties = [sample_novelty(X_synth_pca[:, j], X̂) for j in 1:size(X_synth_pca, 2)]
    div = sample_diversity(samples_vec)

    # feature fidelity
    comparison = feature_summary_comparison(df_clean, df_synth, kept_cols)
    mean_re = mean(comparison.Mean_Rel_Error)
    median_re = median(comparison.Mean_Rel_Error)

    # correlation
    cor_real = column_correlation_matrix(df_clean, kept_cols)
    cor_synth = column_correlation_matrix(df_synth, kept_cols)
    cor_mae = mean(abs.(cor_real .- cor_synth))

    push!(results, (frac, β, noise, pc1_s, pc2_s,
                    pc1_s / real_pc1_std, pc2_s / real_pc2_std,
                    mean_dist, mean_dist / real_mean_dist,
                    mean(novelties), div, mean_re, median_re, cor_mae))

    @info "    PC1 std ratio: $(round(pc1_s / real_pc1_std, digits=3)), " *
          "dist ratio: $(round(mean_dist / real_mean_dist, digits=3)), " *
          "novelty: $(round(mean(novelties), digits=3)), " *
          "mean_rel_err: $(round(mean_re, digits=4))"
end

# save results
CSV.write(joinpath(_PATH_TO_DATA, "beta_sweep_results.csv"), results)
@info "\nβ sweep results:"
show(stdout, results[:, [:β_frac, :β, :noise_scale, :pc1_std_ratio, :dist_ratio,
                          :mean_novelty, :mean_rel_error, :correlation_mae]])

# ══════════════════════════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════════════════════════

# --- Figure 1: Dispersion ratio vs β ---
p1 = plot(results.β_frac, results.dist_ratio,
    xlabel="β / β*", ylabel="Dispersion ratio (synthetic / real)",
    title="Dispersion vs Inverse Temperature",
    lw=2, marker=:circle, markersize=5, label="Mean dist from centroid",
    legend=:topright)
hline!(p1, [1.0], ls=:dash, color=:red, lw=2, label="Target (real = 1.0)")
vline!(p1, [1.0], ls=:dot, color=:gray, lw=1, label="β = β*")
savefig(p1, joinpath(_PATH_TO_FIG, "beta_sweep_dispersion.pdf"))
savefig(p1, joinpath(_PATH_TO_FIG, "beta_sweep_dispersion.png"))

# --- Figure 2: Multi-metric tradeoff ---
p2 = plot(results.β_frac, results.dist_ratio,
    xlabel="β / β*", ylabel="Metric value",
    title="β Sweep: Dispersion–Fidelity Tradeoff",
    lw=2, marker=:circle, label="Dispersion ratio", legend=:right)
plot!(p2, results.β_frac, results.mean_novelty,
    lw=2, marker=:square, label="Novelty")
plot!(p2, results.β_frac, 1.0 .- results.mean_rel_error,
    lw=2, marker=:diamond, label="1 - Mean Rel Error")
plot!(p2, results.β_frac, 1.0 .- results.correlation_mae,
    lw=2, marker=:utriangle, label="1 - Corr MAE")
hline!(p2, [1.0], ls=:dash, color=:red, lw=1, label="")
vline!(p2, [1.0], ls=:dot, color=:gray, lw=1, label="β*")
savefig(p2, joinpath(_PATH_TO_FIG, "beta_sweep_tradeoff.pdf"))
savefig(p2, joinpath(_PATH_TO_FIG, "beta_sweep_tradeoff.png"))

# --- Figure 3: PCA panels at select β values ---
select_fracs = [0.2, 0.5, 0.8, 1.0]
panels = []
for frac in select_fracs
    X_s = synth_pca_collections[frac]
    p = scatter(X_real_pca[1, :], X_real_pca[2, :],
        label="Real", color=:steelblue, alpha=0.6, markersize=4,
        xlabel="PC1", ylabel="PC2",
        title="β = $(round(frac, digits=1))×β*",
        titlefontsize=10)
    scatter!(p, X_s[1, :], X_s[2, :],
        label="Synthetic", color=:coral, alpha=0.4, markersize=3,
        markershape=:diamond)
    push!(panels, p)
end
p3 = plot(panels..., layout=(2, 2), size=(800, 700),
          plot_title="PCA Dispersion at Different β")
savefig(p3, joinpath(_PATH_TO_FIG, "beta_sweep_pca_panels.pdf"))
savefig(p3, joinpath(_PATH_TO_FIG, "beta_sweep_pca_panels.png"))

# --- Figure 4: PC1 std ratio vs β (shows where we match real spread) ---
p4 = plot(results.β_frac, results.pc1_std_ratio,
    xlabel="β / β*", ylabel="σ ratio (synthetic / real)",
    title="PC Standard Deviation Ratio vs β",
    lw=2, marker=:circle, markersize=5, label="PC1 σ ratio")
plot!(p4, results.β_frac, results.pc2_std_ratio,
    lw=2, marker=:square, markersize=5, label="PC2 σ ratio")
hline!(p4, [1.0], ls=:dash, color=:red, lw=2, label="Target (1.0)")
vline!(p4, [1.0], ls=:dot, color=:gray, lw=1, label="β*")
savefig(p4, joinpath(_PATH_TO_FIG, "beta_sweep_pc_std_ratio.pdf"))
savefig(p4, joinpath(_PATH_TO_FIG, "beta_sweep_pc_std_ratio.png"))

@info "\n✓ β sweep complete. Figures saved to figs/"
@info "  Look for the β/β* where dispersion ratio ≈ 1.0 — that's the generative sweet spot."

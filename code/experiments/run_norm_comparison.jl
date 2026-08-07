# ──────────────────────────────────────────────────────────────────────────────
# run_norm_comparison.jl
#
# Compare unit-norm vs raw PCA memory matrices for patient generation.
# The hypothesis: unit-norm projection destroys the anisotropic variance
# structure, causing synthetic patients to cluster in the center of PCA space.
# Raw PCA projections preserve the natural scale of the data.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load data
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading data …"
JLD2.@load joinpath(_PATH_TO_DATA, "patient_memory.jld2") df_clean kept_cols

Z, std_params = standardize(df_clean, kept_cols)

# ══════════════════════════════════════════════════════════════════════════════
# Build both memory matrices
# ══════════════════════════════════════════════════════════════════════════════
@info "\n── Unit-norm memory matrix ──"
X̂_norm, pca_norm, _ = build_memory_matrix(Z; pratio=0.95)
d, K = size(X̂_norm)

@info "\n── Raw PCA memory matrix ──"
X̂_raw, pca_raw, _, norms = build_memory_matrix_raw(Z; pratio=0.95)

# real data reference in PCA space
X_real_pca = MultivariateStats.transform(pca_raw, Z')
real_pc1_std = std(X_real_pca[1, :])
real_pc2_std = std(X_real_pca[2, :])
real_centroid = vec(mean(X_real_pca, dims=2))
real_dists = [norm(X_real_pca[:, k] - real_centroid) for k in 1:K]
real_mean_dist = mean(real_dists)

# ══════════════════════════════════════════════════════════════════════════════
# Find β* for raw memory matrix
# ══════════════════════════════════════════════════════════════════════════════
@info "\n── Attention-entropy transition (raw) ──"
phase_raw = find_entropy_inflection(X̂_raw; α=0.01, n_betas=80,
                                     β_range=(0.001, 100.0))
β_star_raw = phase_raw.β_star

@info "\n── Attention-entropy transition (unit-norm) ──"
phase_norm = find_entropy_inflection(X̂_norm; α=0.01, n_betas=80,
                                      β_range=(0.1, 1000.0))
β_star_norm = phase_norm.β_star

# ══════════════════════════════════════════════════════════════════════════════
# Generate with both methods, sweeping β for raw
# ══════════════════════════════════════════════════════════════════════════════
n_synth = 200
T_steps = 2000
α_step = 0.05

# β sweep for raw memory matrix
β_fractions_raw = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0]

@info "\n── Generating synthetic patients (raw memory, β sweep) ──"
results_raw = DataFrame(
    β_frac = Float64[], β = Float64[],
    pc1_std = Float64[], pc2_std = Float64[],
    pc1_ratio = Float64[], pc2_ratio = Float64[],
    dist_ratio = Float64[], mean_novelty = Float64[],
    mean_rel_error = Float64[], correlation_mae = Float64[]
)

best_raw_pca = nothing
best_raw_df = nothing
best_dist_ratio = Inf

for frac in β_fractions_raw
    global best_raw_pca, best_raw_df, best_dist_ratio
    β = frac * β_star_raw
    @info "  β = $(round(β, digits=4)) ($(frac)×β*_raw)"

    # generate using raw memory
    df_synth = generate_patients(X̂_raw, pca_raw, std_params, n_synth, T_steps;
                                  β=β, α=α_step, seed=42)

    # project to PCA space
    Z_synth = Matrix{Float64}(df_synth[!, kept_cols])
    Z_synth_std = (Z_synth .- std_params.μ') ./ std_params.σ'
    X_synth_pca = MultivariateStats.transform(pca_raw, Z_synth_std')

    pc1_s = std(X_synth_pca[1, :])
    pc2_s = std(X_synth_pca[2, :])
    synth_centroid = vec(mean(X_synth_pca, dims=2))
    synth_dists = [norm(X_synth_pca[:, j] - synth_centroid) for j in 1:size(X_synth_pca, 2)]
    mean_dist = mean(synth_dists)
    dist_ratio = mean_dist / real_mean_dist

    novelties = [sample_novelty(X_synth_pca[:, j], X̂_raw) for j in 1:size(X_synth_pca, 2)]
    comparison = feature_summary_comparison(df_clean, df_synth, kept_cols)
    cor_real = column_correlation_matrix(df_clean, kept_cols)
    cor_synth = column_correlation_matrix(df_synth, kept_cols)
    cor_mae = mean(abs.(cor_real .- cor_synth))

    push!(results_raw, (frac, β, pc1_s, pc2_s,
                        pc1_s/real_pc1_std, pc2_s/real_pc2_std,
                        dist_ratio, mean(novelties),
                        mean(comparison.Mean_Rel_Error), cor_mae))

    @info "    PC1 σ ratio: $(round(pc1_s/real_pc1_std, digits=3)), " *
          "PC2 σ ratio: $(round(pc2_s/real_pc2_std, digits=3)), " *
          "dist ratio: $(round(dist_ratio, digits=3)), " *
          "rel_err: $(round(mean(comparison.Mean_Rel_Error), digits=4))"

    # track best (closest to dist_ratio = 1.0)
    if abs(dist_ratio - 1.0) < abs(best_dist_ratio - 1.0)
        best_dist_ratio = dist_ratio
        best_raw_pca = X_synth_pca
        best_raw_df = df_synth
    end
end

CSV.write(joinpath(_PATH_TO_DATA, "raw_beta_sweep_results.csv"), results_raw)

# unit-norm baseline at β*
@info "\n── Generating synthetic patients (unit-norm, β = β*) ──"
df_synth_norm = generate_patients(X̂_norm, pca_norm, std_params, n_synth, T_steps;
                                   β=β_star_norm, α=α_step, seed=42)
Z_synth_n = Matrix{Float64}(df_synth_norm[!, kept_cols])
Z_synth_n_std = (Z_synth_n .- std_params.μ') ./ std_params.σ'
X_synth_norm_pca = MultivariateStats.transform(pca_norm, Z_synth_n_std')
norm_pc1_ratio = std(X_synth_norm_pca[1, :]) / real_pc1_std
norm_pc2_ratio = std(X_synth_norm_pca[2, :]) / real_pc2_std
@info "  Unit-norm: PC1 σ ratio = $(round(norm_pc1_ratio, digits=3)), PC2 σ ratio = $(round(norm_pc2_ratio, digits=3))"

# ══════════════════════════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════════════════════════
@info "\n── Generating figures ──"

# --- Figure 1: Side-by-side PCA comparison (unit-norm vs raw at best β) ---
p1 = scatter(X_real_pca[1, :], X_real_pca[2, :],
    label="Real", color=:steelblue, alpha=0.7, markersize=5,
    xlabel="PC1", ylabel="PC2",
    title="Unit-Norm (β = β*)")
scatter!(p1, X_synth_norm_pca[1, :], X_synth_norm_pca[2, :],
    label="Synthetic", color=:coral, alpha=0.4, markersize=4,
    markershape=:diamond)

p2 = scatter(X_real_pca[1, :], X_real_pca[2, :],
    label="Real", color=:steelblue, alpha=0.7, markersize=5,
    xlabel="PC1", ylabel="PC2",
    title="Raw PCA (β = $(round(best_dist_ratio, digits=2))×β*)")
scatter!(p2, best_raw_pca[1, :], best_raw_pca[2, :],
    label="Synthetic", color=:coral, alpha=0.4, markersize=4,
    markershape=:diamond)

p_compare = plot(p1, p2, layout=(1, 2), size=(1000, 450),
                 plot_title="Unit-Norm vs Raw PCA Memory")
savefig(p_compare, joinpath(_PATH_TO_FIG, "norm_comparison_pca.pdf"))
savefig(p_compare, joinpath(_PATH_TO_FIG, "norm_comparison_pca.png"))

# --- Figure 2: PC std ratio vs β for raw ---
p3 = plot(results_raw.β_frac, results_raw.pc1_ratio,
    xlabel="β / β*", ylabel="σ ratio (synthetic / real)",
    title="Raw PCA: Dispersion vs β",
    lw=2, marker=:circle, markersize=5, label="PC1 σ ratio")
plot!(p3, results_raw.β_frac, results_raw.pc2_ratio,
    lw=2, marker=:square, markersize=5, label="PC2 σ ratio")
plot!(p3, results_raw.β_frac, results_raw.dist_ratio,
    lw=2, marker=:diamond, markersize=5, label="Overall dist ratio")
hline!(p3, [1.0], ls=:dash, color=:red, lw=2, label="Target (1.0)")
vline!(p3, [1.0], ls=:dot, color=:gray, lw=1, label="β*")
savefig(p3, joinpath(_PATH_TO_FIG, "raw_dispersion_vs_beta.pdf"))
savefig(p3, joinpath(_PATH_TO_FIG, "raw_dispersion_vs_beta.png"))

# --- Figure 3: PCA panels at select β (raw) ---
select_fracs = [1.0, 5.0, 10.0, 30.0]
panels = []
for frac in select_fracs
    β = frac * β_star_raw
    df_s = generate_patients(X̂_raw, pca_raw, std_params, n_synth, T_steps;
                              β=β, α=α_step, seed=42)
    Z_s = Matrix{Float64}(df_s[!, kept_cols])
    Z_s_std = (Z_s .- std_params.μ') ./ std_params.σ'
    X_s_pca = MultivariateStats.transform(pca_raw, Z_s_std')

    p = scatter(X_real_pca[1, :], X_real_pca[2, :],
        label="Real", color=:steelblue, alpha=0.6, markersize=4,
        xlabel="PC1", ylabel="PC2",
        title="Raw, β = $(round(frac, digits=1))×β*",
        titlefontsize=10)
    scatter!(p, X_s_pca[1, :], X_s_pca[2, :],
        label="Synthetic", color=:coral, alpha=0.4, markersize=3,
        markershape=:diamond)
    push!(panels, p)
end
p4 = plot(panels..., layout=(2, 2), size=(800, 700),
          plot_title="Raw PCA Dispersion at Different β")
savefig(p4, joinpath(_PATH_TO_FIG, "raw_pca_panels.pdf"))
savefig(p4, joinpath(_PATH_TO_FIG, "raw_pca_panels.png"))

# --- Figure 4: Entropy curves comparison ---
p5 = plot(phase_norm.βs, phase_norm.Hs_norm,
    xlabel="β", ylabel="H / log(K)",
    title="Entropy Curves",
    xscale=:log10, lw=2, label="Unit-norm", color=:steelblue)
vline!(p5, [β_star_norm], ls=:dash, color=:steelblue, lw=1,
       label="β*_norm = $(round(β_star_norm, digits=2))")
plot!(p5, phase_raw.βs, phase_raw.Hs_norm,
    lw=2, label="Raw PCA", color=:coral)
vline!(p5, [β_star_raw], ls=:dash, color=:coral, lw=1,
       label="β*_raw = $(round(β_star_raw, digits=3))")
savefig(p5, joinpath(_PATH_TO_FIG, "entropy_comparison.pdf"))
savefig(p5, joinpath(_PATH_TO_FIG, "entropy_comparison.png"))

# --- Figure 5: Feature distributions at best raw β ---
n_feat_show = min(6, length(kept_cols))
feat_indices = round.(Int, range(1, length(kept_cols), length=n_feat_show))
p_dists = []
for idx in feat_indices
    col = kept_cols[idx]
    real_vals = collect(skipmissing(df_clean[!, col]))
    synth_vals = best_raw_df[!, col]

    p = histogram(real_vals, normalize=:pdf, alpha=0.5, label="Real",
                  color=:steelblue, bins=20, title=string(col), titlefontsize=9)
    histogram!(p, synth_vals, normalize=:pdf, alpha=0.5, label="Synthetic",
               color=:coral, bins=20)
    push!(p_dists, p)
end
p6 = plot(p_dists..., layout=(2, 3), size=(900, 500),
          plot_title="Feature Distributions (Raw PCA, best β)")
savefig(p6, joinpath(_PATH_TO_FIG, "raw_feature_distributions.pdf"))
savefig(p6, joinpath(_PATH_TO_FIG, "raw_feature_distributions.png"))

@info "\n✓ Done! Key results:"
@info "  Unit-norm β* = $(round(β_star_norm, digits=2)), PC1 σ ratio = $(round(norm_pc1_ratio, digits=3))"
@info "  Raw β* = $(round(β_star_raw, digits=4)), best dist ratio = $(round(best_dist_ratio, digits=3))"
@info "  Figures saved to figs/"

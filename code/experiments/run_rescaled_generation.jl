# ──────────────────────────────────────────────────────────────────────────────
# LEGACY / EXPLORATORY — DO NOT USE FOR PAPER ARTIFACTS.
# This script regenerates synthetic data with non-canonical hyperparameters
# (different α, seed, or N from run_full_longitudinal.jl). Its CSV outputs are
# only consumed by itself and are not used by the paper or supplement.
# Canonical generation: code/experiments/run_full_longitudinal.jl
# ──────────────────────────────────────────────────────────────────────────────
# run_rescaled_generation.jl
#
# Compare three SA generation strategies:
#   A) Unit-norm SA (original — good fidelity, poor dispersion)
#   B) Raw PCA SA (good dispersion, poor fidelity)
#   C) Rescaled SA (unit-norm direction × empirical magnitude)
#
# The rescaled approach decomposes the generation into:
#   - Direction: sampled by SA on the unit sphere (well-behaved Hopfield energy)
#   - Magnitude: drawn from the empirical PCA norm distribution
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load data and build pipelines
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading data …"
JLD2.@load joinpath(_PATH_TO_DATA, "patient_memory.jld2") df_clean kept_cols

Z, std_params = standardize(df_clean, kept_cols)

# unit-norm memory
X̂_unit, pca_model, _ = build_memory_matrix(Z; pratio=0.95)
d, K = size(X̂_unit)

# compute PCA norms (before unit-normalization)
Z_t = Z'
X_pca_raw = MultivariateStats.transform(pca_model, Z_t)  # d × K
pca_norms = [norm(X_pca_raw[:, k]) for k in 1:K]
@info "  PCA norms: range=[$(round(minimum(pca_norms),digits=2)), $(round(maximum(pca_norms),digits=2))], mean=$(round(mean(pca_norms),digits=2)), std=$(round(std(pca_norms),digits=2))"

# real reference
X_real_pca = X_pca_raw  # d × K (raw PCA projections)
real_pc1_std = std(X_real_pca[1, :])
real_pc2_std = std(X_real_pca[2, :])
real_centroid = vec(mean(X_real_pca, dims=2))
real_dists = [norm(X_real_pca[:, k] - real_centroid) for k in 1:K]
real_mean_dist = mean(real_dists)

# find β*
phase = find_entropy_inflection(X̂_unit; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  β* = $(round(β_star, digits=2))"

# ══════════════════════════════════════════════════════════════════════════════
# Generate with all three methods
# ══════════════════════════════════════════════════════════════════════════════
n_synth = 200
T_steps = 2000
α_step = 0.05

# --- Method A: Unit-norm at β* ---
@info "\n── Method A: Unit-norm SA ──"
df_A = generate_patients(X̂_unit, pca_model, std_params, n_synth, T_steps;
                          β=β_star, α=α_step, seed=42)
Z_A = Matrix{Float64}(df_A[!, kept_cols])
Z_A_std = (Z_A .- std_params.μ') ./ std_params.σ'
X_A_pca = MultivariateStats.transform(pca_model, Z_A_std')

# --- Method C: Rescaled SA at β* ---
@info "\n── Method C: Rescaled SA (direction × magnitude) ──"
df_C = generate_patients_rescaled(X̂_unit, pca_model, std_params, pca_norms,
                                   n_synth, T_steps;
                                   β=β_star, α=α_step, seed=42)
Z_C = Matrix{Float64}(df_C[!, kept_cols])
Z_C_std = (Z_C .- std_params.μ') ./ std_params.σ'
X_C_pca = MultivariateStats.transform(pca_model, Z_C_std')

# also try rescaled at lower β (more diffuse directions)
@info "\n── Method C': Rescaled SA at 0.5×β* ──"
df_C2 = generate_patients_rescaled(X̂_unit, pca_model, std_params, pca_norms,
                                    n_synth, T_steps;
                                    β=0.5*β_star, α=α_step, seed=42)
Z_C2 = Matrix{Float64}(df_C2[!, kept_cols])
Z_C2_std = (Z_C2 .- std_params.μ') ./ std_params.σ'
X_C2_pca = MultivariateStats.transform(pca_model, Z_C2_std')

# ══════════════════════════════════════════════════════════════════════════════
# Evaluate all methods
# ══════════════════════════════════════════════════════════════════════════════
@info "\n── Evaluation ──"

function evaluate_method(name, X_synth_pca, df_synth, X_real_pca, df_clean, kept_cols, X̂_unit, real_mean_dist, real_pc1_std, real_pc2_std)
    pc1_s = std(X_synth_pca[1, :])
    pc2_s = std(X_synth_pca[2, :])
    centroid = vec(mean(X_synth_pca, dims=2))
    dists = [norm(X_synth_pca[:, j] - centroid) for j in 1:size(X_synth_pca, 2)]
    dist_ratio = mean(dists) / real_mean_dist

    novelties = [sample_novelty(X_synth_pca[:, j], X̂_unit) for j in 1:size(X_synth_pca, 2)]
    comparison = feature_summary_comparison(df_clean, df_synth, kept_cols)
    cor_real = column_correlation_matrix(df_clean, kept_cols)
    cor_synth = column_correlation_matrix(df_synth, kept_cols)
    cor_mae = mean(abs.(cor_real .- cor_synth))

    @info "  $name:"
    @info "    PC1 σ ratio: $(round(pc1_s/real_pc1_std, digits=3))"
    @info "    PC2 σ ratio: $(round(pc2_s/real_pc2_std, digits=3))"
    @info "    Dist ratio:  $(round(dist_ratio, digits=3))"
    @info "    Novelty:     $(round(mean(novelties), digits=3))"
    @info "    Mean RE:     $(round(mean(comparison.Mean_Rel_Error), digits=4))"
    @info "    Corr MAE:    $(round(cor_mae, digits=4))"

    return (name=name,
            pc1_ratio=pc1_s/real_pc1_std, pc2_ratio=pc2_s/real_pc2_std,
            dist_ratio=dist_ratio, novelty=mean(novelties),
            mean_re=mean(comparison.Mean_Rel_Error), cor_mae=cor_mae)
end

eval_A = evaluate_method("Unit-norm (β*)", X_A_pca, df_A, X_real_pca, df_clean, kept_cols, X̂_unit, real_mean_dist, real_pc1_std, real_pc2_std)
eval_C = evaluate_method("Rescaled (β*)", X_C_pca, df_C, X_real_pca, df_clean, kept_cols, X̂_unit, real_mean_dist, real_pc1_std, real_pc2_std)
eval_C2 = evaluate_method("Rescaled (0.5β*)", X_C2_pca, df_C2, X_real_pca, df_clean, kept_cols, X̂_unit, real_mean_dist, real_pc1_std, real_pc2_std)

# save the best synthetic dataset
CSV.write(joinpath(_PATH_TO_DATA, "synthetic_rescaled_all.csv"), df_C)

# ══════════════════════════════════════════════════════════════════════════════
# Conditioned rescaled generation
# ══════════════════════════════════════════════════════════════════════════════
@info "\n── Conditioned rescaled generation ──"
conditions_to_generate = ["Healthy", "Prior PE", "PCOS"]
n_per_condition = 100

for cond in conditions_to_generate
    mask = BitVector(df_clean.Condition .== cond)
    @info "  $cond (n_real=$(count(mask)))"

    df_cond = generate_conditioned_patients_rescaled(
        X̂_unit, pca_model, std_params, pca_norms, mask,
        n_per_condition, T_steps;
        β=β_star, α=α_step, ρ_ratio=10.0,
        seed=Int(hash(cond) % 2^31))

    fname = replace(lowercase(cond), " " => "_")
    CSV.write(joinpath(_PATH_TO_DATA, "synthetic_rescaled_$(fname).csv"), df_cond)

    comparison = feature_summary_comparison(df_clean[mask, :], df_cond, kept_cols)
    @info "    Mean RE: $(round(mean(comparison.Mean_Rel_Error), digits=4))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════════════════════════
@info "\n── Generating figures ──"

# --- Figure 1: Three-panel PCA comparison ---
p_A = scatter(X_real_pca[1, :], X_real_pca[2, :],
    label="Real", color=:steelblue, alpha=0.7, markersize=5,
    xlabel="PC1", ylabel="PC2", title="A: Unit-Norm", titlefontsize=10)
scatter!(p_A, X_A_pca[1, :], X_A_pca[2, :],
    label="Synthetic", color=:coral, alpha=0.4, markersize=3, markershape=:diamond)

p_C = scatter(X_real_pca[1, :], X_real_pca[2, :],
    label="Real", color=:steelblue, alpha=0.7, markersize=5,
    xlabel="PC1", ylabel="PC2", title="B: Rescaled (β*)", titlefontsize=10)
scatter!(p_C, X_C_pca[1, :], X_C_pca[2, :],
    label="Synthetic", color=:coral, alpha=0.4, markersize=3, markershape=:diamond)

p_C2 = scatter(X_real_pca[1, :], X_real_pca[2, :],
    label="Real", color=:steelblue, alpha=0.7, markersize=5,
    xlabel="PC1", ylabel="PC2", title="C: Rescaled (0.5β*)", titlefontsize=10)
scatter!(p_C2, X_C2_pca[1, :], X_C2_pca[2, :],
    label="Synthetic", color=:coral, alpha=0.4, markersize=3, markershape=:diamond)

p_main = plot(p_A, p_C, p_C2, layout=(1, 3), size=(1200, 400),
              plot_title="SA Generation: Unit-Norm vs Rescaled")
savefig(p_main, joinpath(_PATH_TO_FIG, "rescaled_comparison_pca.pdf"))
savefig(p_main, joinpath(_PATH_TO_FIG, "rescaled_comparison_pca.png"))

# --- Figure 2: Feature distributions (rescaled best) ---
n_feat_show = min(6, length(kept_cols))
feat_indices = round.(Int, range(1, length(kept_cols), length=n_feat_show))
p_dists = []
for idx in feat_indices
    col = kept_cols[idx]
    real_vals = collect(skipmissing(df_clean[!, col]))

    p = histogram(real_vals, normalize=:pdf, alpha=0.5, label="Real",
                  color=:steelblue, bins=20, title=string(col), titlefontsize=9)
    histogram!(p, df_C[!, col], normalize=:pdf, alpha=0.5, label="Rescaled",
               color=:coral, bins=20)
    push!(p_dists, p)
end
p_feat = plot(p_dists..., layout=(2, 3), size=(900, 500),
              plot_title="Feature Distributions: Real vs Rescaled SA")
savefig(p_feat, joinpath(_PATH_TO_FIG, "rescaled_feature_distributions.pdf"))
savefig(p_feat, joinpath(_PATH_TO_FIG, "rescaled_feature_distributions.png"))

# --- Figure 3: Correlation comparison ---
cor_real = column_correlation_matrix(df_clean, kept_cols)
cor_synth_C = column_correlation_matrix(df_C, kept_cols)
p_cr = heatmap(cor_real, title="Real", c=:RdBu, clims=(-1,1), aspect_ratio=1, size=(500,450))
p_cs = heatmap(cor_synth_C, title="Rescaled SA", c=:RdBu, clims=(-1,1), aspect_ratio=1, size=(500,450))
p_cor = plot(p_cr, p_cs, layout=(1,2), size=(1000,450))
savefig(p_cor, joinpath(_PATH_TO_FIG, "rescaled_correlation_comparison.pdf"))
savefig(p_cor, joinpath(_PATH_TO_FIG, "rescaled_correlation_comparison.png"))

# --- Figure 4: Norm distributions (real vs synthetic) ---
synth_norms = [norm(X_C_pca[:, j]) for j in 1:size(X_C_pca, 2)]
p_norms = histogram(pca_norms, normalize=:pdf, alpha=0.5, label="Real",
    color=:steelblue, bins=20, xlabel="PCA norm ‖z‖", ylabel="Density",
    title="PCA Norm Distribution")
histogram!(p_norms, synth_norms, normalize=:pdf, alpha=0.5, label="Synthetic (rescaled)",
    color=:coral, bins=20)
savefig(p_norms, joinpath(_PATH_TO_FIG, "rescaled_norm_distributions.pdf"))
savefig(p_norms, joinpath(_PATH_TO_FIG, "rescaled_norm_distributions.png"))

# --- Figure 5: Conditioned rescaled in PCA space ---
cond_colors = Dict("Healthy" => :steelblue, "Prior PE" => :coral, "PCOS" => :seagreen)
p_cond = scatter(X_real_pca[1, :], X_real_pca[2, :],
    label="Real (all)", color=:lightgray, alpha=0.3, markersize=3,
    xlabel="PC1", ylabel="PC2",
    title="Conditioned Rescaled SA in PCA Space")
for cond in conditions_to_generate
    fname = replace(lowercase(cond), " " => "_")
    df_cs = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_rescaled_$(fname).csv"), DataFrame)
    Z_cs = Matrix{Float64}(df_cs[!, kept_cols])
    Z_cs_std = (Z_cs .- std_params.μ') ./ std_params.σ'
    X_cs_pca = MultivariateStats.transform(pca_model, Z_cs_std')
    color = get(cond_colors, cond, :gray)
    scatter!(p_cond, X_cs_pca[1, :], X_cs_pca[2, :],
        label="Synth: $cond", color=color, alpha=0.5, markersize=4, markershape=:diamond)
end
savefig(p_cond, joinpath(_PATH_TO_FIG, "rescaled_conditioned_pca.pdf"))
savefig(p_cond, joinpath(_PATH_TO_FIG, "rescaled_conditioned_pca.png"))

# --- Figure 6: β sweep for rescaled method ---
@info "\n── β sweep for rescaled method ──"
β_fractions = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
sweep_results = DataFrame(
    β_frac=Float64[], pc1_ratio=Float64[], pc2_ratio=Float64[],
    dist_ratio=Float64[], mean_re=Float64[], cor_mae=Float64[]
)

for frac in β_fractions
    β = frac * β_star
    df_s = generate_patients_rescaled(X̂_unit, pca_model, std_params, pca_norms,
                                       n_synth, T_steps; β=β, α=α_step, seed=42)
    Z_s = Matrix{Float64}(df_s[!, kept_cols])
    Z_s_std = (Z_s .- std_params.μ') ./ std_params.σ'
    X_s_pca = MultivariateStats.transform(pca_model, Z_s_std')

    pc1_r = std(X_s_pca[1, :]) / real_pc1_std
    pc2_r = std(X_s_pca[2, :]) / real_pc2_std
    cent = vec(mean(X_s_pca, dims=2))
    dr = mean([norm(X_s_pca[:, j] - cent) for j in 1:size(X_s_pca, 2)]) / real_mean_dist
    comp = feature_summary_comparison(df_clean, df_s, kept_cols)
    cr = column_correlation_matrix(df_clean, kept_cols)
    cs = column_correlation_matrix(df_s, kept_cols)
    push!(sweep_results, (frac, pc1_r, pc2_r, dr, mean(comp.Mean_Rel_Error), mean(abs.(cr .- cs))))
end

p_sweep = plot(sweep_results.β_frac, sweep_results.pc1_ratio,
    xlabel="β / β*", ylabel="Metric",
    title="Rescaled SA: β Sweep",
    lw=2, marker=:circle, label="PC1 σ ratio")
plot!(p_sweep, sweep_results.β_frac, sweep_results.dist_ratio,
    lw=2, marker=:diamond, label="Dist ratio")
plot!(p_sweep, sweep_results.β_frac, 1.0 .- sweep_results.mean_re,
    lw=2, marker=:square, label="1 - Mean RE")
hline!(p_sweep, [1.0], ls=:dash, color=:red, lw=1, label="Target")
savefig(p_sweep, joinpath(_PATH_TO_FIG, "rescaled_beta_sweep.pdf"))
savefig(p_sweep, joinpath(_PATH_TO_FIG, "rescaled_beta_sweep.png"))

CSV.write(joinpath(_PATH_TO_DATA, "rescaled_beta_sweep.csv"), sweep_results)

@info "\n✓ Done!"

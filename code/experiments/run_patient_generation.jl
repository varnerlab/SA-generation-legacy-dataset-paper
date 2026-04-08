# ──────────────────────────────────────────────────────────────────────────────
# LEGACY / EXPLORATORY — DO NOT USE FOR PAPER ARTIFACTS.
# This script regenerates synthetic data with non-canonical hyperparameters
# (different α, seed, or N from run_full_longitudinal.jl). Its CSV outputs are
# only consumed by itself and are not used by the paper or supplement.
# Canonical generation: code/experiments/run_full_longitudinal.jl
# ──────────────────────────────────────────────────────────────────────────────
# run_patient_generation.jl
#
# Main experiment: generate synthetic patient records from the R61 Legacy
# Biochemical Measurements dataset using Stochastic Attention (SA).
#
# Steps:
#   1. Load and clean the dataset
#   2. Standardize + PCA → memory matrix
#   3. Find β* via entropy inflection
#   4. Generate unconditioned synthetic patients
#   5. Generate condition-specific patients (Healthy, Prior PE, PCOS)
#   6. Evaluate: feature-level statistics, correlation preservation, novelty
#   7. Save results and produce diagnostic figures
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load and clean the dataset
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Loading legacy biochemical dataset …"
xlsx_path = joinpath(_PATH_TO_ORIGINAL_DATA, "20200604_R61_Legacy_Biochemical_Measurements.xlsx")
df_raw = load_legacy_data(xlsx_path)
@info "  Loaded $(nrow(df_raw)) records, $(ncol(df_raw)) columns"
@info "  Conditions: $(sort(collect(countmap(df_raw.Condition))))"

# clean: remove high-missing columns, then incomplete rows
df_clean, kept_cols = clean_assay_data(df_raw; max_missing_frac=0.3)
@info "  Clean dataset: $(nrow(df_clean)) records × $(length(kept_cols)) assays"
@info "  Conditions after cleaning: $(sort(collect(countmap(df_clean.Condition))))"

# save cleaned data
CSV.write(joinpath(_PATH_TO_DATA, "cleaned_legacy_data.csv"), df_clean)
@info "  Saved cleaned data to data/cleaned_legacy_data.csv"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Standardize + PCA → memory matrix
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 2: Building memory matrix …"
Z, std_params = standardize(df_clean, kept_cols)
@info "  Standardized: $(size(Z, 1)) records × $(size(Z, 2)) features"

X̂, pca_model, d_orig = build_memory_matrix(Z; pratio=0.95)
d_pca, K = size(X̂)
@info "  Memory: d=$d_pca, K=$K"

# save memory matrix and pipeline
JLD2.@save joinpath(_PATH_TO_DATA, "patient_memory.jld2") X̂ pca_model std_params df_clean kept_cols

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Find β* via entropy inflection
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 3: Finding β* (entropy inflection) …"
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80,
                                 β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  Using β* = $(round(β_star, digits=2))"

# plot entropy curve
p_entropy = plot(phase.βs, phase.Hs_norm,
    xlabel="β (inverse temperature)", ylabel="H / log(K)",
    title="Attention Entropy vs β  (d=$d_pca, K=$K)",
    xscale=:log10, legend=:topright, lw=2, label="H(β)")
vline!(p_entropy, [β_star], label="β* = $(round(β_star, digits=1))",
       ls=:dash, color=:red, lw=2)
vline!(p_entropy, [phase.β_star_theory], label="√d = $(round(phase.β_star_theory, digits=1))",
       ls=:dot, color=:blue, lw=2)
savefig(p_entropy, joinpath(_PATH_TO_FIG, "entropy_curve.pdf"))
savefig(p_entropy, joinpath(_PATH_TO_FIG, "entropy_curve.png"))
@info "  Saved entropy curve to figs/"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Generate unconditioned synthetic patients
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Generating unconditioned synthetic patients …"
n_synth = 200
T_steps = 2000

df_synth_all = generate_patients(X̂, pca_model, std_params, n_synth, T_steps;
                                  β=β_star, α=0.05, seed=42)
@info "  Generated $(nrow(df_synth_all)) synthetic records"

# evaluate feature-level statistics
comparison_all = feature_summary_comparison(df_clean, df_synth_all, kept_cols)
CSV.write(joinpath(_PATH_TO_DATA, "feature_comparison_all.csv"), comparison_all)
@info "  Mean relative error across features: $(round(mean(comparison_all.Mean_Rel_Error), digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Generate condition-specific patients
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 5: Generating condition-specific synthetic patients …"

conditions_to_generate = ["Healthy", "Prior PE", "PCOS"]
n_per_condition = 100

for cond in conditions_to_generate
    @info "  Condition: $cond"
    mask = BitVector(df_clean.Condition .== cond)
    n_real = count(mask)
    @info "    Real patients in this condition: $n_real / $K"

    df_cond = generate_conditioned_patients(X̂, pca_model, std_params, mask,
                                             n_per_condition, T_steps;
                                             β=β_star, α=0.05, ρ_ratio=10.0,
                                             seed=Int(hash(cond) % 2^31))

    # save
    fname = replace(lowercase(cond), " " => "_")
    CSV.write(joinpath(_PATH_TO_DATA, "synthetic_$(fname).csv"), df_cond)

    # evaluate
    comparison = feature_summary_comparison(df_clean[mask, :], df_cond, kept_cols)
    CSV.write(joinpath(_PATH_TO_DATA, "feature_comparison_$(fname).csv"), comparison)
    @info "    Mean relative error: $(round(mean(comparison.Mean_Rel_Error), digits=4))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Evaluate — novelty, diversity, correlation preservation
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 6: Evaluation diagnostics …"

# project synthetic samples back to PCA space for novelty/diversity
Z_synth = Matrix{Float64}(df_synth_all[!, kept_cols])
Z_synth_std = (Z_synth .- std_params.μ') ./ std_params.σ'
X_synth_pca = MultivariateStats.transform(pca_model, Z_synth_std')  # d_pca × n_synth

novelties = Float64[]
for i in 1:size(X_synth_pca, 2)
    ξ = X_synth_pca[:, i]
    push!(novelties, sample_novelty(ξ, X̂))
end
@info "  Novelty: mean=$(round(mean(novelties), digits=4)), std=$(round(std(novelties), digits=4))"
@info "  Novelty range: [$(round(minimum(novelties), digits=4)), $(round(maximum(novelties), digits=4))]"

# diversity
synth_samples = [X_synth_pca[:, i] for i in 1:size(X_synth_pca, 2)]
div = sample_diversity(synth_samples)
@info "  Diversity (mean pairwise cosine distance): $(round(div, digits=4))"

# correlation preservation
cor_real = column_correlation_matrix(df_clean, kept_cols)
cor_synth = column_correlation_matrix(df_synth_all, kept_cols)
cor_diff = abs.(cor_real .- cor_synth)
@info "  Correlation MAE: $(round(mean(cor_diff), digits=4))"
@info "  Correlation max absolute error: $(round(maximum(cor_diff), digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Diagnostic figures
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 7: Generating diagnostic figures …"

# --- Figure 1: PCA projection (real vs synthetic) ---
# project real data to first 2 PCA components
X_real_pca = MultivariateStats.transform(pca_model, Z')  # d_pca × K

p_pca = scatter(X_real_pca[1, :], X_real_pca[2, :],
    label="Real", color=:steelblue, alpha=0.7, markersize=5,
    xlabel="PC1", ylabel="PC2",
    title="PCA Space: Real vs Synthetic Patients")
scatter!(p_pca, X_synth_pca[1, :], X_synth_pca[2, :],
    label="Synthetic", color=:coral, alpha=0.5, markersize=4, markershape=:diamond)
savefig(p_pca, joinpath(_PATH_TO_FIG, "pca_real_vs_synthetic.pdf"))
savefig(p_pca, joinpath(_PATH_TO_FIG, "pca_real_vs_synthetic.png"))

# --- Figure 2: PCA by condition (real data colored by condition) ---
conditions_unique = unique(df_clean.Condition)
cond_colors = Dict("Healthy" => :steelblue, "Prior PE" => :coral,
                    "PCOS" => :seagreen, "Multi" => :purple)

p_cond = plot(xlabel="PC1", ylabel="PC2",
              title="PCA Space by Clinical Condition")
for cond in conditions_unique
    mask = df_clean.Condition .== cond
    color = get(cond_colors, cond, :gray)
    scatter!(p_cond, X_real_pca[1, mask], X_real_pca[2, mask],
        label=cond, color=color, alpha=0.7, markersize=5)
end
savefig(p_cond, joinpath(_PATH_TO_FIG, "pca_by_condition.pdf"))
savefig(p_cond, joinpath(_PATH_TO_FIG, "pca_by_condition.png"))

# --- Figure 3: Feature distribution comparison (select top features) ---
# pick 6 representative features
n_feat_show = min(6, length(kept_cols))
feat_indices = round.(Int, range(1, length(kept_cols), length=n_feat_show))
p_dists = []
for idx in feat_indices
    col = kept_cols[idx]
    real_vals = collect(skipmissing(df_clean[!, col]))
    synth_vals = df_synth_all[!, col]

    p = histogram(real_vals, normalize=:pdf, alpha=0.5, label="Real",
                  color=:steelblue, bins=20, title=string(col), titlefontsize=9)
    histogram!(p, synth_vals, normalize=:pdf, alpha=0.5, label="Synthetic",
               color=:coral, bins=20)
    push!(p_dists, p)
end
p_combined = plot(p_dists..., layout=(2, 3), size=(900, 500),
                  plot_title="Feature Distributions: Real vs Synthetic")
savefig(p_combined, joinpath(_PATH_TO_FIG, "feature_distributions.pdf"))
savefig(p_combined, joinpath(_PATH_TO_FIG, "feature_distributions.png"))

# --- Figure 4: Correlation heatmap comparison ---
p_cor_real = heatmap(cor_real, title="Correlation (Real)", c=:RdBu,
                     clims=(-1, 1), aspect_ratio=1, size=(500, 450))
p_cor_synth = heatmap(cor_synth, title="Correlation (Synthetic)", c=:RdBu,
                      clims=(-1, 1), aspect_ratio=1, size=(500, 450))
p_cor = plot(p_cor_real, p_cor_synth, layout=(1, 2), size=(1000, 450))
savefig(p_cor, joinpath(_PATH_TO_FIG, "correlation_comparison.pdf"))
savefig(p_cor, joinpath(_PATH_TO_FIG, "correlation_comparison.png"))

# --- Figure 5: Conditioned generation PCA overlay ---
p_condgen = plot(xlabel="PC1", ylabel="PC2",
                 title="Condition-Specific Synthetic Patients in PCA Space")

# plot real data lightly
scatter!(p_condgen, X_real_pca[1, :], X_real_pca[2, :],
    label="Real (all)", color=:lightgray, alpha=0.3, markersize=3)

for cond in conditions_to_generate
    fname = replace(lowercase(cond), " " => "_")
    df_csynth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_$(fname).csv"), DataFrame)

    # project to PCA space
    Z_c = Matrix{Float64}(df_csynth[!, kept_cols])
    Z_c_std = (Z_c .- std_params.μ') ./ std_params.σ'
    X_c_pca = MultivariateStats.transform(pca_model, Z_c_std')

    color = get(cond_colors, cond, :gray)
    scatter!(p_condgen, X_c_pca[1, :], X_c_pca[2, :],
        label="Synth: $cond", color=color, alpha=0.6, markersize=4,
        markershape=:diamond)
end
savefig(p_condgen, joinpath(_PATH_TO_FIG, "conditioned_generation_pca.pdf"))
savefig(p_condgen, joinpath(_PATH_TO_FIG, "conditioned_generation_pca.png"))

# --- Figure 6: Novelty histogram ---
p_novelty = histogram(novelties, bins=30, normalize=:pdf,
    xlabel="Novelty (1 - max cosine similarity)",
    ylabel="Density", title="Novelty of Synthetic Patients",
    color=:steelblue, alpha=0.7, label="")
vline!(p_novelty, [mean(novelties)], label="mean = $(round(mean(novelties), digits=3))",
       ls=:dash, color=:red, lw=2)
savefig(p_novelty, joinpath(_PATH_TO_FIG, "novelty_histogram.pdf"))
savefig(p_novelty, joinpath(_PATH_TO_FIG, "novelty_histogram.png"))

@info "\n✓ All done! Results saved to data/ and figs/"
@info "  Synthetic datasets:"
@info "    data/cleaned_legacy_data.csv"
@info "    data/feature_comparison_all.csv"
for cond in conditions_to_generate
    fname = replace(lowercase(cond), " " => "_")
    @info "    data/synthetic_$(fname).csv"
    @info "    data/feature_comparison_$(fname).csv"
end
@info "  Figures: figs/*.pdf"

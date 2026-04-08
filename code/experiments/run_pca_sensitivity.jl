# PCA sensitivity analysis: vary PCA variance threshold and measure SA quality
include(joinpath(@__DIR__, "..", "Include.jl"))

# ── Load and prepare data ────────────────────────────────────────────────────
@info "Loading data …"
df_real = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)

# Complete subjects
all_s = unique(df_real.SubjectID)
complete_subjects = [s for s in all_s
    if sum(df_real.SubjectID .== s) == 3 &&
       all(v -> any((df_real.SubjectID .== s) .& (df_real.Visit .== v)), 1:3)]

# Feature columns (same as synthetic CSV)
df_synth_ref = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
synth_meta = [:SyntheticID, :Visit]
kept_cols = [c for c in Symbol.(names(df_synth_ref)) if c ∉ synth_meta]
n_assays = length(kept_cols)
K = length(complete_subjects)
n_visits = 3

@info "  K=$K patients, $n_assays features/visit"

# Build concatenated DataFrame
concat_cols = Symbol[]
for v in 1:n_visits
    for col in kept_cols
        push!(concat_cols, Symbol("V$(v)_$(col)"))
    end
end

concat_data = Matrix{Float64}(undef, K, length(concat_cols))
for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_real.SubjectID .== s) .& (df_real.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            concat_data[i, offset + j] = Float64(df_real[row, col])
        end
    end
end

df_concat = DataFrame(concat_data, concat_cols)
@info "  Concatenated: $(length(concat_cols)) total features"

# ── Run SA at different PCA thresholds ───────────────────────────────────────
thresholds = [0.85, 0.90, 0.95, 0.975, 0.99]
n_synth = 100
T_steps = 2000
α_step = 0.01  # canonical (matches run_full_longitudinal.jl)

results_all = DataFrame(
    Threshold = Float64[],
    d_PCA = Int[],
    K_over_d = Float64[],
    β_star = Float64[],
    Novelty = Float64[],
    Median_MRE = Float64[],
    Mean_MRE = Float64[],
    Max_MRE = Float64[],
    Corr_MAE = Float64[],
)

for thresh in thresholds
    @info "\n═══ PCA threshold = $(round(Int, thresh*100))% ═══"

    # Standardize
    Z_concat, std_params = standardize(df_concat, concat_cols)

    # PCA + norms (raw, then normalize separately)
    X_pca_raw, pca_model, d_orig, pca_norms = build_memory_matrix_raw(Z_concat; pratio=thresh)
    d_pca = size(X_pca_raw, 1)

    # Unit-normalize
    X̂ = copy(X_pca_raw)
    for k in 1:K
        X̂[:, k] ./= (norm(X̂[:, k]) + 1e-12)
    end

    @info "  d_PCA = $d_pca, K/d = $(round(K/d_pca, digits=2))"

    # Find β*
    phase = find_entropy_inflection(X̂; α=α_step, n_betas=80, β_range=(0.1, 1000.0))
    β_star = phase.β_star
    @info "  β* = $(round(β_star, digits=2))"

    # Generate patients
    @info "  Generating $n_synth patients …"
    df_synth_gen = generate_patients_rescaled(X̂, pca_model, std_params, pca_norms,
                                              n_synth, T_steps; β=β_star, α=α_step, seed=42)

    # Compute MREs (compare concat columns)
    mres = Float64[]
    for col in concat_cols
        real_mean = mean(df_concat[!, col])
        synth_mean = mean(df_synth_gen[!, col])
        if abs(real_mean) > 1e-12
            push!(mres, abs(synth_mean - real_mean) / abs(real_mean))
        end
    end

    # Correlation MAE
    X_real_mat = Matrix{Float64}(df_concat[!, concat_cols])
    X_synth_mat = Matrix{Float64}(df_synth_gen)
    C_real = cor(X_real_mat)
    C_synth = cor(X_synth_mat)
    C_real[isnan.(C_real)] .= 0.0
    C_synth[isnan.(C_synth)] .= 0.0
    upper = [i < j for i in 1:size(C_real,1), j in 1:size(C_real,2)]
    corr_mae = mean(abs.(C_real[upper] - C_synth[upper]))

    # Novelty
    novelties = Float64[]
    for i in 1:n_synth
        x_synth = Vector{Float64}(df_synth_gen[i, :])
        z_synth = (x_synth .- std_params.μ) ./ std_params.σ
        ξ_synth = MultivariateStats.transform(pca_model, z_synth)
        ξ_dir = ξ_synth ./ (norm(ξ_synth) + 1e-12)
        max_cos = maximum(dot(ξ_dir, X̂[:, k]) for k in 1:K)
        push!(novelties, 1.0 - max_cos)
    end
    mean_novelty = mean(novelties)

    @info "  Results: Med MRE=$(round(median(mres), digits=4)), Corr MAE=$(round(corr_mae, digits=4)), Novelty=$(round(mean_novelty, digits=3))"

    push!(results_all, (thresh, d_pca, round(K/d_pca, digits=2), round(β_star, digits=2),
                         round(mean_novelty, digits=3),
                         round(median(mres), digits=4), round(mean(mres), digits=4),
                         round(maximum(mres), digits=4), round(corr_mae, digits=4)))
end

# ── Print summary ────────────────────────────────────────────────────────────
@info "\n═══ PCA Sensitivity Summary ═══"
for row in eachrow(results_all)
    @info "  $(round(row.Threshold*100, digits=1))%: d=$(row.d_PCA), MRE=$(row.Median_MRE), Corr=$(row.Corr_MAE), Nov=$(row.Novelty)"
end

CSV.write(joinpath(_PATH_TO_DATA, "pca_sensitivity_results.csv"), results_all)
@info "  Saved pca_sensitivity_results.csv"

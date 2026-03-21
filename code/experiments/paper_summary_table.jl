# ──────────────────────────────────────────────────────────────────────────────
# paper_summary_table.jl
#
# Per-feature, per-visit summary statistics: real vs synthetic.
# Outputs CSV + printed table.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load data
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading data …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms n_assays df_clean kept_cols

# ══════════════════════════════════════════════════════════════════════════════
# Generate 500 synthetic patients
# ══════════════════════════════════════════════════════════════════════════════
d_pca, K_mem = size(X̂)
n_synth = 500
T_steps = 2000

phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star

Random.seed!(123)
synth_concat = Matrix{Float64}(undef, n_synth, length(std_params_concat.col_names))
for i in 1:n_synth
    k0 = rand(1:K_mem)
    ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d_pca)
    ξ₀ ./= norm(ξ₀)
    result = sample(X̂, ξ₀, T_steps; β=β_star, α=0.05)
    ξ_dir = result.Ξ[end, :]
    ξ_dir ./= (norm(ξ_dir) + 1e-12)
    r = rand(pca_norms)
    ξ_rescaled = r .* ξ_dir
    z_std = vec(MultivariateStats.reconstruct(pca_concat, ξ_rescaled))
    synth_concat[i, :] = z_std .* std_params_concat.σ .+ std_params_concat.μ
end

# split into visits
df_synth_visits = DataFrame[]
for v in 1:3
    offset = (v - 1) * n_assays
    cols_data = synth_concat[:, (offset+1):(offset+n_assays)]
    df_v = DataFrame(cols_data, kept_cols)
    df_v.Visit = fill(v, n_synth)
    push!(df_synth_visits, df_v)
end
df_synth = vcat(df_synth_visits...)

# ══════════════════════════════════════════════════════════════════════════════
# Build summary table
# ══════════════════════════════════════════════════════════════════════════════
@info "Building summary table …"

rows = []
for v in 1:3
    real_v = df_clean[df_clean.Visit .== v, :]
    synth_v = df_synth[df_synth.Visit .== v, :]
    for col in kept_cols
        real_vals = collect(skipmissing(real_v[!, col]))
        synth_vals = synth_v[!, col]

        r_mean = mean(real_vals)
        r_std = std(real_vals)
        s_mean = mean(synth_vals)
        s_std = std(synth_vals)

        mean_rel_err = abs(r_mean) > 1e-12 ? abs(s_mean - r_mean) / abs(r_mean) : NaN
        std_rel_err = abs(r_std) > 1e-12 ? abs(s_std - r_std) / abs(r_std) : NaN

        push!(rows, (
            Visit = v,
            Feature = string(col),
            Real_Mean = round(r_mean, digits=2),
            Real_Std = round(r_std, digits=2),
            Synth_Mean = round(s_mean, digits=2),
            Synth_Std = round(s_std, digits=2),
            Mean_Rel_Err = round(mean_rel_err, digits=4),
            Std_Rel_Err = round(std_rel_err, digits=4),
        ))
    end
end

df_table = DataFrame(rows)

# save CSV
CSV.write(joinpath(_PATH_TO_DATA, "paper_summary_statistics.csv"), df_table)

# print per-visit aggregates
println("\n" * "="^70)
println("Per-Visit Aggregate Summary")
println("="^70)
for v in 1:3
    vt = df_table[df_table.Visit .== v, :]
    avg_mean_err = mean(filter(!isnan, vt.Mean_Rel_Err))
    avg_std_err = mean(filter(!isnan, vt.Std_Rel_Err))
    med_mean_err = median(filter(!isnan, vt.Mean_Rel_Err))
    med_std_err = median(filter(!isnan, vt.Std_Rel_Err))
    println("  Visit $v:  Mean RE = $(round(avg_mean_err*100, digits=2))% (median $(round(med_mean_err*100, digits=2))%)  |  Std RE = $(round(avg_std_err*100, digits=2))% (median $(round(med_std_err*100, digits=2))%)")
end

# print full table
println("\n" * "="^70)
println("Full Feature Table")
println("="^70)
show(stdout, df_table; allrows=true, allcols=true, truncate=0)
println()

@info "✓ Saved to data/paper_summary_statistics.csv"

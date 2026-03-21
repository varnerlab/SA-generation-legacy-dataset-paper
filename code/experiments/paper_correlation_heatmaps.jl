# ──────────────────────────────────────────────────────────────────────────────
# paper_correlation_heatmaps.jl
#
# Side-by-side correlation heatmaps: real vs synthetic, per visit.
# Shows that SA preserves cross-feature correlation structure.
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
# Compute correlation matrices
# ══════════════════════════════════════════════════════════════════════════════
@info "Computing correlations …"

# short labels for readability
short_labels = String[]
for col in kept_cols
    s = string(col)
    # truncate long names
    if length(s) > 15
        push!(short_labels, s[1:13] * "…")
    else
        push!(short_labels, s)
    end
end

function corr_matrix(df, cols)
    n = length(cols)
    C = Matrix{Float64}(undef, n, n)
    for i in 1:n
        vi = collect(skipmissing(df[!, cols[i]]))
        for j in 1:n
            if i == j
                C[i,j] = 1.0
            else
                vj = collect(skipmissing(df[!, cols[j]]))
                # need paired complete cases
                if length(vi) == length(vj)
                    C[i,j] = cor(vi, vj)
                else
                    # row-wise pairing
                    xi = Float64[]
                    xj = Float64[]
                    for row in 1:nrow(df)
                        a = df[row, cols[i]]
                        b = df[row, cols[j]]
                        if !ismissing(a) && !ismissing(b)
                            push!(xi, a)
                            push!(xj, b)
                        end
                    end
                    C[i,j] = length(xi) > 2 ? cor(xi, xj) : NaN
                end
            end
        end
    end
    return C
end

# ══════════════════════════════════════════════════════════════════════════════
# Build 2 × 3 figure: [Real V1, Real V2, Real V3; Synth V1, Synth V2, Synth V3]
# ══════════════════════════════════════════════════════════════════════════════
@info "Plotting …"

visit_titles = ["Visit 1 (baseline)", "Visit 2 (1st trimester)", "Visit 3 (3rd trimester)"]
panels = []

for v in 1:3
    real_v = df_clean[df_clean.Visit .== v, :]
    C_real = corr_matrix(real_v, kept_cols)
    p = heatmap(C_real,
        title="Real — $(visit_titles[v])",
        titlefontsize=9,
        xticks=(1:n_assays, short_labels), xrotation=90,
        yticks=(1:n_assays, short_labels),
        tickfontsize=5,
        clims=(-1, 1), color=:RdBu,
        aspect_ratio=1,
        background_color=:white,
        foreground_color=:black,
        size=(500, 500))
    push!(panels, p)
end

for v in 1:3
    synth_v = df_synth[df_synth.Visit .== v, :]
    C_synth = corr_matrix(synth_v, kept_cols)
    p = heatmap(C_synth,
        title="Synthetic — $(visit_titles[v])",
        titlefontsize=9,
        xticks=(1:n_assays, short_labels), xrotation=90,
        yticks=(1:n_assays, short_labels),
        tickfontsize=5,
        clims=(-1, 1), color=:RdBu,
        aspect_ratio=1,
        background_color=:white,
        foreground_color=:black,
        size=(500, 500))
    push!(panels, p)
end

p_corr = plot(panels..., layout=(2, 3), size=(1800, 1200),
    plot_title="Feature Correlation Structure: Real vs Synthetic",
    plot_titlefontsize=14,
    background_color=:white,
    margin=8Plots.mm)

savefig(p_corr, joinpath(_PATH_TO_FIG, "correlation_heatmaps_paper.pdf"))
savefig(p_corr, joinpath(_PATH_TO_FIG, "correlation_heatmaps_paper.png"))

# ══════════════════════════════════════════════════════════════════════════════
# Also compute element-wise correlation error
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^60)
println("Correlation Matrix Agreement (Frobenius norm of difference)")
println("="^60)
for v in 1:3
    real_v = df_clean[df_clean.Visit .== v, :]
    synth_v = df_synth[df_synth.Visit .== v, :]
    C_real = corr_matrix(real_v, kept_cols)
    C_synth = corr_matrix(synth_v, kept_cols)
    # replace NaN with 0 for comparison
    C_real[isnan.(C_real)] .= 0
    C_synth[isnan.(C_synth)] .= 0
    diff = C_real - C_synth
    frob = norm(diff) / n_assays  # normalized
    mae = mean(abs.(diff))
    println("  Visit $v:  MAE = $(round(mae, digits=4))  |  Frobenius/d = $(round(frob, digits=4))")
end

@info "✓ Saved to figs/correlation_heatmaps_paper.pdf"

# ──────────────────────────────────────────────────────────────────────────────
# paper_pca_by_visit.jl
#
# PCA projection of individual records (not concatenated visits).
# Real and synthetic colored by visit, showing pregnancy trajectory shift.
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
# Per-record PCA (fit on real data, transform both)
# ══════════════════════════════════════════════════════════════════════════════
@info "Fitting per-record PCA …"

# extract real data matrix (complete rows only)
real_rows = []
real_visits = Int[]
for row in 1:nrow(df_clean)
    vals = Float64[]
    complete = true
    for col in kept_cols
        v = df_clean[row, col]
        if ismissing(v)
            complete = false
            break
        end
        push!(vals, v)
    end
    if complete
        push!(real_rows, vals)
        push!(real_visits, df_clean[row, :Visit])
    end
end
X_real = reduce(hcat, real_rows)  # n_assays × n_real

# standardize (fit on real)
μ_r = vec(mean(X_real, dims=2))
σ_r = vec(std(X_real, dims=2))
for i in eachindex(σ_r)
    σ_r[i] < 1e-12 && (σ_r[i] = 1.0)
end
Z_real = (X_real .- μ_r) ./ σ_r

# fit PCA on real
pca_rec = MultivariateStats.fit(PCA, Z_real; maxoutdim=2)
P_real = MultivariateStats.transform(pca_rec, Z_real)  # 2 × n_real

# transform synthetic
X_synth = Matrix{Float64}(undef, n_assays, nrow(df_synth))
for i in 1:nrow(df_synth)
    for (j, col) in enumerate(kept_cols)
        X_synth[j, i] = df_synth[i, col]
    end
end
Z_synth = (X_synth .- μ_r) ./ σ_r
P_synth = MultivariateStats.transform(pca_rec, Z_synth)  # 2 × n_synth_total

synth_visits = df_synth.Visit

var1 = round(100 * MultivariateStats.principalvars(pca_rec)[1] / MultivariateStats.tvar(pca_rec), digits=1)
var2 = round(100 * MultivariateStats.principalvars(pca_rec)[2] / MultivariateStats.tvar(pca_rec), digits=1)

@info "  PC1: $(var1)% variance, PC2: $(var2)% variance"

# ══════════════════════════════════════════════════════════════════════════════
# Colors (matching paper plot scheme)
# ══════════════════════════════════════════════════════════════════════════════
synth_colors = [
    RGB(0.55, 0.82, 0.82),  # V1: teal
    RGB(0.55, 0.78, 0.55),  # V2: green
    RGB(0.90, 0.65, 0.50),  # V3: salmon
]
real_colors = [
    RGB(0.05, 0.55, 0.55),  # V1: dark teal
    RGB(0.10, 0.52, 0.10),  # V2: dark green
    RGB(0.80, 0.30, 0.10),  # V3: dark coral
]
visit_labels = ["V1: not pregnant", "V2: first trimester", "V3: third trimester"]

# ══════════════════════════════════════════════════════════════════════════════
# Plot: single panel, all visits overlaid
# ══════════════════════════════════════════════════════════════════════════════
@info "Plotting …"

p = plot(
    xlabel="PC1 ($(var1)%)", ylabel="PC2 ($(var2)%)",
    title="PCA Projection: Real vs Synthetic by Visit",
    titlefontsize=13,
    background_color=:white,
    foreground_color=:black,
    grid=false,
    framestyle=:box,
    legend=:outerright,
    legendfontsize=8,
    size=(900, 550),
    margin=5Plots.mm,
)

# layer 1: synthetic (background, smaller)
for v in 1:3
    mask = synth_visits .== v
    scatter!(p, P_synth[1, mask], P_synth[2, mask],
        label="$(visit_labels[v]) (synth)",
        color=synth_colors[v],
        alpha=0.4, markersize=3, markerstrokewidth=0,
        markershape=:circle)
end

# layer 2: real (foreground, larger, with stroke)
for v in 1:3
    mask = real_visits .== v
    scatter!(p, P_real[1, mask], P_real[2, mask],
        label="$(visit_labels[v]) (real)",
        color=real_colors[v],
        alpha=0.9, markersize=7,
        markerstrokecolor=:white, markerstrokewidth=0.8,
        markershape=:circle)
end

# add visit centroid arrows (real data)
real_centroids = [mean(P_real[:, real_visits .== v], dims=2) for v in 1:3]
for i in 1:2
    c1 = real_centroids[i]
    c2 = real_centroids[i+1]
    plot!(p, [c1[1], c2[1]], [c1[2], c2[2]],
        arrow=true, color=:black, lw=2.5, label=(i==1 ? "Pregnancy trajectory" : ""),
        linestyle=:solid)
end

savefig(p, joinpath(_PATH_TO_FIG, "pca_by_visit_paper.pdf"))
savefig(p, joinpath(_PATH_TO_FIG, "pca_by_visit_paper.png"))

@info "✓ Saved to figs/pca_by_visit_paper.pdf"

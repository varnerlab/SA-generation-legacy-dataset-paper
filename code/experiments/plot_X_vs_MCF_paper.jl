# ──────────────────────────────────────────────────────────────────────────────
# plot_X_vs_MCF_paper.jl
#
# Paper version: white background, colored plot areas per panel.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load data
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading data …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms n_assays df_clean kept_cols

x_col = Symbol("X")
mcf_col = Symbol("TF Only MCF (mm)")

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

synth_all_x = df_synth[!, x_col]
synth_all_mcf = df_synth[!, mcf_col]

# ══════════════════════════════════════════════════════════════════════════════
# Axis limits
# ══════════════════════════════════════════════════════════════════════════════
all_x = vcat(collect(skipmissing(df_clean[!, x_col])), synth_all_x)
all_mcf = vcat(collect(skipmissing(df_clean[!, mcf_col])), synth_all_mcf)
x_lim = (floor(minimum(all_x)/10)*10 - 5, ceil(maximum(all_x)/10)*10 + 5)
mcf_lim = (floor(minimum(all_mcf)/5)*5 - 2, ceil(maximum(all_mcf)/5)*5 + 2)

# ══════════════════════════════════════════════════════════════════════════════
# Colors
# ══════════════════════════════════════════════════════════════════════════════
# soft tinted backgrounds for each panel
panel_bg = [
    RGB(0.92, 0.97, 0.97),  # V1: ice blue
    RGB(0.93, 0.97, 0.93),  # V2: pale green
    RGB(0.98, 0.94, 0.92),  # V3: pale peach
]

# synthetic visit-specific (lighter, more transparent)
synth_colors = [
    RGB(0.55, 0.82, 0.82),  # V1: teal
    RGB(0.55, 0.78, 0.55),  # V2: green
    RGB(0.90, 0.65, 0.50),  # V3: salmon
]

# real/actual (darker, bolder)
real_colors = [
    RGB(0.05, 0.55, 0.55),  # V1: dark teal
    RGB(0.10, 0.52, 0.10),  # V2: dark green
    RGB(0.80, 0.30, 0.10),  # V3: dark coral
]

# gray for background cloud
gray_cloud = RGB(0.35, 0.35, 0.35)

visit_labels = ["V1: not pregnant", "V2: first trimester", "V3: third trimester"]

# ══════════════════════════════════════════════════════════════════════════════
# Build figure
# ══════════════════════════════════════════════════════════════════════════════
panels = []
for v in 1:3

    # layer 1: ALL synthetic (gray cloud)
    p = scatter(synth_all_x, synth_all_mcf,
        label="Synthetic (all)",
        color=gray_cloud,
        alpha=0.45, markersize=2, markerstrokewidth=0,
        xlabel="Factor X (% nominal)", ylabel="MCF (mm)",
        title=visit_labels[v], titlefontsize=11,
        xlims=x_lim, ylims=mcf_lim,
        background_color_inside=panel_bg[v],
        background_color_outside=:white,
        foreground_color=:black,
        grid=false,
        legend=:topleft, legendfontsize=7,
        framestyle=:box)

    # layer 2: visit-specific synthetic
    synth_v = df_synth[df_synth.Visit .== v, :]
    scatter!(p, synth_v[!, x_col], synth_v[!, mcf_col],
        label="Visit $v (synthetic)",
        color=synth_colors[v],
        alpha=0.55, markersize=3.5, markerstrokewidth=0)

    # layer 3: visit-specific real (on top)
    real_v = df_clean[df_clean.Visit .== v, :]
    real_x = Float64[]
    real_mcf = Float64[]
    for row in 1:nrow(real_v)
        xv = real_v[row, x_col]
        mv = real_v[row, mcf_col]
        if !ismissing(xv) && !ismissing(mv)
            push!(real_x, xv)
            push!(real_mcf, mv)
        end
    end
    scatter!(p, real_x, real_mcf,
        label="Visit $v (actual)",
        color=real_colors[v],
        alpha=0.9, markersize=6,
        markerstrokecolor=:white, markerstrokewidth=0.8)

    push!(panels, p)
end

p_paper = plot(panels..., layout=(1, 3), size=(1500, 500),
               plot_title="Factor X vs Maximum Clot Firmness (MCF)",
               plot_titlefontsize=14,
               background_color=:white,
               margin=5Plots.mm)

savefig(p_paper, joinpath(_PATH_TO_FIG, "X_vs_MCF_paper.pdf"))
savefig(p_paper, joinpath(_PATH_TO_FIG, "X_vs_MCF_paper.png"))

@info "✓ Saved to figs/X_vs_MCF_paper.pdf"

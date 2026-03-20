# ──────────────────────────────────────────────────────────────────────────────
# plot_X_vs_MCF.jl
#
# Reproduce the JuliaCon slide 8 style: Factor X vs MCF scatter plot
# showing the pregnancy trajectory (lower-left → upper-right) across visits.
#
# Three layers per panel:
#   1. Gray cloud: ALL synthetic patients (all visits) — context
#   2. Colored: Visit-specific synthetic patients
#   3. Darker/larger: Visit-specific real/actual patients
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load data
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading data …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") df_clean kept_cols

df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

x_col = Symbol("X")
mcf_col = Symbol("TF Only MCF (mm)")

@info "  Real data: $(nrow(df_clean)) records"
@info "  Synthetic data: $(nrow(df_synth)) records"

# ══════════════════════════════════════════════════════════════════════════════
# Collect all synthetic data (all visits pooled) for background cloud
# ══════════════════════════════════════════════════════════════════════════════
synth_all_x = df_synth[!, x_col]
synth_all_mcf = df_synth[!, mcf_col]

# ══════════════════════════════════════════════════════════════════════════════
# Build the 1×3 figure
# ══════════════════════════════════════════════════════════════════════════════
visit_labels = ["V1: not pregnant", "V2: first trimester", "V3: third trimester"]
visit_colors_synth = [
    RGB(0.4, 0.85, 0.85),   # teal/cyan for V1
    RGB(0.6, 0.9, 0.6),     # light green for V2
    RGB(1.0, 0.75, 0.6),    # salmon/peach for V3
]
visit_colors_real = [
    RGB(0.1, 0.7, 0.7),     # darker teal for V1
    RGB(0.2, 0.7, 0.2),     # darker green for V2
    RGB(0.9, 0.4, 0.2),     # darker coral for V3
]

# compute common axis limits
all_x_vals = vcat(
    collect(skipmissing(df_clean[!, x_col])),
    synth_all_x
)
all_mcf_vals = vcat(
    collect(skipmissing(df_clean[!, mcf_col])),
    synth_all_mcf
)
x_lim = (floor(minimum(all_x_vals)/10)*10 - 5, ceil(maximum(all_x_vals)/10)*10 + 5)
mcf_lim = (floor(minimum(all_mcf_vals)/5)*5 - 2, ceil(maximum(all_mcf_vals)/5)*5 + 2)

panels = []
for v in 1:3
    # layer 1: ALL synthetic (gray background cloud)
    p = scatter(synth_all_x, synth_all_mcf,
        label="Synthetic (all)",
        color=RGB(0.65, 0.65, 0.65),
        alpha=0.3, markersize=2.5, markerstrokewidth=0,
        xlabel="X (% nominal)", ylabel="MCF (mm)",
        title=visit_labels[v], titlefontsize=12,
        xlims=x_lim, ylims=mcf_lim,
        background_color=RGB(0.12, 0.12, 0.12),
        foreground_color=:white,
        grid=false, legend=:topleft, legendfontsize=7,
        legendbgcolor=RGBA(0,0,0,0.5))

    # layer 2: visit-specific synthetic
    synth_v = df_synth[df_synth.Visit .== v, :]
    scatter!(p, synth_v[!, x_col], synth_v[!, mcf_col],
        label="Visit $v (synthetic)",
        color=visit_colors_synth[v],
        alpha=0.6, markersize=4, markerstrokewidth=0)

    # layer 3: visit-specific real (on top, larger, with edge)
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
        color=visit_colors_real[v],
        alpha=0.9, markersize=6,
        markerstrokecolor=:white, markerstrokewidth=0.5)

    push!(panels, p)
end

p_final = plot(panels..., layout=(1, 3), size=(1500, 500),
               plot_title="Factor X vs Maximum Clot Firmness (MCF) — SA Generated",
               plot_titlefontsize=14,
               plot_titlecolor=:white,
               background_color=RGB(0.12, 0.12, 0.12),
               margin=5Plots.mm)

savefig(p_final, joinpath(_PATH_TO_FIG, "X_vs_MCF_juliacon_style.pdf"))
savefig(p_final, joinpath(_PATH_TO_FIG, "X_vs_MCF_juliacon_style.png"))
@info "  Saved to figs/X_vs_MCF_juliacon_style.pdf"

# also make a higher-density version with more synthetic patients
@info "\nGenerating higher-density version (500 patients) …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms n_assays

d_pca, K_mem = size(X̂)
n_synth_hd = 500
T_steps = 2000

phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star

Random.seed!(123)
synth_hd = Matrix{Float64}(undef, n_synth_hd, length(std_params_concat.col_names))
for i in 1:n_synth_hd
    k0 = rand(1:K_mem)
    ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d_pca)
    ξ₀ ./= norm(ξ₀)
    result = sample(X̂, ξ₀, T_steps; β=β_star, α=0.05)
    ξ_dir = result.Ξ[end, :]
    ξ_dir ./= (norm(ξ_dir) + 1e-12)
    r = rand(pca_norms)
    ξ_rescaled = r .* ξ_dir
    z_std = vec(MultivariateStats.reconstruct(pca_concat, ξ_rescaled))
    synth_hd[i, :] = z_std .* std_params_concat.σ .+ std_params_concat.μ
end

# split into visits
df_hd_visits = DataFrame[]
for v in 1:3
    offset = (v - 1) * n_assays
    x_idx = findfirst(kept_cols .== x_col)
    mcf_idx = findfirst(kept_cols .== mcf_col)
    cols_data = synth_hd[:, (offset+1):(offset+n_assays)]
    df_v = DataFrame(cols_data, kept_cols)
    df_v.Visit = fill(v, n_synth_hd)
    push!(df_hd_visits, df_v)
end
df_hd = vcat(df_hd_visits...)

hd_all_x = df_hd[!, x_col]
hd_all_mcf = df_hd[!, mcf_col]

# recompute limits including HD data
all_x2 = vcat(all_x_vals, hd_all_x)
all_mcf2 = vcat(all_mcf_vals, hd_all_mcf)
x_lim2 = (floor(minimum(all_x2)/10)*10 - 5, ceil(maximum(all_x2)/10)*10 + 5)
mcf_lim2 = (floor(minimum(all_mcf2)/5)*5 - 2, ceil(maximum(all_mcf2)/5)*5 + 2)

panels_hd = []
for v in 1:3
    p = scatter(hd_all_x, hd_all_mcf,
        label="Synthetic (all)",
        color=RGB(0.65, 0.65, 0.65),
        alpha=0.25, markersize=2, markerstrokewidth=0,
        xlabel="X (% nominal)", ylabel="MCF (mm)",
        title=visit_labels[v], titlefontsize=12,
        xlims=x_lim2, ylims=mcf_lim2,
        background_color=RGB(0.12, 0.12, 0.12),
        foreground_color=:white,
        grid=false, legend=:topleft, legendfontsize=7,
        legendbgcolor=RGBA(0,0,0,0.5))

    hd_v = df_hd[df_hd.Visit .== v, :]
    scatter!(p, hd_v[!, x_col], hd_v[!, mcf_col],
        label="Visit $v (synthetic)",
        color=visit_colors_synth[v],
        alpha=0.5, markersize=3.5, markerstrokewidth=0)

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
        color=visit_colors_real[v],
        alpha=0.9, markersize=7,
        markerstrokecolor=:white, markerstrokewidth=0.5)

    push!(panels_hd, p)
end

p_hd = plot(panels_hd..., layout=(1, 3), size=(1500, 500),
            plot_title="Factor X vs Maximum Clot Firmness (MCF) — SA Generated (N=500)",
            plot_titlefontsize=14,
            plot_titlecolor=:white,
            background_color=RGB(0.12, 0.12, 0.12),
            margin=5Plots.mm)

savefig(p_hd, joinpath(_PATH_TO_FIG, "X_vs_MCF_juliacon_style_hd.pdf"))
savefig(p_hd, joinpath(_PATH_TO_FIG, "X_vs_MCF_juliacon_style_hd.png"))
@info "  Saved HD version to figs/X_vs_MCF_juliacon_style_hd.pdf"

@info "\n✓ Done!"

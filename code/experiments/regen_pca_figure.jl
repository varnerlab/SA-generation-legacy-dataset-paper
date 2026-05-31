# Standalone script to regenerate the SA vs MVN PCA figure
# Matches the original paper_sa_vs_mvn.jl approach: per-visit rows, standardized, single PCA
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots, LinearAlgebra, Random, Distributions, MultivariateStats

const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ── Load data ────────────────────────────────────────────────────────────────
@info "Loading data …"
df_real  = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

# Identify assay columns from synthetic CSV (no metadata)
synth_meta = [:SyntheticID, :Visit]
kept_cols = [c for c in Symbol.(names(df_synth)) if c ∉ synth_meta]
n_assays = length(kept_cols)

# ── Build per-visit row matrices (n_assays × n_rows) ─────────────────────────
# Real: each row in df_real with complete data becomes one column
real_rows = Vector{Float64}[]
real_visits = Int[]
for row in 1:nrow(df_real)
    vals = Float64[]
    complete = true
    for col in kept_cols
        v = df_real[row, col]
        if ismissing(v)
            complete = false; break
        end
        push!(vals, Float64(v))
    end
    if complete
        push!(real_rows, vals)
        push!(real_visits, Int(df_real[row, :Visit]))
    end
end
X_real = reduce(hcat, real_rows)  # n_assays × n_real_rows

# Standardize (fit on real)
μ_r = vec(mean(X_real, dims=2))
σ_r = vec(std(X_real, dims=2))
for i in eachindex(σ_r); σ_r[i] < 1e-12 && (σ_r[i] = 1.0); end
Z_real = (X_real .- μ_r) ./ σ_r

# Fit PCA on standardized real data
pca_rec = fit(PCA, Z_real; maxoutdim=2)
P_real = MultivariateStats.transform(pca_rec, Z_real)
var1 = round(100 * principalvars(pca_rec)[1] / tvar(pca_rec), digits=1)
var2 = round(100 * principalvars(pca_rec)[2] / tvar(pca_rec), digits=1)

# SA synthetic
X_sa = Matrix{Float64}(undef, n_assays, nrow(df_synth))
for i in 1:nrow(df_synth)
    for (j, col) in enumerate(kept_cols)
        X_sa[j, i] = Float64(df_synth[i, col])
    end
end
Z_sa = (X_sa .- μ_r) ./ σ_r
P_sa = MultivariateStats.transform(pca_rec, Z_sa)
sa_visits = Int.(df_synth.Visit)

# MVN synthetic — concatenated MVN (single 216-dim distribution, matching methods section)
complete_subjects = let
    all_s = unique(df_real.SubjectID)
    [s for s in all_s if sum(df_real.SubjectID .== s) == 3 &&
       all(v -> any((df_real.SubjectID .== s) .& (df_real.Visit .== v)), 1:3)]
end
K = length(complete_subjects)
n_visits = 3
d_concat = n_assays * n_visits
n_mvn = 100

# Build concatenated real matrix (K × d_concat) for MVN fitting
X_real_concat = Matrix{Float64}(undef, K, d_concat)
for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        row = findfirst((df_real.SubjectID .== s) .& (df_real.Visit .== v))
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_real_concat[i, offset + j] = Float64(df_real[row, col])
        end
    end
end
μ_mvn = vec(mean(X_real_concat, dims=1))
Σ_mvn = cov(X_real_concat); Σ_mvn = (Σ_mvn + Σ_mvn') / 2; Σ_mvn += 1e-6 * I

Random.seed!(42)
X_mvn_concat = rand(MvNormal(μ_mvn, Σ_mvn), n_mvn)'  # n_mvn × d_concat

# Split into per-visit rows for PCA projection (same as real and SA)
mvn_rows = Vector{Float64}[]
mvn_visits_vec = Int[]
for i in 1:n_mvn
    for v in 1:n_visits
        offset = (v - 1) * n_assays
        push!(mvn_rows, X_mvn_concat[i, offset+1:offset+n_assays])
        push!(mvn_visits_vec, v)
    end
end
X_mvn = reduce(hcat, mvn_rows)
Z_mvn = (X_mvn .- μ_r) ./ σ_r
P_mvn = MultivariateStats.transform(pca_rec, Z_mvn)

@info "  Real rows: $(size(X_real, 2)), SA rows: $(nrow(df_synth)), MVN rows: $(length(mvn_visits_vec))"

# ── Per-column (per-visit) shared axis limits ────────────────────────────────
# Each visit's SA (top) and MVN (bottom) panel share xlims/ylims so the two
# synthetic methods are compared on an identical frame, with the real cloud
# anchored at the same scale in both. Limits span real+SA+MVN for that visit.
function visit_limits(v)
    rm = real_visits .== v
    sm = sa_visits .== v
    mm = mvn_visits_vec .== v
    xs = vcat(P_real[1, rm], P_sa[1, sm], P_mvn[1, mm])
    ys = vcat(P_real[2, rm], P_sa[2, sm], P_mvn[2, mm])
    padx = 0.05 * (maximum(xs) - minimum(xs))
    pady = 0.05 * (maximum(ys) - minimum(ys))
    return ((minimum(xs) - padx, maximum(xs) + padx),
            (minimum(ys) - pady, maximum(ys) + pady))
end
col_lims = [visit_limits(v) for v in 1:3]

# ── Colors ───────────────────────────────────────────────────────────────────
synth_colors = [RGB(0.45, 0.70, 0.88), RGB(0.55, 0.80, 0.55), RGB(0.92, 0.60, 0.40)]
real_colors  = [RGB(0.10, 0.45, 0.70), RGB(0.25, 0.65, 0.25), RGB(0.85, 0.35, 0.10)]
visit_labels = ["V1: not pregnant", "V2: first trimester", "V3: third trimester"]
bg_color = RGB(0.97, 0.97, 0.98)

# ── Build 2×3 panel figure ───────────────────────────────────────────────────
panels = []

for v in 1:3
    rm = real_visits .== v
    sm = sa_visits .== v
    p = scatter(P_real[1, rm], P_real[2, rm];
        label="Real", color=real_colors[v], alpha=0.85, markersize=7,
        markerstrokecolor=:white, markerstrokewidth=0.5,
        xlabel="PC1 ($var1%)", ylabel="PC2 ($var2%)",
        title="$(visit_labels[v])", titlefontsize=15,
        xguidefontsize=13, yguidefontsize=13, tickfontsize=11,
        xlims=col_lims[v][1], ylims=col_lims[v][2],
        background_color_inside=bg_color,
        grid=false, framestyle=:box, legendfontsize=9, margin=3Plots.mm)
    scatter!(p, P_sa[1, sm], P_sa[2, sm];
        label="SA synth", color=synth_colors[v], alpha=0.7, markersize=6,
        markerstrokecolor=:white, markerstrokewidth=0.5)
    push!(panels, p)
end

for v in 1:3
    rm = real_visits .== v
    mm = mvn_visits_vec .== v
    p = scatter(P_real[1, rm], P_real[2, rm];
        label="Real", color=real_colors[v], alpha=0.85, markersize=7,
        markerstrokecolor=:white, markerstrokewidth=0.5,
        xlabel="PC1 ($var1%)", ylabel="PC2 ($var2%)",
        title=" ", titlefontsize=15,  # blank title reserves height to match SA row
        xguidefontsize=13, yguidefontsize=13, tickfontsize=11,
        xlims=col_lims[v][1], ylims=col_lims[v][2],
        background_color_inside=bg_color,
        grid=false, framestyle=:box, legendfontsize=9, margin=3Plots.mm)
    scatter!(p, P_mvn[1, mm], P_mvn[2, mm];
        label="MVN synth", color=synth_colors[v], alpha=0.7, markersize=6,
        markerstrokecolor=:white, markerstrokewidth=0.5, markershape=:diamond)
    push!(panels, p)
end

p_pca = plot(panels...; layout=(2, 3), size=(1500, 900),
    margin=5Plots.mm)
savefig(p_pca, joinpath(_PATH_TO_FIG, "sa_vs_mvn_pca_by_visit_v2.pdf"))
savefig(p_pca, joinpath(_PATH_TO_FIG, "sa_vs_mvn_pca_by_visit_v2.png"))
@info "  Saved sa_vs_mvn_pca_by_visit_v2.pdf"

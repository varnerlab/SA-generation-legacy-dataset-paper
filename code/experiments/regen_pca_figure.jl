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

# MVN synthetic — generate per-visit (independent MVN per visit, like the JuliaCon baseline)
complete_subjects = let
    all_s = unique(df_real.SubjectID)
    [s for s in all_s if sum(df_real.SubjectID .== s) == 3 &&
       all(v -> any((df_real.SubjectID .== s) .& (df_real.Visit .== v)), 1:3)]
end
n_mvn = 100

Random.seed!(42)
mvn_rows = Vector{Float64}[]
mvn_visits_vec = Int[]
for v in 1:3
    # get real data for this visit
    visit_data = Float64[]
    for s in complete_subjects
        row = findfirst((df_real.SubjectID .== s) .& (df_real.Visit .== v))
        for col in kept_cols
            push!(visit_data, Float64(df_real[row, col]))
        end
    end
    X_v = reshape(visit_data, n_assays, length(complete_subjects))
    μ_v = vec(mean(X_v, dims=2))
    Σ_v = cov(X_v'); Σ_v = (Σ_v + Σ_v') / 2; Σ_v += 1e-6 * I
    samples = rand(MvNormal(μ_v, Σ_v), n_mvn)
    for i in 1:n_mvn
        push!(mvn_rows, samples[:, i])
        push!(mvn_visits_vec, v)
    end
end
X_mvn = reduce(hcat, mvn_rows)
Z_mvn = (X_mvn .- μ_r) ./ σ_r
P_mvn = MultivariateStats.transform(pca_rec, Z_mvn)

@info "  Real rows: $(size(X_real, 2)), SA rows: $(nrow(df_synth)), MVN rows: $(length(mvn_visits_vec))"

# ── Colors ───────────────────────────────────────────────────────────────────
synth_colors = [RGB(0.35, 0.70, 0.70), RGB(0.40, 0.68, 0.40), RGB(0.85, 0.50, 0.35)]
real_colors  = [RGB(0.05, 0.40, 0.40), RGB(0.05, 0.40, 0.05), RGB(0.65, 0.20, 0.05)]
visit_labels = ["V1: not pregnant", "V2: first trimester", "V3: third trimester"]
bg_color = RGB(0.97, 0.97, 0.98)

# ── Build 2×3 panel figure ───────────────────────────────────────────────────
panels = []

for v in 1:3
    rm = real_visits .== v
    sm = sa_visits .== v
    p = scatter(P_real[1, rm], P_real[2, rm];
        label="Real", color=real_colors[v], alpha=0.8, markersize=5,
        markerstrokecolor=:white, markerstrokewidth=0.5,
        xlabel="PC1 ($var1%)", ylabel="PC2 ($var2%)",
        title="SA — $(visit_labels[v])", titlefontsize=9,
        background_color_inside=bg_color,
        grid=false, framestyle=:box, legendfontsize=7, margin=3Plots.mm)
    scatter!(p, P_sa[1, sm], P_sa[2, sm];
        label="SA synth", color=synth_colors[v], alpha=0.7, markersize=4.5,
        markerstrokecolor=:white, markerstrokewidth=0.5)
    push!(panels, p)
end

for v in 1:3
    rm = real_visits .== v
    mm = mvn_visits_vec .== v
    p = scatter(P_real[1, rm], P_real[2, rm];
        label="Real", color=real_colors[v], alpha=0.8, markersize=5,
        markerstrokecolor=:white, markerstrokewidth=0.5,
        xlabel="PC1 ($var1%)", ylabel="PC2 ($var2%)",
        title="MVN — $(visit_labels[v])", titlefontsize=9,
        background_color_inside=bg_color,
        grid=false, framestyle=:box, legendfontsize=7, margin=3Plots.mm)
    scatter!(p, P_mvn[1, mm], P_mvn[2, mm];
        label="MVN synth", color=synth_colors[v], alpha=0.7, markersize=4.5,
        markerstrokecolor=:white, markerstrokewidth=0.5, markershape=:diamond)
    push!(panels, p)
end

p_pca = plot(panels...; layout=(2, 3), size=(1500, 900),
    plot_title="PCA Projection: SA vs MVN by Visit",
    plot_titlefontsize=14, margin=5Plots.mm)
savefig(p_pca, joinpath(_PATH_TO_FIG, "sa_vs_mvn_pca_by_visit.pdf"))
savefig(p_pca, joinpath(_PATH_TO_FIG, "sa_vs_mvn_pca_by_visit.png"))
@info "  Saved sa_vs_mvn_pca_by_visit.pdf"

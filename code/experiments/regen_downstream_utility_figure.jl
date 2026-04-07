# Standalone script: regenerate downstream utility figure from cached results CSV
# No calibration — just reads downstream_utility_results.csv and plots
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots

const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

const FEATURE_LABELS = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

# ── Load cached results ──────────────────────────────────────────────────────
@info "Loading cached downstream utility results …"
results = CSV.read(joinpath(_PATH_TO_DATA, "downstream_utility_results.csv"), DataFrame)
@info "  $(nrow(results)) rows loaded"

# ── Scatter plot: synth-calibrated vs real-calibrated predictions ────────────
bg_color = RGB(0.97, 0.97, 0.98)
panels = []
for feat in FEATURE_LABELS
    fsub = filter(r -> r.Feature == feat, results)

    all_vals = vcat(fsub.Pred_Real, fsub.Pred_Synth)
    lo, hi = minimum(all_vals), maximum(all_vals)
    margin = 0.1 * (hi - lo)
    lo -= margin; hi += margin

    p = plot(; xlabel="Real-calibrated", ylabel="Synth-calibrated",
             title=feat, titlefontsize=18, guidefontsize=16, tickfontsize=13,
             legendfontsize=13,
             xlim=(lo, hi), ylim=(lo, hi), aspect_ratio=1,
             background_color_inside=bg_color, grid=false, framestyle=:box,
             legend= feat == FEATURE_LABELS[end] ? :topleft : false,
             margin=4Plots.mm)
    plot!(p, [lo, hi], [lo, hi]; ls=:dash, lc=:gray50, lw=1.5, label="y=x")

    v2 = filter(r -> r.Visit == 2, fsub)
    v3 = filter(r -> r.Visit == 3, fsub)

    scatter!(p, v2.Pred_Real, v2.Pred_Synth;
             label="V2 (1st tri)", color=RGB(0.25, 0.65, 0.25), ms=6, ma=0.8,
             markerstrokecolor=:white, markerstrokewidth=0.8)
    scatter!(p, v3.Pred_Real, v3.Pred_Synth;
             label="V3 (3rd tri)", color=RGB(0.85, 0.35, 0.10), ms=6, ma=0.8,
             markerstrokecolor=:white, markerstrokewidth=0.8)

    push!(panels, p)
end

p_util = plot(panels...; layout=(1, 5), size=(2200, 500), margin=6Plots.mm,
              left_margin=10Plots.mm)
savefig(p_util, joinpath(_PATH_TO_FIG, "downstream_utility_scatter.pdf"))
savefig(p_util, joinpath(_PATH_TO_FIG, "downstream_utility_scatter.png"))
@info "  Saved downstream_utility_scatter.pdf"

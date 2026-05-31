# Standalone script: regenerate downstream utility figure from cached results CSV
# No calibration — just reads downstream_utility_results.csv and plots
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots

const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

const FEATURE_LABELS = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

# ── Pick ~n round tick values spanning [lo, hi] (avoids crowded 4-digit ticks) ─
function nice_ticks(lo, hi; n=4)
    raw = (hi - lo) / n
    mag = 10.0^floor(log10(raw))
    r = raw / mag
    step = mag * (r <= 1.5 ? 1.0 : r <= 3.0 ? 2.0 : r <= 7.0 ? 5.0 : 10.0)
    t0 = ceil(lo / step) * step
    return collect(t0:step:hi)
end

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
    tks = nice_ticks(lo, hi)

    p = plot(; xlabel="Real-calibrated", ylabel="Synth-calibrated",
             title=feat, titlefontsize=18, guidefontsize=16, tickfontsize=13,
             legendfontsize=13,
             xlim=(lo, hi), ylim=(lo, hi), xticks=tks, yticks=tks,
             background_color_inside=bg_color, grid=false, framestyle=:box,
             legend= feat == FEATURE_LABELS[end] ? :topleft : false,
             top_margin=2Plots.mm, bottom_margin=15Plots.mm,
             left_margin=(feat == FEATURE_LABELS[1] ? 16Plots.mm : 9Plots.mm),
             right_margin=3Plots.mm)
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

p_util = plot(panels...; layout=(1, 5), size=(2400, 600))
savefig(p_util, joinpath(_PATH_TO_FIG, "downstream_utility_scatter_v2.pdf"))
savefig(p_util, joinpath(_PATH_TO_FIG, "downstream_utility_scatter_v2.png"))
@info "  Saved downstream_utility_scatter_v2.pdf"

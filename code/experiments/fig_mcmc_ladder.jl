# ──────────────────────────────────────────────────────────────────────────────
# fig_mcmc_ladder.jl
#
# Supplement figure for the SA sampling-correctness investigation (npj major
# revision). Plots (from already-computed CSVs — no new sampling is run here):
#   Panel A — V0→V3 MCMC ladder: DCR(synth→real) per rung + MIA AUC (V0, V3)
#             on a twin axis, against the real→real DCR baseline.
#   Panel B — V0→V3 ladder fidelity cost: Median_MRE per rung.
#   Panel C — β-privacy frontier: DCR(synth→real) vs β, against the same
#             real→real baseline, with Median_MRE on a twin axis and the
#             operating point β* marked.
#
# Reads:  data/mcmc_ladder_results.csv, data/mcmc_ladder_mia.csv,
#         data/beta_privacy_sweep.csv
# Writes: figs/mcmc_ladder.pdf, figs/mcmc_ladder.png
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
gr()

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load data
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Loading CSVs …"

df_ladder = CSV.read(joinpath(_PATH_TO_DATA, "mcmc_ladder_results.csv"), DataFrame)
df_mia    = CSV.read(joinpath(_PATH_TO_DATA, "mcmc_ladder_mia.csv"), DataFrame)
df_beta   = CSV.read(joinpath(_PATH_TO_DATA, "beta_privacy_sweep.csv"), DataFrame)

# ladder is ordered V0..V3 in the CSV; keep that order explicitly
rung_order = ["V0", "V1", "V2", "V3"]
df_ladder = df_ladder[indexin(rung_order, df_ladder.Rung), :]

# real→real DCR baseline — identical on every row by construction
baseline_dcr = df_ladder.DCR_real_to_real_median[1]
@assert all(isapprox.(df_ladder.DCR_real_to_real_median, baseline_dcr; atol=1e-8))
@assert all(isapprox.(df_beta.DCR_real_to_real_median, baseline_dcr; atol=1e-8))

dcr_ladder = df_ladder.DCR_synth_to_real_median
mre_ladder = 100 .* df_ladder.Median_MRE

# self-checks: the whole point of the figure is that DCR never reaches baseline
@assert all(dcr_ladder .< baseline_dcr) "Ladder DCR must stay below the real→real baseline"
@assert all(df_beta.DCR_synth_to_real_median .< baseline_dcr) "β-sweep DCR must stay below the real→real baseline"

# MIA AUC summary rows (only V0 and V3 have MIA, by design — ladder endpoints)
get_stat(scheme, stat) = only(
    df_mia[(df_mia.Scheme .== scheme) .& (df_mia.Kind .== "summary") .& (df_mia.Rep .== stat), :AUC]
)
v0_auc_mean, v0_auc_sd = get_stat("V0", "mean_auc"), get_stat("V0", "sd_auc")
v3_auc_mean, v3_auc_sd = get_stat("V3", "mean_auc"), get_stat("V3", "sd_auc")

@info "  Ladder DCR (synth→real): $(round.(dcr_ladder, digits=3))"
@info "  Ladder baseline (real→real): $(round(baseline_dcr, digits=3))"
@info "  MIA AUC — V0: $(round(v0_auc_mean,digits=4)) ± $(round(v0_auc_sd,digits=4)),  V3: $(round(v3_auc_mean,digits=4)) ± $(round(v3_auc_sd,digits=4))"
@info "  Ladder Median_MRE (%): $(round.(mre_ladder, digits=2))"

# β-sweep
β_frac   = df_beta.beta_frac
dcr_beta = df_beta.DCR_synth_to_real_median
mre_beta = 100 .* df_beta.Median_MRE

# closest-approach gap (largest DCR in-sweep is the closest approach to baseline)
gap_closest, idx_closest = findmin(baseline_dcr .- dcr_beta)
@info "  β-sweep closest-approach gap: $(round(gap_closest,digits=2)) at β_frac=$(β_frac[idx_closest])"

β_star_frac = 1.0
@assert β_star_frac in β_frac "β_frac=1.0 (β*) must be one of the swept points"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Shared style (author convention — gray panel bg, no grid, boxed frame,
# legends on every panel, no overall figure title)
# ══════════════════════════════════════════════════════════════════════════════
bg_color       = RGB(0.96, 0.96, 0.96)
privacy_blue   = RGB(0.20, 0.50, 0.72)   # DCR / privacy series
fidelity_orange = RGB(0.85, 0.45, 0.25)  # MRE / fidelity series
auc_color      = RGB(0.55, 0.12, 0.15)   # MIA AUC markers (maroon — distinct from the black dashed baseline)

common_attrs = (background_color=bg_color, background_color_subplot=bg_color,
                foreground_color=:black,
                grid=false, framestyle=:box,
                guidefontsize=10, tickfontsize=8, legendfontsize=7,
                margin=7Plots.mm)

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Panel A — ladder privacy: DCR per rung + MIA AUC (twin axis)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 3: Building Panel A (ladder privacy) …"

x_rung = 1:4

pA = plot(x_rung, dcr_ladder;
    xticks=(x_rung, rung_order), xlabel="Ladder rung (V0 → V3: growing pooled memory)",
    ylabel="DCR (synthetic → real)", ylims=(13.5, 16.6),
    label="DCR synth→real", color=privacy_blue, lw=2.5,
    marker=:circle, markersize=7, markerstrokecolor=:white, markerstrokewidth=0.8,
    legend=:bottomleft,
    title="A", titlelocation=:left, titlefontsize=12,
    common_attrs...)

hline!(pA, [baseline_dcr]; color=:black, ls=:dash, lw=1.5, label="real→real baseline")
annotate!(pA, [(1.9, baseline_dcr - 0.30,
    text("real→real = $(round(baseline_dcr, digits=2))", 8, :black, :center))])

pA2 = twinx(pA)
plot!(pA2, [1, 4], [0.5, 0.5]; color=:gray50, ls=:dot, lw=1.2,
    ylabel="MIA AUC", ylims=(0.4, 1.08), right_margin=14Plots.mm,
    label="chance (AUC = 0.5)", legend=:topright,
    background_color_subplot=bg_color,
    guidefontsize=10, tickfontsize=8, legendfontsize=7)
scatter!(pA2, [1, 4], [v0_auc_mean, v3_auc_mean];
    yerror=[v0_auc_sd, v3_auc_sd],
    color=auc_color, markerstrokecolor=:black, markershape=:diamond, markersize=7,
    markerstrokewidth=1.8, linecolor=:black, msc=:black,
    label="MIA AUC (mean±SD)")
annotate!(pA2, [(1.15, v0_auc_mean - 0.10, text("$(round(v0_auc_mean,digits=3))±$(round(v0_auc_sd,digits=3))", 7, :black, :left))])
annotate!(pA2, [(3.5, v3_auc_mean - 0.13, text("$(round(v3_auc_mean,digits=3))±$(round(v3_auc_sd,digits=3))", 7, :black, :left))])

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Panel B — ladder fidelity cost: Median_MRE per rung
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 4: Building Panel B (ladder fidelity) …"

pB = bar(x_rung, mre_ladder;
    xticks=(x_rung, rung_order), xlabel="Ladder rung (V0 → V3: growing pooled memory)",
    ylabel="Median MRE (%)", ylims=(0, 2.6),
    label="Median MRE", color=fidelity_orange, alpha=0.85, bar_width=0.6,
    legend=:topleft,
    title="B", titlelocation=:left, titlefontsize=12,
    common_attrs...)

for (xi, yi) in zip(x_rung, mre_ladder)
    annotate!(pB, [(xi, yi + 0.10, text("$(round(yi, digits=2))%", 8, :black, :center))])
end
annotate!(pB, [(2.5, 2.45, text("pooling (V2/V3) raises MRE\nvs V0/V1", 7, :black, :center))])

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Panel C — β-privacy frontier: DCR vs β, MRE twin axis, β* marker
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 5: Building Panel C (β-privacy frontier) …"

pC = plot(β_frac, dcr_beta;
    xlabel="β / β*  (β_frac)", ylabel="DCR (synthetic → real)", ylims=(13.0, 16.6),
    label="DCR synth→real", color=privacy_blue, lw=2.5,
    marker=:circle, markersize=7, markerstrokecolor=:white, markerstrokewidth=0.8,
    legend=:bottomleft,
    title="C", titlelocation=:left, titlefontsize=12,
    common_attrs...)

hline!(pC, [baseline_dcr]; color=:black, ls=:dash, lw=1.5, label="real→real baseline")
vline!(pC, [β_star_frac]; color=:gray20, ls=:dashdot, lw=1.5, label="operating point β* (β≈2.94)")
annotate!(pC, [(β_frac[idx_closest] + 0.05, (baseline_dcr + dcr_beta[idx_closest]) / 2,
    text("gap ≈ $(round(gap_closest, digits=2))", 8, :black, :left))])
plot!(pC, [β_frac[idx_closest], β_frac[idx_closest]], [dcr_beta[idx_closest], baseline_dcr];
    color=:gray30, lw=1, ls=:dot, label=nothing, arrow=false)

pC2 = twinx(pC)
plot!(pC2, β_frac, mre_beta;
    ylabel="Median MRE (%)", ylims=(1.0, 1.6), right_margin=10Plots.mm,
    color=fidelity_orange, lw=2.5, marker=:utriangle, markersize=6,
    markerstrokecolor=:white, markerstrokewidth=0.5,
    label="Median MRE", legend=:topright,
    background_color_subplot=bg_color,
    guidefontsize=10, tickfontsize=8, legendfontsize=7)
annotate!(pC2, [(β_star_frac + 0.06, mre_beta[findfirst(==(β_star_frac), β_frac)] + 0.03,
    text("MRE min. near β*", 7.5, :black, :left))])

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Compose and save (no overall title — author convention)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 6: Composing and saving figure …"

p_all = plot(pA, pB, pC; layout=(1, 3), size=(2100, 580),
    background_color=:white, left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=4Plots.mm)

savefig(p_all, joinpath(_PATH_TO_FIG, "mcmc_ladder.pdf"))
savefig(p_all, joinpath(_PATH_TO_FIG, "mcmc_ladder.png"))
@info "  Saved figs/mcmc_ladder.pdf and figs/mcmc_ladder.png"

# ══════════════════════════════════════════════════════════════════════════════
# Final report
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^78)
println("MCMC Ladder + β-Privacy Frontier Figure — Summary")
println("="^78)
println("  Panel A: DCR per rung = $(round.(dcr_ladder, digits=3))  (baseline = $(round(baseline_dcr,digits=3)))")
println("           MIA AUC:  V0 = $(round(v0_auc_mean,digits=4)) ± $(round(v0_auc_sd,digits=4))" *
        "   V3 = $(round(v3_auc_mean,digits=4)) ± $(round(v3_auc_sd,digits=4))")
println("  Panel B: Median MRE per rung (%) = $(round.(mre_ladder, digits=2))")
println("  Panel C: β_frac sweep = $(β_frac)")
println("           DCR(β)      = $(round.(dcr_beta, digits=3))")
println("           closest-approach gap = $(round(gap_closest,digits=3)) at β_frac=$(β_frac[idx_closest])")
println("           β* = $(β_star_frac) (β ≈ $(round(df_beta.beta[findfirst(==(β_star_frac), β_frac)], digits=2)))")
println("="^78)

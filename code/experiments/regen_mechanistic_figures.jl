# Standalone script: regenerate Level 4 mechanistic figures from cached results CSV
# No ODE solving — just reads validate_mechanistic_results.csv and plots
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots, StatsPlots, StatsBase, HypothesisTests, Printf

const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ── Load cached results ──────────────────────────────────────────────────────
@info "Loading cached mechanistic results …"
results = CSV.read(joinpath(_PATH_TO_DATA, "validate_mechanistic_results.csv"), DataFrame)
@info "  $(nrow(results)) rows loaded"

const FEATURE_LABELS = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

# ── Colors ───────────────────────────────────────────────────────────────────
visit_colors = [RGB(0.10, 0.45, 0.70), RGB(0.25, 0.65, 0.25), RGB(0.85, 0.35, 0.10)]
visit_labels_short = ["V1 (BL)", "V2 (1st tri)", "V3 (3rd tri)"]
bg_color = RGB(0.97, 0.97, 0.98)

# ── Plot 1: Predicted vs Measured scatter ────────────────────────────────────
for cond in ["TF-only", "TF+TM"]
    csub = filter(r -> r.Condition == cond, results)

    plots_list = []
    for (fi, feat) in enumerate(FEATURE_LABELS)
        fsub = filter(r -> r.Feature == feat, csub)

        all_vals = vcat(fsub.Measured, fsub.Predicted)
        lo, hi = minimum(all_vals), maximum(all_vals)
        margin = 0.1 * max(hi - lo, 1e-6)
        lo -= margin; hi += margin

        p = plot(; xlabel="Dataset value", ylabel= fi == 1 ? "ODE predicted" : "",
                 title=feat,
                 legend= fi == 5 ? :topleft : false,
                 xlim=(lo, hi), ylim=(lo, hi),
                 guidefontsize=9, titlefontsize=10, tickfontsize=7,
                 background_color_inside=bg_color,
                 grid=false, framestyle=:box,
                 margin=4Plots.mm, aspect_ratio=1)
        plot!(p, [lo, hi], [lo, hi]; ls=:dash, lc=:gray50, lw=1.5, label="y=x")

        # Real patients: filled circles
        for v in 1:3
            vsub = filter(r -> r.Source == "Real" && r.Visit == v, fsub)
            if nrow(vsub) > 0
                scatter!(p, vsub.Measured, vsub.Predicted;
                         label="Real $(visit_labels_short[v])",
                         color=visit_colors[v], ms=6, ma=0.85,
                         markerstrokecolor=:white, markerstrokewidth=0.8,
                         markershape=:circle)
            end
        end

        # Synthetic patients: open rings
        for v in 1:3
            vsub = filter(r -> r.Source == "Synthetic" && r.Visit == v, fsub)
            if nrow(vsub) > 0
                scatter!(p, vsub.Measured, vsub.Predicted;
                         label="Synth $(visit_labels_short[v])",
                         color=visit_colors[v], ms=4, ma=0.5,
                         markerstrokecolor=visit_colors[v], markerstrokewidth=1.2,
                         markershape=:circle,
                         markerstrokealpha=0.6,
                         fillalpha=0.0)
            end
        end

        push!(plots_list, p)
    end

    # Save standalone scatter (for supplement TF+TM)
    cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
    p_scatter_only = plot(plots_list...; layout=(1, 5), size=(2200, 500),
                 margin=6Plots.mm, left_margin=10Plots.mm)
    savefig(p_scatter_only, joinpath(_PATH_TO_FIG, "validate_mechanistic_$(cond_tag).pdf"))
    savefig(p_scatter_only, joinpath(_PATH_TO_FIG, "validate_mechanistic_$(cond_tag).png"))
    @info "  Saved validate_mechanistic_$(cond_tag).pdf"

    # Build ratio histograms (bottom row)
    ratio_panels = []
    for (fi, feat) in enumerate(FEATURE_LABELS)
        fsub = filter(r -> r.Feature == feat, csub)
        valid = filter(r -> abs(r.Measured) > 1e-12, fsub)

        real_ratios  = filter(r -> r.Source == "Real", valid).Predicted ./
                       filter(r -> r.Source == "Real", valid).Measured
        synth_ratios = filter(r -> r.Source == "Synthetic", valid).Predicted ./
                       filter(r -> r.Source == "Synthetic", valid).Measured

        real_ratios  = clamp.(real_ratios, 0.0, 5.0)
        synth_ratios = clamp.(synth_ratios, 0.0, 5.0)

        # Compute cloud overlap and KS test
        ks_D = NaN
        ks_p = NaN
        overlap = NaN
        if length(real_ratios) > 2 && length(synth_ratios) > 2
            real_lo, real_hi = quantile(real_ratios, 0.05), quantile(real_ratios, 0.95)
            overlap = mean((synth_ratios .>= real_lo) .& (synth_ratios .<= real_hi))

            ks_test = ApproximateTwoSampleKSTest(real_ratios, synth_ratios)
            ks_D = ks_test.δ
            ks_p = pvalue(ks_test)
        end

        p = plot(; xlabel="ODE pred / Dataset", ylabel= fi == 1 ? "Density" : "",
                 title="",
                 legend= fi == 5 ? :topright : false,
                 guidefontsize=9, titlefontsize=10, tickfontsize=7,
                 background_color_inside=bg_color,
                 grid=false, framestyle=:box,
                 margin=4Plots.mm)

        if length(real_ratios) > 2
            histogram!(p, real_ratios; normalize=:pdf, bins=20,
                       fillalpha=0.5, fc=RGB(0.2, 0.4, 0.65), lc=RGB(0.2, 0.4, 0.65),
                       label="Real")
        end
        if length(synth_ratios) > 2
            histogram!(p, synth_ratios; normalize=:pdf, bins=20,
                       fillalpha=0.4, fc=RGB(0.85, 0.45, 0.25), lc=RGB(0.85, 0.45, 0.25),
                       label="Synthetic")
        end
        vline!(p, [1.0]; lc=:black, ls=:dash, lw=1.5, label="Perfect")

        # Annotate: cloud overlap + KS test
        if !isnan(overlap)
            # Get y-axis top for positioning
            yl = Plots.ylims(p)
            y_top = yl[2]
            x_range = Plots.xlims(p)
            x_right = x_range[2] - 0.05 * (x_range[2] - x_range[1])

            ks_p_str = ks_p < 0.001 ? "p<0.001" : @sprintf("p=%.3f", ks_p)
            p_color = ks_p >= 0.05 ? :forestgreen : :red

            annotate!(p, x_right, y_top * 0.95,
                      text("OL: $(round(Int, overlap*100))%", 8, :right, :forestgreen))
            annotate!(p, x_right, y_top * 0.82,
                      text("D=$(round(ks_D, digits=3))", 8, :right, :gray30))
            annotate!(p, x_right, y_top * 0.69,
                      text(ks_p_str, 8, :right, p_color))

            @info "    $feat: KS D=$(round(ks_D, digits=3)), $ks_p_str, overlap=$(round(Int, overlap*100))%"
        end

        push!(ratio_panels, p)
    end

    # Save standalone ratio histograms
    p_ratio_only = plot(ratio_panels...; layout=(1, 5), size=(2200, 500),
                   margin=6Plots.mm, left_margin=10Plots.mm)
    savefig(p_ratio_only, joinpath(_PATH_TO_FIG, "validate_mechanistic_ratios_$(cond_tag).pdf"))
    savefig(p_ratio_only, joinpath(_PATH_TO_FIG, "validate_mechanistic_ratios_$(cond_tag).png"))
    @info "  Saved validate_mechanistic_ratios_$(cond_tag).pdf"

    # Combined 2×5 figure (scatter on top, ratios on bottom)
    all_panels = vcat(plots_list, ratio_panels)
    p_combined = plot(all_panels...; layout=(2, 5), size=(2200, 900),
                 margin=5Plots.mm, left_margin=10Plots.mm)
    savefig(p_combined, joinpath(_PATH_TO_FIG, "validate_mechanistic_combined_$(cond_tag).pdf"))
    savefig(p_combined, joinpath(_PATH_TO_FIG, "validate_mechanistic_combined_$(cond_tag).png"))
    @info "  Saved validate_mechanistic_combined_$(cond_tag).pdf"
end

# ── Plot 3: Rank correlation bar charts ──────────────────────────────────────
rank_results = CSV.read(joinpath(_PATH_TO_DATA, "validate_mechanistic_rankcorr.csv"), DataFrame)

if nrow(rank_results) > 0
    for cond in ["TF-only", "TF+TM"]
        cond_ranks = filter(r -> r.Condition == cond, rank_results)
        nrow(cond_ranks) == 0 && continue

        agg = combine(groupby(cond_ranks, [:Source, :Feature]),
                      :Spearman => mean => :MeanRho)

        feats = FEATURE_LABELS
        real_rhos = Float64[]
        synth_rhos = Float64[]
        for feat in feats
            r_real  = filter(r -> r.Source == "Real" && r.Feature == feat, agg)
            r_synth = filter(r -> r.Source == "Synthetic" && r.Feature == feat, agg)
            push!(real_rhos, nrow(r_real) > 0 ? r_real.MeanRho[1] : 0.0)
            push!(synth_rhos, nrow(r_synth) > 0 ? r_synth.MeanRho[1] : 0.0)
        end

        short_labels = [split(f)[1] for f in feats]
        x = 1:length(feats)

        p = groupedbar(
            [real_rhos synth_rhos];
            bar_position=:dodge, bar_width=0.35,
            xticks=(x, short_labels),
            ylabel="Mean Spearman ρ",
            title="Rank Correlations — $cond",
            label=["Real" "Synthetic"],
            color=[RGB(0.2, 0.4, 0.65) RGB(0.85, 0.45, 0.25)],
            fillalpha=0.8,
            legend=:topright,
            size=(600, 400),
            ylim=(-0.5, 1.0),
            background_color_inside=bg_color,
            grid=false, framestyle=:box,
        )
        hline!(p, [0.0]; lc=:gray, ls=:dash, lw=1, label="")

        cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
        savefig(p, joinpath(_PATH_TO_FIG, "validate_mechanistic_rankcorr_$(cond_tag).pdf"))
        savefig(p, joinpath(_PATH_TO_FIG, "validate_mechanistic_rankcorr_$(cond_tag).png"))
        @info "  Saved validate_mechanistic_rankcorr_$(cond_tag).pdf"
    end
end

# ── Plot 4: Visit generalization boxplots ────────────────────────────────────
for cond in ["TF-only", "TF+TM"]
    csub = filter(r -> r.Condition == cond && r.Source == "Real", results)

    plots_list = []
    for (fi, feat) in enumerate(FEATURE_LABELS)
        fsub = filter(r -> r.Feature == feat, csub)
        valid = filter(r -> abs(r.Measured) > 1e-12, fsub)

        ratios_by_visit = [
            filter(r -> r.Visit == v, valid).Predicted ./ filter(r -> r.Visit == v, valid).Measured
            for v in 1:3
        ]

        p = boxplot(["V1\n(train)" "V2\n(test)" "V3\n(test)"],
                    hcat([length(r) > 0 ? r : [NaN] for r in ratios_by_visit]...);
                    legend=false, ylabel= fi == 1 ? "Pred / Meas" : "",
                    title=split(feat)[1],
                    fillalpha=0.6, fc=RGB(0.2, 0.4, 0.65),
                    guidefontsize=9, titlefontsize=10, tickfontsize=7,
                    background_color_inside=bg_color,
                    grid=false, framestyle=:box,
                    margin=4Plots.mm)
        hline!(p, [1.0]; lc=:red, ls=:dash, lw=1.5, label="")

        push!(plots_list, p)
    end

    cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
    p_gen = plot(plots_list...; layout=(1, 5), size=(1600, 400),
                 margin=6Plots.mm, left_margin=10Plots.mm)
    savefig(p_gen, joinpath(_PATH_TO_FIG, "validate_mechanistic_generalization_$(cond_tag).pdf"))
    savefig(p_gen, joinpath(_PATH_TO_FIG, "validate_mechanistic_generalization_$(cond_tag).png"))
    @info "  Saved validate_mechanistic_generalization_$(cond_tag).pdf"
end

@info "\nDone! All mechanistic figures regenerated from cached data."

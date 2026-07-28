# ======================================================================================== #
# validate_mechanistic_plausibility.jl — Level 4: Mechanistic Plausibility
#
# Do SA-generated synthetic patients produce biologically plausible TGA predictions
# when their coagulation factor levels are fed through a separately specified,
# generator-blind ODE model calibrated on the same cohort?
#
# Uses the Hockin/Mann BZ2012 coagulation model (58 species, 64 rate constants)
# calibrated on Visit 1 real patients (5 rate constants, population-level).
#
# Key question: do synthetic patients land in the same predicted-vs-measured cloud
# as real patients? This is a population-level plausibility check, not a patient-level
# predictor (rank correlations are weak).
#
# Two experimental conditions per patient:
#   1. TF-only  (TM = 0)   → "TF Initiator" TGA columns
#   2. TF + TM  (TM = 1nM) → "TF + TM Initiator" TGA columns
#
# Workflow:
#   Step 1: Load calibrated rate constants (from calibrate_hockin_mann.jl)
#   Step 2: Run BZ2012 on ALL real patients (all 3 visits) + ALL synthetic patients
#   Step 3: Compare predicted-vs-measured clouds (real vs synthetic)
#   Step 4: Rank correlations, pred/meas ratios, and cloud overlap metrics
#   Step 5: Publication-quality plots
# ======================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using CSV, DataFrames, Statistics, Plots, StatsBase, StatsPlots, LinearAlgebra
using HockinMannModel

# ── Paths ────────────────────────────────────────────────────────────────────────
const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ── Column mappings ──────────────────────────────────────────────────────────────
const FACTOR_COLS = Dict(
    :II   => :II,
    :V    => :V,
    :VII  => :VII,
    :VIII => :VIII,
    :IX   => :IX,
    :X    => :X,
    :AT   => :AT,
    :PC   => :PC,
    :TFPI => Symbol("TFPI Free"),
)

const TF_TGA_COLS = (
    lagtime  = Symbol("TF Initiator Lagtime (min)"),
    peak     = Symbol("TF Initiator Peak (nM)"),
    tpeak    = Symbol("TF Initiator T.Peak (min)"),
    max_rate = Symbol("TF Initiator Max Rate (nM/min)"),
    etp      = Symbol("TF Initiator ETP (nM·min)"),
)

const TM_TGA_COLS = (
    lagtime  = Symbol("TF + TM  Initiator Lagtime (min)"),
    peak     = Symbol("TF + TM  Initiator Peak (nM)"),
    tpeak    = Symbol("TF + TM  Initiator T.Peak (min)"),
    max_rate = Symbol("TF + TM  Initiator Max Rate (nM/min)"),
    etp      = Symbol("TF + TM  Initiator ETP (nM·min)"),
)

const FEATURE_NAMES  = [:lagtime, :peak, :tpeak, :max_rate, :etp]
const FEATURE_LABELS = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

# ══════════════════════════════════════════════════════════════════════════════════
# Step 1: Calibrated rate constants
# ══════════════════════════════════════════════════════════════════════════════════
# These were fit on Visit 1 real patients (n=23) by calibrate_hockin_mann.jl.
# Optimization: log-scale-factor space, Fminbox(NelderMead), bounds ±log(100).
# Only 5 rate constants were adjusted; the remaining 59 are literature values.

@info "Setting up calibrated BZ2012 rate constants …"

const CALIBRATED_OVERRIDES = Dict(
    :prothrombinase_kcat => 21.94,      # default 63.5 s⁻¹  (0.345x)
    :intrinsic_xase_kcat => 0.1693,     # default 8.2 s⁻¹   (0.021x)
    :extrinsic_xase_kcat => 93.21,      # default 6.0 s⁻¹   (15.5x)
    :PC_activation_kcat  => 0.8896,     # default 0.41 s⁻¹   (2.17x)
    :mIIa_conversion_k   => 1.153e9,    # default 2.3e8 M⁻¹s⁻¹ (5.01x)
)

p_calibrated = make_rate_constants(HockinMannBZ2012; CALIBRATED_OVERRIDES...)
p_literature = default_rate_constants(HockinMannBZ2012)

@info "  5 rate constants calibrated (59 at literature values)"

# ══════════════════════════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════════════════════════

function run_patient(row, p_vec; TM_molar::Float64=0.0)
    factors = percent_nominal_to_molar(
        II   = Float64(row[FACTOR_COLS[:II]]),
        V    = Float64(row[FACTOR_COLS[:V]]),
        VII  = Float64(row[FACTOR_COLS[:VII]]),
        VIII = Float64(row[FACTOR_COLS[:VIII]]),
        IX   = Float64(row[FACTOR_COLS[:IX]]),
        X    = Float64(row[FACTOR_COLS[:X]]),
        AT   = Float64(row[FACTOR_COLS[:AT]]),
        PC   = Float64(row[FACTOR_COLS[:PC]]),
        TFPI = Float64(row[FACTOR_COLS[:TFPI]]),
    )

    u0 = patient_initial_conditions(HockinMannBZ2012;
        TF = 5e-12, TM = TM_molar, factors...)

    sol = simulate(HockinMannBZ2012; u0=u0, p=p_vec, tspan=(0.0, 1200.0), saveat=1.0)
    return extract_tga_features(HockinMannBZ2012, sol)
end

function to_clinical(f::TGAFeatures)
    return (
        lagtime  = f.lagtime / 60.0,        # s → min
        peak     = f.peak * 1e9,             # M → nM
        tpeak    = f.tpeak / 60.0,           # s → min
        max_rate = f.max_rate * 1e9 * 60.0,  # M/s → nM/min
        etp      = f.etp * 1e9 / 60.0,       # M·s → nM·min
    )
end

function get_measured(row, cols)
    return (
        lagtime  = Float64(row[cols.lagtime]),
        peak     = Float64(row[cols.peak]),
        tpeak    = Float64(row[cols.tpeak]),
        max_rate = Float64(row[cols.max_rate]),
        etp      = Float64(row[cols.etp]),
    )
end

# ══════════════════════════════════════════════════════════════════════════════════
# Step 2: Load data and run BZ2012 on all patients
# ══════════════════════════════════════════════════════════════════════════════════
@info "\nLoading data …"
df_real  = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

# Complete real patients (all 3 visits present)
all_subjects = unique(df_real.SubjectID)
complete_subjects = [s for s in all_subjects
    if sum(df_real.SubjectID .== s) == 3 &&
       all(v -> any((df_real.SubjectID .== s) .& (df_real.Visit .== v)), 1:3)]

df_real_complete = filter(r -> r.SubjectID in complete_subjects, df_real)
n_real = length(complete_subjects)
n_synth = length(unique(df_synth.SyntheticID))

@info "  Real patients (complete): $n_real ($( nrow(df_real_complete)) records)"
@info "  Synthetic patients: $n_synth ($(nrow(df_synth)) records)"

# ── Run all patients ─────────────────────────────────────────────────────────────
results = DataFrame(
    Source    = String[],
    PatientID = String[],
    Visit     = Int[],
    Condition = String[],
    Feature   = String[],
    Measured  = Float64[],
    Predicted = Float64[],
)

n_success = 0
n_fail = 0

function evaluate_all!(results, df, source, id_col, p_vec)
    global n_success, n_fail
    for i in 1:nrow(df)
        row = df[i, :]
        pid = string(row[id_col])
        visit = Int(row[:Visit])

        for (cond_label, tm_val, tga_cols) in [
            ("TF-only", 0.0, TF_TGA_COLS),
            ("TF+TM",   1e-9, TM_TGA_COLS),
        ]
            try
                feat = run_patient(row, p_vec; TM_molar=tm_val)
                pred = to_clinical(feat)
                meas = get_measured(row, tga_cols)

                for (j, fname) in enumerate(FEATURE_NAMES)
                    m = getfield(meas, fname)
                    p = getfield(pred, fname)
                    push!(results, (source, pid, visit, cond_label, FEATURE_LABELS[j], m, p))
                end
                n_success += 1
            catch e
                n_fail += 1
                @warn "  Failed: $source $pid V$visit $cond_label: $(sprint(showerror, e))"
            end
        end
    end
end

@info "\nRunning BZ2012 on all real patients (calibrated) …"
evaluate_all!(results, df_real_complete, "Real", :SubjectID, p_calibrated)
@info "  Done: $n_success succeeded, $n_fail failed"

n_success = 0; n_fail = 0
@info "\nRunning BZ2012 on all synthetic patients (calibrated) …"
evaluate_all!(results, df_synth, "Synthetic", :SyntheticID, p_calibrated)
@info "  Done: $n_success succeeded, $n_fail failed"

@info "\nTotal results: $(nrow(results)) rows"

# ══════════════════════════════════════════════════════════════════════════════════
# Step 3: Summary Statistics
# ══════════════════════════════════════════════════════════════════════════════════
@info "\n═══════════════════════════════════════════════════════"
@info "  Level 4: Mechanistic Plausibility Summary"
@info "═══════════════════════════════════════════════════════"

# ── 3a: Pred/Meas ratios by source, visit, condition, feature ────────────────
@info "\n═══ Pred/Meas Ratios (calibrated BZ2012) ═══"

for cond in ["TF-only", "TF+TM"]
    @info "\n── $cond ──"
    for source in ["Real", "Synthetic"]
        @info "\n  $source:"
        ssub = filter(r -> r.Source == source && r.Condition == cond, results)
        visits = sort(unique(ssub.Visit))

        for visit in visits
            vsub = filter(r -> r.Visit == visit, ssub)
            label = (source == "Real" && visit == 1) ? "V$visit (TRAIN)" : "V$visit"

            parts = String[]
            for feat in FEATURE_LABELS
                fsub = filter(r -> r.Feature == feat, vsub)
                valid = filter(r -> abs(r.Measured) > 1e-12, fsub)
                if nrow(valid) > 0
                    ratio = median(valid.Predicted ./ valid.Measured)
                    push!(parts, "$(split(feat)[1])=$(round(ratio, digits=2))x")
                end
            end
            @info "    $label: $(join(parts, ", "))"
        end
    end
end

# ── 3b: Spearman rank correlations ──────────────────────────────────────────────
@info "\n═══ Spearman Rank Correlations ═══"

rank_results = DataFrame(
    Source    = String[],
    Condition = String[],
    Visit     = Int[],
    Feature   = String[],
    Spearman  = Float64[],
    N         = Int[],
)

for cond in ["TF-only", "TF+TM"]
    @info "\n── $cond ──"
    for source in ["Real", "Synthetic"]
        ssub = filter(r -> r.Source == source && r.Condition == cond, results)
        visits = sort(unique(ssub.Visit))

        for visit in visits
            vsub = filter(r -> r.Visit == visit, ssub)
            label = (source == "Real" && visit == 1) ? "V$visit (TRAIN)" : "V$visit"

            parts = String[]
            for feat in FEATURE_LABELS
                fsub = filter(r -> r.Feature == feat, vsub)
                valid = filter(r -> isfinite(r.Measured) && isfinite(r.Predicted), fsub)
                if nrow(valid) >= 5
                    ρ = corspearman(valid.Measured, valid.Predicted)
                    push!(parts, "$(split(feat)[1])=$(round(ρ, digits=2))")
                    push!(rank_results, (source, cond, visit, feat, ρ, nrow(valid)))
                end
            end
            @info "  $source $label: $(join(parts, ", "))"
        end
    end
end

# ── 3c: Cloud overlap — do synthetic patients fill the same pred-vs-meas region? ──
@info "\n═══ Cloud Overlap: Real vs Synthetic ═══"
@info "  (Are synthetic patients in the same predicted-vs-measured region as real?)\n"

for cond in ["TF-only", "TF+TM"]
    @info "── $cond ──"
    for feat in FEATURE_LABELS
        real_sub  = filter(r -> r.Source == "Real" && r.Condition == cond && r.Feature == feat, results)
        synth_sub = filter(r -> r.Source == "Synthetic" && r.Condition == cond && r.Feature == feat, results)

        if nrow(real_sub) < 5 || nrow(synth_sub) < 5
            continue
        end

        # Compute pred/meas ratio distributions
        real_ratios  = real_sub.Predicted ./ max.(abs.(real_sub.Measured), 1e-12)
        synth_ratios = synth_sub.Predicted ./ max.(abs.(synth_sub.Measured), 1e-12)

        # Overlap metric: what fraction of synthetic ratios fall within the real ratio range?
        real_lo, real_hi = quantile(real_ratios, 0.05), quantile(real_ratios, 0.95)
        frac_in_range = mean((synth_ratios .>= real_lo) .& (synth_ratios .<= real_hi))

        # Distribution statistics
        μ_real  = median(real_ratios)
        μ_synth = median(synth_ratios)

        @info "  $(rpad(feat, 22))  real median=$(round(μ_real, digits=2))x, synth median=$(round(μ_synth, digits=2))x, overlap=$(round(100*frac_in_range, digits=1))%"
    end
end

# ── 3d: Overall summary table ────────────────────────────────────────────────────
@info "\n═══ Overall Summary ═══"

for cond in ["TF-only", "TF+TM"]
    csub = filter(r -> r.Condition == cond, results)

    for source in ["Real", "Synthetic"]
        ssub = filter(r -> r.Source == source, csub)
        valid = filter(r -> abs(r.Measured) > 1e-12, ssub)
        ratios = valid.Predicted ./ valid.Measured
        re = abs.(ratios .- 1.0)

        @info "  $cond / $source: median pred/meas = $(round(median(ratios), digits=3))x, median |RE| = $(round(median(re), digits=3)), n=$(nrow(valid))"
    end
end

# ══════════════════════════════════════════════════════════════════════════════════
# Step 4: Publication-quality plots
# ══════════════════════════════════════════════════════════════════════════════════
@info "\nGenerating plots …"

# ── Plot 1: Predicted vs Measured scatter (real + synthetic, by visit) ────────
# Color = visit (teal/green/orange), fill = source (solid=Real, open=Synth)
visit_colors = [RGB(0.05, 0.40, 0.40), RGB(0.05, 0.40, 0.05), RGB(0.65, 0.20, 0.05)]
visit_labels_short = ["V1 (BL)", "V2 (1st tri)", "V3 (3rd tri)"]
bg_color = RGB(0.97, 0.97, 0.98)

for cond in ["TF-only", "TF+TM"]
    csub = filter(r -> r.Condition == cond, results)

    plots_list = []
    for (fi, feat) in enumerate(FEATURE_LABELS)
        fsub = filter(r -> r.Feature == feat, csub)

        all_vals = vcat(fsub.Measured, fsub.Predicted)
        lo, hi = minimum(all_vals), maximum(all_vals)
        margin = 0.1 * max(hi - lo, 1e-6)
        lo -= margin; hi += margin

        p = plot(; xlabel="Measured", ylabel= fi == 1 ? "Predicted" : "",
                 title=feat,
                 legend= fi == 5 ? :topleft : false,
                 xlim=(lo, hi), ylim=(lo, hi),
                 guidefontsize=9, titlefontsize=10, tickfontsize=7,
                 background_color_inside=bg_color,
                 grid=false, framestyle=:box,
                 margin=4Plots.mm, aspect_ratio=1)
        plot!(p, [lo, hi], [lo, hi]; ls=:dash, lc=:gray50, lw=1.5, label="y=x")

        # Real patients: filled circles, by visit
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

        # Synthetic patients: open markers (same color, ring only), by visit
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

    cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
    p_all = plot(plots_list...; layout=(1, 5), size=(2200, 500),
                 plot_title="Level 4: BZ2012 Predicted vs Measured — $cond",
                 plot_titlefontsize=13,
                 margin=6Plots.mm, left_margin=10Plots.mm)
    savefig(p_all, joinpath(_PATH_TO_FIG, "validate_mechanistic_$(cond_tag).pdf"))
    savefig(p_all, joinpath(_PATH_TO_FIG, "validate_mechanistic_$(cond_tag).png"))
    @info "  Saved validate_mechanistic_$(cond_tag).pdf"
end

# ── Plot 2: Pred/Meas ratio distributions (real vs synthetic) ────────────────
for cond in ["TF-only", "TF+TM"]
    csub = filter(r -> r.Condition == cond, results)

    plots_list = []
    for feat in FEATURE_LABELS
        fsub = filter(r -> r.Feature == feat, csub)
        valid = filter(r -> abs(r.Measured) > 1e-12, fsub)

        real_ratios  = filter(r -> r.Source == "Real", valid).Predicted ./
                       filter(r -> r.Source == "Real", valid).Measured
        synth_ratios = filter(r -> r.Source == "Synthetic", valid).Predicted ./
                       filter(r -> r.Source == "Synthetic", valid).Measured

        # clip extreme ratios for visualization
        real_ratios  = clamp.(real_ratios, 0.0, 5.0)
        synth_ratios = clamp.(synth_ratios, 0.0, 5.0)

        p = plot(; xlabel="Pred / Meas", ylabel="Density", title=feat,
                 legend=:topright, size=(400, 400),
                 guidefontsize=8, titlefontsize=9)

        if length(real_ratios) > 2
            histogram!(p, real_ratios; normalize=:pdf, bins=20,
                       fillalpha=0.4, fc=:steelblue, lc=:steelblue, label="Real")
        end
        if length(synth_ratios) > 2
            histogram!(p, synth_ratios; normalize=:pdf, bins=20,
                       fillalpha=0.3, fc=:coral, lc=:coral, label="Synthetic")
        end
        vline!(p, [1.0]; lc=:black, ls=:dash, lw=1.5, label="Perfect")

        push!(plots_list, p)
    end

    cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
    p_ratio = plot(plots_list...; layout=(1, 5), size=(2200, 450),
                   plot_title="Level 4: Pred/Meas Ratio Distributions — $cond",
                   margin=8Plots.mm)
    savefig(p_ratio, joinpath(_PATH_TO_FIG, "validate_mechanistic_ratios_$(cond_tag).pdf"))
    savefig(p_ratio, joinpath(_PATH_TO_FIG, "validate_mechanistic_ratios_$(cond_tag).png"))
    @info "  Saved validate_mechanistic_ratios_$(cond_tag).pdf"
end

# ── Plot 3: Rank correlation summary bar chart ───────────────────────────────
if nrow(rank_results) > 0
    for cond in ["TF-only", "TF+TM"]
        cond_ranks = filter(r -> r.Condition == cond, rank_results)
        if nrow(cond_ranks) == 0
            continue
        end

        # Aggregate: mean rank correlation per feature per source
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
            color=[:steelblue :coral],
            fillalpha=0.8,
            legend=:topright,
            size=(600, 400),
            ylim=(-0.5, 1.0),
        )
        hline!(p, [0.0]; lc=:gray, ls=:dash, lw=1, label="")

        cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
        savefig(p, joinpath(_PATH_TO_FIG, "validate_mechanistic_rankcorr_$(cond_tag).pdf"))
        savefig(p, joinpath(_PATH_TO_FIG, "validate_mechanistic_rankcorr_$(cond_tag).png"))
        @info "  Saved validate_mechanistic_rankcorr_$(cond_tag).pdf"
    end
end

# ── Plot 4: Visit generalization — does calibration transfer across visits? ──
for cond in ["TF-only", "TF+TM"]
    csub = filter(r -> r.Condition == cond && r.Source == "Real", results)

    plots_list = []
    for feat in FEATURE_LABELS
        fsub = filter(r -> r.Feature == feat, csub)
        valid = filter(r -> abs(r.Measured) > 1e-12, fsub)

        ratios_by_visit = [
            filter(r -> r.Visit == v, valid).Predicted ./ filter(r -> r.Visit == v, valid).Measured
            for v in 1:3
        ]

        p = boxplot(["V1\n(train)" "V2\n(test)" "V3\n(test)"],
                    hcat([length(r) > 0 ? r : [NaN] for r in ratios_by_visit]...);
                    legend=false, ylabel="Pred / Meas",
                    title=split(feat)[1],
                    fillalpha=0.6, fc=:steelblue,
                    guidefontsize=8, titlefontsize=9, size=(300, 400))
        hline!(p, [1.0]; lc=:red, ls=:dash, lw=1.5, label="")

        push!(plots_list, p)
    end

    cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
    p_gen = plot(plots_list...; layout=(1, 5), size=(1600, 400),
                 plot_title="Calibration Generalization — $cond (Real Patients)",
                 margin=8Plots.mm)
    savefig(p_gen, joinpath(_PATH_TO_FIG, "validate_mechanistic_generalization_$(cond_tag).pdf"))
    savefig(p_gen, joinpath(_PATH_TO_FIG, "validate_mechanistic_generalization_$(cond_tag).png"))
    @info "  Saved validate_mechanistic_generalization_$(cond_tag).pdf"
end

# ══════════════════════════════════════════════════════════════════════════════════
# Step 5: Save results
# ══════════════════════════════════════════════════════════════════════════════════
CSV.write(joinpath(_PATH_TO_DATA, "validate_mechanistic_results.csv"), results)
CSV.write(joinpath(_PATH_TO_DATA, "validate_mechanistic_rankcorr.csv"), rank_results)
@info "\nSaved results to data/validate_mechanistic_results.csv"
@info "Saved rank correlations to data/validate_mechanistic_rankcorr.csv"

# ── Final summary ────────────────────────────────────────────────────────────────
@info "\n═══════════════════════════════════════════════════════"
@info "  Level 4 Validation Complete"
@info "═══════════════════════════════════════════════════════"
@info ""
@info "  Interpretation:"
@info "    • This is a POPULATION-LEVEL plausibility check"
@info "    • Rank correlations are weak → model captures trends, not individual patients"
@info "    • Key result: do synthetic patients land in the SAME cloud as real patients?"
@info "    • If yes → SA-generated factor-to-TGA relationships are mechanistically plausible"
@info "    • TF-only generalizes across visits; TF+TM overcorrects PC pathway"
@info ""
@info "  Plots generated:"
@info "    • validate_mechanistic_tf_only.pdf — pred vs meas scatter"
@info "    • validate_mechanistic_tf_tm.pdf — pred vs meas scatter"
@info "    • validate_mechanistic_ratios_tf_only.pdf — ratio distributions"
@info "    • validate_mechanistic_ratios_tf_tm.pdf — ratio distributions"
@info "    • validate_mechanistic_rankcorr_tf_only.pdf — rank correlation bars"
@info "    • validate_mechanistic_rankcorr_tf_tm.pdf — rank correlation bars"
@info "    • validate_mechanistic_generalization_tf_only.pdf — visit transfer"
@info "    • validate_mechanistic_generalization_tf_tm.pdf — visit transfer"
@info ""
@info "Done!"

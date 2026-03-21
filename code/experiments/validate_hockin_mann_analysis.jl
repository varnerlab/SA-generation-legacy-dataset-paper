# ======================================================================================== #
# validate_hockin_mann_analysis.jl — Rank correlations + IIa-only thrombin readout
#
# Quick follow-up to validate_hockin_mann.jl:
#   1. Spearman rank correlations on existing results
#   2. Re-run with IIa-only (no mIIa) to test if meizothrombin explains overprediction
# ======================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots, StatsBase
using HockinMannModel

# ── paths ──────────────────────────────────────────────────────────────────────
const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: Rank correlations from existing results
# ══════════════════════════════════════════════════════════════════════════════
@info "═══ Part 1: Rank Correlations (IIa + 1.2·mIIa) ═══"

results = CSV.read(joinpath(_PATH_TO_DATA, "hm_validation_sample.csv"), DataFrame)

for cond in ["TF-only", "TF+TM"]
    @info "\n── $cond ──"
    sub = filter(r -> r.Condition == cond, results)
    for feat in unique(sub.Feature)
        fsub = filter(r -> r.Feature == feat, sub)
        if nrow(fsub) >= 3
            ρ = corspearman(fsub.Measured, fsub.Predicted)
            @info "  $feat: Spearman ρ = $(round(ρ, digits=3)) (n=$(nrow(fsub)))"
        end
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: IIa-only thrombin — does removing mIIa fix the overprediction?
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Part 2: IIa-only vs IIa+1.2·mIIa ═══"

# Column mappings (same as main script)
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

"""
Extract TGA features using IIa-only (species index 7, no mIIa contribution).
"""
function extract_tga_iia_only(sol; thrombin_threshold::Float64=10e-9)
    T = sol.t
    # IIa only — species 7 for BZ2012
    thrombin = sol[7, :]
    n = length(T)

    peak_idx = argmax(thrombin)
    peak = thrombin[peak_idx]
    tpeak = T[peak_idx]

    lag_idx = findfirst(x -> x > thrombin_threshold, thrombin)
    lagtime = isnothing(lag_idx) ? T[end] : T[lag_idx]

    max_rate = 0.0
    for i in 2:n
        dt = T[i] - T[i-1]
        if dt > 0
            rate = (thrombin[i] - thrombin[i-1]) / dt
            rate > max_rate && (max_rate = rate)
        end
    end

    etp = 0.0
    for i in 2:n
        dt = T[i] - T[i-1]
        etp += 0.5 * (thrombin[i] + thrombin[i-1]) * dt
    end

    return (lagtime=lagtime, peak=peak, tpeak=tpeak, max_rate=max_rate, etp=etp)
end

function run_patient(row; TM_molar::Float64=0.0)
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

    sol = simulate(HockinMannBZ2012; u0=u0, tspan=(0.0, 1200.0), saveat=1.0)

    # both readouts
    feat_total = extract_tga_features(HockinMannBZ2012, sol)
    feat_iia   = extract_tga_iia_only(sol)

    return (sol=sol, total=feat_total, iia_only=feat_iia)
end

function to_clinical(f)
    return (
        lagtime  = (hasproperty(f, :lagtime) ? f.lagtime : f[:lagtime]) / 60.0,
        peak     = (hasproperty(f, :peak) ? f.peak : f[:peak]) * 1e9,
        tpeak    = (hasproperty(f, :tpeak) ? f.tpeak : f[:tpeak]) / 60.0,
        max_rate = (hasproperty(f, :max_rate) ? f.max_rate : f[:max_rate]) * 1e9 * 60.0,
        etp      = (hasproperty(f, :etp) ? f.etp : f[:etp]) * 1e9 / 60.0,
    )
end

# Load data
df_real  = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
df_real_v1  = filter(r -> r.Visit == 1, df_real)
df_synth_v1 = filter(r -> r.Visit == 1, df_synth)

n_test = min(5, nrow(df_real_v1), nrow(df_synth_v1))

# Collect results for both readouts
feature_names = [:lagtime, :peak, :tpeak, :max_rate, :etp]
feature_labels = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

comparison = DataFrame(
    Source     = String[],
    PatientID  = String[],
    Condition  = String[],
    Feature    = String[],
    Measured   = Float64[],
    Pred_Total = Float64[],  # IIa + 1.2·mIIa
    Pred_IIa   = Float64[],  # IIa only
)

for (label, df_subset, id_col) in [
    ("Real", df_real_v1, :SubjectID),
    ("Synthetic", df_synth_v1, :SyntheticID),
]
    for i in 1:n_test
        row = df_subset[i, :]
        pid = string(row[id_col])

        for (cond_label, tm_val, tga_cols) in [
            ("TF-only", 0.0, TF_TGA_COLS),
            ("TF+TM",   1e-9, TM_TGA_COLS),
        ]
            @info "  $label patient $pid, $cond_label …"
            result = run_patient(row; TM_molar=tm_val)

            pred_total = to_clinical(result.total)
            pred_iia   = to_clinical(result.iia_only)

            for (j, fname) in enumerate(feature_names)
                m  = Float64(row[getfield(tga_cols, fname)])
                pt = getfield(pred_total, fname)
                pi = getfield(pred_iia, fname)
                push!(comparison, (label, pid, cond_label, feature_labels[j], m, pt, pi))
            end
        end
    end
end

# ── Print comparison table ─────────────────────────────────────────────────────
@info "\n═══ IIa+1.2·mIIa vs IIa-only: Side-by-side ═══"
for cond in ["TF-only", "TF+TM"]
    @info "\n── $cond ──"
    sub = filter(r -> r.Condition == cond, comparison)

    for feat in feature_labels
        fsub = filter(r -> r.Feature == feat, sub)
        re_total = mean(abs.(fsub.Pred_Total .- fsub.Measured) ./ max.(abs.(fsub.Measured), 1e-12))
        re_iia   = mean(abs.(fsub.Pred_IIa .- fsub.Measured) ./ max.(abs.(fsub.Measured), 1e-12))

        ratio_total = mean(fsub.Pred_Total ./ max.(fsub.Measured, 1e-12))
        ratio_iia   = mean(fsub.Pred_IIa ./ max.(fsub.Measured, 1e-12))

        @info "  $feat:"
        @info "    IIa+1.2·mIIa: mean pred/meas = $(round(ratio_total, digits=2))x, mean RE = $(round(re_total, digits=3))"
        @info "    IIa-only:     mean pred/meas = $(round(ratio_iia, digits=2))x, mean RE = $(round(re_iia, digits=3))"
    end

    # Rank correlations for both
    @info "\n  Spearman rank correlations:"
    for feat in feature_labels
        fsub = filter(r -> r.Feature == feat, sub)
        if nrow(fsub) >= 3
            ρ_total = corspearman(fsub.Measured, fsub.Pred_Total)
            ρ_iia   = corspearman(fsub.Measured, fsub.Pred_IIa)
            @info "    $feat: ρ(total)=$(round(ρ_total, digits=3)), ρ(IIa-only)=$(round(ρ_iia, digits=3))"
        end
    end
end

# ── Scatter: IIa-only predicted vs measured ────────────────────────────────────
@info "\nGenerating IIa-only scatter plots …"
for cond in ["TF-only", "TF+TM"]
    sub = filter(r -> r.Condition == cond, comparison)

    plots_list = []
    for fl in feature_labels
        fsub = filter(r -> r.Feature == fl, sub)
        real_sub = filter(r -> r.Source == "Real", fsub)
        synth_sub = filter(r -> r.Source == "Synthetic", fsub)

        all_vals = vcat(fsub.Measured, fsub.Pred_IIa)
        lo, hi = minimum(all_vals), maximum(all_vals)
        margin = 0.1 * (hi - lo)
        lo -= margin; hi += margin

        p = plot(; xlabel="Measured", ylabel="Predicted (IIa only)", title=fl,
                 legend=:topleft, size=(400, 400), xlim=(lo, hi), ylim=(lo, hi))
        plot!(p, [lo, hi], [lo, hi]; ls=:dash, lc=:gray, label="y=x")
        scatter!(p, real_sub.Measured, real_sub.Pred_IIa;
                 label="Real", mc=:steelblue, ms=6, ma=0.8)
        scatter!(p, synth_sub.Measured, synth_sub.Pred_IIa;
                 label="Synthetic", mc=:coral, ms=6, ma=0.8, shape=:diamond)
        push!(plots_list, p)
    end

    p_all = plot(plots_list...; layout=(1, 5), size=(2000, 400),
                 plot_title="BZ2012 IIa-only — $cond", margin=5Plots.mm)

    cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
    savefig(p_all, joinpath(_PATH_TO_FIG, "hm_validation_iia_only_$(cond_tag).pdf"))
    savefig(p_all, joinpath(_PATH_TO_FIG, "hm_validation_iia_only_$(cond_tag).png"))
    @info "  Saved hm_validation_iia_only_$(cond_tag).pdf"
end

# ── Save ───────────────────────────────────────────────────────────────────────
CSV.write(joinpath(_PATH_TO_DATA, "hm_validation_iia_comparison.csv"), comparison)
@info "\nSaved to data/hm_validation_iia_comparison.csv"
@info "Done!"

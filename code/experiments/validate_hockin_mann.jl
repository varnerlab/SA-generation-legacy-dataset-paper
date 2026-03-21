# ======================================================================================== #
# validate_hockin_mann.jl — Independent mechanistic validation of synthetic patients
#
# Uses the Hockin/Mann BZ2012 model (58 species, 64 frozen literature rate constants)
# to predict TGA features from patient coagulation factor levels. Compares predicted
# features to the SA-generated TGA features as an independent plausibility check.
#
# Two conditions per patient:
#   1. TF-only  (TM = 0)   → compare to "TF Initiator" TGA columns
#   2. TF + TM  (TM > 0)   → compare to "TF + TM Initiator" TGA columns
# ======================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots, PrettyTables
using HockinMannModel

# ── paths ──────────────────────────────────────────────────────────────────────
const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ── Load data ──────────────────────────────────────────────────────────────────
@info "Loading real and synthetic data …"
df_real  = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

# ── Column mappings ────────────────────────────────────────────────────────────
# Coagulation factor columns (% nominal in the data)
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

# TGA feature columns — TF-only condition
const TF_TGA_COLS = (
    lagtime  = Symbol("TF Initiator Lagtime (min)"),
    peak     = Symbol("TF Initiator Peak (nM)"),
    tpeak    = Symbol("TF Initiator T.Peak (min)"),
    max_rate = Symbol("TF Initiator Max Rate (nM/min)"),
    etp      = Symbol("TF Initiator ETP (nM·min)"),
)

# TGA feature columns — TF + TM condition
const TM_TGA_COLS = (
    lagtime  = Symbol("TF + TM  Initiator Lagtime (min)"),
    peak     = Symbol("TF + TM  Initiator Peak (nM)"),
    tpeak    = Symbol("TF + TM  Initiator T.Peak (min)"),
    max_rate = Symbol("TF + TM  Initiator Max Rate (nM/min)"),
    etp      = Symbol("TF + TM  Initiator ETP (nM·min)"),
)

# ── Helper: run BZ2012 for one patient row ─────────────────────────────────────
"""
    predict_tga(row; TM_molar=0.0) → TGAFeatures

Given a DataFrame row with coagulation factors in % nominal, convert to molar,
run BZ2012, and return predicted TGA features.
"""
function predict_tga(row; TM_molar::Float64=0.0)

    # convert % nominal → molar
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

    # build patient-specific initial conditions
    u0 = patient_initial_conditions(HockinMannBZ2012;
        TF   = 5e-12,      # 5 pM TF trigger (standard TGA)
        TM   = TM_molar,
        factors...,
    )

    # simulate
    sol = simulate(HockinMannBZ2012; u0=u0, tspan=(0.0, 1200.0), saveat=1.0)

    # extract features
    return extract_tga_features(HockinMannBZ2012, sol)
end

"""
    tga_to_clinical_units(f::TGAFeatures) → NamedTuple

Convert TGA features from model units (M, seconds) to clinical units (nM, minutes).
"""
function tga_to_clinical_units(f::TGAFeatures)
    return (
        lagtime  = f.lagtime / 60.0,        # s → min
        peak     = f.peak * 1e9,             # M → nM
        tpeak    = f.tpeak / 60.0,           # s → min
        max_rate = f.max_rate * 1e9 * 60.0,  # M/s → nM/min
        etp      = f.etp * 1e9 / 60.0,       # M·s → nM·min
    )
end

"""
    get_measured_tga(row, cols) → NamedTuple

Extract measured TGA features from a DataFrame row.
"""
function get_measured_tga(row, cols)
    return (
        lagtime  = Float64(row[cols.lagtime]),
        peak     = Float64(row[cols.peak]),
        tpeak    = Float64(row[cols.tpeak]),
        max_rate = Float64(row[cols.max_rate]),
        etp      = Float64(row[cols.etp]),
    )
end

# ── Run validation on a few patients ───────────────────────────────────────────
@info "Running BZ2012 validation on a sample of patients …"

# Pick first 5 real patients at Visit 1, and first 5 synthetic at Visit 1
df_real_v1 = filter(r -> r.Visit == 1, df_real)
df_synth_v1 = filter(r -> r.Visit == 1, df_synth)

n_test = min(5, nrow(df_real_v1), nrow(df_synth_v1))

results = DataFrame(
    Source      = String[],
    PatientID   = String[],
    Condition   = String[],  # TF-only or TF+TM
    Feature     = String[],
    Measured    = Float64[],
    Predicted   = Float64[],
    RelError    = Float64[],
)

feature_names = [:lagtime, :peak, :tpeak, :max_rate, :etp]
feature_labels = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

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

            # predict
            pred_raw = predict_tga(row; TM_molar=tm_val)
            pred = tga_to_clinical_units(pred_raw)

            # measured
            meas = get_measured_tga(row, tga_cols)

            # store per-feature comparison
            for (j, fname) in enumerate(feature_names)
                m = getfield(meas, fname)
                p = getfield(pred, fname)
                re = abs(m) > 1e-12 ? abs(p - m) / abs(m) : NaN
                push!(results, (label, pid, cond_label, feature_labels[j], m, p, re))
            end
        end
    end
end

# ── Display results ────────────────────────────────────────────────────────────
@info "\n═══ Validation Results ═══"

for cond in ["TF-only", "TF+TM"]
    @info "\n── $cond condition ──"
    sub = filter(r -> r.Condition == cond, results)

    # pivot: one row per patient, columns = features
    for src in ["Real", "Synthetic"]
        @info "\n  $src patients:"
        ssub = filter(r -> r.Source == src, sub)
        pids = unique(ssub.PatientID)

        header = ["Patient", feature_labels..., "Mean RE"]
        data = Matrix{Any}(undef, length(pids), length(feature_labels) + 2)
        for (i, pid) in enumerate(pids)
            data[i, 1] = pid
            psub = filter(r -> r.PatientID == pid, ssub)
            res = Float64[]
            for (j, fl) in enumerate(feature_labels)
                row = filter(r -> r.Feature == fl, psub)
                m = row.Measured[1]
                p = row.Predicted[1]
                re = row.RelError[1]
                push!(res, re)
                data[i, j+1] = "$(round(m, digits=2)) → $(round(p, digits=2))"
            end
            data[i, end] = round(mean(filter(!isnan, res)), digits=3)
        end
        pretty_table(data; column_labels=header)
    end
end

# ── Summary statistics ─────────────────────────────────────────────────────────
@info "\n═══ Summary ═══"
for cond in ["TF-only", "TF+TM"]
    sub = filter(r -> r.Condition == cond, results)
    for src in ["Real", "Synthetic"]
        ssub = filter(r -> r.Source == src, sub)
        valid_re = filter(!isnan, ssub.RelError)
        @info "  $cond / $src: median RE = $(round(median(valid_re), digits=3)), mean RE = $(round(mean(valid_re), digits=3))"
    end
end

# ── Scatter plots: predicted vs measured ───────────────────────────────────────
@info "\nGenerating scatter plots …"
for cond in ["TF-only", "TF+TM"]
    sub = filter(r -> r.Condition == cond, results)

    plots_list = []
    for fl in feature_labels
        fsub = filter(r -> r.Feature == fl, sub)
        real_sub = filter(r -> r.Source == "Real", fsub)
        synth_sub = filter(r -> r.Source == "Synthetic", fsub)

        all_vals = vcat(fsub.Measured, fsub.Predicted)
        lo, hi = minimum(all_vals), maximum(all_vals)
        margin = 0.1 * (hi - lo)
        lo -= margin; hi += margin

        p = plot(; xlabel="Measured", ylabel="Predicted", title=fl,
                 legend=:topleft, size=(400, 400), xlim=(lo, hi), ylim=(lo, hi))
        plot!(p, [lo, hi], [lo, hi]; ls=:dash, lc=:gray, label="y=x")
        scatter!(p, real_sub.Measured, real_sub.Predicted;
                 label="Real", mc=:steelblue, ms=6, ma=0.8)
        scatter!(p, synth_sub.Measured, synth_sub.Predicted;
                 label="Synthetic", mc=:coral, ms=6, ma=0.8, shape=:diamond)
        push!(plots_list, p)
    end

    p_all = plot(plots_list...; layout=(1, 5), size=(2000, 400),
                 plot_title="BZ2012 Validation — $cond", margin=5Plots.mm)

    cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
    savefig(p_all, joinpath(_PATH_TO_FIG, "hm_validation_$(cond_tag).pdf"))
    savefig(p_all, joinpath(_PATH_TO_FIG, "hm_validation_$(cond_tag).png"))
    @info "  Saved hm_validation_$(cond_tag).pdf"
end

# ── Save results table ─────────────────────────────────────────────────────────
CSV.write(joinpath(_PATH_TO_DATA, "hm_validation_sample.csv"), results)
@info "\nSaved results to data/hm_validation_sample.csv"
@info "Done!"

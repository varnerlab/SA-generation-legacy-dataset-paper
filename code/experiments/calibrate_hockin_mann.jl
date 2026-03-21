# ======================================================================================== #
# calibrate_hockin_mann.jl — Global calibration of BZ2012 on Visit 1 real patients
#
# Fits a small number of rate constants (population-level, NOT per-patient) to minimize
# TGA prediction error across all 23 complete real patients at Visit 1.
#
# Train: Visit 1 real patients (n=23)
# Test:  Visit 2+3 real patients (held-out timepoints)
#        All visits synthetic patients
# ======================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots, Optim, PrettyTables
using HockinMannModel

# ── paths ──────────────────────────────────────────────────────────────────────
const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ── Column mappings ────────────────────────────────────────────────────────────
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

# ── Parameters to calibrate ───────────────────────────────────────────────────
# These are the rate constants most likely responsible for the systematic
# overprediction of peak and max_rate:
#   - prothrombinase_kcat (p[31]): THE key catalytic rate for thrombin burst
#   - intrinsic_xase_kcat (p[22]): amplification-phase Xa generation
#   - extrinsic_xase_kcat (p[10]): initiation-phase Xa generation
#   - PC_activation_kcat (p[47]): protein C activation (TF+TM condition)
#   - mIIa_conversion_k (p[32]): meizothrombin → thrombin conversion

const CALIBRATION_PARAMS = [
    :prothrombinase_kcat,    # default 63.5 s⁻¹
    :intrinsic_xase_kcat,    # default 8.2 s⁻¹
    :extrinsic_xase_kcat,    # default 6.0 s⁻¹
    :PC_activation_kcat,     # default 0.41 s⁻¹
    :mIIa_conversion_k,      # default 2.3e8 M⁻¹s⁻¹
]

# Get default values and indices
const PARAM_INDICES = [rate_index(HockinMannBZ2012, name) for name in CALIBRATION_PARAMS]
const PARAM_DEFAULTS = [default_rate_constants(HockinMannBZ2012)[i] for i in PARAM_INDICES]

@info "Calibrating $(length(CALIBRATION_PARAMS)) rate constants:"
for (name, idx, val) in zip(CALIBRATION_PARAMS, PARAM_INDICES, PARAM_DEFAULTS)
    @info "  $name (p[$idx]) = $val"
end

# ══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════════════════════

function run_patient_with_params(row, p_vec; TM_molar::Float64=0.0)
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
        lagtime  = f.lagtime / 60.0,
        peak     = f.peak * 1e9,
        tpeak    = f.tpeak / 60.0,
        max_rate = f.max_rate * 1e9 * 60.0,
        etp      = f.etp * 1e9 / 60.0,
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

"""
Compute normalized relative squared error between predicted and measured TGA features.
Each feature is weighted equally (relative error).
"""
function patient_error(pred, meas)
    err = 0.0
    n = 0
    for fname in (:lagtime, :peak, :tpeak, :max_rate, :etp)
        m = getfield(meas, fname)
        p = getfield(pred, fname)
        if abs(m) > 1e-12 && isfinite(p)
            err += ((p - m) / m)^2
            n += 1
        end
    end
    return n > 0 ? err / n : Inf
end

# ══════════════════════════════════════════════════════════════════════════════
# Load data and identify training set
# ══════════════════════════════════════════════════════════════════════════════
@info "\nLoading data …"
df_real  = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

# Complete patients (all 3 visits)
all_subjects = unique(df_real.SubjectID)
complete_subjects = [s for s in all_subjects
    if sum(df_real.SubjectID .== s) == 3 &&
       all(v -> any((df_real.SubjectID .== s) .& (df_real.Visit .== v)), 1:3)]

@info "  Complete patients: $(length(complete_subjects))"

# Training set: Visit 1 of complete patients
df_train = filter(r -> r.Visit == 1 && r.SubjectID in complete_subjects, df_real)
@info "  Training rows (Visit 1): $(nrow(df_train))"

# ══════════════════════════════════════════════════════════════════════════════
# Loss function: average error across all training patients, both conditions
# ══════════════════════════════════════════════════════════════════════════════

function total_loss(log_scale_factors::Vector{Float64})
    # Parameters are optimized in log-scale-factor space:
    #   p[i] = default[i] * exp(log_scale_factors[i])
    # This ensures positivity and normalizes the search space.
    p_vec = default_rate_constants(HockinMannBZ2012)
    for (j, idx) in enumerate(PARAM_INDICES)
        p_vec[idx] = PARAM_DEFAULTS[j] * exp(log_scale_factors[j])
    end

    total_err = 0.0
    n_patients = 0

    for i in 1:nrow(df_train)
        row = df_train[i, :]

        for (tm_val, tga_cols) in [(0.0, TF_TGA_COLS), (1e-9, TM_TGA_COLS)]
            try
                feat = run_patient_with_params(row, p_vec; TM_molar=tm_val)
                pred = to_clinical(feat)
                meas = get_measured(row, tga_cols)
                err = patient_error(pred, meas)
                if isfinite(err)
                    total_err += err
                    n_patients += 1
                end
            catch e
                # ODE solver failure → penalize
                total_err += 100.0
                n_patients += 1
            end
        end
    end

    return n_patients > 0 ? total_err / n_patients : Inf
end

# ══════════════════════════════════════════════════════════════════════════════
# Optimization
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Starting calibration ═══"
@info "  Initial loss (literature defaults): $(round(total_loss(zeros(length(CALIBRATION_PARAMS))), digits=4))"

# Optimize in log-scale-factor space (start at 0 = default values)
x0 = zeros(length(CALIBRATION_PARAMS))

# Bounds: allow each parameter to vary by up to 100x in either direction
lower = fill(-log(100.0), length(CALIBRATION_PARAMS))
upper = fill(log(100.0), length(CALIBRATION_PARAMS))

result = optimize(
    total_loss,
    lower, upper, x0,
    Fminbox(NelderMead()),
    Optim.Options(
        iterations = 500,
        show_trace = true,
        show_every = 50,
        f_tol = 1e-6,
    ),
)

@info "\n═══ Calibration Results ═══"
@info "  Converged: $(Optim.converged(result))"
@info "  Final loss: $(round(Optim.minimum(result), digits=4))"
@info "  Iterations: $(Optim.iterations(result))"

# Extract calibrated parameters
log_sf = Optim.minimizer(result)
p_calibrated = default_rate_constants(HockinMannBZ2012)
@info "\n  Calibrated rate constants:"
for (j, (name, idx)) in enumerate(zip(CALIBRATION_PARAMS, PARAM_INDICES))
    scale = exp(log_sf[j])
    p_calibrated[idx] = PARAM_DEFAULTS[j] * scale
    @info "    $name: $(PARAM_DEFAULTS[j]) → $(round(p_calibrated[idx], sigdigits=4)) ($(round(scale, digits=3))x)"
end

# ══════════════════════════════════════════════════════════════════════════════
# Evaluate on ALL data: train (V1) + test (V2, V3) real + synthetic
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Evaluation ═══"

feature_names = [:lagtime, :peak, :tpeak, :max_rate, :etp]
feature_labels = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

results_all = DataFrame(
    Source    = String[],
    PatientID = String[],
    Visit     = Int[],
    Condition = String[],
    Feature   = String[],
    Measured  = Float64[],
    Pred_Lit  = Float64[],   # literature (frozen) parameters
    Pred_Cal  = Float64[],   # calibrated parameters
)

p_literature = default_rate_constants(HockinMannBZ2012)

function evaluate_dataset!(results_all, df, source, id_col, p_lit, p_cal)
    for i in 1:nrow(df)
        row = df[i, :]
        pid = string(row[id_col])
        visit = Int(row[:Visit])

        for (cond_label, tm_val, tga_cols) in [
            ("TF-only", 0.0, TF_TGA_COLS),
            ("TF+TM",   1e-9, TM_TGA_COLS),
        ]
            try
                feat_lit = run_patient_with_params(row, p_lit; TM_molar=tm_val)
                feat_cal = run_patient_with_params(row, p_cal; TM_molar=tm_val)
                pred_lit = to_clinical(feat_lit)
                pred_cal = to_clinical(feat_cal)
                meas = get_measured(row, tga_cols)

                for (j, fname) in enumerate(feature_names)
                    m  = getfield(meas, fname)
                    pl = getfield(pred_lit, fname)
                    pc = getfield(pred_cal, fname)
                    push!(results_all, (source, pid, visit, cond_label, feature_labels[j], m, pl, pc))
                end
            catch e
                @warn "  Failed: $source patient $pid visit $visit $cond_label: $e"
            end
        end
    end
end

# Real patients (complete only)
df_real_complete = filter(r -> r.SubjectID in complete_subjects, df_real)
@info "  Evaluating $(nrow(df_real_complete)) real patient records …"
evaluate_dataset!(results_all, df_real_complete, "Real", :SubjectID, p_literature, p_calibrated)

# Synthetic patients (all visits)
@info "  Evaluating $(nrow(df_synth)) synthetic patient records …"
evaluate_dataset!(results_all, df_synth, "Synthetic", :SyntheticID, p_literature, p_calibrated)

# ══════════════════════════════════════════════════════════════════════════════
# Summary statistics by visit and condition
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Summary: Median Relative Error ═══"
@info "  (Literature → Calibrated)"

for source in ["Real", "Synthetic"]
    @info "\n  ── $source ──"
    ssub = filter(r -> r.Source == source, results_all)

    for cond in ["TF-only", "TF+TM"]
        csub = filter(r -> r.Condition == cond, ssub)

        for visit in sort(unique(csub.Visit))
            vsub = filter(r -> r.Visit == visit, csub)
            label = visit == 1 ? "V$visit (TRAIN)" : "V$visit (TEST)"

            re_lit = [abs(r.Pred_Lit - r.Measured) / max(abs(r.Measured), 1e-12) for r in eachrow(vsub)]
            re_cal = [abs(r.Pred_Cal - r.Measured) / max(abs(r.Measured), 1e-12) for r in eachrow(vsub)]

            @info "    $cond $label: median RE lit=$(round(median(re_lit), digits=3)) → cal=$(round(median(re_cal), digits=3))"
        end
    end
end

# ── Per-feature breakdown for real patients ────────────────────────────────────
@info "\n═══ Per-Feature: Real Patients (pred/meas ratio) ═══"
for cond in ["TF-only", "TF+TM"]
    @info "\n  ── $cond ──"
    csub = filter(r -> r.Source == "Real" && r.Condition == cond, results_all)

    for visit in 1:3
        vsub = filter(r -> r.Visit == visit, csub)
        label = visit == 1 ? "V$visit (TRAIN)" : "V$visit (TEST)"

        @info "    $label:"
        for feat in feature_labels
            fsub = filter(r -> r.Feature == feat, vsub)
            ratio_lit = mean(fsub.Pred_Lit ./ max.(abs.(fsub.Measured), 1e-12))
            ratio_cal = mean(fsub.Pred_Cal ./ max.(abs.(fsub.Measured), 1e-12))
            @info "      $feat: lit=$(round(ratio_lit, digits=2))x → cal=$(round(ratio_cal, digits=2))x"
        end
    end
end

# ── Rank correlations ──────────────────────────────────────────────────────────
@info "\n═══ Spearman Rank Correlations (Calibrated) ═══"
using StatsBase

for source in ["Real", "Synthetic"]
    @info "\n  ── $source ──"
    ssub = filter(r -> r.Source == source, results_all)

    for cond in ["TF-only", "TF+TM"]
        csub = filter(r -> r.Condition == cond, ssub)

        for visit in sort(unique(csub.Visit))
            vsub = filter(r -> r.Visit == visit, csub)
            label = visit == 1 ? "V$visit (TRAIN)" : "V$visit (TEST)"

            parts = String[]
            for feat in feature_labels
                fsub = filter(r -> r.Feature == feat, vsub)
                if nrow(fsub) >= 3
                    ρ = corspearman(fsub.Measured, fsub.Pred_Cal)
                    push!(parts, "$(split(feat)[1])=$(round(ρ, digits=2))")
                end
            end
            @info "    $cond $label: $(join(parts, ", "))"
        end
    end
end

# ── Save ───────────────────────────────────────────────────────────────────────
CSV.write(joinpath(_PATH_TO_DATA, "hm_calibration_results.csv"), results_all)
@info "\nSaved results to data/hm_calibration_results.csv"

# ── Scatter plots: calibrated predicted vs measured ────────────────────────────
@info "\nGenerating scatter plots …"
for cond in ["TF-only", "TF+TM"]
    for source in ["Real"]
        sub = filter(r -> r.Condition == cond && r.Source == source, results_all)

        plots_list = []
        for fl in feature_labels
            fsub = filter(r -> r.Feature == fl, sub)

            v1 = filter(r -> r.Visit == 1, fsub)
            v2 = filter(r -> r.Visit == 2, fsub)
            v3 = filter(r -> r.Visit == 3, fsub)

            all_vals = vcat(fsub.Measured, fsub.Pred_Cal)
            lo, hi = minimum(all_vals), maximum(all_vals)
            margin = 0.1 * (hi - lo)
            lo -= margin; hi += margin

            p = plot(; xlabel="Measured", ylabel="Predicted (cal)", title=fl,
                     legend=:topleft, size=(400, 400), xlim=(lo, hi), ylim=(lo, hi))
            plot!(p, [lo, hi], [lo, hi]; ls=:dash, lc=:gray, label="y=x")
            scatter!(p, v1.Measured, v1.Pred_Cal; label="V1 (train)", mc=:steelblue, ms=5, ma=0.7)
            scatter!(p, v2.Measured, v2.Pred_Cal; label="V2 (test)", mc=:coral, ms=5, ma=0.7, shape=:diamond)
            scatter!(p, v3.Measured, v3.Pred_Cal; label="V3 (test)", mc=:green, ms=5, ma=0.7, shape=:utriangle)
            push!(plots_list, p)
        end

        cond_tag = cond == "TF-only" ? "tf_only" : "tf_tm"
        p_all = plot(plots_list...; layout=(1, 5), size=(2000, 400),
                     plot_title="BZ2012 Calibrated — $cond ($source)", margin=5Plots.mm)
        savefig(p_all, joinpath(_PATH_TO_FIG, "hm_calibrated_$(cond_tag)_$(lowercase(source)).pdf"))
        savefig(p_all, joinpath(_PATH_TO_FIG, "hm_calibrated_$(cond_tag)_$(lowercase(source)).png"))
        @info "  Saved hm_calibrated_$(cond_tag)_$(lowercase(source)).pdf"
    end
end

@info "\nDone!"

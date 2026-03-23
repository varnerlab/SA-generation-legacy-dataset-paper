# ======================================================================================== #
# downstream_utility_calibration.jl — R1.4: Downstream Utility Validation
#
# Can synthetic data substitute for real data when calibrating a mechanistic model?
#
# Experiment:
#   1. Calibrate BZ2012 on SYNTHETIC V1 patients (N=100)
#   2. Calibrate BZ2012 on REAL V1 patients (K=23) [already done, use saved values]
#   3. Evaluate BOTH calibrations on REAL V2 + V3 patients (held-out)
#   4. Compare: does the synth-calibrated model predict real patients comparably?
# ======================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using CSV, DataFrames, Statistics, Plots, Optim, LinearAlgebra, StatsBase, Random
using HockinMannModel

const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ── Column mappings ──────────────────────────────────────────────────────────
const FACTOR_COLS = Dict(
    :II   => :II, :V => :V, :VII => :VII, :VIII => :VIII,
    :IX   => :IX, :X => :X, :AT  => :AT, :PC => :PC,
    :TFPI => Symbol("TFPI Free"),
)

const TF_TGA_COLS = (
    lagtime  = Symbol("TF Initiator Lagtime (min)"),
    peak     = Symbol("TF Initiator Peak (nM)"),
    tpeak    = Symbol("TF Initiator T.Peak (min)"),
    max_rate = Symbol("TF Initiator Max Rate (nM/min)"),
    etp      = Symbol("TF Initiator ETP (nM·min)"),
)

const FEATURE_NAMES = [:lagtime, :peak, :tpeak, :max_rate, :etp]
const FEATURE_LABELS = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

# ── Helpers ──────────────────────────────────────────────────────────────────
function run_patient(row, p_vec; TM_molar::Float64=0.0)
    factors = percent_nominal_to_molar(
        II=Float64(row[FACTOR_COLS[:II]]), V=Float64(row[FACTOR_COLS[:V]]),
        VII=Float64(row[FACTOR_COLS[:VII]]), VIII=Float64(row[FACTOR_COLS[:VIII]]),
        IX=Float64(row[FACTOR_COLS[:IX]]), X=Float64(row[FACTOR_COLS[:X]]),
        AT=Float64(row[FACTOR_COLS[:AT]]), PC=Float64(row[FACTOR_COLS[:PC]]),
        TFPI=Float64(row[FACTOR_COLS[:TFPI]]),
    )
    u0 = patient_initial_conditions(HockinMannBZ2012; TF=5e-12, TM=TM_molar, factors...)
    sol = simulate(HockinMannBZ2012; u0=u0, p=p_vec, tspan=(0.0, 1200.0), saveat=1.0)
    return extract_tga_features(HockinMannBZ2012, sol)
end

function to_clinical(f::TGAFeatures)
    return (lagtime=f.lagtime/60.0, peak=f.peak*1e9, tpeak=f.tpeak/60.0,
            max_rate=f.max_rate*1e9*60.0, etp=f.etp*1e9/60.0)
end

function get_measured(row, cols)
    return (lagtime=Float64(row[cols.lagtime]), peak=Float64(row[cols.peak]),
            tpeak=Float64(row[cols.tpeak]), max_rate=Float64(row[cols.max_rate]),
            etp=Float64(row[cols.etp]))
end

function patient_error(pred, meas)
    err = 0.0; n = 0
    for fname in FEATURE_NAMES
        m = getfield(meas, fname); p = getfield(pred, fname)
        if abs(m) > 1e-12 && isfinite(p)
            err += ((p - m) / m)^2; n += 1
        end
    end
    return n > 0 ? err / n : Inf
end

# ── Calibration parameters ───────────────────────────────────────────────────
const CALIBRATION_PARAMS = [
    :prothrombinase_kcat, :intrinsic_xase_kcat, :extrinsic_xase_kcat,
    :PC_activation_kcat, :mIIa_conversion_k,
]
const PARAM_INDICES = [rate_index(HockinMannBZ2012, name) for name in CALIBRATION_PARAMS]
const PARAM_DEFAULTS = [default_rate_constants(HockinMannBZ2012)[i] for i in PARAM_INDICES]

# ── Load data ────────────────────────────────────────────────────────────────
@info "Loading data …"
df_real = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

# Complete real patients
all_s = unique(df_real.SubjectID)
complete_subjects = [s for s in all_s
    if sum(df_real.SubjectID .== s) == 3 &&
       all(v -> any((df_real.SubjectID .== s) .& (df_real.Visit .== v)), 1:3)]

df_real_complete = filter(r -> r.SubjectID in complete_subjects, df_real)

# Training sets
df_real_v1 = filter(r -> r.Visit == 1, df_real_complete)
df_synth_v1 = filter(r -> r.Visit == 1, df_synth)

# Test set: real V2 + V3
df_real_test = filter(r -> r.Visit in [2, 3], df_real_complete)

@info "  Real V1 (train): $(nrow(df_real_v1)) patients"
@info "  Synth V1 (train): $(nrow(df_synth_v1)) patients"
@info "  Real V2+V3 (test): $(nrow(df_real_test)) records"

# ══════════════════════════════════════════════════════════════════════════════
# Calibration on REAL V1 (reference — already known)
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Calibrating on REAL V1 patients ═══"

function make_loss(df_train)
    return function(log_sf)
        p_vec = default_rate_constants(HockinMannBZ2012)
        for (j, idx) in enumerate(PARAM_INDICES)
            p_vec[idx] = PARAM_DEFAULTS[j] * exp(log_sf[j])
        end
        total_err = 0.0; n = 0
        for i in 1:nrow(df_train)
            row = df_train[i, :]
            for (tm_val, tga_cols) in [(0.0, TF_TGA_COLS)]  # TF-only
                try
                    feat = run_patient(row, p_vec; TM_molar=tm_val)
                    pred = to_clinical(feat)
                    meas = get_measured(row, tga_cols)
                    err = patient_error(pred, meas)
                    if isfinite(err); total_err += err; n += 1; end
                catch; total_err += 100.0; n += 1; end
            end
        end
        return n > 0 ? total_err / n : Inf
    end
end

lower = fill(-log(100.0), length(CALIBRATION_PARAMS))
upper = fill(log(100.0), length(CALIBRATION_PARAMS))
n_dim = length(CALIBRATION_PARAMS)

# Multi-start optimization: run from N_STARTS random initial points + literature default
N_STARTS = 10
Random.seed!(42)

function multi_start_calibrate(loss_fn, label)
    @info "\n═══ Calibrating on $label ═══"

    # Generate starting points: x0 = 0 (literature) + N_STARTS random
    starts = [zeros(n_dim)]  # literature default
    for _ in 1:N_STARTS
        x0_rand = lower .+ (upper .- lower) .* rand(n_dim)
        push!(starts, x0_rand)
    end

    best_loss = Inf
    best_result = nothing

    for (si, x0) in enumerate(starts)
        label_s = si == 1 ? "literature" : "random-$si"
        try
            result = optimize(loss_fn, lower, upper, x0,
                Fminbox(NelderMead()), Optim.Options(iterations=500, f_tol=1e-6, show_trace=false))
            fl = Optim.minimum(result)
            conv = Optim.converged(result)
            @info "  Start $label_s: loss=$(round(fl, digits=4)), converged=$conv"
            if fl < best_loss
                best_loss = fl
                best_result = result
            end
        catch e
            @warn "  Start $label_s: FAILED ($e)"
        end
    end

    @info "  Best loss: $(round(best_loss, digits=4))"

    p_cal = default_rate_constants(HockinMannBZ2012)
    for (j, idx) in enumerate(PARAM_INDICES)
        p_cal[idx] = PARAM_DEFAULTS[j] * exp(Optim.minimizer(best_result)[j])
    end

    for (j, name) in enumerate(CALIBRATION_PARAMS)
        @info "    $name: $(PARAM_DEFAULTS[j]) → $(round(p_cal[PARAM_INDICES[j]], sigdigits=4))"
    end

    return p_cal, best_loss
end

loss_real = make_loss(df_real_v1)
p_cal_real, loss_val_real = multi_start_calibrate(loss_real, "REAL V1 ($( nrow(df_real_v1)) patients)")

loss_synth = make_loss(df_synth_v1)
p_cal_synth, loss_val_synth = multi_start_calibrate(loss_synth, "SYNTHETIC V1 ($(nrow(df_synth_v1)) patients)")

# ══════════════════════════════════════════════════════════════════════════════
# Evaluate BOTH calibrations on REAL V2+V3 (held-out test set)
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Evaluating on REAL V2+V3 (held-out) ═══"

results = DataFrame(
    Visit=Int[], Feature=String[],
    Measured=Float64[], Pred_Real=Float64[], Pred_Synth=Float64[],
)

for i in 1:nrow(df_real_test)
    row = df_real_test[i, :]
    visit = Int(row[:Visit])
    try
        feat_real = run_patient(row, p_cal_real; TM_molar=0.0)
        feat_synth = run_patient(row, p_cal_synth; TM_molar=0.0)
        pred_r = to_clinical(feat_real)
        pred_s = to_clinical(feat_synth)
        meas = get_measured(row, TF_TGA_COLS)
        for (j, fname) in enumerate(FEATURE_NAMES)
            m = getfield(meas, fname)
            pr = getfield(pred_r, fname)
            ps = getfield(pred_s, fname)
            push!(results, (visit, FEATURE_LABELS[j], m, pr, ps))
        end
    catch e
        @warn "  Failed: visit $visit: $e"
    end
end

@info "  $(nrow(results)) predictions on held-out data"

# ── Compare calibrations ─────────────────────────────────────────────────────
@info "\n═══ Head-to-head: Real-calibrated vs Synth-calibrated on held-out V2+V3 ═══"
@info "  $(rpad("Feature", 22))  Real-cal MdRE  Synth-cal MdRE  Ratio"

for feat in FEATURE_LABELS
    fsub = filter(r -> r.Feature == feat, results)
    valid = filter(r -> abs(r.Measured) > 1e-12, fsub)

    re_real = abs.(valid.Pred_Real .- valid.Measured) ./ abs.(valid.Measured)
    re_synth = abs.(valid.Pred_Synth .- valid.Measured) ./ abs.(valid.Measured)

    md_real = median(re_real)
    md_synth = median(re_synth)
    ratio = md_synth / max(md_real, 1e-12)

    @info "  $(rpad(feat, 22))  $(round(md_real, digits=3))          $(round(md_synth, digits=3))           $(round(ratio, digits=2))x"
end

# Overall
all_re_real = Float64[]
all_re_synth = Float64[]
for row in eachrow(results)
    if abs(row.Measured) > 1e-12
        push!(all_re_real, abs(row.Pred_Real - row.Measured) / abs(row.Measured))
        push!(all_re_synth, abs(row.Pred_Synth - row.Measured) / abs(row.Measured))
    end
end

@info "\n  Overall median RE on held-out V2+V3:"
@info "    Real-calibrated:  $(round(median(all_re_real), digits=3))"
@info "    Synth-calibrated: $(round(median(all_re_synth), digits=3))"
@info "    Ratio: $(round(median(all_re_synth) / median(all_re_real), digits=2))x"

# ── Save results ─────────────────────────────────────────────────────────────
CSV.write(joinpath(_PATH_TO_DATA, "downstream_utility_results.csv"), results)
@info "\n  Saved downstream_utility_results.csv"

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
             title=feat, titlefontsize=10, guidefontsize=9, tickfontsize=7,
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

@info "\nDone!"

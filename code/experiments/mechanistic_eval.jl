# ──────────────────────────────────────────────────────────────────────────────
# mechanistic_eval.jl — includable BZ2012 mechanistic-evaluation helper
#
# Faithful extraction (not a rewrite) of the calibrated Hockin/Mann BZ2012
# "run patient through the ODE, compare to measured TGA, compute cloud
# overlap + KS" harness that previously lived only as un-importable
# script-level code, so that E1 (PCA ablation), E2 (out-of-hull plausibility),
# and SIM arm 2 do not each re-implement it.
#
# Lifted from the FROZEN reference scripts (neither is modified, and this file
# must reproduce their published TF-only numbers: cloud overlap 0.86-0.93,
# KS D 0.079-0.127, all p > 0.30):
#   - validate_mechanistic_plausibility.jl
#       CALIBRATED_OVERRIDES ............ lines ~79-85 (copied verbatim)
#       FACTOR_COLS / TF_TGA_COLS / TM_TGA_COLS ... lines ~39-68 (copied verbatim)
#       run_patient / to_clinical / get_measured ... lines ~96-134
#       evaluate_all! per-row/per-feature loop body ... lines ~170-198
#   - regen_mechanistic_figures.jl
#       cloud overlap + KS computation ... lines ~91-109
#
# NOT run standalone. The including script is expected to have already run
#   using HockinMannModel, HypothesisTests, DataFrames, Statistics
# (and CSV, if it loads its own data) before
#   include(joinpath(@__DIR__, "mechanistic_eval.jl"))
# ──────────────────────────────────────────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════════
# Calibrated rate constants — copied verbatim from
# validate_mechanistic_plausibility.jl:79-85 (fit on Visit-1 real patients by
# calibrate_hockin_mann.jl; 5 of 64 rate constants adjusted, 59 at literature values)
# ══════════════════════════════════════════════════════════════════════════════════
const CALIBRATED_OVERRIDES = Dict(
    :prothrombinase_kcat => 21.94,      # default 63.5 s⁻¹  (0.345x)
    :intrinsic_xase_kcat => 0.1693,     # default 8.2 s⁻¹   (0.021x)
    :extrinsic_xase_kcat => 93.21,      # default 6.0 s⁻¹   (15.5x)
    :PC_activation_kcat  => 0.8896,     # default 0.41 s⁻¹   (2.17x)
    :mIIa_conversion_k   => 1.153e9,    # default 2.3e8 M⁻¹s⁻¹ (5.01x)
)

const _P_CALIBRATED = make_rate_constants(HockinMannBZ2012; CALIBRATED_OVERRIDES...)

# ── Column mappings — copied verbatim from validate_mechanistic_plausibility.jl:39-68 ──
const _FACTOR_COLS = Dict(
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

const _FEATURE_NAMES  = [:lagtime, :peak, :tpeak, :max_rate, :etp]
const _FEATURE_LABELS = ["Lagtime (min)", "Peak (nM)", "T.Peak (min)", "Max Rate (nM/min)", "ETP (nM·min)"]

# ══════════════════════════════════════════════════════════════════════════════════
# run_patient / to_clinical / get_measured — copied verbatim (logic-for-logic)
# from validate_mechanistic_plausibility.jl:96-134
# ══════════════════════════════════════════════════════════════════════════════════

function _run_patient(row, p_vec; TM_molar::Float64=0.0)
    factors = percent_nominal_to_molar(
        II   = Float64(row[_FACTOR_COLS[:II]]),
        V    = Float64(row[_FACTOR_COLS[:V]]),
        VII  = Float64(row[_FACTOR_COLS[:VII]]),
        VIII = Float64(row[_FACTOR_COLS[:VIII]]),
        IX   = Float64(row[_FACTOR_COLS[:IX]]),
        X    = Float64(row[_FACTOR_COLS[:X]]),
        AT   = Float64(row[_FACTOR_COLS[:AT]]),
        PC   = Float64(row[_FACTOR_COLS[:PC]]),
        TFPI = Float64(row[_FACTOR_COLS[:TFPI]]),
    )

    u0 = patient_initial_conditions(HockinMannBZ2012;
        TF = 5e-12, TM = TM_molar, factors...)

    sol = simulate(HockinMannBZ2012; u0=u0, p=p_vec, tspan=(0.0, 1200.0), saveat=1.0)
    return extract_tga_features(HockinMannBZ2012, sol)
end

function _to_clinical(f::TGAFeatures)
    return (
        lagtime  = f.lagtime / 60.0,        # s → min
        peak     = f.peak * 1e9,             # M → nM
        tpeak    = f.tpeak / 60.0,           # s → min
        max_rate = f.max_rate * 1e9 * 60.0,  # M/s → nM/min
        etp      = f.etp * 1e9 / 60.0,       # M·s → nM·min
    )
end

function _get_measured(row, cols)
    return (
        lagtime  = Float64(row[cols.lagtime]),
        peak     = Float64(row[cols.peak]),
        tpeak    = Float64(row[cols.tpeak]),
        max_rate = Float64(row[cols.max_rate]),
        etp      = Float64(row[cols.etp]),
    )
end

# ══════════════════════════════════════════════════════════════════════════════════
# Public helper 1: bz2012_ratios
# ══════════════════════════════════════════════════════════════════════════════════

"""
    bz2012_ratios(df, kept_cols; TM=0.0) -> DataFrame

Run the calibrated BZ2012 ODE (`CALIBRATED_OVERRIDES`) on every row of `df`
and pair the resulting predicted TGA features against the measured features
named by `kept_cols` (e.g. `TF_TGA_COLS` for TM=0.0, or `TM_TGA_COLS` for
TM=1e-9), returning a long DataFrame with columns
`PatientID, Feature, Predicted, Measured, Ratio`
(`Ratio = Predicted / Measured`), 5 rows per patient row of `df` (one per TGA
feature, in `_FEATURE_NAMES` order: lagtime, peak, tpeak, max_rate, etp).

`df` must have either a `SubjectID` column (real patients) or a `SyntheticID`
column (synthetic patients) — auto-detected — plus the 9 coagulation-factor
columns (`II,V,VII,VIII,IX,X,AT,PC,"TFPI Free"`) and the TGA columns named by
`kept_cols`.

Faithful extraction of `run_patient`/`to_clinical`/`get_measured` plus the
per-row/per-feature loop body of `evaluate_all!` in
`validate_mechanistic_plausibility.jl` (lines 96-134, 170-198), restricted to
a single TM condition per call — the frozen script's two-condition sweep is
reproduced by calling this twice: once with `(TF_TGA_COLS; TM=0.0)` and once
with `(TM_TGA_COLS; TM=1e-9)`.

Rows for which the ODE integration or feature extraction fails are skipped
(with a `@warn`), matching the frozen script's try/catch + warn-and-continue
behaviour.
"""
function bz2012_ratios(df::DataFrame, kept_cols; TM::Float64=0.0)
    id_col = hasproperty(df, :SubjectID) ? :SubjectID : :SyntheticID

    out = DataFrame(PatientID=String[], Feature=String[],
                     Predicted=Float64[], Measured=Float64[], Ratio=Float64[])

    for i in 1:nrow(df)
        row = df[i, :]
        pid = string(row[id_col])
        try
            feat = _run_patient(row, _P_CALIBRATED; TM_molar=TM)
            pred = _to_clinical(feat)
            meas = _get_measured(row, kept_cols)

            for (j, fname) in enumerate(_FEATURE_NAMES)
                m = getfield(meas, fname)
                p = getfield(pred, fname)
                push!(out, (pid, _FEATURE_LABELS[j], p, m, p / m))
            end
        catch e
            @warn "  bz2012_ratios: failed on $pid (TM=$TM): $(sprint(showerror, e))"
        end
    end

    return out
end

# ══════════════════════════════════════════════════════════════════════════════════
# Public helper 2: overlap_ks
# ══════════════════════════════════════════════════════════════════════════════════

"""
    overlap_ks(real_ratios::AbstractVector{<:Real}, synth_ratios::AbstractVector{<:Real})
        -> (frac_in_range, ks_D, ks_p)

Cloud-overlap + two-sample KS comparison of a real vs. synthetic Pred/Meas
ratio distribution. Mirrors `regen_mechanistic_figures.jl:96-109`:

1. clamp both vectors to `[0, 5]`,
2. band = 5th/95th percentile of the (clamped) real ratios,
3. `frac_in_range` = fraction of (clamped) synthetic ratios inside that band,
4. `ApproximateTwoSampleKSTest(real, synth)` on the clamped vectors → `(δ, pvalue)`.
"""
function overlap_ks(real_ratios::AbstractVector{<:Real}, synth_ratios::AbstractVector{<:Real})
    real_c  = clamp.(Float64.(real_ratios), 0.0, 5.0)
    synth_c = clamp.(Float64.(synth_ratios), 0.0, 5.0)

    real_lo, real_hi = quantile(real_c, 0.05), quantile(real_c, 0.95)
    frac_in_range = mean((synth_c .>= real_lo) .& (synth_c .<= real_hi))

    ks_test = ApproximateTwoSampleKSTest(real_c, synth_c)
    ks_D = ks_test.δ
    ks_p = pvalue(ks_test)

    return (frac_in_range, ks_D, ks_p)
end

# ──────────────────────────────────────────────────────────────────────────────
# test_mechanistic_eval.jl
#
# Equivalence check (not a package test — this repo has no test suite):
# confirms that mechanistic_eval.jl's `bz2012_ratios`/`overlap_ks` faithfully
# reproduce the published TF-only mechanistic-plausibility numbers computed by
# the FROZEN reference scripts `validate_mechanistic_plausibility.jl` /
# `regen_mechanistic_figures.jl` (neither is modified or re-run here), by:
#
#   (a) re-running the calibrated BZ2012 ODE fresh on the same 23 real
#       (complete-subject) + 100 synthetic patients (TF-only, TM=0) and
#       comparing the fresh Predicted/Measured values, row-for-row, against
#       the cached `data/validate_mechanistic_results.csv` (produced by
#       validate_mechanistic_plausibility.jl and never regenerated here), and
#
#   (b) checking `overlap_ks` on the fresh ratios lands in the published
#       TF-only band: cloud overlap 0.86-0.93, KS D 0.079-0.127, all p > 0.30.
#
# Running BZ2012 on ~123 patients solves ~123 stiff ODEs — expect a few
# minutes, plus HockinMannModel precompilation on first load.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests

include(joinpath(@__DIR__, "mechanistic_eval.jl"))

# ══════════════════════════════════════════════════════════════════════════════════
# Load data — same real-cohort selection as validate_mechanistic_plausibility.jl
# Step 2 (lines ~140-154): complete real subjects (all 3 visits present).
# ══════════════════════════════════════════════════════════════════════════════════
@info "Loading data …"
df_real  = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

all_subjects = unique(df_real.SubjectID)
complete_subjects = [s for s in all_subjects
    if sum(df_real.SubjectID .== s) == 3 &&
       all(v -> any((df_real.SubjectID .== s) .& (df_real.Visit .== v)), 1:3)]

df_real_complete = filter(r -> r.SubjectID in complete_subjects, df_real)
n_real  = length(complete_subjects)
n_synth = length(unique(df_synth.SyntheticID))

@info "  Real patients (complete): $n_real ($(nrow(df_real_complete)) records)"
@info "  Synthetic patients: $n_synth ($(nrow(df_synth)) records)"

# ══════════════════════════════════════════════════════════════════════════════════
# (a) Fresh bz2012_ratios (TF-only, TM=0.0) vs cached validate_mechanistic_results.csv
# ══════════════════════════════════════════════════════════════════════════════════
@info "\nRunning bz2012_ratios on real cohort (TF-only, calibrated) …"
real_r = bz2012_ratios(df_real_complete, TF_TGA_COLS; TM=0.0)
@info "  real_r: $(nrow(real_r)) rows"

@info "Running bz2012_ratios on synthetic cohort (TF-only, calibrated) …"
synth_r = bz2012_ratios(df_synth, TF_TGA_COLS; TM=0.0)
@info "  synth_r: $(nrow(synth_r)) rows"

@info "\nLoading cached validate_mechanistic_results.csv …"
cached = CSV.read(joinpath(_PATH_TO_DATA, "validate_mechanistic_results.csv"), DataFrame)
cached_real  = filter(r -> r.Source == "Real"      && r.Condition == "TF-only", cached)
cached_synth = filter(r -> r.Source == "Synthetic" && r.Condition == "TF-only", cached)

# evaluate_all! (frozen script) pushes 5 features per row, in _FEATURE_NAMES
# order, for each row of the SAME filtered df in the SAME order — so fresh
# and cached line up positionally without needing to sort (sorting on
# PatientID/Feature would be ambiguous: a real SubjectID repeats across its
# 3 visits with the same Feature labels).
function compare_to_cache(fresh::DataFrame, cached_sub::DataFrame; label::String, rtol::Float64)
    @assert nrow(fresh) == nrow(cached_sub) "$label: row count mismatch ($(nrow(fresh)) vs $(nrow(cached_sub)))"
    @assert fresh.PatientID == string.(cached_sub.PatientID) "$label: PatientID sequence mismatch"
    @assert fresh.Feature == cached_sub.Feature "$label: Feature sequence mismatch"

    rel_pred = abs.(fresh.Predicted .- cached_sub.Predicted) ./ (abs.(cached_sub.Predicted) .+ 1e-8)
    rel_meas = abs.(fresh.Measured  .- cached_sub.Measured)  ./ (abs.(cached_sub.Measured)  .+ 1e-8)
    @info "  $label: max rel err  Predicted=$(maximum(rel_pred)),  Measured=$(maximum(rel_meas))"

    @assert isapprox(fresh.Predicted, cached_sub.Predicted; rtol=rtol, atol=1e-6) "$label: Predicted diverges from cached CSV"
    @assert isapprox(fresh.Measured,  cached_sub.Measured;  rtol=1e-8, atol=1e-8) "$label: Measured diverges from cached CSV (should be an exact column re-read)"
end

const _RTOL = 1e-4   # ODE integration is deterministic given fixed p, u0, tspan, saveat — small headroom for solver/BLAS version drift
compare_to_cache(real_r,  cached_real;  label="Real TF-only",      rtol=_RTOL)
compare_to_cache(synth_r, cached_synth; label="Synthetic TF-only", rtol=_RTOL)

@info "\n✓ fresh bz2012_ratios reproduces cached validate_mechanistic_results.csv at rtol=$_RTOL"

# ══════════════════════════════════════════════════════════════════════════════════
# (b) overlap_ks vs published TF-only band
# ══════════════════════════════════════════════════════════════════════════════════
frac, ksD, ksp = overlap_ks(real_r.Ratio, synth_r.Ratio)   # pooled across features
@info "\n  Pooled TF-only: overlap=$(round(frac, digits=3)), KS D=$(round(ksD, digits=3)), p=$(round(ksp, digits=3))"

# Cross-check per-feature against the cached CSV's own per-feature overlap/KS
# (same computation regen_mechanistic_figures.jl performs and logs), in case
# the pooled number is ambiguous relative to the published per-feature range.
@info "\n  Per-feature (fresh):"
for feat in _FEATURE_LABELS
    r_sub = filter(r -> r.Feature == feat, real_r)
    s_sub = filter(r -> r.Feature == feat, synth_r)
    f, d, p = overlap_ks(r_sub.Ratio, s_sub.Ratio)
    @info "    $(rpad(feat, 22)) overlap=$(round(f, digits=3)), KS D=$(round(d, digits=3)), p=$(round(p, digits=3))"
end

# Published canonical (paper): TF-only cloud overlap 0.86-0.93, KS D 0.079-0.127, all p > 0.30
@assert 0.80 <= frac <= 0.95   "overlap $frac outside published TF-only band"
@assert ksD <= 0.20            "KS D $ksD larger than published"

println("mechanistic helper reproduces published TF-only envelope ✓ (overlap=$(round(frac,digits=3)), KS D=$(round(ksD,digits=3)))")

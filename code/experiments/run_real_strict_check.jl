# ──────────────────────────────────────────────────────────────────────────────
# run_real_strict_check.jl   (npj round 2, Reviewer 3 comment 3; spec W2)
#
# The strict mechanistic check requires all 15 TF-only BZ2012 ratios of a
# profile (5 TGA features x 3 visits) to fall inside each feature's real-cohort
# 5th-95th percentile band. This script applies the same check, and the
# biological check, to the 23 real patients, so the synthetic pass rates
# (54, 41, and 24 of 100) can be read against the real ones.
#
# Rows:
#   Real (in-sample bands)       the bands used for the synthetic profiles
#   Real (leave-one-out bands)   bands rebuilt without the patient being tested
#   Published synthetic          recomputed here; must equal (54, 41, 24)
#
# Reads:  data/full_longitudinal_memory.jld2, data/cleaned_full_data.csv,
#         data/synthetic_full_longitudinal.csv
# Writes: data/real_strict_check.csv
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))
include(joinpath(@__DIR__, "revision2_helpers.jl"))

ref = load_reference()
bands = strict_bands(ref)

# ── Real patients ──────────────────────────────────────────────────────────────
n_real = length(ref.complete_subjects)
real_bio = 0
real_strict = 0
real_both = 0
real_strict_loo = 0
real_both_loo = 0
for sid in ref.complete_subjects
    rows = ref.df_real_complete[ref.df_real_complete.SubjectID .== sid, :]
    ratios = ref.real_ratios[ref.real_ratios.PatientID .== string(sid), :]
    bio = bio_plausible_patient(rows, ref)
    strict = strict_mech_pass(ratios, bands, R2_N_VISITS)

    # leave-one-out bands: rebuild each feature's band without this patient
    others = ref.real_ratios[ref.real_ratios.PatientID .!= string(sid), :]
    bands_loo = Dict{String,Tuple{Float64,Float64}}()
    for f in unique(others.Feature)
        r = others.Ratio[others.Feature .== f]
        bands_loo[f] = (quantile(r, 0.05), quantile(r, 0.95))
    end
    strict_loo = strict_mech_pass(ratios, bands_loo, R2_N_VISITS)

    global real_bio += bio
    global real_strict += strict
    global real_both += (bio && strict)
    global real_strict_loo += strict_loo
    global real_both_loo += (bio && strict_loo)
end
@info "Real (in-sample bands): bio=$real_bio strict=$real_strict both=$real_both of $n_real"
@info "Real (leave-one-out bands): bio=$real_bio strict=$real_strict_loo both=$real_both_loo of $n_real"

# ── Published synthetic cohort (recomputed with the same functions) ────────────
df_pub = published_cohort_long(ref)
pub_ids = unique(df_pub.SyntheticID)
pub_bio = 0
pub_strict = 0
pub_both = 0
for id in pub_ids
    rows = df_pub[df_pub.SyntheticID .== id, :]
    bio = bio_plausible_patient(rows, ref)
    ratios = bz2012_ratios(rows, TF_TGA_COLS; TM=0.0)
    strict = strict_mech_pass(ratios, bands, R2_N_VISITS)
    global pub_bio += bio
    global pub_strict += strict
    global pub_both += (bio && strict)
end
@info "Published synthetic: bio=$pub_bio strict=$pub_strict both=$pub_both of $(length(pub_ids))"
@assert (pub_bio, pub_strict, pub_both) == (54, 41, 24) "published synthetic counts changed: ($pub_bio, $pub_strict, $pub_both)"
@info "  ✓ published synthetic counts reproduce (54, 41, 24)"

df_out = DataFrame(Cohort=["Real (in-sample bands)", "Real (leave-one-out bands)", "Published synthetic"],
                   N=[n_real, n_real, length(pub_ids)],
                   Bio_Pass=[real_bio, real_bio, pub_bio],
                   Strict_Mech_Pass=[real_strict, real_strict_loo, pub_strict],
                   Both_Pass=[real_both, real_both_loo, pub_both])
CSV.write(joinpath(_PATH_TO_DATA, "real_strict_check.csv"), df_out)
@info "Wrote data/real_strict_check.csv"

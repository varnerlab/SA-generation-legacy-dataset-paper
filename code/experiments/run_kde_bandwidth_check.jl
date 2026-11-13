# ──────────────────────────────────────────────────────────────────────────────
# run_kde_bandwidth_check.jl   (npj round 2, Reviewer 3 comment 1; spec W1)
#
# With equal multiplicity weights and unit-normalized memories, the SA sampling
# distribution is an isotropic Gaussian KDE on the 23 unit-normalized PCA
# profiles with bandwidth h = beta^(-1/2). This script computes the bandwidths
# that standard rules (Scott, Silverman, leave-one-out likelihood) would choose
# on those profiles, converts each to beta = 1/h^2, and generates 100 profiles
# at each beta with the exact sampler and the canonical decode. The beta* row
# must reproduce the Exact row of exact_sampler_reference.csv.
#
# Reads:  data/full_longitudinal_memory.jld2, data/exact_sampler_reference.csv
# Writes: data/kde_bandwidth_check.csv
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))
include(joinpath(@__DIR__, "revision2_helpers.jl"))

ref = load_reference()
bw = kde_bandwidths(ref.X̂)
@info "Bandwidths: sigma=$(round(bw.sigma, digits=4)) scott=$(round(bw.scott, digits=4)) " *
      "silverman=$(round(bw.silverman, digits=4)) loo=$(round(bw.loo, digits=4))"

rules = [("EntropyTransition", 1 / sqrt(ref.β_star)),
         ("Scott", bw.scott),
         ("Silverman", bw.silverman),
         ("LeaveOneOut", bw.loo)]

rows = NamedTuple[]
for (name, h) in rules
    β = 1 / h^2
    dirs, _ = exact_directions(ref, β, MersenneTwister(42), R2_N_SYNTH)
    df = decode_cohort(ref, dirs, ref.shared_magnitudes)
    m = evaluate_cohort(ref, dirs, df; with_mech=true)
    push!(rows, (Rule=name, h=h, Beta=β, Beta_over_BetaStar=β / ref.β_star,
                 Median_MRE=m.Median_MRE, CrossVisit_Frob=m.CrossVisit_Frob,
                 Mech_Overlap_TFonly=m.Mech_Overlap_TFonly, Median_Novelty=m.Median_Novelty,
                 Frac_Novelty_gt_thresh=m.Frac_Novelty_gt_thresh,
                 DCR_synth_to_real_median=m.DCR_synth_to_real_median))
    @info "  $name: h=$(round(h, digits=4)) β=$(round(β, digits=3)) β/β*=$(round(β / ref.β_star, digits=3)) " *
          "MRE=$(round(100 * m.Median_MRE, digits=2))% Frob=$(round(m.CrossVisit_Frob, digits=3)) " *
          "overlap=$(round(m.Mech_Overlap_TFonly, digits=3)) nov=$(round(m.Median_Novelty, digits=3)) " *
          "DCR=$(round(m.DCR_synth_to_real_median, digits=3))"
end
df_out = DataFrame(rows)

# Reproducibility check: the beta* row equals the published exact-sampler row.
exact_row = CSV.read(joinpath(_PATH_TO_DATA, "exact_sampler_reference.csv"), DataFrame)[1, :]
for f in (:Median_MRE, :CrossVisit_Frob, :Mech_Overlap_TFonly, :Median_Novelty, :DCR_synth_to_real_median)
    @assert isapprox(df_out[1, f], exact_row[f]; rtol=1e-10) "EntropyTransition $f does not match exact_sampler_reference.csv"
end
@info "  ✓ beta* row reproduces exact_sampler_reference.csv"

CSV.write(joinpath(_PATH_TO_DATA, "kde_bandwidth_check.csv"), df_out)
@info "Wrote data/kde_bandwidth_check.csv ($(nrow(df_out)) rows)"

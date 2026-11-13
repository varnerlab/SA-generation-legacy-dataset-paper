# ──────────────────────────────────────────────────────────────────────────────
# run_log_transform_check.jl   (npj round 2, Reviewer 3 comment 3; spec W3b)
#
# Tests a positive representation. Each representation rebuilds the memory
# matrix from the 23 real profiles with the natural log applied to selected
# columns, then follows the canonical pipeline (z-score, PCA at 95%, unit
# memories, its own beta*), draws 100 profiles with the exact sampler
# (latent seed 42; magnitudes from its own PCA norms, seed 424242), decodes,
# exponentiates the logged columns, and is scored in original units.
#
#   Untransformed    no log (must reproduce exact_sampler_reference.csv)
#   Log hormones     log of the 6 estradiol and progesterone columns
#   Log all assays   log of all 216 concatenated columns (all are positive)
#
# Reads:  data/full_longitudinal_memory.jld2, data/cleaned_full_data.csv,
#         data/exact_sampler_reference.csv
# Writes: data/log_transform_check.csv
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))
include(joinpath(@__DIR__, "revision2_helpers.jl"))

ref = load_reference()
A = hormone_affine(ref)
d_concat = ref.n_assays * R2_N_VISITS

representations = [("Untransformed", Int[]),
                   ("Log hormones", A.idx),
                   ("Log all assays", collect(1:d_concat))]

rows = NamedTuple[]
for (name, log_idx) in representations
    rep = build_representation(ref, log_idx)
    d = size(rep.X̂, 1)
    out = exact_sample(rep.X̂, R2_N_SYNTH; β=rep.β_star, rng=MersenneTwister(42))
    dirs = Matrix{Float64}(undef, d, R2_N_SYNTH)
    for i in 1:R2_N_SYNTH
        dirs[:, i] = r2_unitize(out.Ξ[i, :])
    end
    Random.seed!(R2_MAGNITUDE_SEED)
    mags = [rand(rep.pca_norms) for _ in 1:R2_N_SYNTH]
    df = decode_cohort_rep(ref, rep, dirs, mags)

    # Fidelity in original units; novelty in this representation's own PCA space.
    m = evaluate_cohort(merge(ref, (X̂=rep.X̂,)), dirs, df; with_mech=true)
    p = plausibility_summary(ref, df)
    push!(rows, (Representation=name, d=d, Beta_star=rep.β_star,
                 Median_MRE=m.Median_MRE, Hormone_MRE=p.Hormone_MRE,
                 CrossVisit_Frob=m.CrossVisit_Frob, Mech_Overlap_TFonly=m.Mech_Overlap_TFonly,
                 Bio_Pass=p.Bio_Pass, Hormone_Neg_Profiles=p.Hormone_Neg_Profiles,
                 Strict_Mech_Pass=p.Strict_Mech_Pass, Both_Pass=p.Both_Pass,
                 Median_Novelty=m.Median_Novelty, DCR_synth_to_real_median=m.DCR_synth_to_real_median))
    @info "$name: d=$d β*=$(round(rep.β_star, digits=3)) MRE=$(round(100 * m.Median_MRE, digits=2))% " *
          "hMRE=$(round(100 * p.Hormone_MRE, digits=1))% Frob=$(round(m.CrossVisit_Frob, digits=3)) " *
          "overlap=$(round(m.Mech_Overlap_TFonly, digits=3)) bio=$(p.Bio_Pass) neg=$(p.Hormone_Neg_Profiles) " *
          "strict=$(p.Strict_Mech_Pass) both=$(p.Both_Pass) nov=$(round(m.Median_Novelty, digits=3)) " *
          "DCR=$(round(m.DCR_synth_to_real_median, digits=3))"
end
df_out = DataFrame(rows)

# Reproducibility check: the untransformed row equals the published exact-sampler row.
exact_row = CSV.read(joinpath(_PATH_TO_DATA, "exact_sampler_reference.csv"), DataFrame)[1, :]
for f in (:Median_MRE, :CrossVisit_Frob, :Mech_Overlap_TFonly, :Median_Novelty, :DCR_synth_to_real_median)
    @assert isapprox(df_out[1, f], exact_row[f]; rtol=1e-10) "Untransformed $f does not match exact_sampler_reference.csv"
end
@info "  ✓ untransformed row reproduces exact_sampler_reference.csv"

CSV.write(joinpath(_PATH_TO_DATA, "log_transform_check.csv"), df_out)
@info "Wrote data/log_transform_check.csv"

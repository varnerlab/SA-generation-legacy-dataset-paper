# ──────────────────────────────────────────────────────────────────────────────
# run_convex_hull_analysis.jl
#
# E2: Convex-hull / extrapolation analysis. Quantifies whether SA's Langevin
# sampler interpolates within the convex hull of the 23 real patients (in raw
# 18-D PCA coordinates) or extrapolates beyond it, and — for points that do
# extrapolate — whether they remain biologically/mechanistically plausible
# ("controlled novel extrapolation") or drift into implausible territory.
#
# Method: for each of the 100 canonical synthetic patients' raw PCA coords
# (Ysyn, 18×100), test hull membership by solving the QP
# min ‖Yw − p‖² s.t. sum(w)=1, w≥0 (a convex combination of the 23 real
# columns of Y that best reproduces the point); the optimal value is 0 (point
# is a convex combination, i.e. inside the hull) or > 0, in which case its
# square root is the signed Euclidean distance from the point to the hull.
# For the out-of-hull subset, decode back to measurement space and check:
#   (i)  BioPlausible   — non-negativity + <5σ outlier predicates from
#                         validate_biological_constraints.jl's Check 3
#   (ii) MechInEnvelope — per-patient BZ2012 TF-only Predicted/Measured ratio
#                         falls inside the real cohort's [5th,95th] percentile
#                         band, via the Task 1c mechanistic_eval.jl helper
#
# A leave-one-out context check is also printed: with only K=23 real patients
# in d=18 PCA dimensions, the hull is a very sparse polytope, so even a
# genuine held-out real patient typically falls outside the hull of the other
# 22 by a wide margin — this calibrates how "extreme" a given out-of-hull
# distance for a synthetic point actually is.
#
# Reads:  data/full_longitudinal_memory.jld2   (X̂, pca_norms, pca_concat,
#                                                std_params_concat, kept_cols,
#                                                n_assays, complete_subjects,
#                                                df_clean)
#         data/synthetic_full_longitudinal.csv (canonical N=100 SA cohort)
# Writes: data/convex_hull_results.csv
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using JuMP, HiGHS
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load canonical SA memory + synthetic cohort
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading canonical SA memory + synthetic cohort …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_norms pca_concat std_params_concat kept_cols n_assays complete_subjects df_clean

d_pca, K = size(X̂)
n_visits = length(std_params_concat.col_names) ÷ n_assays
@assert length(kept_cols) == n_assays
@info "  X̂: $(size(X̂))  (d=$d_pca, K=$K real patients)   n_assays=$n_assays  n_visits=$n_visits"

Y = X̂ .* pca_norms'   # d×23 raw (un-normalized) PCA coordinates of the real cohort
@info "  Raw PCA coords Y: $(size(Y))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 1 (moved here from "Add the LP/QP solver dependency"): hull membership
# via LP feasibility (Yw = p, sum(w)=1, w≥0), with a QP objective (min ‖s‖²,
# s = Yw − p slack) so that infeasible points still yield a signed distance to
# the hull. `dist <= tol` ⟺ p is (numerically) inside conv(cols of Y).
# ══════════════════════════════════════════════════════════════════════════════
"""
    hull_membership(Y, p; tol=1e-7) -> (inhull::Bool, dist::Float64)

Is `p` (length-`d`) in the convex hull of the `K` columns of `Y` (`d×K`)?
Solves the QP `min ‖Yw − p‖²  s.t.  sum(w) = 1, w ≥ 0` via HiGHS/JuMP. The
optimal objective value is the squared Euclidean distance from `p` to the
hull (0 iff `p` is inside). Returns `(dist <= tol, dist)`.

Implementation note: the objective is built directly in terms of `w`
(`sum((Y*w - p)^2)`) rather than via an auxiliary slack vector `s` with
equality constraints `Yw - p == s` (the brief's originally-sketched
formulation). Both are mathematically identical QPs, but the slack-variable
version intermittently made HiGHS's QP crossover report a spurious
`OTHER_ERROR` ("claims optimality, but with ... primal infeasibilities" of
~1e-4) for a handful of the 100 synthetic points sitting close to a hull
facet — a solver robustness quirk, not a modeling ambiguity. The no-slack
formulation resolved all 100 points to `OPTIMAL` in testing and reproduces
the slack version's objective value exactly wherever the slack version does
converge.
"""
function hull_membership(Y::AbstractMatrix{<:Real}, p::AbstractVector{<:Real}; tol::Float64=1e-7)
    d, K = size(Y)
    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, w[1:K] >= 0)
    @constraint(m, sum(w) == 1)
    @objective(m, Min, sum((sum(Y[i, k] * w[k] for k in 1:K) - p[i])^2 for i in 1:d))
    optimize!(m)
    @assert termination_status(m) == MOI.OPTIMAL "hull_membership: HiGHS did not reach OPTIMAL (status=$(termination_status(m)))"
    dist = sqrt(max(objective_value(m), 0.0))
    return (dist <= tol, dist)
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: self-check (TDD RED→GREEN gate). Before hull_membership existed,
# running this script failed with `UndefVarError: hull_membership not
# defined` (confirmed) — now that it is implemented immediately above, the
# same block must pass.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nSelf-check: hull_membership on a real vertex (should be in-hull) and a far point (should be out) …"
let (inb, db) = hull_membership(Y, Y[:, 1])
    @assert inb && db < 1e-5 "self-check failed: real vertex should be in-hull with dist≈0 (got inhull=$inb, dist=$db)"
end
let (ino, dout) = hull_membership(Y, Y[:, 1] .+ 100)
    @assert !ino "self-check failed: far point should be out-of-hull (dist=$dout)"
end
@info "  ✓ hull_membership self-checks pass"

# ══════════════════════════════════════════════════════════════════════════════
# Build Ysyn (18×100): rebuild the 100×216 concatenated matrix from the
# canonical synthetic CSV (SAME block layout as run_full_longitudinal.jl /
# run_pca_ablation.jl's longdf_to_concat: V1 block, V2 block, V3 block,
# offset=(v-1)*n_assays), standardize with std_params_concat, then transform
# through the canonical pca_concat model.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nLoading synthetic cohort …"
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
synth_ids = sort(unique(df_synth.SyntheticID))
n_synth = length(synth_ids)
@info "  Synthetic patients: $n_synth"

"""
    longdf_to_concat(df, id_col, ids, kept_cols, n_assays, n_visits) -> Matrix (N × n_assays*n_visits)

Reshape a long-format (one row per patient×visit) DataFrame into a
concatenated-visit matrix, columns ordered V1 block, V2 block, V3 block —
same layout as `run_full_longitudinal.jl`'s `X_concat` / `run_pca_ablation.jl`'s
`longdf_to_concat`.
"""
function longdf_to_concat(df::DataFrame, id_col::Symbol, ids::Vector, kept_cols::Vector{Symbol},
                           n_assays::Int, n_visits::Int)
    N = length(ids)
    concat = Matrix{Float64}(undef, N, n_assays * n_visits)
    for (i, id) in enumerate(ids)
        for v in 1:n_visits
            row = findfirst((df[!, id_col] .== id) .& (df.Visit .== v))
            row === nothing && error("longdf_to_concat: missing record for id=$id, visit=$v")
            offset = (v - 1) * n_assays
            for (j, col) in enumerate(kept_cols)
                concat[i, offset + j] = Float64(df[row, col])
            end
        end
    end
    return concat
end

concat_syn = longdf_to_concat(df_synth, :SyntheticID, synth_ids, kept_cols, n_assays, n_visits)
Zsyn = (concat_syn .- std_params_concat.μ') ./ std_params_concat.σ'
Ysyn = MultivariateStats.transform(pca_concat, Zsyn')   # d_pca × n_synth
@assert size(Ysyn) == (d_pca, n_synth) "Ysyn shape mismatch: got $(size(Ysyn)), expected ($d_pca, $n_synth)"
@info "  Ysyn (synthetic, raw PCA coords): $(size(Ysyn))"

# ══════════════════════════════════════════════════════════════════════════════
# Round-trip check (auditable version of a claim the report otherwise only
# asserts in prose): rebuild the REAL cohort's concatenated matrix — same
# `complete_subjects` patient order and `kept_cols`/`n_assays`/`n_visits`
# block layout used to originally build X̂/Y in run_full_longitudinal.jl —
# via `longdf_to_concat`, push it through the IDENTICAL
# standardize(std_params_concat) → transform(pca_concat, ·) pipeline used
# above for Ysyn, and confirm it reproduces Y = X̂ .* pca_norms'. If this
# ever failed, Y and Ysyn would silently NOT be in the same coordinate
# system and every hull-membership distance computed below would be
# meaningless.
# ══════════════════════════════════════════════════════════════════════════════
concat_real_roundtrip = longdf_to_concat(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays, n_visits)
Z_real_roundtrip = (concat_real_roundtrip .- std_params_concat.μ') ./ std_params_concat.σ'
Y_roundtrip = MultivariateStats.transform(pca_concat, Z_real_roundtrip')
@assert isapprox(Y_roundtrip, Y; atol=1e-8) "round-trip check FAILED: rebuilding the real cohort through standardize(std_params_concat)+transform(pca_concat,·) does not reproduce Y = X̂ .* pca_norms' — Y and Ysyn are NOT in the same coordinate system"
@info "  ✓ Round-trip check: real-cohort rebuild via standardize+transform(pca_concat) reproduces Y (atol=1e-8)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3a: Classify all 100 synthetic points for hull membership.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nClassifying $n_synth synthetic points for hull membership …"
inhull_flags = falses(n_synth)
dists = zeros(n_synth)
for k in 1:n_synth
    inb, d = hull_membership(Y, Ysyn[:, k])
    inhull_flags[k] = inb
    dists[k] = d
end

n_in = sum(inhull_flags)
n_out = n_synth - n_in
frac_in_hull = n_in / n_synth
frac_out_hull = n_out / n_synth
@info "  In-hull:  $n_in/$n_synth ($(round(100*frac_in_hull, digits=1))%)"
@info "  Out-of-hull: $n_out/$n_synth ($(round(100*frac_out_hull, digits=1))%)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3b: Out-of-hull plausibility — (i) BioPlausible, (ii) MechInEnvelope.
# ══════════════════════════════════════════════════════════════════════════════
out_idx = findall(.!inhull_flags)
@info "\nAssessing plausibility of $(length(out_idx)) out-of-hull point(s) …"

BioPlausible = Vector{Union{Bool,Missing}}(missing, n_synth)
MechInEnvelope = Vector{Union{Bool,Missing}}(missing, n_synth)

if !isempty(out_idx)

    # ── (i) Bio-constraint predicates — reused from validate_biological_constraints.jl
    #    Check 3 (non-negativity of concentration-type features; no >5σ outliers
    #    vs the real cohort). Real cohort = cleaned_full_data.csv (SAME source
    #    validate_biological_constraints.jl uses for df_real, i.e. all subjects
    #    with any visit present, not just complete-case), pooled across visits.
    positive_cols = [:II, :V, :VII, :VIII, :IX, :X, :AT, :PC, :Fbgn, :Estradiol,
                     Symbol("TF Initiator Peak (nM)"), Symbol("TF Initiator ETP (nM·min)"),
                     Symbol("TF Initiator Lagtime (min)")]
    outlier_features = [:II, :V, :VII, :VIII, :IX, :X, :AT, :PC, :Fbgn,
                        Symbol("TF Initiator Peak (nM)"), Symbol("TF Initiator ETP (nM·min)")]
    positive_cols = filter(c -> c in kept_cols, positive_cols)
    outlier_features = filter(c -> c in kept_cols, outlier_features)

    df_real_full = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
    μ_real = Dict(col => mean(collect(skipmissing(df_real_full[:, col]))) for col in outlier_features)
    σ_real = Dict(col => std(collect(skipmissing(df_real_full[:, col]))) for col in outlier_features)

    """
        bio_plausible_patient(rows) -> Bool

    True iff, across the patient's visit rows, no `positive_cols` value is
    negative and no `outlier_features` value is >5σ from the real-cohort mean
    (Checks 3 of validate_biological_constraints.jl).
    """
    function bio_plausible_patient(rows::AbstractDataFrame)
        for col in positive_cols
            any(v -> v < 0, rows[:, col]) && return false
        end
        for col in outlier_features
            σ_real[col] > 1e-12 || continue
            any(v -> abs(v - μ_real[col]) > 5 * σ_real[col], rows[:, col]) && return false
        end
        return true
    end

    # ── (ii) Mechanistic-envelope check — Task 1c helper (mechanistic_eval.jl).
    #    Real-cohort ratio band computed ONCE, on complete-subject rows, mirroring
    #    test_mechanistic_eval.jl / run_pca_ablation.jl's shared cohort prep.
    @info "  Building real-cohort BZ2012 TF-only ratio band (complete-subject rows) …"
    df_real_complete = df_clean[in.(df_clean.SubjectID, Ref(complete_subjects)), :]
    real_ratios = bz2012_ratios(df_real_complete, TF_TGA_COLS; TM=0.0)
    @info "  real_ratios: $(nrow(real_ratios)) rows"

    # NOTE — banding convention (deliberate, documented here; not a bug):
    # `band` below is computed PER TGA FEATURE — one separate [5th,95th]
    # percentile band each for lagtime, peak, tpeak, max_rate, and ETP,
    # built only from that feature's own real-cohort Predicted/Measured
    # ratios. This is a deliberate departure from the sibling E1 helper
    # `overlap_ks` (mechanistic_eval.jl), which instead POOLS all features'
    # ratios into a single vector before taking one 5th/95th band. Pooling
    # here would mix heterogeneous quantities (a lagtime ratio, a peak
    # ratio, an ETP ratio, …) into one distribution, which is less
    # defensible than judging each feature against its own real-cohort
    # spread — hence the per-feature choice. The consequence is that
    # `mech_in_envelope_patient`'s "ALL features, at ALL 3 visits, must
    # fall in-band" criterion (below) is a strict, per-feature-and-per-visit
    # test — this per-feature banding is precisely what makes that
    # all-features-and-all-visits criterion so demanding relative to E1's
    # single pooled overlap_ks pass rate.
    feature_labels = unique(real_ratios.Feature)
    band = Dict(f => (quantile(filter(r -> r.Feature == f, real_ratios).Ratio, 0.05),
                       quantile(filter(r -> r.Feature == f, real_ratios).Ratio, 0.95))
                for f in feature_labels)
    for f in feature_labels
        lo, hi = band[f]
        @info "    band[$f] = [$(round(lo,digits=3)), $(round(hi,digits=3))]"
    end

    """
        mech_in_envelope_patient(rows) -> Bool

    Runs `bz2012_ratios` (TF-only, calibrated) on one out-of-hull synthetic
    patient's decoded visit rows and checks every resulting feature ratio
    against the real-cohort [5th,95th] band. Conservative on PARTIAL ODE
    failure, not just total failure: the full expected ratio-row count is
    `nrow(rows)` visits × `length(TF_TGA_COLS)` features (one row per
    feature per visit); if the ODE integration/feature-extraction fails on
    ANY visit — not only if it fails on every visit — the patient's ratio
    set is incomplete and envelope membership cannot be confirmed on the
    missing rows, so the patient is conservatively marked false rather than
    evaluated on whatever partial subset of rows happened to succeed.
    """
    function mech_in_envelope_patient(rows::AbstractDataFrame)
        r = bz2012_ratios(rows, TF_TGA_COLS; TM=0.0)
        expected_n = nrow(rows) * length(TF_TGA_COLS)
        nrow(r) == expected_n || return false
        for row in eachrow(r)
            lo, hi = band[row.Feature]
            (row.Ratio < lo || row.Ratio > hi) && return false
        end
        return true
    end

    # ── Decode the out-of-hull subset to measurement space and evaluate both checks.
    out_ids = synth_ids[out_idx]

    concat_out = Matrix{Float64}(undef, length(out_idx), n_assays * n_visits)
    for (i, k) in enumerate(out_idx)
        concat_out[i, :] = decode_sample(Ysyn[:, k], pca_concat, std_params_concat)
    end

    df_out_visits = DataFrame[]
    for v in 1:n_visits
        offset = (v - 1) * n_assays
        cols = concat_out[:, (offset + 1):(offset + n_assays)]
        df_v = DataFrame(cols, kept_cols)
        df_v.SyntheticID = out_ids
        df_v.Visit = fill(v, length(out_idx))
        push!(df_out_visits, df_v)
    end
    df_out_long = vcat(df_out_visits...)
    sort!(df_out_long, [:SyntheticID, :Visit])

    for (i, k) in enumerate(out_idx)
        sid = out_ids[i]
        rows = df_out_long[df_out_long.SyntheticID .== sid, :]
        bp = bio_plausible_patient(rows)
        BioPlausible[k] = bp
        me = mech_in_envelope_patient(rows)
        MechInEnvelope[k] = me
        @info "  [out-of-hull SyntheticID=$sid] dist=$(round(dists[k],digits=3))  BioPlausible=$bp  MechInEnvelope=$me"
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Assemble + write results
# ══════════════════════════════════════════════════════════════════════════════
df_results = DataFrame(
    SyntheticID    = synth_ids,
    InHull         = inhull_flags,
    DistToHull     = dists,
    BioPlausible   = BioPlausible,
    MechInEnvelope = MechInEnvelope,
)
CSV.write(joinpath(_PATH_TO_DATA, "convex_hull_results.csv"), df_results)
@info "\n✓ Wrote data/convex_hull_results.csv"

# ══════════════════════════════════════════════════════════════════════════════
# Context check: with only K=23 real patients in d=18 PCA dimensions, how far
# outside the hull does a genuine, held-out REAL patient fall relative to the
# other 22? (Leave-one-out.) This calibrates whether the synthetic
# dist-to-hull values above are "extreme" or in line with what the extreme
# sparsity of a 23-point hull in 18-D implies for ANY new point, real or not.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nContext: leave-one-out real-patient distances to the hull of the other $(K-1) …"
loo_dists = Vector{Float64}(undef, K)
for i in 1:K
    others = setdiff(1:K, i)
    _, dd = hull_membership(Y[:, others], Y[:, i])
    loo_dists[i] = dd
end
@info "  Real-cohort LOO dist-to-hull: median=$(round(median(loo_dists),digits=2))  " *
      "range=[$(round(minimum(loo_dists),digits=2)), $(round(maximum(loo_dists),digits=2))]"
if isempty(out_idx)
    @info "  Synthetic out-of-hull dist:   n/a (no out-of-hull points)"
else
    @info "  Synthetic out-of-hull dist:   median=$(round(median(dists[out_idx]),digits=2))  " *
          "range=[$(round(minimum(dists[out_idx]),digits=2)), $(round(maximum(dists[out_idx]),digits=2))]"
end

# ── Persist the LOO comparison as a committed artifact (this comparison is
#    load-bearing for a manuscript claim — R2.4 — so it must not live only
#    in console @info output). Rows: real_loo (K=23 leave-one-out real
#    patients vs. the hull of the other 22) and synth_outhull (the SA
#    synthetic points classified out-of-hull above vs. the full K=23 hull).
loo_summary_df = DataFrame(Group=String[], Median=Float64[], Min=Float64[], Max=Float64[], N=Int[])
push!(loo_summary_df, ("real_loo", median(loo_dists), minimum(loo_dists), maximum(loo_dists), length(loo_dists)))
if isempty(out_idx)
    push!(loo_summary_df, ("synth_outhull", NaN, NaN, NaN, 0))
else
    push!(loo_summary_df, ("synth_outhull", median(dists[out_idx]), minimum(dists[out_idx]), maximum(dists[out_idx]), length(out_idx)))
end
CSV.write(joinpath(_PATH_TO_DATA, "convex_hull_loo_summary.csv"), loo_summary_df)
@info "  ✓ Wrote data/convex_hull_loo_summary.csv"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: headline numbers + invariant checks
# ══════════════════════════════════════════════════════════════════════════════
plausible_mask = out_idx .|> k -> (BioPlausible[k] === true) && (MechInEnvelope[k] === true)
frac_plausible = isempty(out_idx) ? NaN : mean(plausible_mask)

bio_only_mask  = out_idx .|> k -> (BioPlausible[k] === true)
mech_only_mask = out_idx .|> k -> (MechInEnvelope[k] === true)
frac_bio_only  = isempty(out_idx) ? NaN : mean(bio_only_mask)
frac_mech_only = isempty(out_idx) ? NaN : mean(mech_only_mask)

println("\n" * "="^90)
println("E2: Convex-hull / extrapolation analysis — SA in raw 18-D PCA coordinates")
println("="^90)
println("  Inside hull  (interpolation):  $(round(100*frac_in_hull, digits=1))%  ($n_in/$n_synth)")
println("  Outside hull (extrapolation):  $(round(100*frac_out_hull, digits=1))%  ($n_out/$n_synth)")
println("  Outside-but-plausible (of the outside subset, Bio AND Mech): " *
        (isempty(out_idx) ? "n/a (no out-of-hull points)" :
         "$(round(100*frac_plausible, digits=1))%  ($(sum(plausible_mask))/$n_out)"))
if !isempty(out_idx)
    # NOTE: these are MARGINAL pass rates (each computed independently over
    # the out-of-hull subset), not mutually exclusive categories — they can
    # (and do) overlap with the joint Bio-AND-Mech rate above and with each
    # other, so they do not sum to 100% with the joint/fail-both rates.
    println("    (breakdown — BioPlausible (marginal): $(round(100*frac_bio_only,digits=1))%; " *
            "MechInEnvelope (considered alone): $(round(100*frac_mech_only,digits=1))%)")
end
println("="^90)

@info "\nInvariant checks …"
@assert isapprox(frac_in_hull + frac_out_hull, 1.0; atol=1e-12) "frac_in_hull + frac_out_hull must equal 1"
@assert isempty(out_idx) || (0.0 <= frac_plausible <= 1.0) "frac_plausible must be in [0,1]"
@info "  ✓ frac_in_hull + frac_out_hull == 1"
@info "  ✓ frac_plausible ∈ [0,1]"

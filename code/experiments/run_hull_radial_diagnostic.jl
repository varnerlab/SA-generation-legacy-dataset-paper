# ──────────────────────────────────────────────────────────────────────────────
# run_hull_radial_diagnostic.jl
#
# Diagnostic for E2: WHY is every SA synthetic outside the convex hull of the
# 23 real patients? Two candidate mechanisms:
#
#   (M1) Langevin noise nudges points just past a hull facet. Predicts the
#        radial overshoot 1/t* is barely above 1 (points sit just outside).
#   (M2) The generator samples a DIRECTION on the unit sphere (normalized
#        memories, renormalization at run_full_longitudinal.jl:129) and then
#        attaches an INDEPENDENTLY resampled radius from another patient
#        (line 130). Predicts a large overshoot, because the hull's radial
#        extent in a generic direction is far below the vertex-scale radius
#        that gets attached.
#
# FINDINGS (this script):
#   - M1 is refuted: median overshoot 6.75, not the ≈1.05 a facet nudge implies.
#   - M2's radius half is refuted too: median t* is 0.1481 in raw PCA space vs
#     0.1496 for the same directions on the unit sphere (ratio 1.01), so the
#     resampled radius is NOT what puts points outside. The real patients'
#     norms are tightly clustered (IQR 12.0-15.4), so the drawn radius lands
#     where the patient's own radius would have.
#   - What remains is M2's normalization half: renormalizing onto the sphere
#     (all 100 directions have t* < 1 in conv(X̂), max 0.231, as the convexity
#     argument requires) combined with hull thinness (the hull reaches the
#     patients' typical radius 13.8 only along the 23 vertex directions and
#     extends ~2.0 in a generic direction).
#   - SEVERITY MATTERS for the real-vs-synthetic comparison: the LOO real
#     loses its OWN vertex, so synthetics must lose their nearest vertex for a
#     fair test (section E/F). Matched, synthetics are NOT meaningfully less
#     extrapolative than held-out reals — anchor-free hull distance 11.0 vs
#     11.9, Mann-Whitney p=0.071. The unmatched comparison (8.04 vs 15.19)
#     overstates the gap and should not be used. Treat the two as
#     indistinguishable; the 100 synthetics also share one memory set, so they
#     are not independent and the p-values are optimistic.
#
# Method: for a point p, solve the LP  max t  s.t.  t·p = Y w, Σw = 1, w ≥ 0.
# t* is where the ray from the origin through p exits the hull (the origin is
# inside conv(Y) because PCA coordinates are centered). t* ≥ 1 ⟺ p in hull;
# 1/t* is the radial overshoot factor.
#
# The decisive control is the same statistic for a held-out REAL patient
# (LOO hull of the other 22). If real patients overshoot as much as
# synthetics, hull membership is not measuring extrapolation at this K.
#
# A third check tests the geometric claim directly in the NORMALIZED sampling
# space: any unit vector in the convex hull of unit-norm memories must equal a
# memory, so every renormalized sample should have t* < 1 in conv(X̂).
#
# Reads:  data/full_longitudinal_memory.jld2, data/synthetic_full_longitudinal.csv
# Writes: data/hull_radial_diagnostic.csv
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using JuMP, HiGHS
using Statistics, LinearAlgebra, HypothesisTests

@info "Loading canonical SA memory + synthetic cohort …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_norms pca_concat std_params_concat kept_cols n_assays complete_subjects df_clean

d_pca, K = size(X̂)
n_visits = length(std_params_concat.col_names) ÷ n_assays
Y = X̂ .* pca_norms'
@info "  d=$d_pca, K=$K real patients"

# ── radial hull exit factor ───────────────────────────────────────────────────
"""
    radial_hull_factor(Y, p; c=nothing) -> t*

Largest `t ≥ 0` with `c + t·(p − c) ∈ conv(columns of Y)`, i.e. where the ray
from the interior point `c` through `p` exits the hull. `c` defaults to the
centroid of `Y`'s columns, which is always in the hull, so `t = 0` is always
feasible. `t* ≥ 1` iff `p` itself is in the hull.

The centroid anchor matters for the leave-one-out control: dropping a patient
shifts the centroid off the origin, and the origin is then generally NOT in
the hull of the remaining 22, so an origin-anchored ray can miss the hull
entirely (LP infeasible) rather than returning a meaningful exit factor.
"""
function radial_hull_factor(Y::AbstractMatrix{<:Real}, p::AbstractVector{<:Real};
                            c::Union{Nothing,AbstractVector{<:Real}} = nothing)
    d, K = size(Y)
    cc = c === nothing ? vec(mean(Y, dims = 2)) : c
    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, w[1:K] >= 0)
    @variable(m, t >= 0)
    @constraint(m, sum(w) == 1)
    @constraint(m, [i = 1:d], sum(Y[i, k] * w[k] for k in 1:K) == cc[i] + t * (p[i] - cc[i]))
    @objective(m, Max, t)
    optimize!(m)
    @assert termination_status(m) == MOI.OPTIMAL "radial_hull_factor: status=$(termination_status(m))"
    return value(t)
end

# self-check: a real vertex must have t* ≈ 1 (it is on the hull boundary)
let t_v = radial_hull_factor(Y, Y[:, 1])
    @assert isapprox(t_v, 1.0; atol = 1e-4) "self-check failed: real vertex t*=$t_v, expected ≈1"
    @info "  ✓ self-check: real vertex t* = $(round(t_v, digits=6))"
end
# self-check: a point at half the radius of a vertex must be strictly inside
let t_h = radial_hull_factor(Y, 0.5 .* Y[:, 1])
    @assert t_h > 1.5 "self-check failed: half-radius vertex t*=$t_h, expected >1.5"
    @info "  ✓ self-check: half-radius vertex t* = $(round(t_h, digits=4))"
end

# ── rebuild synthetic PCA coords (same path as run_convex_hull_analysis.jl) ───
function longdf_to_concat(df::DataFrame, id_col::Symbol, ids::Vector, kept_cols::Vector{Symbol},
                          n_assays::Int, n_visits::Int)
    N = length(ids)
    concat = Matrix{Float64}(undef, N, n_assays * n_visits)
    for (i, id) in enumerate(ids)
        for v in 1:n_visits
            row = findfirst((df[!, id_col] .== id) .& (df.Visit .== v))
            row === nothing && error("missing record for id=$id, visit=$v")
            offset = (v - 1) * n_assays
            for (j, col) in enumerate(kept_cols)
                concat[i, offset+j] = Float64(df[row, col])
            end
        end
    end
    return concat
end

df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
synth_ids = sort(unique(df_synth.SyntheticID))
n_synth = length(synth_ids)
concat_syn = longdf_to_concat(df_synth, :SyntheticID, synth_ids, kept_cols, n_assays, n_visits)
Zsyn = (concat_syn .- std_params_concat.μ') ./ std_params_concat.σ'
Ysyn = MultivariateStats.transform(pca_concat, Zsyn')
@info "  Ysyn: $(size(Ysyn))"

# ══════════════════════════════════════════════════════════════════════════════
# A. Synthetics: radial overshoot vs the full 23-patient hull
# ══════════════════════════════════════════════════════════════════════════════
@info "\nA. Synthetic radial hull factors …"
t_syn      = [radial_hull_factor(Y, Ysyn[:, i]) for i in 1:n_synth]
norm_syn   = [norm(Ysyn[:, i]) for i in 1:n_synth]
rhull_syn  = t_syn .* norm_syn          # hull's radial extent in each synthetic's direction
overshoot_syn = 1.0 ./ t_syn

# ══════════════════════════════════════════════════════════════════════════════
# B. CONTROL — held-out real patients vs the hull of the other 22
# ══════════════════════════════════════════════════════════════════════════════
@info "B. Leave-one-out real radial hull factors …"
t_loo = Float64[]
for k in 1:K
    others = setdiff(1:K, k)
    push!(t_loo, radial_hull_factor(Y[:, others], Y[:, k]))
end
overshoot_loo = 1.0 ./ t_loo

# ══════════════════════════════════════════════════════════════════════════════
# C. Normalized sampling space — is the renormalized direction outside conv(X̂)?
# ══════════════════════════════════════════════════════════════════════════════
@info "C. Unit-sphere directions vs the hull of the unit-norm memories …"
t_unit = [radial_hull_factor(X̂, Ysyn[:, i] ./ norm_syn[i]) for i in 1:n_synth]

# ══════════════════════════════════════════════════════════════════════════════
# D. Vertex-count control. B tests each real against 22 vertices while A tests
# each synthetic against 23, and a hull with fewer vertices is thinner, which
# would inflate B's overshoot on its own. Re-test the synthetics against a
# 22-vertex hull (one real dropped at random, seeded) so both sides face the
# same vertex count.
# ══════════════════════════════════════════════════════════════════════════════
@info "D. Vertex-count control — synthetics vs 22-vertex hulls …"
Random.seed!(42)
t_syn22 = Float64[]
for i in 1:n_synth
    drop = rand(1:K)
    others = setdiff(1:K, drop)
    push!(t_syn22, radial_hull_factor(Y[:, others], Ysyn[:, i]))
end
overshoot_syn22 = 1.0 ./ t_syn22

# ══════════════════════════════════════════════════════════════════════════════
# E. HARSH control. D drops a RANDOM vertex from the synthetics' hull, but the
# LOO real loses its OWN vertex — the one best placed to contain it. That
# asymmetry inflates B and biases the comparison toward "reals extrapolate
# more". Re-test each synthetic with its NEAREST real vertex removed, which
# is the matched-severity version of the LOO operation.
# ══════════════════════════════════════════════════════════════════════════════
@info "E. Harsh control — synthetics with their NEAREST real vertex removed …"
nearest_real = [argmin([norm(Ysyn[:, i] .- Y[:, k]) for k in 1:K]) for i in 1:n_synth]
t_syn_harsh = Float64[]
for i in 1:n_synth
    others = setdiff(1:K, nearest_real[i])
    push!(t_syn_harsh, radial_hull_factor(Y[:, others], Ysyn[:, i]))
end
overshoot_syn_harsh = 1.0 ./ t_syn_harsh

# ══════════════════════════════════════════════════════════════════════════════
# F. ANCHOR-FREE control. The radial factor depends on where the ray is
# anchored. Euclidean distance to the hull does not. Same nearest-vertex-
# removed severity on both sides.
# ══════════════════════════════════════════════════════════════════════════════
"""
    hull_distance(Y, p) -> Euclidean distance from `p` to conv(columns of Y)

Solves `min ‖Yw − p‖² s.t. Σw = 1, w ≥ 0` (same formulation as
run_convex_hull_analysis.jl's `hull_membership`). Zero iff `p` is inside.
"""
function hull_distance(Y::AbstractMatrix{<:Real}, p::AbstractVector{<:Real})
    d, K = size(Y)
    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, w[1:K] >= 0)
    @constraint(m, sum(w) == 1)
    @objective(m, Min, sum((sum(Y[i, k] * w[k] for k in 1:K) - p[i])^2 for i in 1:d))
    optimize!(m)
    @assert termination_status(m) == MOI.OPTIMAL "hull_distance: status=$(termination_status(m))"
    return sqrt(max(objective_value(m), 0.0))
end

@info "F. Anchor-free hull distances (nearest vertex removed on both sides) …"
dist_syn_harsh = [hull_distance(Y[:, setdiff(1:K, nearest_real[i])], Ysyn[:, i]) for i in 1:n_synth]
dist_loo       = [hull_distance(Y[:, setdiff(1:K, k)], Y[:, k]) for k in 1:K]
p_overshoot = pvalue(MannWhitneyUTest(overshoot_syn_harsh, overshoot_loo))
p_distance = pvalue(MannWhitneyUTest(dist_syn_harsh, dist_loo))

# ══════════════════════════════════════════════════════════════════════════════
# Report
# ══════════════════════════════════════════════════════════════════════════════
@assert maximum(abs.(vec(mean(Y, dims = 2)))) < 1e-8 "Y is not centered; the origin-anchored ray in A is not the centroid ray"
q(v) = (round(quantile(v, 0.25), digits=3), round(median(v), digits=3), round(quantile(v, 0.75), digits=3))

@info "\n══════ RESULTS ══════"
@info "A. Synthetics (n=$n_synth) vs hull of all 23 reals:"
@info "   t* (exit factor)      IQR/median: $(q(t_syn))     frac in hull (t*≥1): $(round(mean(t_syn .>= 1.0), digits=3))"
@info "   overshoot 1/t*        IQR/median: $(q(overshoot_syn))"
@info "   ‖p‖ (drawn radius)    IQR/median: $(q(norm_syn))"
@info "   hull extent in dir p  IQR/median: $(q(rhull_syn))"
@info "   real ‖y‖ (pca_norms)  IQR/median: $(q(pca_norms))"
@info ""
@info "B. CONTROL — held-out real patients (n=$K) vs hull of other 22:"
@info "   t*                    IQR/median: $(q(t_loo))     frac in hull: $(round(mean(t_loo .>= 1.0), digits=3))"
@info "   overshoot 1/t*        IQR/median: $(q(overshoot_loo))"
@info ""
@info "C. Renormalized directions vs conv(X̂) (unit-norm memories):"
@info "   t*                    IQR/median: $(q(t_unit))     frac with t*<1: $(round(mean(t_unit .< 1.0), digits=3))"
@info "   max t* over all $n_synth: $(round(maximum(t_unit), digits=6))  (proof says < 1 unless the sample IS a memory)"

# M1 vs M2 discriminator
@info "\n══════ MECHANISM ══════"
@info "M1 (noise nudges past a facet) predicts median overshoot ≈ 1.0-1.1."
@info "M2 (direction/radius decoupling) predicts median overshoot >> 1."
@info "   observed median overshoot (synthetics, 23 vertices): $(round(median(overshoot_syn), digits=3))"
@info "   observed median overshoot (synthetics, 22 vertices): $(round(median(overshoot_syn22), digits=3))"
@info "   observed median overshoot (real LOO,   22 vertices): $(round(median(overshoot_loo), digits=3))"
@info "   vertex-count-matched ratio synth/real: $(round(median(overshoot_syn22)/median(overshoot_loo), digits=3))"
@info ""
@info "══════ SEVERITY-MATCHED CONTROLS (nearest vertex removed on both sides) ══════"
@info "   E. radial overshoot   synth (harsh): $(q(overshoot_syn_harsh))   real LOO: $(q(overshoot_loo))"
@info "      frac of synthetics below the real median: $(round(mean(overshoot_syn_harsh .< median(overshoot_loo)), digits=3))"
@info "      Mann-Whitney p: $(round(p_overshoot, sigdigits=3))"
@info "   F. hull DISTANCE (anchor-free)  synth (harsh): $(q(dist_syn_harsh))   real LOO: $(q(dist_loo))"
@info "      frac of synthetics below the real median: $(round(mean(dist_syn_harsh .< median(dist_loo)), digits=3))"
@info "      Mann-Whitney p: $(round(p_distance, sigdigits=3))"
@info ""
@info "   Radius resampling contribution: median t* in raw space $(round(median(t_syn), digits=4))"
@info "   vs median t* for the same directions on the unit sphere $(round(median(t_unit), digits=4))"
@info "   ratio $(round(median(t_unit)/median(t_syn), digits=3)) (≈1 ⟹ the drawn radius is NOT what puts points outside)"

out = DataFrame(
    Kind = vcat(fill("synthetic", n_synth), fill("real_loo", K), fill("unit_dir", n_synth)),
    Index = vcat(1:n_synth, 1:K, 1:n_synth),
    TStar = vcat(t_syn, t_loo, t_unit),
    Overshoot = vcat(overshoot_syn, overshoot_loo, 1.0 ./ t_unit),
    Radius = vcat(norm_syn, pca_norms, fill(1.0, n_synth)),
    HullExtent = vcat(rhull_syn, t_loo .* pca_norms, t_unit),
)
CSV.write(joinpath(_PATH_TO_DATA, "hull_radial_diagnostic.csv"), out)
@info "\nWrote data/hull_radial_diagnostic.csv"

# Persist the severity-matched comparisons used in the manuscript separately
# from convex_hull_loo_summary.csv, whose synthetic row is intentionally
# unmatched (all 23 real vertices). Repeating the comparison p-value on the two
# group rows keeps the artifact tidy and directly machine-readable.
matched_summary = DataFrame(
    Metric = ["radial_overshoot", "radial_overshoot", "hull_distance", "hull_distance"],
    Group = ["synthetic_nearest_vertex_removed", "real_leave_one_out",
             "synthetic_nearest_vertex_removed", "real_leave_one_out"],
    N = [length(overshoot_syn_harsh), length(overshoot_loo),
         length(dist_syn_harsh), length(dist_loo)],
    Q1 = [quantile(overshoot_syn_harsh, 0.25), quantile(overshoot_loo, 0.25),
          quantile(dist_syn_harsh, 0.25), quantile(dist_loo, 0.25)],
    Median = [median(overshoot_syn_harsh), median(overshoot_loo),
              median(dist_syn_harsh), median(dist_loo)],
    Q3 = [quantile(overshoot_syn_harsh, 0.75), quantile(overshoot_loo, 0.75),
          quantile(dist_syn_harsh, 0.75), quantile(dist_loo, 0.75)],
    Min = [minimum(overshoot_syn_harsh), minimum(overshoot_loo),
           minimum(dist_syn_harsh), minimum(dist_loo)],
    Max = [maximum(overshoot_syn_harsh), maximum(overshoot_loo),
           maximum(dist_syn_harsh), maximum(dist_loo)],
    MannWhitneyP = [p_overshoot, p_overshoot, p_distance, p_distance],
)
CSV.write(joinpath(_PATH_TO_DATA, "hull_severity_matched_summary.csv"), matched_summary)
@info "Wrote data/hull_severity_matched_summary.csv"

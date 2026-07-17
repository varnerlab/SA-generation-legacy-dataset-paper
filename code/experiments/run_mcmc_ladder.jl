# ──────────────────────────────────────────────────────────────────────────────
# run_mcmc_ladder.jl
#
# Task S3 — the V0→V3 MCMC ladder. The published generation pipeline
# (run_full_longitudinal.jl, "V0" below) draws each synthetic patient from a
# single anchor-initialized ULA chain run for T=2000 steps and keeps only the
# endpoint. Task S2 found that at β*≈2.94 the target π is MULTIMODAL (23 sharp
# memory basins): individual ULA chains lock into one basin and do not cross-
# mix within any tested horizon, so a mechanical R̂-based burn-in is an
# artifact (R̂ RISES with window start, B=0 mechanically but not trustworthy).
#
# This script isolates which piece of V0's design — anchor init, absence of
# burn-in/pooling, or ULA discretization bias — drives the SA memorization
# finding (E3: DCR synth→real < DCR real→real), by holding β*, α, the 18-D
# memory X̂, N=100, and the direction→magnitude decode fixed and varying ONLY
# the direction sampler across four rungs:
#
#   V0  anchor-init  + ULA    + endpoint            100 chains × 1
#   V1  sphere-init  + ULA    + endpoint            100 chains × 1
#   V2  sphere-init  + ULA    + burn-in B, thin τ, pool   20 chains × 5
#   V3  sphere-init  + MALA   + burn-in B, thin τ, pool   20 chains × 5
#
# B and τ are NOT hardcoded. Step 1 derives them from PER-CHAIN settling
# times (a chain enters and stabilizes within its own attractor basin — this
# is a single-chain stationarity property, distinct from the cross-chain R̂
# agreement that S2 showed fails here) and validates with a Geweke-style
# first-10%-vs-last-50% z-score check on each chain's post-B energy segment.
#
# Reads:  data/full_longitudinal_memory.jld2   (X̂, pca_concat, std_params_concat,
#                                                pca_norms, complete_subjects,
#                                                kept_cols, n_assays, df_clean)
#         data/mcmc_mixing_diagnostic.csv       (S2's B=0/τ=107 — printed only,
#                                                 as the starting/transient-
#                                                 inflated estimate; NOT reused)
# Writes: data/mcmc_ladder_burnin.csv           (B, τ, per-chain settling times,
#                                                Geweke z-scores)
#         data/mcmc_ladder_results.csv          (one row per rung, fidelity + DCR)
#
# Real/synthetic concat-matrix reconstruction uses `build_concat_matrix`
# (code/src/Patient.jl, Task S1) — no re-copied inline loop. The mechanistic
# TF-only axis uses `mechanistic_eval.jl`'s `bz2012_ratios`/`overlap_ks`
# (Task S1c), NOT a re-implementation. Mechanistic re-evaluation runs BZ2012
# across 23 (real, once) + 100×4 (rungs) patients — expect several minutes.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))

using Random, LinearAlgebra, Statistics

# ══════════════════════════════════════════════════════════════════════════════
# Step 0: Load canonical memory + real reference (held constant across rungs)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 0: Loading canonical pipeline memory …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
JLD2.@load mem_path X̂ pca_concat std_params_concat pca_norms complete_subjects kept_cols n_assays df_clean

d, K = size(X̂)
n_visits = 3
@info "  X̂: $(size(X̂))  (d=$d, K=$K real patients),  n_assays=$n_assays"

@info "Step 0b: β* via find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1,1000.0)) …"
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  β* = $(round(β_star, digits=4))"

# S2's mechanical (transient-inflated) B/τ — read for context/printing only,
# per the brief: "for the starting τ estimate only". NOT used downstream.
s2_path = joinpath(_PATH_TO_DATA, "mcmc_mixing_diagnostic.csv")
if isfile(s2_path)
    df_s2 = CSV.read(s2_path, DataFrame)
    s2_summary = df_s2[df_s2.row_type .== "summary", :]
    if nrow(s2_summary) == 1
        @info "  S2 mechanical estimate (context only, NOT reused): B=$(s2_summary.B[1]), τ=$(s2_summary.tau[1]) " *
              "(R̂-increasing artifact per task-S2-report.md — this script re-derives B/τ from settling times)"
    end
end

# Real reference matrices — via the Task-S1 helper, not a re-copied inline loop.
X_real = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)
@info "  X_real (build_concat_matrix): $(size(X_real))"

function safe_cor(X)
    # Copied per validate_cross_visit_covariance.jl:93 (local, as directed by the brief).
    C = cor(X)
    C[isnan.(C)] .= 0.0
    return C
end

C_real = safe_cor(X_real)

Zreal = (X_real .- std_params_concat.μ') ./ std_params_concat.σ'
rr = [minimum(norm(Zreal[i, :] .- Zreal[j, :]) for j in 1:K if j != i) for i in 1:K]
median_rr = median(rr)
@info "  DCR real→real median (reference, all rungs): $(round(median_rr, digits=4))"

df_real_complete = df_clean[in.(df_clean.SubjectID, Ref(complete_subjects)), :]
@info "  Running BZ2012 (TF-only, calibrated) on the $K real complete-case patients …"
real_ratios = bz2012_ratios(df_real_complete, TF_TGA_COLS; TM=0.0)
@info "  real_ratios: $(nrow(real_ratios)) rows"

# Shared decode channel: ONE seeded draw of 100 magnitudes from pca_norms,
# reused identically for every rung so fidelity/DCR differences are
# attributable only to the direction sampler (brief: "held constant across
# ALL rungs ... magnitude drawn from the 23 pca_norms (seeded)").
const MAGNITUDE_SEED = 424_242
Random.seed!(MAGNITUDE_SEED)
shared_magnitudes = [rand(pca_norms) for _ in 1:100]
@info "  Shared magnitude draw (seed=$MAGNITUDE_SEED): range [$(round(minimum(shared_magnitudes),digits=3)), $(round(maximum(shared_magnitudes),digits=3))]"

# ══════════════════════════════════════════════════════════════════════════════
# Helper: integrated_acf — copied verbatim from run_sampling_diagnostics.jl /
# run_mcmc_mixing_diagnostic.jl (standalone scripts, not include-able here).
# ══════════════════════════════════════════════════════════════════════════════
function integrated_acf(e::Vector{Float64}; max_lag::Int=500, threshold::Float64=0.05)
    n = length(e)
    μ = mean(e)
    var_e = var(e)
    var_e < 1e-30 && return (τ_int=1.0, acf=Float64[1.0])

    acf_vals = Float64[]
    τ = 0.5
    for lag in 1:min(max_lag, n - 1)
        c = mean((e[1:n-lag] .- μ) .* (e[lag+1:n] .- μ)) / var_e
        push!(acf_vals, c)
        c < threshold && break
        τ += c
    end
    return (τ_int=2.0 * τ, acf=acf_vals)
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Derive + validate burn-in B and thinning τ from settling times.
#
# S2 found that the mechanical split-R̂ burn-in rule is an artifact here (R̂
# RISES with window start — persistent, non-shrinking between-chain
# heterogeneity from basin-locking, not evidence of shared convergence).
# Instead: run M=20 sphere-init ULA chains long (T=6000), and for each chain
# find its own SETTLING STEP — the first t after which the trailing W=500
# window mean energy stays within ε of that chain's OWN final-basin mean
# energy for the rest of the run (a single-chain stationarity property,
# unrelated to whether chains agree with each other).
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 1: Deriving burn-in B + thinning τ from settling times …"

const M = 20
const T_SETTLE = 6000
const α_STEP = 0.01
const W_SETTLE = 500
const SPHERE_SEED_BASE = 9000
const MALA_SEED_BASE = 9500

@info "  Running M=$M sphere-init ULA chains (T=$T_SETTLE, α=$α_STEP) …"
sphere_energies = Matrix{Float64}(undef, M, T_SETTLE + 1)
sphere_traj = Vector{Matrix{Float64}}(undef, M)

for c in 1:M
    seed = SPHERE_SEED_BASE + c
    Random.seed!(seed)
    ξ0 = randn(d)
    ξ0 ./= norm(ξ0)
    res = sample(X̂, ξ0, T_SETTLE; β=β_star, α=α_STEP, seed=seed)
    sphere_traj[c] = res.Ξ
    for t in 0:T_SETTLE
        sphere_energies[c, t+1] = hopfield_energy(res.Ξ[t+1, :], X̂, β_star)
    end
end
@info "  Done."

"""
    settling_step(E; W=500) -> Int

First 0-based iteration `t` (of an energy trace `E` with `E[t+1]` the value
at time `t`) after which the trailing `W`-window mean stays within
`ε = max(2·std(final W-window), 1e-6)` of the chain's own final-basin mean
energy (mean of its last `W` samples) for the remainder of the trace — the
chain has entered and stabilized in its attractor basin. Returns
`length(E)-1` (never settles by this criterion) if no such `t` exists.

NOTE (documentation only, not a behavior change): the trailing-window width
`W=500` mechanically FLOORS the earliest reportable settling time at
`t=W-1=499` — no chain can ever be reported as settling before the window
itself has filled once. That is why 18/20 chains here report exactly t=499
(see task-S3-report.md). Consequently the derived B should not be trusted on
the settling statistic alone; it is cross-validated downstream by the Geweke
stationarity gate (`derive_tau_and_geweke`/`geweke_z`), which is the actual
adequacy check for whether B is large enough.
"""
function settling_step(E::Vector{Float64}; W::Int=500)
    T1 = length(E)
    W >= T1 && return T1 - 1

    final_window = @view E[end-W+1:end]
    final_mean = mean(final_window)
    ε = max(2 * std(final_window), 1e-6)

    n_win = T1 - W + 1
    tw_means = Vector{Float64}(undef, n_win)
    for ii in 1:n_win
        idx = W + ii - 1
        tw_means[ii] = mean(@view E[idx-W+1:idx])
    end
    within = abs.(tw_means .- final_mean) .< ε
    last_bad = findlast(!, within)
    settle_ii = last_bad === nothing ? 1 : min(last_bad + 1, n_win)
    settle_t = W - 2 + settle_ii
    return max(settle_t, 0)
end

settling_times = [settling_step(sphere_energies[c, :]; W=W_SETTLE) for c in 1:M]
@info "  Settling times (0-based t): min=$(minimum(settling_times)), median=$(median(settling_times)), " *
      "p90=$(round(quantile(settling_times,0.9),digits=1)), max=$(maximum(settling_times))"

B = ceil(Int, 1.5 * quantile(settling_times, 0.90))
@assert B < T_SETTLE "Derived B=$B must be less than T_SETTLE=$T_SETTLE"
@info "  B (burn-in) = ceil(1.5 × p90 settling time) = $B"

"""
    geweke_z(seg; frac1=0.1, frac2=0.5) -> Float64

Geweke-style stationarity z-score: compares the mean of the first `frac1`
fraction of `seg` against the mean of the last `frac2` fraction.

The energy trace here has substantial autocorrelation (τ≈90-190, verified
directly — NOT an artifact of the `integrated_acf` max_lag=500 cutoff: a
max_lag=3000/threshold=0.01 rerun reproduces the same τ within a few percent
on multiple chains). A first implementation used the raw per-segment
variance/n (i.i.d. assumption) and found the Geweke check FAILED at ~15-20%
pass fraction no matter how much B was escalated (0→749→867→…→1346, τ stable
at ≈118-120 throughout) — a naive i.i.d. z-score is inflated by ≈√τ for an
autocorrelated series and will reject stationarity for almost any chain with
τ>1, regardless of true burn-in. The standard fix (this is what Geweke's
actual diagnostic does, via a spectral-density/long-run-variance estimator)
is to divide each segment's variance by its OWN effective sample size
`n/τ_int(segment)` rather than by `n`. With that correction, the same B=749
(no escalation) already gives a 90% pass fraction (verified in a standalone
diagnostic against these exact chains) — the per-chain process is legitimately
stationary from early on; the earlier failure was the uncorrected variance
estimator, not a genuine finding about the chains.
"""
function geweke_z(seg::Vector{Float64}; frac1::Float64=0.1, frac2::Float64=0.5)
    n = length(seg)
    n1 = max(2, round(Int, frac1 * n))
    n2 = max(2, round(Int, frac2 * n))
    seg1 = seg[1:n1]
    seg2 = seg[end-n2+1:end]
    τ1 = integrated_acf(seg1; max_lag=min(500, n1 - 1)).τ_int
    τ2 = integrated_acf(seg2; max_lag=min(500, n2 - 1)).τ_int
    m1, v1 = mean(seg1), var(seg1)
    m2, v2 = mean(seg2), var(seg2)
    se = sqrt(v1 * τ1 / n1 + v2 * τ2 / n2)
    se < 1e-12 && return 0.0
    return (m1 - m2) / se
end

"""
    derive_tau_and_geweke(energies, B) -> (τ, τ_per_chain, z_vals, pass_frac)

Recompute τ as the ceil of the mean per-chain integrated ACF on the POST-B
energy segments (the settled-segment autocorrelation, which may differ
sharply from S2's transient-inflated τ=107), then Geweke-check each chain's
post-B segment.
"""
function derive_tau_and_geweke(energies::Matrix{Float64}, B::Int)
    Mloc = size(energies, 1)
    post = energies[:, B+1:end]
    τ_per_chain = [integrated_acf(post[c, :]).τ_int for c in 1:Mloc]
    τloc = ceil(Int, mean(τ_per_chain))
    z_vals = [geweke_z(post[c, :]) for c in 1:Mloc]
    pass_frac = mean(abs.(z_vals) .< 2)
    return τloc, τ_per_chain, z_vals, pass_frac
end

τ, τ_per_chain, geweke_z_vals, geweke_frac = derive_tau_and_geweke(sphere_energies, B)
@info "  τ (settled-segment, post-B) = $τ   Geweke pass fraction = $(round(geweke_frac,digits=3))"

n_escalations = 0
while geweke_frac < 0.8 && n_escalations < 5
    global B, τ, τ_per_chain, geweke_z_vals, geweke_frac, n_escalations
    n_escalations += 1
    B_old = B
    B = B + max(τ, 100)
    @assert B < T_SETTLE "Escalated B=$B must remain < T_SETTLE=$T_SETTLE"
    τ, τ_per_chain, geweke_z_vals, geweke_frac = derive_tau_and_geweke(sphere_energies, B)
    @warn "  Geweke pass fraction was < 0.8 at B=$B_old; escalated to B=$B (τ=$τ, pass=$(round(geweke_frac,digits=3))) [escalation #$n_escalations]"
end
@assert geweke_frac >= 0.8 "Geweke stationarity check failed after $n_escalations escalation(s): pass fraction=$geweke_frac"
@info "  FINAL: B=$B, τ=$τ, Geweke pass fraction=$(round(geweke_frac,digits=3)) (≥0.8 ✓, $n_escalations escalation(s))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Four direction samplers (shared decode)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 2: Drawing directions for all 4 rungs …"

unitize(v::Vector{Float64}) = v ./ (norm(v) + 1e-12)

"""
    pool_from_trajectories(trajs, B, τ; n_draws=5) -> d×(M·n_draws) directions

Pool `n_draws` post-burn-in, τ-thinned states per chain across all chains in
`trajs` (each a (T+1)×d trajectory matrix), normalized to unit directions.
"""
function pool_from_trajectories(trajs::Vector{Matrix{Float64}}, B::Int, τ::Int; n_draws::Int=5)
    Mloc = length(trajs)
    dloc = size(trajs[1], 2)
    dirs = Matrix{Float64}(undef, dloc, Mloc * n_draws)
    idx = 1
    for c in 1:Mloc
        for j in 1:n_draws
            t = B + j * τ
            dirs[:, idx] = unitize(trajs[c][t+1, :])
            idx += 1
        end
    end
    return dirs
end

# Extend the sphere-init ULA chains (reused for V2 pooling — S2's framing:
# "the 20 sphere-init chains act as a multi-start covering the basins" — is
# the intended design here) if B+5τ exceeds the T=6000 already run. Re-running
# with the SAME seed reproduces the existing prefix bit-for-bit and simply
# continues further (sample() reseeds deterministically at seed).
needed_T = B + 5 * τ
if needed_T > T_SETTLE
    @info "  Extending sphere-init ULA chains from T=$T_SETTLE to T=$needed_T for V2 pooling …"
    for c in 1:M
        seed = SPHERE_SEED_BASE + c
        Random.seed!(seed)
        ξ0 = randn(d)
        ξ0 ./= norm(ξ0)
        res = sample(X̂, ξ0, needed_T; β=β_star, α=α_STEP, seed=seed)
        sphere_traj[c] = res.Ξ
    end
else
    @info "  B+5τ=$needed_T ≤ T_SETTLE=$T_SETTLE: reusing the Step-1 sphere-init trajectories directly."
end

@info "  Running M=$M fresh sphere-init MALA chains (T=$needed_T, α=$α_STEP) for V3 …"
# Rigor fix (Task S3 revision): V3's per-chain initial ξ₀ is drawn from the
# SAME seed sequence V2 uses (SPHERE_SEED_BASE), not MALA_SEED_BASE, so chain
# c starts at the IDENTICAL sphere point in both V2 and V3 (common random
# numbers for initialization). MALA_SEED_BASE still seeds the propagation
# noise inside `mala_sample` (its own RNG stream, independent of the ξ0 draw)
# — the ONLY thing that differs between V2 and V3 is now the sampler itself
# (`sample` [ULA] vs `mala_sample` [MALA]), not the 20 starting points.
mala_traj = Vector{Matrix{Float64}}(undef, M)
mala_accept = Vector{Float64}(undef, M)
for c in 1:M
    init_seed = SPHERE_SEED_BASE + c
    Random.seed!(init_seed)
    ξ0 = randn(d)
    ξ0 ./= norm(ξ0)
    seed = MALA_SEED_BASE + c
    res = mala_sample(X̂, ξ0, needed_T; β=β_star, α=α_STEP, seed=seed)
    mala_traj[c] = res.Ξ
    mala_accept[c] = res.accept_rate
end
@info "  MALA acceptance (V3 chains): $(round(mean(mala_accept),digits=4)) ± $(round(std(mala_accept),digits=4))"

const SEED_V0 = 42   # matches run_full_longitudinal.jl's canonical Random.seed!(42)
const SEED_V1 = 43

"""
    draw_directions(rung, X̂, β_star; kwargs...) -> d×100 unit-direction matrix

rung ∈ {:V0, :V1, :V2, :V3}. V0/V1 draw 100 independent chains × 1 endpoint
each (T=2000, matching the published horizon); V2/V3 pool 5 post-burn-in,
τ-thinned draws from each of the 20 (`sphere_traj`/`mala_traj`) chains.
"""
function draw_directions(rung::Symbol, X̂::Matrix{Float64}, β_star::Float64;
    α::Float64=0.01, T0::Int=2000,
    B::Int=0, τ::Int=1,
    sphere_traj::Union{Vector{Matrix{Float64}},Nothing}=nothing,
    mala_traj::Union{Vector{Matrix{Float64}},Nothing}=nothing)

    dloc, Kloc = size(X̂)

    if rung == :V0
        Random.seed!(SEED_V0)
        dirs = Matrix{Float64}(undef, dloc, 100)
        for i in 1:100
            k0 = rand(1:Kloc)
            ξ0 = X̂[:, k0] .+ 0.01 .* randn(dloc)
            ξ0 ./= norm(ξ0)
            res = sample(X̂, ξ0, T0; β=β_star, α=α)
            dirs[:, i] = unitize(res.Ξ[end, :])
        end
        return dirs
    elseif rung == :V1
        Random.seed!(SEED_V1)
        dirs = Matrix{Float64}(undef, dloc, 100)
        for i in 1:100
            ξ0 = randn(dloc)
            ξ0 ./= norm(ξ0)
            res = sample(X̂, ξ0, T0; β=β_star, α=α)
            dirs[:, i] = unitize(res.Ξ[end, :])
        end
        return dirs
    elseif rung == :V2
        sphere_traj === nothing && error("V2 requires sphere_traj")
        return pool_from_trajectories(sphere_traj, B, τ)
    elseif rung == :V3
        mala_traj === nothing && error("V3 requires mala_traj")
        return pool_from_trajectories(mala_traj, B, τ)
    else
        error("unknown rung $rung")
    end
end

dirs_v0 = draw_directions(:V0, X̂, β_star)
dirs_v1 = draw_directions(:V1, X̂, β_star)
dirs_v2 = draw_directions(:V2, X̂, β_star; B=B, τ=τ, sphere_traj=sphere_traj)
dirs_v3 = draw_directions(:V3, X̂, β_star; B=B, τ=τ, mala_traj=mala_traj)
@info "  Drew directions: V0 $(size(dirs_v0)), V1 $(size(dirs_v1)), V2 $(size(dirs_v2)), V3 $(size(dirs_v3))"

"""
    decode_cohort(dirs, magnitudes, pca_model, std_params, kept_cols, n_assays, n_visits)
        -> DataFrame

Shared decode: for each direction column i, `decode_sample(dirs[:,i] .*
magnitudes[i], pca_model, std_params)` → 216-vector, split into the 3-visit
long format (kept_cols + SyntheticID + Visit) matching
`synthetic_full_longitudinal.csv`'s layout.
"""
function decode_cohort(dirs::Matrix{Float64}, magnitudes::Vector{Float64},
    pca_model, std_params, kept_cols::Vector{Symbol}, n_assays::Int, n_visits::Int)

    dloc, n = size(dirs)
    concat = Matrix{Float64}(undef, n, n_assays * n_visits)
    for i in 1:n
        concat[i, :] = decode_sample(dirs[:, i] .* magnitudes[i], pca_model, std_params)
    end
    df_visits = DataFrame[]
    for v in 1:n_visits
        offset = (v - 1) * n_assays
        cols = concat[:, (offset+1):(offset+n_assays)]
        df_v = DataFrame(cols, kept_cols)
        df_v.SyntheticID = 1:n
        df_v.Visit = fill(v, n)
        push!(df_visits, df_v)
    end
    df = vcat(df_visits...)
    sort!(df, [:SyntheticID, :Visit])
    return df
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Self-check (the failing test) — each rung's decoded cohort must be
# finite 100×216, and V0 must reproduce the published fidelity.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 3: Self-check …"

rung_dirs = Dict(:V0 => dirs_v0, :V1 => dirs_v1, :V2 => dirs_v2, :V3 => dirs_v3)
rung_order = [:V0, :V1, :V2, :V3]
rung_cohorts = Dict{Symbol,DataFrame}()

for rung in rung_order
    dirs = rung_dirs[rung]
    @assert size(dirs) == (d, 100) "$rung: expected $d×100 directions, got $(size(dirs))"
    @assert all(isfinite, dirs) "$rung: non-finite directions"

    df_rung = decode_cohort(dirs, shared_magnitudes, pca_concat, std_params_concat, kept_cols, n_assays, n_visits)
    X_check = build_concat_matrix(df_rung, :SyntheticID, collect(1:100), kept_cols, n_assays)
    @assert size(X_check) == (100, n_assays * n_visits) "$rung: decoded cohort wrong shape $(size(X_check))"
    @assert all(isfinite, X_check) "$rung: decoded cohort has non-finite values"

    rung_cohorts[rung] = df_rung
    @info "  $rung: decoded cohort finite 100×$(n_assays*n_visits) ✓"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Fidelity + DCR per rung
# ══════════════════════════════════════════════════════════════════════════════
"""
    evaluate_rung(name, dirs, df_rung) -> NamedTuple

Pooled median MRE (per-visit `feature_summary_comparison` on the K=23
complete-case real cohort → pooled median); cross-visit Frobenius (`safe_cor`,
`build_concat_matrix`); BZ2012 TF-only mechanistic overlap + KS
(`bz2012_ratios`/`overlap_ks`); median novelty (`sample_novelty` on the unit
directions, pre-decode); DCR synth→real median (Euclidean, canonical
standardized space, `std_params_concat`) alongside the fixed DCR real→real
reference.
"""
function evaluate_rung(name::String, dirs::Matrix{Float64}, df_rung::DataFrame)
    n = size(dirs, 2)

    all_mre = Float64[]
    for v in 1:n_visits
        real_v = df_real_complete[df_real_complete.Visit .== v, :]
        synth_v = df_rung[df_rung.Visit .== v, :]
        comp = feature_summary_comparison(real_v, synth_v, kept_cols)
        append!(all_mre, comp.Mean_Rel_Error)
    end
    med_mre = median(all_mre)

    X_rung = build_concat_matrix(df_rung, :SyntheticID, collect(1:n), kept_cols, n_assays)
    C_rung = safe_cor(X_rung)
    frob = norm(C_rung - C_real) / norm(C_real)

    @info "  [$name] BZ2012 mechanistic eval (TF-only) on $n synthetic patients …"
    synth_ratios = bz2012_ratios(df_rung, TF_TGA_COLS; TM=0.0)
    frac, ksD, ksp = overlap_ks(real_ratios.Ratio, synth_ratios.Ratio)

    nov = median([sample_novelty(dirs[:, i], X̂) for i in 1:n])

    Zrung = (X_rung .- std_params_concat.μ') ./ std_params_concat.σ'
    dcr_sr = [minimum(norm(Zrung[i, :] .- Zreal[j, :]) for j in 1:K) for i in 1:n]
    med_dcr_sr = median(dcr_sr)

    @info "  [$name] MRE=$(round(100*med_mre,digits=2))%  Frob=$(round(frob,digits=4))  " *
          "Overlap=$(round(frac,digits=3))  KS D=$(round(ksD,digits=3))  Nov=$(round(nov,digits=3))  " *
          "DCR(s→r)=$(round(med_dcr_sr,digits=3))  [DCR(r→r)=$(round(median_rr,digits=3))]"

    return (Rung=name, Median_MRE=med_mre, CrossVisit_Frob=frob,
        Mech_Overlap_TFonly=frac, Mech_KS_D=ksD, Mech_KS_p=ksp,
        Median_Novelty=nov, DCR_synth_to_real_median=med_dcr_sr,
        DCR_real_to_real_median=median_rr)
end

@info "\nStep 4: Fidelity + DCR per rung …"
results = NamedTuple[]
v0_median_mre = NaN

for rung in rung_order
    dirs = rung_dirs[rung]
    df_rung = rung_cohorts[rung]
    r = evaluate_rung(string(rung), dirs, df_rung)
    push!(results, r)

    if rung == :V0
        global v0_median_mre = r.Median_MRE
        @info "  >>> V0 cross-check: pooled median MRE = $(round(100*v0_median_mre, digits=3))% (target ≈1.2% ± 0.4pp)"
        @assert isapprox(v0_median_mre, 0.012; atol=0.004) "V0 does not reproduce published fidelity — harness mis-wired (got $(100*v0_median_mre)%)"
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Assemble results, mark best_rung, write CSVs, verify invariants
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 5: Assembling results …"

dcr_list = [r.DCR_synth_to_real_median for r in results]
best_idx = argmax(dcr_list)

df_results = DataFrame(
    Rung=[r.Rung for r in results],
    Median_MRE=[r.Median_MRE for r in results],
    CrossVisit_Frob=[r.CrossVisit_Frob for r in results],
    Mech_Overlap_TFonly=[r.Mech_Overlap_TFonly for r in results],
    Mech_KS_D=[r.Mech_KS_D for r in results],
    Mech_KS_p=[r.Mech_KS_p for r in results],
    Median_Novelty=[r.Median_Novelty for r in results],
    DCR_synth_to_real_median=[r.DCR_synth_to_real_median for r in results],
    DCR_real_to_real_median=[r.DCR_real_to_real_median for r in results],
    best_rung=[i == best_idx for i in 1:length(results)],
)
CSV.write(joinpath(_PATH_TO_DATA, "mcmc_ladder_results.csv"), df_results)

df_burnin_summary = DataFrame(
    row_type=["summary"], M=[M], T_settle=[needed_T > T_SETTLE ? needed_T : T_SETTLE],
    alpha=[α_STEP], beta_star=[β_star], W=[W_SETTLE],
    B=[B], tau=[τ], n_escalations=[n_escalations],
    geweke_pass_frac=[geweke_frac], ESS_post_B=[(T_SETTLE - B) / τ],
    mala_accept_mean=[mean(mala_accept)], mala_accept_std=[std(mala_accept)],
)
df_burnin_chains = DataFrame(
    row_type=fill("chain", M),
    chain=1:M, seed=[SPHERE_SEED_BASE + c for c in 1:M],
    settling_time=settling_times,
    tau_per_chain=τ_per_chain,
    geweke_z=geweke_z_vals,
    geweke_pass=abs.(geweke_z_vals) .< 2,
)
df_burnin = vcat(df_burnin_summary, df_burnin_chains; cols=:union)
CSV.write(joinpath(_PATH_TO_DATA, "mcmc_ladder_burnin.csv"), df_burnin)

@info "  Wrote data/mcmc_ladder_results.csv ($(nrow(df_results)) rows) and data/mcmc_ladder_burnin.csv ($(nrow(df_burnin)) rows)"

# ── Invariant checks ──
@assert nrow(df_results) == 4 "expected 4 rung rows, got $(nrow(df_results))"
@assert all(df_results.Median_MRE .> 0) "every Median_MRE must be > 0"
@assert all(0 .<= df_results.Mech_Overlap_TFonly .<= 1) "Mech_Overlap_TFonly must be in [0,1]"
@assert all(0 .<= df_results.Median_Novelty .<= 1) "Median_Novelty must be in [0,1]"
@assert sum(df_results.best_rung) == 1 "expected exactly one best_rung=true, got $(sum(df_results.best_rung))"
@assert isapprox(v0_median_mre, 0.012; atol=0.004) "V0 MRE cross-check failed at final assembly"
@assert geweke_frac >= 0.8 "Geweke check must pass (≥0.8) at final assembly"
@info "\nAll self-checks passed."

# ══════════════════════════════════════════════════════════════════════════════
# Final report
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^100)
println("MCMC Ladder (Task S3) — Burn-in derivation")
println("="^100)
println("  M=$M sphere-init ULA chains, T_settle=$T_SETTLE, W=$W_SETTLE, α=$α_STEP, β*=$(round(β_star,digits=4))")
println("  Settling times: min=$(minimum(settling_times)), median=$(median(settling_times)), " *
        "p90=$(round(quantile(settling_times,0.9),digits=1)), max=$(maximum(settling_times))")
println("  B (burn-in)  = $B   [= ceil(1.5 × p90 settling time), $n_escalations escalation(s)]")
println("  τ (thinning) = $τ   [ceil(mean integrated ACF, post-B settled segments)]")
println("  Geweke pass fraction (|z|<2, post-B, first-10%-vs-last-50%) = $(round(geweke_frac,digits=3))  (≥0.8 required)")
println("  needed_T for V2/V3 pooling = B+5τ = $needed_T")
println("-"^100)
println("MCMC Ladder — Rung × Metric Table")
println("="^100)
pretty_table(df_results)

println("\n" * "="^100)
println("Honest read-out")
println("="^100)
v0_row = df_results[findfirst(==("V0"), df_results.Rung), :]
v1_row = df_results[findfirst(==("V1"), df_results.Rung), :]
v2_row = df_results[findfirst(==("V2"), df_results.Rung), :]
v3_row = df_results[findfirst(==("V3"), df_results.Rung), :]
println("  V1 vs V0 (init only):     DCR(s→r) $(round(v0_row.DCR_synth_to_real_median,digits=3)) → $(round(v1_row.DCR_synth_to_real_median,digits=3))" *
        "   MRE $(round(100*v0_row.Median_MRE,digits=2))% → $(round(100*v1_row.Median_MRE,digits=2))%")
println("  V3 vs V2 (ULA→MALA only): DCR(s→r) $(round(v2_row.DCR_synth_to_real_median,digits=3)) → $(round(v3_row.DCR_synth_to_real_median,digits=3))" *
        "   MRE $(round(100*v2_row.Median_MRE,digits=2))% → $(round(100*v3_row.Median_MRE,digits=2))%")
println("  V2 vs V0 (pooling/burn-in): DCR(s→r) $(round(v0_row.DCR_synth_to_real_median,digits=3)) → $(round(v2_row.DCR_synth_to_real_median,digits=3))" *
        "   MRE $(round(100*v0_row.Median_MRE,digits=2))% → $(round(100*v2_row.Median_MRE,digits=2))%")
for row in eachrow(df_results)
    moved = row.DCR_synth_to_real_median > median_rr ? "ABOVE real→real ($(round(median_rr,digits=3))) — privacy improves" :
            "still below real→real ($(round(median_rr,digits=3)))"
    println("  $(row.Rung): DCR(s→r)=$(round(row.DCR_synth_to_real_median,digits=3)) — $moved")
end
best_row = df_results[df_results.best_rung, :][1, :]
println("-"^100)
println("  best_rung (max DCR synth→real) = $(best_row.Rung)  (DCR=$(round(best_row.DCR_synth_to_real_median,digits=3)), " *
        "MRE=$(round(100*best_row.Median_MRE,digits=2))%)")
println("="^100)

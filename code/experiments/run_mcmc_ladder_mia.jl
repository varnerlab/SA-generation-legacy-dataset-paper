# ──────────────────────────────────────────────────────────────────────────────
# run_mcmc_ladder_mia.jl — Task S4: repeated-split membership-inference (MIA)
# confirmation on V0 and the best-DCR rung (V3)
#
# The S3 MCMC ladder (run_mcmc_ladder.jl) found DCR synth→real never rises
# above the DCR real→real baseline across any of its four direction-sampler
# rungs (V0 anchor-init/ULA/endpoint … V3 sphere-init/MALA/burn-in+pooled) —
# the privacy leak persists even under a "proper" MCMC scheme (MALA, near-zero
# discretization bias, multi-chain pooling past a derived+validated burn-in).
# This script asks the same question with the STRONGER membership-inference
# metric: run E3b's adaptive repeated-split MIA (`run_privacy_split_sensitivity.jl`)
# on V0 (the published scheme, reproducing E3b's mean AUC ≈0.975 as a wiring
# cross-check) AND on V3 (the best-DCR rung), holding everything EXCEPT the
# direction sampler fixed for the 15-subject-subset holdout retrain at each
# repeated split — so a difference in AUC(V0) vs AUC(V3) is attributable only
# to the direction-sampler change (anchor-init ULA endpoint vs sphere-init
# MALA burn-in+thinned pooling), the same controlled-comparison logic S3 used
# for fidelity/DCR, now applied to the privacy metric itself.
#
# `split_mia_auc` (E3b, run_privacy_split_sensitivity.jl) is generalized here
# into a `scheme`-parameterized retrain: :V0 calls the frozen
# `sa_generate_from_matrix` helper directly (Patient.jl, Task S1) — this IS
# E3b's V0 retrain, unchanged; :V3 is a new `v3_retrain` that mirrors
# `sa_generate_from_matrix`'s standardize→PCA(0.95)→unit-norm→pca_norms→β*
# prefix on the TRAIN SUBSET's own memory, then swaps only the direction step
# for V3's scheme lifted from `run_mcmc_ladder.jl` (20 sphere-init MALA
# chains, run B+5τ steps, pool 5 τ-thinned post-burn-in draws per chain — B/τ
# read from `mcmc_ladder_burnin.csv`'s derived+validated summary row, NOT
# hardcoded). E3's frozen scripts (`run_privacy_split_sensitivity.jl`,
# `run_privacy_analysis.jl`, `run_mcmc_ladder.jl`, `src/Patient.jl`) are READ
# and adapted from, never modified.
#
# Reads:  data/full_longitudinal_memory.jld2   (df_clean, kept_cols, n_assays,
#                                                complete_subjects,
#                                                std_params_concat — same
#                                                subset E3b reads; MIA works in
#                                                standardized measurement space
#                                                and needs neither the full
#                                                cohort X̂ nor the mechanistic
#                                                model)
#         data/mcmc_ladder_results.csv          (confirms best_rung == "V3")
#         data/mcmc_ladder_burnin.csv           (B, τ — derived+validated in S3)
# Writes: data/mcmc_ladder_mia.csv              (per-rep + summary AUC rows,
#                                                 Scheme ∈ {V0, V3})
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

using LinearAlgebra, Statistics, Random

# ══════════════════════════════════════════════════════════════════════════════
# Step 0: Load canonical pipeline memory + cleaned data (same subset E3b
# reads), confirm best_rung == "V3" from the S3 ladder, and read B/τ from the
# S3-derived burn-in summary row (not hardcoded).
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 0: Loading canonical pipeline memory + cleaned data …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
d = JLD2.load(mem_path)
df_clean          = d["df_clean"]
kept_cols         = d["kept_cols"]
n_assays          = d["n_assays"]
complete_subjects = d["complete_subjects"]
std_params_concat = d["std_params_concat"]

K = length(complete_subjects)
n_visits = 3
d_concat = n_assays * n_visits
@info "  Complete subjects K=$K, n_assays=$n_assays, d_concat=$d_concat"

@info "Step 0b: Confirming best_rung from data/mcmc_ladder_results.csv …"
ladder_results_path = joinpath(_PATH_TO_DATA, "mcmc_ladder_results.csv")
isfile(ladder_results_path) || error("data/mcmc_ladder_results.csv not found — run run_mcmc_ladder.jl (S3) first.")
df_ladder_results = CSV.read(ladder_results_path, DataFrame)
best_rows = df_ladder_results[df_ladder_results.best_rung .== true, :]
@assert nrow(best_rows) == 1 "Expected exactly one best_rung==true row in mcmc_ladder_results.csv, got $(nrow(best_rows))"
const BEST_RUNG = String(best_rows.Rung[1])
@info "  best_rung (from CSV) = $BEST_RUNG"
@assert BEST_RUNG == "V3" "Task brief's cross-task context asserts best_rung == V3 (max DCR synth→real); CSV reports \"$BEST_RUNG\" instead. STOP — re-check before proceeding."

@info "Step 0c: Reading B, τ from data/mcmc_ladder_burnin.csv summary row …"
burnin_path = joinpath(_PATH_TO_DATA, "mcmc_ladder_burnin.csv")
isfile(burnin_path) || error("data/mcmc_ladder_burnin.csv not found — run run_mcmc_ladder.jl (S3) first.")
df_burnin = CSV.read(burnin_path, DataFrame)
burnin_summary = df_burnin[df_burnin.row_type .== "summary", :]
@assert nrow(burnin_summary) == 1 "Expected exactly one burn-in summary row, got $(nrow(burnin_summary))"
const B_V3 = Int(burnin_summary.B[1])
const τ_V3 = Int(burnin_summary.tau[1])
@info "  B (burn-in) = $B_V3, τ (thinning) = $τ_V3  (S3-derived+validated; read, not hardcoded)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Rebuild the 23×216 real concatenated matrix + standardize → Zreal
# (identical construction to run_privacy_split_sensitivity.jl Step 1, via the
# Task-S1 helper `build_concat_matrix` instead of a re-copied inline loop).
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Rebuilding real concatenated matrix + standardizing …"
X_concat = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)
@info "  X_concat (real): $(size(X_concat))"

Zreal = (X_concat .- std_params_concat.μ') ./ std_params_concat.σ'
@info "  Zreal: $(size(Zreal))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Generalized split-MIA — a `scheme`-parameterized retrain.
#
# :V0 = published anchor-init single-ULA-endpoint scheme, via the FROZEN
#       `sa_generate_from_matrix` helper (Patient.jl, Task S1) — this is
#       exactly E3b's retrain call, unmodified.
# :V3 = sphere-init MALA, pooled 20 chains × 5 thinned post-burn-in draws,
#       lifted/adapted from `run_mcmc_ladder.jl`'s `draw_directions(:V3,...)`
#       + `pool_from_trajectories`, applied to the TRAIN SUBSET's own memory.
# ══════════════════════════════════════════════════════════════════════════════
const ALPHA_STEP = 0.01
const T0_V0      = 2000   # V0's published ULA horizon (sa_generate_from_matrix default)
const M_CHAINS   = 20     # V3 pooling: 20 chains × 5 draws = 100
const N_DRAWS    = 5

unitize(v::AbstractVector{Float64}) = v ./ (norm(v) + 1e-12)

"""
    v3_retrain(X_train, kept_cols, n_assays, B, τ; gen_seed=42, pratio=0.95)
        -> synth_df::DataFrame

Regenerate the N=100 holdout-retrain synthetics using V3's direction-sampler
scheme (sphere-init MALA, pooled 20 chains × 5 τ-thinned post-burn-in draws),
applied to `X_train`'s OWN memory: standardize → PCA(`pratio`) → unit-norm →
`pca_norms` → β* via `find_entropy_inflection` on the subset — this prefix
mirrors `sa_generate_from_matrix` (Patient.jl:639-671) exactly, so the V0 and
V3 retrain paths share an IDENTICAL subset-memory construction and differ
ONLY in the direction sampler (the same controlled-comparison design S3 used
across its four rungs). The direction step itself is lifted from
`run_mcmc_ladder.jl`'s `draw_directions(:V3,...)`/`pool_from_trajectories`,
parameterized by `gen_seed` (in place of the ladder's hardcoded
SPHERE_SEED_BASE/MALA_SEED_BASE) so it is deterministic given `gen_seed` and
the ONLY thing that varies across repeated-split reps is `X_train` (the
partition) — matching V0's `gen_seed=42`-fixed convention.
"""
function v3_retrain(X_train::Matrix{Float64}, kept_cols::Vector{Symbol}, n_assays::Int,
                     B::Int, τ::Int; gen_seed::Int=42, pratio::Float64=0.95)

    K_sub, d_concat_sub = size(X_train)
    n_visits_sub = d_concat_sub ÷ n_assays

    concat_col_names = Symbol[]
    for v in 1:n_visits_sub, col in kept_cols
        push!(concat_col_names, Symbol("V$(v)_$(col)"))
    end

    μ_concat = vec(mean(X_train, dims=1))
    σ_concat = vec(std(X_train, dims=1))
    for i in eachindex(σ_concat)
        σ_concat[i] < 1e-12 && (σ_concat[i] = 1.0)
    end
    Z_concat = (X_train .- μ_concat') ./ σ_concat'
    std_params_sub = StandardizationParams(μ_concat, σ_concat, concat_col_names)

    pca_sub = MultivariateStats.fit(PCA, Z_concat'; pratio=pratio)
    d_pca_sub = MultivariateStats.outdim(pca_sub)
    X_pca_sub = MultivariateStats.transform(pca_sub, Z_concat')

    pca_norms_sub = [norm(X_pca_sub[:, k]) for k in 1:K_sub]
    X̂_sub = copy(X_pca_sub)
    for k in 1:K_sub
        X̂_sub[:, k] ./= (norm(X̂_sub[:, k]) + 1e-12)
    end

    phase = find_entropy_inflection(X̂_sub; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
    β_star_sub = phase.β_star

    # ── V3 direction sampler: 20 sphere-init MALA chains, B+5τ steps, pool 5
    # τ-thinned post-burn-in draws each (lifted from run_mcmc_ladder.jl's
    # draw_directions(:V3,...) / pool_from_trajectories) ──
    needed_T = B + N_DRAWS * τ
    sphere_base = gen_seed * 1000
    mala_base   = gen_seed * 1000 + 500

    dirs = Matrix{Float64}(undef, d_pca_sub, M_CHAINS * N_DRAWS)
    idx = 1
    for c in 1:M_CHAINS
        init_seed = sphere_base + c
        Random.seed!(init_seed)
        ξ0 = randn(d_pca_sub)
        ξ0 ./= norm(ξ0)
        mseed = mala_base + c
        res = mala_sample(X̂_sub, ξ0, needed_T; β=β_star_sub, α=ALPHA_STEP, seed=mseed)
        for j in 1:N_DRAWS
            t = B + j * τ
            dirs[:, idx] = unitize(res.Ξ[t+1, :])
            idx += 1
        end
    end

    # Magnitude draws: reseeded independently (deterministic given gen_seed,
    # not dependent on the number of RNG draws consumed by the chains above).
    Random.seed!(gen_seed)
    n = size(dirs, 2)
    magnitudes = [rand(pca_norms_sub) for _ in 1:n]

    concat = Matrix{Float64}(undef, n, n_assays * n_visits_sub)
    for i in 1:n
        concat[i, :] = decode_sample(dirs[:, i] .* magnitudes[i], pca_sub, std_params_sub)
    end

    df_visits = DataFrame[]
    for v in 1:n_visits_sub
        offset = (v - 1) * n_assays
        cols = concat[:, (offset+1):(offset+n_assays)]
        df_v = DataFrame(cols, kept_cols)
        df_v.SyntheticID = 1:n
        df_v.Visit = fill(v, n)
        push!(df_visits, df_v)
    end
    synth_df = vcat(df_visits...)
    sort!(synth_df, [:SyntheticID, :Visit])
    return synth_df
end

"""
    split_mia_auc(scheme, Zreal, X_concat, train_idx, holdout_idx, kept_cols,
                  n_assays; gen_seed=42, B=B_V3, τ=τ_V3) -> Float64

Generalized E3b split-MIA (factored from `run_privacy_split_sensitivity.jl`'s
`split_mia_auc`): retrain SA on `X_concat[train_idx,:]` using the `scheme`
direction sampler (`:V0` → `sa_generate_from_matrix`, unmodified; `:V3` →
`v3_retrain` above), standardize the resulting 100 synthetics with the SAME
canonical `std_params_concat` (baked into `Zreal`), score real records in the
standardized space by −(min distance to a holdout-retrain synthetic), and
return the Mann–Whitney-style MIA AUC (members = `train_idx`, non-members =
`holdout_idx`). Identical scoring/tie-handling to E3b.
"""
function split_mia_auc(scheme::Symbol, Zreal::Matrix{Float64}, X_concat::Matrix{Float64},
                        train_idx::Vector{Int}, holdout_idx::Vector{Int},
                        kept_cols::Vector{Symbol}, n_assays::Int;
                        gen_seed::Int=42, B::Int=B_V3, τ::Int=τ_V3)

    X_train = X_concat[train_idx, :]

    if scheme == :V0
        synth_df_ho, _mem_ho = sa_generate_from_matrix(X_train, kept_cols, n_assays, 100, T0_V0; seed=gen_seed)
    elseif scheme == :V3
        synth_df_ho = v3_retrain(X_train, kept_cols, n_assays, B, τ; gen_seed=gen_seed)
    else
        error("unknown scheme $scheme (expected :V0 or :V3)")
    end

    ho_synth_ids = sort(unique(synth_df_ho.SyntheticID))
    X_synth_ho = build_concat_matrix(synth_df_ho, :SyntheticID, ho_synth_ids, kept_cols, n_assays)

    Z_synth_ho = (X_synth_ho .- std_params_concat.μ') ./ std_params_concat.σ'

    Zreal_train   = Zreal[train_idx, :]
    Zreal_holdout = Zreal[holdout_idx, :]

    nn_dist_to_synth_ho(z) = minimum(norm(z - Z_synth_ho[k, :]) for k in 1:size(Z_synth_ho, 1))

    member_nn_dist    = [nn_dist_to_synth_ho(Zreal_train[i, :])   for i in 1:length(train_idx)]
    nonmember_nn_dist = [nn_dist_to_synth_ho(Zreal_holdout[i, :]) for i in 1:length(holdout_idx)]

    member_scores    = -member_nn_dist
    nonmember_scores = -nonmember_nn_dist

    auc_pairs = [m > n ? 1.0 : (m == n ? 0.5 : 0.0) for m in member_scores, n in nonmember_scores]
    return mean(auc_pairs)
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Fast sanity check — split_mia_auc(:V0, ...) on E3's FIXED
# 15/8 split must reproduce E3's AUC = 1.0 exactly (mirrors E3b's own Step-3
# self-check). Cheap (one retrain) — catches wiring bugs before the expensive
# repeated-split loop below.
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 3: Fast sanity check — split_mia_auc(:V0, E3 fixed split) …"
sorted_subjects = sort(complete_subjects)
e3_train_ids   = sorted_subjects[1:15]
e3_holdout_ids = sorted_subjects[16:23]
e3_train_idx   = [findfirst(==(s), complete_subjects) for s in e3_train_ids]
e3_holdout_idx = [findfirst(==(s), complete_subjects) for s in e3_holdout_ids]

e3_auc = split_mia_auc(:V0, Zreal, X_concat, e3_train_idx, e3_holdout_idx, kept_cols, n_assays; gen_seed=42)
@info "  split_mia_auc(:V0, E3 fixed split) = $(round(e3_auc, digits=4))  (expect 1.0)"
@assert isapprox(e3_auc, 1.0; atol=1e-9) "Sanity check FAILED: generalized split_mia_auc(:V0,...) does not reproduce E3's fixed-split AUC=1.0 (got $e3_auc). Aborting before trusting the repeated-split loop."
@info "  ✓ Fast sanity check PASSED — :V0 scheme matches E3/E3b exactly on the fixed split."

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Adaptive repeated-split MIA loop (factored, generalized over
# `scheme`) — same adaptive stopping rule as E3b: floor R=40, cap R=200, stop
# once SD/sqrt(R) < 0.01 after the floor. A FRESH MersenneTwister(rng_seed) is
# constructed per call so V0 and V3 draw the IDENTICAL partition sequence
# (common random numbers) — isolating the direction sampler as the only
# varying factor between the two AUC distributions, mirroring the S3 ladder's
# controlled-comparison design (and its V2/V3 common-random-numbers fix).
# ══════════════════════════════════════════════════════════════════════════════
function run_repeated_split_mia(scheme::Symbol, Zreal::Matrix{Float64}, X_concat::Matrix{Float64},
                                 complete_subjects::Vector, kept_cols::Vector{Symbol}, n_assays::Int;
                                 gen_seed::Int=42, B::Int=B_V3, τ::Int=τ_V3,
                                 rng_seed::Int=2026, floor_r::Int=40, cap_r::Int=200, se_tol::Float64=0.01)

    Kloc = length(complete_subjects)
    rng = MersenneTwister(rng_seed)

    rep_nums     = Int[]
    aucs         = Float64[]
    train_sets   = Vector{Vector{Int}}()
    holdout_sets = Vector{Vector{Int}}()

    stop_reason = "cap"
    final_R = cap_r

    R = 0
    while R < cap_r
        R += 1
        perm = randperm(rng, Kloc)
        train_idx   = sort(perm[1:15])
        holdout_idx = sort(perm[16:23])

        auc = split_mia_auc(scheme, Zreal, X_concat, train_idx, holdout_idx, kept_cols, n_assays;
                             gen_seed=gen_seed, B=B, τ=τ)
        @assert 0.0 <= auc <= 1.0 "[$scheme] AUC out of range at rep $R: $auc"

        push!(rep_nums, R)
        push!(aucs, auc)
        push!(train_sets, complete_subjects[train_idx])
        push!(holdout_sets, complete_subjects[holdout_idx])

        running_mean = mean(aucs)
        running_sd   = R > 1 ? std(aucs) : 0.0
        running_se   = running_sd / sqrt(R)

        @info "  [$scheme] Rep $R: AUC=$(round(auc, digits=4))  running mean=$(round(running_mean, digits=4))  SD=$(round(running_sd, digits=4))  SE=$(round(running_se, digits=5))"

        if R >= floor_r && running_se < se_tol
            stop_reason = "SE<$(se_tol)_after_floor"
            final_R = R
            break
        end
        if R == cap_r
            stop_reason = "cap"
            final_R = R
            break
        end
    end

    mean_auc = mean(aucs)
    sd_auc   = std(aucs)
    se_final = sd_auc / sqrt(final_R)

    return (scheme=scheme, final_R=final_R, stop_reason=stop_reason,
            aucs=aucs, rep_nums=rep_nums, train_sets=train_sets, holdout_sets=holdout_sets,
            mean_auc=mean_auc, sd_auc=sd_auc, se_final=se_final,
            min_auc=minimum(aucs), max_auc=maximum(aucs),
            p05=quantile(aucs, 0.05), p50=quantile(aucs, 0.50), p95=quantile(aucs, 0.95),
            frac_eq_1=mean(aucs .== 1.0), frac_gt_0p9=mean(aucs .> 0.9))
end

@info "Step 4: Repeated-split MIA — V0 (published scheme, the Step-2 self-check / wiring anchor) …"
v0_result = run_repeated_split_mia(:V0, Zreal, X_concat, complete_subjects, kept_cols, n_assays; gen_seed=42)
@info "  V0 done: final_R=$(v0_result.final_R) stop_reason=$(v0_result.stop_reason) mean=$(round(v0_result.mean_auc,digits=4)) sd=$(round(v0_result.sd_auc,digits=4))"

@info "Step 4b: SELF-CHECK — V0 repeated-split mean AUC must reproduce E3b ≈ 0.975 …"
@assert isapprox(v0_result.mean_auc, 0.975; atol=0.03) "V0 MIA does not reproduce E3b — sampler wiring wrong (got mean AUC=$(v0_result.mean_auc))"
@info "  ✓ V0 cross-check PASSED — mean AUC=$(round(v0_result.mean_auc,digits=4)) ≈ 0.975 (E3b), sampler wiring confirmed correct."

@info "Step 5: Repeated-split MIA — V3 (best-DCR rung: sphere-init MALA, pooled, burn-in B=$B_V3 τ=$τ_V3) …"
v3_result = run_repeated_split_mia(:V3, Zreal, X_concat, complete_subjects, kept_cols, n_assays; gen_seed=42)
@info "  V3 done: final_R=$(v3_result.final_R) stop_reason=$(v3_result.stop_reason) mean=$(round(v3_result.mean_auc,digits=4)) sd=$(round(v3_result.sd_auc,digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Assemble + write data/mcmc_ladder_mia.csv (per-rep + summary rows,
# Scheme ∈ {V0, V3})
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 6: Writing data/mcmc_ladder_mia.csv …"

function scheme_rows(res)
    kind  = String[]
    rep   = String[]
    auc   = Union{Float64,Missing}[]
    ntr   = Union{Int,Missing}[]
    nho   = Union{Int,Missing}[]
    extra = String[]

    for i in 1:res.final_R
        push!(kind, "rep")
        push!(rep, string(res.rep_nums[i]))
        push!(auc, res.aucs[i])
        push!(ntr, length(res.train_sets[i]))
        push!(nho, length(res.holdout_sets[i]))
        push!(extra, "train=" * join(res.train_sets[i], ";") * "|holdout=" * join(res.holdout_sets[i], ";"))
    end

    summary_pairs = [
        ("mean_auc", res.mean_auc), ("sd_auc", res.sd_auc),
        ("p05", res.p05), ("p50", res.p50), ("p95", res.p95),
        ("frac_auc_eq_1", res.frac_eq_1), ("frac_auc_gt_0.9", res.frac_gt_0p9),
        ("final_R", Float64(res.final_R)), ("se_final", res.se_final),
        ("min_auc", res.min_auc), ("max_auc", res.max_auc),
    ]
    for (label, val) in summary_pairs
        push!(kind, "summary"); push!(rep, label); push!(auc, val)
        push!(ntr, missing); push!(nho, missing); push!(extra, "")
    end
    push!(kind, "summary"); push!(rep, "stop_reason"); push!(auc, missing)
    push!(ntr, missing); push!(nho, missing); push!(extra, res.stop_reason)

    return DataFrame(Scheme=fill(string(res.scheme), length(kind)), Kind=kind, Rep=rep, AUC=auc,
                      N_train=ntr, N_holdout=nho, Extra=extra)
end

df_out = vcat(scheme_rows(v0_result), scheme_rows(v3_result))
CSV.write(joinpath(_PATH_TO_DATA, "mcmc_ladder_mia.csv"), df_out)
@info "  Wrote data/mcmc_ladder_mia.csv ($(nrow(df_out)) rows)"

# ── Invariant checks ──
@assert Set(unique(df_out.Scheme)) == Set(["V0", "V3"]) "expected Scheme ∈ {V0,V3} only"
rep_aucs = df_out[df_out.Kind .== "rep", :AUC]
@assert all(a -> 0.0 <= a <= 1.0, skipmissing(rep_aucs)) "all rep AUCs must be in [0,1]"
@assert isapprox(v0_result.mean_auc, 0.975; atol=0.03) "V0 cross-check must still hold at final assembly"
@info "All self-checks passed."

# ══════════════════════════════════════════════════════════════════════════════
# Final report — honest read-out (no "right answer" except the V0 anchor)
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^100)
println("Task S4 — Repeated-split MIA: V0 (published) vs V3 (best-DCR rung)")
println("="^100)
df_summary = DataFrame(
    Scheme=["V0", "V3"],
    final_R=[v0_result.final_R, v3_result.final_R],
    stop_reason=[v0_result.stop_reason, v3_result.stop_reason],
    mean_AUC=[v0_result.mean_auc, v3_result.mean_auc],
    sd_AUC=[v0_result.sd_auc, v3_result.sd_auc],
    p05=[v0_result.p05, v3_result.p05],
    p50=[v0_result.p50, v3_result.p50],
    p95=[v0_result.p95, v3_result.p95],
    frac_eq_1=[v0_result.frac_eq_1, v3_result.frac_eq_1],
    frac_gt_0p9=[v0_result.frac_gt_0p9, v3_result.frac_gt_0p9],
)
pretty_table(df_summary)

println("-"^100)
println("  V0 cross-check (wiring anchor): mean AUC = $(round(v0_result.mean_auc,digits=4)) " *
        "(target ≈0.975, E3b) → $(isapprox(v0_result.mean_auc, 0.975; atol=0.03) ? "PASSED" : "FAILED")")
println("-"^100)
verdict = v3_result.mean_auc < 0.6 ?
    "V3 mean AUC ≈ $(round(v3_result.mean_auc,digits=3)) (≈0.5) — the sphere-init/MALA/burn-in+pooled " *
    "sampler change FIXES the membership-inference leak relative to V0." :
    "V3 mean AUC ≈ $(round(v3_result.mean_auc,digits=3)) (still ≈0.9x) — the leak is INTRINSIC to π at β*; " *
    "proper MCMC mixing (near-zero discretization bias, multi-chain pooling past a validated burn-in) does " *
    "NOT rescue privacy."
println("  HONEST VERDICT: $verdict")
println("="^100)
println("  Output: data/mcmc_ladder_mia.csv")

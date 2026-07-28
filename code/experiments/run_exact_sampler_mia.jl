# ──────────────────────────────────────────────────────────────────────────────
# run_exact_sampler_mia.jl
#
# Task 2.3 of the exact-sampler update. E3b/S4 established that the
# membership-inference signal survives every sampler change tested (V0 mean AUC
# 0.975, V3 mean AUC 0.963 — proper MALA mixing past a validated burn-in does
# not rescue privacy). That left one alternative explanation open: both rungs
# are finite-time Markov chains, so the leak could still be a dynamics artifact
# rather than a property of the target.
#
# The analytic sampler closes that gap. Drawing ancestrally from
# pi = sum_k q_k N(m_k, beta^-1 I) has no chain, no burn-in, and no
# discretization bias. If the AUC survives here too, the signal is a property of
# the TARGET at K=23, which is what the manuscript already claims.
#
# Protocol: the SAME 40 train/holdout partitions used by V0 and V3. The split
# stream (MersenneTwister(2026), randperm) is regenerated and then asserted
# against the subject IDs frozen in mcmc_ladder_mia.csv, so an identical
# partition sequence is verified rather than assumed. Scoring is byte-identical
# to E3b: score = -(min distance to a holdout-retrain synthetic) in the canonical
# standardized space, AUC = Mann-Whitney over member/non-member pairs.
#
# Reads:  data/full_longitudinal_memory.jld2
#         data/mcmc_ladder_mia.csv        (V0/V3 splits + summaries; NOT modified)
# Writes: data/exact_sampler_mia.csv      (per-rep + summary rows, Scheme = "Exact")
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using Random, LinearAlgebra, Statistics

const N_SYNTH_MIA = 100
const GEN_SEED    = 42
const SPLIT_RNG   = 2026     # identical to run_mcmc_ladder_mia.jl's rng_seed
const N_TRAIN     = 15

unitize(v::Vector{Float64}) = v ./ (norm(v) + 1e-12)

# ══════════════════════════════════════════════════════════════════════════════
# Step 0/1: canonical memory, real matrix, standardization (mirrors S4 exactly)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 0: Loading canonical pipeline memory …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
dmem = JLD2.load(mem_path)
df_clean          = dmem["df_clean"]
kept_cols         = dmem["kept_cols"]
n_assays          = dmem["n_assays"]
complete_subjects = dmem["complete_subjects"]
std_params_concat = dmem["std_params_concat"]

K = length(complete_subjects)
n_visits = 3
@info "  K=$K, n_assays=$n_assays"

X_concat = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)
Zreal = (X_concat .- std_params_concat.μ') ./ std_params_concat.σ'
@info "  X_concat $(size(X_concat)), Zreal $(size(Zreal))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: exact-sampler retrain — mirrors v3_retrain, ancestral directions
# ══════════════════════════════════════════════════════════════════════════════
"""
    exact_retrain(X_train, kept_cols, n_assays; gen_seed) -> DataFrame

Refit standardization, PCA, unit memories and β* on the training split exactly
as V0/V3 do, then draw `N_SYNTH_MIA` directions ancestrally from the exact
target instead of running a chain. Decode is the published pushforward.
"""
function exact_retrain(X_train::Matrix{Float64}, kept_cols::Vector{Symbol}, n_assays::Int;
                       gen_seed::Int=GEN_SEED, pratio::Float64=0.95)

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

    β_star_sub = find_entropy_inflection(X̂_sub; α=0.01, n_betas=80, β_range=(0.1, 1000.0)).β_star

    # ── ancestral direction draw (no chain) ──
    out = exact_sample(X̂_sub, N_SYNTH_MIA; β=β_star_sub,
                       rng=MersenneTwister(gen_seed * 1000 + 777))
    dirs = Matrix{Float64}(undef, d_pca_sub, N_SYNTH_MIA)
    for i in 1:N_SYNTH_MIA
        dirs[:, i] = unitize(out.Ξ[i, :])
    end

    # Magnitudes reseeded independently — identical convention to v3_retrain.
    Random.seed!(gen_seed)
    magnitudes = [rand(pca_norms_sub) for _ in 1:N_SYNTH_MIA]

    concat = Matrix{Float64}(undef, N_SYNTH_MIA, n_assays * n_visits_sub)
    for i in 1:N_SYNTH_MIA
        concat[i, :] = decode_sample(dirs[:, i] .* magnitudes[i], pca_sub, std_params_sub)
    end

    df_visits = DataFrame[]
    for v in 1:n_visits_sub
        offset = (v - 1) * n_assays
        df_v = DataFrame(concat[:, (offset+1):(offset+n_assays)], kept_cols)
        df_v.SyntheticID = 1:N_SYNTH_MIA
        df_v.Visit = fill(v, N_SYNTH_MIA)
        push!(df_visits, df_v)
    end
    synth_df = vcat(df_visits...)
    sort!(synth_df, [:SyntheticID, :Visit])
    return synth_df
end

"""E3b scoring, unchanged: −(NN distance to a holdout-retrain synthetic), MW AUC."""
function exact_split_mia_auc(train_idx::Vector{Int}, holdout_idx::Vector{Int})
    synth_df_ho = exact_retrain(X_concat[train_idx, :], kept_cols, n_assays; gen_seed=GEN_SEED)
    ids = sort(unique(synth_df_ho.SyntheticID))
    X_synth_ho = build_concat_matrix(synth_df_ho, :SyntheticID, ids, kept_cols, n_assays)
    Z_synth_ho = (X_synth_ho .- std_params_concat.μ') ./ std_params_concat.σ'

    nn(z) = minimum(norm(z - Z_synth_ho[k, :]) for k in 1:size(Z_synth_ho, 1))
    member_scores    = [-nn(Zreal[i, :]) for i in train_idx]
    nonmember_scores = [-nn(Zreal[i, :]) for i in holdout_idx]

    pairs = [m > n ? 1.0 : (m == n ? 0.5 : 0.0) for m in member_scores, n in nonmember_scores]
    return mean(pairs)
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: regenerate the V0/V3 split sequence and VERIFY against the frozen CSV
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 3: Regenerating the S4 split sequence and checking it against mcmc_ladder_mia.csv …"
mia_path = joinpath(_PATH_TO_DATA, "mcmc_ladder_mia.csv")
isfile(mia_path) || error("data/mcmc_ladder_mia.csv not found — run run_mcmc_ladder_mia.jl (S4) first.")
df_mia = CSV.read(mia_path, DataFrame)

v0_reps = df_mia[(df_mia.Scheme .== "V0") .& (df_mia.Kind .== "rep"), :]
R_TARGET = nrow(v0_reps)
@info "  V0 contributed $R_TARGET repeated splits; matching that count."

rng = MersenneTwister(SPLIT_RNG)
splits = Tuple{Vector{Int},Vector{Int}}[]
for _ in 1:R_TARGET
    perm = randperm(rng, K)
    push!(splits, (sort(perm[1:N_TRAIN]), sort(perm[(N_TRAIN+1):K])))
end

for r in 1:R_TARGET
    recorded = v0_reps.Extra[r]
    tr_str = split(split(recorded, "|")[1], "=")[2]
    tr_recorded = parse.(Int, split(tr_str, ";"))
    tr_regenerated = complete_subjects[splits[r][1]]
    @assert tr_recorded == tr_regenerated "split $r differs from the frozen V0 partition:\n  recorded    = $tr_recorded\n  regenerated = $tr_regenerated"
end
@info "  ✓ All $R_TARGET partitions match the frozen V0/V3 splits exactly."

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: repeated-split MIA under the exact sampler
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Repeated-split MIA — Exact (analytic ancestral) over $R_TARGET splits …"
aucs = Float64[]
for r in 1:R_TARGET
    tr, ho = splits[r]
    auc = exact_split_mia_auc(tr, ho)
    @assert 0.0 <= auc <= 1.0 "AUC out of range at rep $r: $auc"
    push!(aucs, auc)
    r % 10 == 0 && @info "  … rep $r/$R_TARGET  running mean=$(round(mean(aucs),digits=4))"
end

mean_auc = mean(aucs)
sd_auc   = std(aucs)
se_auc   = sd_auc / sqrt(R_TARGET)
frac_1   = mean(aucs .== 1.0)
frac_09  = mean(aucs .> 0.9)

@info "\n  Exact: mean AUC=$(round(mean_auc,digits=4))  sd=$(round(sd_auc,digits=4))  " *
      "range [$(round(minimum(aucs),digits=4)), $(round(maximum(aucs),digits=4))]  " *
      "frac perfect=$(round(frac_1,digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: write exact_sampler_mia.csv (mcmc_ladder_mia.csv left untouched)
# ══════════════════════════════════════════════════════════════════════════════
kind  = String[]; rep = String[]; auc_col = Union{Float64,Missing}[]
ntr   = Union{Int,Missing}[]; nho = Union{Int,Missing}[]; extra = String[]

for r in 1:R_TARGET
    push!(kind, "rep"); push!(rep, string(r)); push!(auc_col, aucs[r])
    push!(ntr, N_TRAIN); push!(nho, K - N_TRAIN)
    push!(extra, "train=" * join(complete_subjects[splits[r][1]], ";") *
                 "|holdout=" * join(complete_subjects[splits[r][2]], ";"))
end
for (label, val) in [("mean_auc", mean_auc), ("sd_auc", sd_auc), ("se_final", se_auc),
                     ("p05", quantile(aucs, 0.05)), ("p50", quantile(aucs, 0.50)),
                     ("p95", quantile(aucs, 0.95)),
                     ("frac_auc_eq_1", frac_1), ("frac_auc_gt_0.9", frac_09),
                     ("min_auc", minimum(aucs)), ("max_auc", maximum(aucs)),
                     ("final_R", Float64(R_TARGET))]
    push!(kind, "summary"); push!(rep, label); push!(auc_col, val)
    push!(ntr, missing); push!(nho, missing); push!(extra, "")
end

df_out = DataFrame(Scheme=fill("Exact", length(kind)), Kind=kind, Rep=rep, AUC=auc_col,
                   N_train=ntr, N_holdout=nho, Extra=extra)
CSV.write(joinpath(_PATH_TO_DATA, "exact_sampler_mia.csv"), df_out)
@info "  Wrote data/exact_sampler_mia.csv ($(nrow(df_out)) rows)"

# ══════════════════════════════════════════════════════════════════════════════
# Report
# ══════════════════════════════════════════════════════════════════════════════
get_summary(scheme, label) = begin
    rows = df_mia[(df_mia.Scheme .== scheme) .& (df_mia.Kind .== "summary") .& (df_mia.Rep .== label), :]
    nrow(rows) == 1 ? rows.AUC[1] : missing
end

println("\n" * "="^88)
println("REPEATED-SPLIT MEMBERSHIP INFERENCE — chains vs the analytic reference")
println("="^88)
println(rpad("Scheme", 10), rpad("Sampler", 34), lpad("mean AUC", 10), lpad("sd", 9), lpad("frac=1.0", 11))
println("-"^88)
for (s, desc) in [("V0", "anchor+ULA+endpoint (published)"), ("V3", "sphere+MALA+pool")]
    println(rpad(s, 10), rpad(desc, 34),
            lpad(round(get_summary(s, "mean_auc"), digits=4), 10),
            lpad(round(get_summary(s, "sd_auc"), digits=4), 9),
            lpad(round(get_summary(s, "frac_auc_eq_1"), digits=4), 11))
end
println(rpad("Exact", 10), rpad("analytic ancestral (no chain)", 34),
        lpad(round(mean_auc, digits=4), 10), lpad(round(sd_auc, digits=4), 9),
        lpad(round(frac_1, digits=4), 11))
println("-"^88)
Δ = abs(mean_auc - get_summary("V0", "mean_auc"))
println("  |mean AUC(Exact) − mean AUC(V0)| = ", round(Δ, digits=4),
        "  (predeclared materially-close band: 0.05)")
println("  VERDICT: ", Δ <= 0.05 ?
    "the membership signal survives exact ancestral sampling — it is a property of the\n" *
    "           target at K=23, not of ULA discretization, anchor initialization, or finite-time mixing." :
    "the exact reference DIVERGES from V0 — do not claim sampler equivalence; report the discrepancy.")
println("="^88)

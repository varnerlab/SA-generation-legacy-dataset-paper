# ──────────────────────────────────────────────────────────────────────────────
# run_privacy_split_sensitivity.jl — Task 4b (E3b): Adaptive repeated-split
# membership-inference sensitivity analysis
#
# E3 (run_privacy_analysis.jl) found MIA_AUC = 1.0 on ONE deterministic
# 15-train/8-holdout split of the 23 complete subjects. This script asks:
# is that a robust property of the SA generator at this (K, β*, seed)
# operating point, or a one-split artifact?
#
# Method: factor the E3 MIA computation into `split_mia_auc(...)`, then repeat
# it over many RANDOM 15/8 partitions of the 23 complete subjects (single
# seeded RNG, MersenneTwister(2026), via `randperm`), holding the SA
# generation seed FIXED at 42 across every repetition — so partition identity
# is the ONLY thing that varies from rep to rep. An adaptive stopping rule
# halts the loop once the running mean AUC is estimated to Monte-Carlo
# precision SE = SD/sqrt(R) < 0.01 (with a floor of R=40 reps so the rule
# can't trivially stop after 1-2 reps on zero variance), or at a hard cap of
# R=200 reps.
#
# DCR (distance-to-closest-record) is NOT recomputed here — it is
# split-independent (computed once, canonically, in E3) and its two medians
# (synth→real 13.78, real→real 15.86) are only reprinted for context.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

using LinearAlgebra, Statistics, Random

# ══════════════════════════════════════════════════════════════════════════════
# Step 0: Load canonical pipeline memory + cleaned data (same as E3)
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

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Rebuild the 23×216 real concatenated matrix (canonical layout,
# identical to run_privacy_analysis.jl Step 1), then standardize → Zreal
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Rebuilding real concatenated matrix + standardizing …"
X_concat = Matrix{Float64}(undef, K, d_concat)
for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_clean.SubjectID .== s) .& (df_clean.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_concat[i, offset + j] = df_clean[row, col]
        end
    end
end
@info "  X_concat (real): $(size(X_concat))"

Zreal = (X_concat .- std_params_concat.μ') ./ std_params_concat.σ'
@info "  Zreal: $(size(Zreal))"

# small helper, identical convention to run_privacy_analysis.jl (not part of
# the SA generation pipeline itself, so it is not "copying the pipeline" —
# it just reshapes a long-format synthetic DataFrame into the K×216 layout)
function rebuild_concat_matrix(df_long::DataFrame, ids::AbstractVector, id_col::Symbol,
                                kept_cols::Vector{Symbol}, n_assays::Int, n_visits::Int)
    n = length(ids)
    X = Matrix{Float64}(undef, n, n_assays * n_visits)
    for (i, sid) in enumerate(ids)
        for v in 1:n_visits
            mask = (df_long[!, id_col] .== sid) .& (df_long.Visit .== v)
            row = findfirst(mask)
            offset = (v - 1) * n_assays
            for (j, col) in enumerate(kept_cols)
                X[i, offset + j] = df_long[row, col]
            end
        end
    end
    return X
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Factored MIA computation — split_mia_auc(...)
#
# Given a train/holdout partition (as index vectors into `complete_subjects`),
# retrain SA on the train subjects only, standardize the resulting synthetics
# with the SAME canonical std_params_concat, score every real record by
# −(min Euclidean distance to a holdout-retrain synthetic), and return the
# rank-statistic MIA AUC (members = train reals, non-members = holdout reals).
# This is the E3 Step 5-7 machinery, factored so it can be called in a loop
# with only the partition (and, incidentally, the gen_seed) varying.
# ══════════════════════════════════════════════════════════════════════════════
"""
    split_mia_auc(Zreal, X_concat, train_idx, holdout_idx, kept_cols, n_assays;
                  gen_seed=42) -> Float64

Retrain SA on `X_concat[train_idx, :]` via `sa_generate_from_matrix(...;
seed=gen_seed)`, standardize the resulting 100 synthetics with the SAME
`std_params_concat` implicitly baked into `Zreal` (i.e. `Zreal` must already
be standardized with the canonical params), score real records in the
standardized space by −(min distance to a holdout-retrain synthetic), and
return the Mann–Whitney-style MIA AUC (members = `train_idx`, non-members =
`holdout_idx`).
"""
function split_mia_auc(Zreal::Matrix{Float64}, X_concat::Matrix{Float64},
                        train_idx::Vector{Int}, holdout_idx::Vector{Int},
                        kept_cols::Vector{Symbol}, n_assays::Int;
                        gen_seed::Int=42)

    X_train = X_concat[train_idx, :]

    synth_df_ho, _mem_ho = sa_generate_from_matrix(X_train, kept_cols, n_assays, 100, 2000; seed=gen_seed)

    ho_synth_ids = sort(unique(synth_df_ho.SyntheticID))
    X_synth_ho = rebuild_concat_matrix(synth_df_ho, ho_synth_ids, :SyntheticID,
                                        kept_cols, n_assays, 3)

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
# Step 3: TDD self-check — split_mia_auc on the E3 FIXED split must reproduce
# AUC = 1.0 before the repeated-split loop is trusted.
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 3: Self-check — split_mia_auc on E3's fixed split reproduces AUC=1.0 …"
sorted_subjects = sort(complete_subjects)
e3_train_ids   = sorted_subjects[1:15]
e3_holdout_ids = sorted_subjects[16:23]
e3_train_idx   = [findfirst(==(s), complete_subjects) for s in e3_train_ids]
e3_holdout_idx = [findfirst(==(s), complete_subjects) for s in e3_holdout_ids]

e3_auc = split_mia_auc(Zreal, X_concat, e3_train_idx, e3_holdout_idx, kept_cols, n_assays; gen_seed=42)
@info "  split_mia_auc(E3 fixed split) = $(round(e3_auc, digits=4))  (expect 1.0)"
@assert isapprox(e3_auc, 1.0; atol=1e-9) "Self-check FAILED: factored split_mia_auc does not reproduce E3's AUC=1.0 on the fixed split (got $e3_auc). Aborting before trusting the repeated-split loop."
@info "  ✓ Self-check PASSED — factored function matches E3 exactly."

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Adaptive repeated-split loop
#
# rng = MersenneTwister(2026); each rep draws a fresh random 15/8 partition
# via randperm(rng, 23). gen_seed FIXED at 42 across every rep (only the
# partition varies). Stop when SD/sqrt(R) < 0.01 AND R >= 40 (floor), or at
# R == 200 (cap).
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 4: Adaptive repeated-split MIA loop (rng=MersenneTwister(2026), gen_seed=42 fixed) …"

rng = MersenneTwister(2026)
const FLOOR_R = 40
const CAP_R   = 200
const SE_TOL  = 0.01

rep_nums    = Int[]
aucs        = Float64[]
train_sets  = Vector{Vector{Int}}()
holdout_sets = Vector{Vector{Int}}()

stop_reason = "cap"  # default if we run all the way to CAP_R
final_R = CAP_R

R = 0
while R < CAP_R
    global R += 1
    perm = randperm(rng, K)                 # K == 23
    train_idx   = sort(perm[1:15])
    holdout_idx = sort(perm[16:23])

    auc = split_mia_auc(Zreal, X_concat, train_idx, holdout_idx, kept_cols, n_assays; gen_seed=42)
    @assert 0.0 <= auc <= 1.0 "AUC out of range at rep $R: $auc"

    push!(rep_nums, R)
    push!(aucs, auc)
    push!(train_sets, complete_subjects[train_idx])
    push!(holdout_sets, complete_subjects[holdout_idx])

    running_mean = mean(aucs)
    running_sd   = R > 1 ? std(aucs) : 0.0
    running_se   = running_sd / sqrt(R)

    @info "  Rep $R: AUC=$(round(auc, digits=4))  running mean=$(round(running_mean, digits=4))  SD=$(round(running_sd, digits=4))  SE=$(round(running_se, digits=5))"

    if R >= FLOOR_R && running_se < SE_TOL
        global stop_reason = "SE<$(SE_TOL)_after_floor"
        global final_R = R
        break
    end
    if R == CAP_R
        global stop_reason = "cap"
        global final_R = R
        break
    end
end

@info "Step 4 complete: final_R=$final_R, stop_reason=$stop_reason"

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Summary statistics over the R AUC values
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 5: Summary statistics …"

mean_auc = mean(aucs)
sd_auc   = std(aucs)
se_final = sd_auc / sqrt(final_R)
min_auc  = minimum(aucs)
max_auc  = maximum(aucs)
p05      = quantile(aucs, 0.05)
p50      = quantile(aucs, 0.50)
p95      = quantile(aucs, 0.95)
frac_eq_1   = mean(aucs .== 1.0)
frac_gt_0p9 = mean(aucs .> 0.9)

@info "  final_R                = $final_R"
@info "  stop_reason            = $stop_reason"
@info "  mean(AUC)              = $(round(mean_auc, digits=4))"
@info "  SD(AUC)                = $(round(sd_auc, digits=4))"
@info "  SE (final)             = $(round(se_final, digits=5))"
@info "  min / max               = $(round(min_auc, digits=4)) / $(round(max_auc, digits=4))"
@info "  p05 / p50 / p95         = $(round(p05, digits=4)) / $(round(p50, digits=4)) / $(round(p95, digits=4))"
@info "  frac(AUC == 1.0)        = $(round(frac_eq_1, digits=4))"
@info "  frac(AUC > 0.9)         = $(round(frac_gt_0p9, digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Context — reprint the canonical (split-independent) DCR medians
# from E3's data/privacy_results.csv summary rows (not recomputed here).
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 6: Reprinting canonical DCR medians from privacy_results.csv (context only) …"
median_dcr_sr = NaN
median_rr = NaN
privacy_csv_path = joinpath(_PATH_TO_DATA, "privacy_results.csv")
if isfile(privacy_csv_path)
    df_privacy = CSV.read(privacy_csv_path, DataFrame)
    summ = df_privacy[df_privacy.Kind .== "summary", :]
    row_sr = summ[summ.RecordID .== "median_dcr_synth_to_real", :Value]
    row_rr = summ[summ.RecordID .== "median_dcr_real_to_real", :Value]
    if !isempty(row_sr); global median_dcr_sr = row_sr[1]; end
    if !isempty(row_rr); global median_rr = row_rr[1]; end
    @info "  median(DCR synth→real) = $(round(median_dcr_sr, digits=4))  (from E3, split-independent)"
    @info "  median(DCR real→real)  = $(round(median_rr, digits=4))  (from E3, split-independent)"
else
    @warn "  data/privacy_results.csv not found — skipping DCR context (run E3 first for full context)."
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Write data/privacy_split_sensitivity.csv
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 7: Writing data/privacy_split_sensitivity.csv …"

rows_kind    = String[]
rows_field1  = String[]  # Rep (as string) or summary label
rows_auc     = Union{Float64,Missing}[]
rows_ntrain  = Union{Int,Missing}[]
rows_nholdout = Union{Int,Missing}[]
rows_extra   = String[]  # train/holdout subject IDs (rep rows) or blank

for i in 1:final_R
    push!(rows_kind, "rep")
    push!(rows_field1, string(rep_nums[i]))
    push!(rows_auc, aucs[i])
    push!(rows_ntrain, length(train_sets[i]))
    push!(rows_nholdout, length(holdout_sets[i]))
    push!(rows_extra, "train=" * join(train_sets[i], ";") * "|holdout=" * join(holdout_sets[i], ";"))
end

summary_pairs = [
    ("mean_auc", mean_auc),
    ("sd_auc", sd_auc),
    ("p05", p05),
    ("p50", p50),
    ("p95", p95),
    ("frac_auc_eq_1", frac_eq_1),
    ("frac_auc_gt_0.9", frac_gt_0p9),
    ("final_R", Float64(final_R)),
    ("se_final", se_final),
    ("min_auc", min_auc),
    ("max_auc", max_auc),
]
for (label, val) in summary_pairs
    push!(rows_kind, "summary")
    push!(rows_field1, label)
    push!(rows_auc, val)
    push!(rows_ntrain, missing)
    push!(rows_nholdout, missing)
    push!(rows_extra, "")
end
# stop_reason is a string, not a Float64 — its own row, AUC column left missing
push!(rows_kind, "summary")
push!(rows_field1, "stop_reason")
push!(rows_auc, missing)
push!(rows_ntrain, missing)
push!(rows_nholdout, missing)
push!(rows_extra, stop_reason)

df_out = DataFrame(Kind=rows_kind, Rep=rows_field1, AUC=rows_auc,
                    N_train=rows_ntrain, N_holdout=rows_nholdout, Extra=rows_extra)
CSV.write(joinpath(_PATH_TO_DATA, "privacy_split_sensitivity.csv"), df_out)

@info "\n✓ Done!"
@info "  final_R      = $final_R"
@info "  stop_reason  = $stop_reason"
@info "  mean ± SD    = $(round(mean_auc, digits=4)) ± $(round(sd_auc, digits=4))"
@info "  [p05, p95]   = [$(round(p05, digits=4)), $(round(p95, digits=4))]"
@info "  frac == 1.0  = $(round(frac_eq_1, digits=4))"
@info "  frac > 0.9   = $(round(frac_gt_0p9, digits=4))"
@info "  (context) median DCR synth→real=$(round(median_dcr_sr, digits=4)) vs real→real=$(round(median_rr, digits=4))"
@info "  Output: data/privacy_split_sensitivity.csv"

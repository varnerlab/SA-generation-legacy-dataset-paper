# ──────────────────────────────────────────────────────────────────────────────
# run_privacy_analysis.jl — E3: Membership-inference / privacy analysis
#
# Two diagnostics, both in the standardized [V1|V2|V3] concatenated space
# (23 × 216, standardize via the canonical `std_params_concat`):
#
#   1. DCR (distance to closest record), Euclidean:
#        synth_to_real   — each of the 100 canonical synthetics' min distance
#                           to the 23 reals
#        real_to_real    — each real's leave-one-out min distance to the
#                           other 22 reals (the "how close are reals to each
#                           other anyway" reference distribution)
#        holdout_to_train — each of the 8 holdout reals' min distance to the
#                           15 train reals (natural train/holdout gap, for
#                           context when interpreting the MIA AUC below)
#
#   2. NN-MIA AUC: deterministic 15-train / 8-holdout split of the 23
#      complete subjects; regenerate N=100 SA synthetics from the 15 TRAIN
#      subjects only (via the Task-1b helper `sa_generate_from_matrix`, NOT
#      the inline pipeline). An attacker scores each real record by
#      −(distance to nearest holdout-retrain synthetic); members = the 15
#      train reals, non-members = the 8 holdout reals. AUC = P(score_member
#      > score_nonmember) via the rank statistic (ties count as 1/2, the
#      standard Mann–Whitney U convention — matters only if two distances
#      tie exactly, which is measure-zero for continuous Euclidean distances
#      but is handled for correctness).
#
# Target/expected: median(dcr_sr) ≳ median(rr) (synthetics no closer to
# reals than reals are to each other) and MIA_AUC ≈ 0.5 (no membership
# leakage). Reported honestly either way — see task-4-report.md.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

using LinearAlgebra, Statistics, Random

# ══════════════════════════════════════════════════════════════════════════════
# Step 0: Load canonical pipeline memory + cleaned data
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 0: Loading canonical pipeline memory + cleaned data …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
d = JLD2.load(mem_path)
df_clean            = d["df_clean"]
kept_cols            = d["kept_cols"]
n_assays             = d["n_assays"]
complete_subjects    = d["complete_subjects"]
std_params_concat    = d["std_params_concat"]

K = length(complete_subjects)
n_visits = 3
d_concat = n_assays * n_visits
@info "  Complete subjects K=$K, n_assays=$n_assays, d_concat=$d_concat"

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Rebuild the 23×216 real concatenated matrix (canonical layout:
# offset=(v-1)*n_assays, kept_cols order — same as run_full_longitudinal.jl
# and test_sa_helper.jl), then standardize with std_params_concat → Zreal
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

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Load the canonical 100 synthetics, rebuild 100×216, standardize
# with the SAME std_params_concat → Zsyn
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 2: Rebuilding canonical synthetic matrix + standardizing …"
df_synth_canon = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
n_synth_canon = length(unique(df_synth_canon.SyntheticID))

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

synth_canon_ids = sort(unique(df_synth_canon.SyntheticID))
X_synth_canon = rebuild_concat_matrix(df_synth_canon, synth_canon_ids, :SyntheticID,
                                       kept_cols, n_assays, n_visits)
@info "  X_synth_canon: $(size(X_synth_canon))"

Zsyn = (X_synth_canon .- std_params_concat.μ') ./ std_params_concat.σ'
@info "  Zsyn: $(size(Zsyn))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: DCR distributions, Euclidean in standardized space
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 3: DCR distributions …"

# each column of A → Euclidean distance to its closest column in B
dcr(A, B) = [minimum(norm(a - b) for b in eachcol(B)) for a in eachcol(A)]

dcr_sr = dcr(permutedims(Zsyn), permutedims(Zreal))              # 100 values
rr     = [minimum(norm(Zreal[i, :] - Zreal[j, :]) for j in 1:K if j != i) for i in 1:K]  # 23 values

@assert length(dcr_sr) == n_synth_canon "expected $(n_synth_canon) synth_to_real records, got $(length(dcr_sr))"
@assert length(rr) == K "expected $K real_to_real records, got $(length(rr))"

median_dcr_sr = median(dcr_sr)
median_rr = median(rr)
@info "  median(DCR synth→real) = $(round(median_dcr_sr, digits=4))"
@info "  median(DCR real→real)  = $(round(median_rr, digits=4))"
sanity_ok = median_dcr_sr >= median_rr
@info "  Sanity: median(dcr_sr) ≳ median(rr)? $(sanity_ok)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Deterministic 15-train / 8-holdout split of the 23 complete subjects
#
# Fixed split (no unseeded `rand`): first 15 / last 8 of sort(complete_subjects).
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 4: Deterministic train/holdout split …"
sorted_subjects = sort(complete_subjects)
train_ids   = sorted_subjects[1:15]
holdout_ids = sorted_subjects[16:23]
@info "  train_ids   (n=$(length(train_ids)))   = $train_ids"
@info "  holdout_ids (n=$(length(holdout_ids))) = $holdout_ids"

train_idx   = [findfirst(==(s), complete_subjects) for s in train_ids]
holdout_idx = [findfirst(==(s), complete_subjects) for s in holdout_ids]

X_train = X_concat[train_idx, :]      # 15×216, RAW measurement space (helper standardizes internally)
@info "  X_train: $(size(X_train))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Holdout-retrain synthetics via the Task-1b helper (NOT the inline
# pipeline) — 100 synthetics generated from the 15 TRAIN subjects only.
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 5: Generating 100 holdout-retrain synthetics via sa_generate_from_matrix(...; seed=42) …"
synth_df_ho, mem_ho = sa_generate_from_matrix(X_train, kept_cols, n_assays, 100, 2000; seed=42)
@info "  Holdout-retrain: d_pca=$(mem_ho.d_pca), β*=$(round(mem_ho.β_star, digits=3))"

ho_synth_ids = sort(unique(synth_df_ho.SyntheticID))
X_synth_ho = rebuild_concat_matrix(synth_df_ho, ho_synth_ids, :SyntheticID,
                                    kept_cols, n_assays, n_visits)
@info "  X_synth_ho: $(size(X_synth_ho))"

# standardize the holdout-retrain synthetics with the SAME canonical
# std_params_concat used for Zreal/Zsyn, so member/non-member distances are
# computed in a common, fixed (non-holdout-dependent) space.
Z_synth_ho = (X_synth_ho .- std_params_concat.μ') ./ std_params_concat.σ'

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: holdout_to_train DCR — natural train/holdout gap (context for the
# MIA AUC below): each holdout real's min distance to the 15 train reals.
# ══════════════════════════════════════════════════════════════════════════════
Zreal_train   = Zreal[train_idx, :]
Zreal_holdout = Zreal[holdout_idx, :]
holdout_to_train = [minimum(norm(Zreal_holdout[i, :] - Zreal_train[j, :]) for j in 1:length(train_idx))
                     for i in 1:length(holdout_idx)]
@info "  median(DCR holdout→train) = $(round(median(holdout_to_train), digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: NN-MIA — score = −(distance to nearest holdout-retrain synthetic);
# members = 15 train reals, non-members = 8 holdout reals.
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 7: NN-MIA AUC …"

nn_dist_to_synth_ho(z) = minimum(norm(z - Z_synth_ho[k, :]) for k in 1:size(Z_synth_ho, 1))

member_nn_dist    = [nn_dist_to_synth_ho(Zreal_train[i, :])   for i in 1:length(train_idx)]
nonmember_nn_dist = [nn_dist_to_synth_ho(Zreal_holdout[i, :]) for i in 1:length(holdout_idx)]

member_scores    = -member_nn_dist
nonmember_scores = -nonmember_nn_dist

# AUC = P(score_member > score_nonmember), rank statistic with tie-handling
# (Mann-Whitney convention: ties count as 1/2; matters only for exact ties,
# measure-zero for continuous Euclidean distances but handled for correctness).
auc_pairs = [m > n ? 1.0 : (m == n ? 0.5 : 0.0) for m in member_scores, n in nonmember_scores]
mia_auc = mean(auc_pairs)

@assert 0.0 <= mia_auc <= 1.0 "MIA AUC out of range: $mia_auc"
@assert length(member_scores) == 15 "expected 15 member scores, got $(length(member_scores))"
@assert length(nonmember_scores) == 8 "expected 8 non-member scores, got $(length(nonmember_scores))"

@info "  MIA AUC = $(round(mia_auc, digits=4))"
mia_ok = isapprox(mia_auc, 0.5; atol=0.15)
@info "  Target ≈0.5 (no membership leakage)? within ±0.15: $(mia_ok)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 8: Write data/privacy_results.csv
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 8: Writing data/privacy_results.csv …"

rows_kind   = String[]
rows_record = String[]
rows_value  = Float64[]

for (i, sid) in enumerate(synth_canon_ids)
    push!(rows_kind, "synth_to_real"); push!(rows_record, "Syn$(sid)"); push!(rows_value, dcr_sr[i])
end
for (i, s) in enumerate(complete_subjects)
    push!(rows_kind, "real_to_real"); push!(rows_record, "R$(s)"); push!(rows_value, rr[i])
end
for (i, s) in enumerate(holdout_ids)
    push!(rows_kind, "holdout_to_train"); push!(rows_record, "H$(s)"); push!(rows_value, holdout_to_train[i])
end
for (i, s) in enumerate(train_ids)
    push!(rows_kind, "mia_member_nn_dist"); push!(rows_record, "Train$(s)"); push!(rows_value, member_nn_dist[i])
end
for (i, s) in enumerate(holdout_ids)
    push!(rows_kind, "mia_nonmember_nn_dist"); push!(rows_record, "Holdout$(s)"); push!(rows_value, nonmember_nn_dist[i])
end

# summary rows
summary_pairs = [
    ("MIA_AUC", mia_auc),
    ("median_dcr_synth_to_real", median_dcr_sr),
    ("median_dcr_real_to_real", median_rr),
    ("median_dcr_holdout_to_train", median(holdout_to_train)),
    ("n_train", Float64(length(train_ids))),
    ("n_holdout", Float64(length(holdout_ids))),
    ("n_synth_canon", Float64(n_synth_canon)),
]
for (label, val) in summary_pairs
    push!(rows_kind, "summary"); push!(rows_record, label); push!(rows_value, val)
end

df_out = DataFrame(Kind=rows_kind, RecordID=rows_record, Value=rows_value)
CSV.write(joinpath(_PATH_TO_DATA, "privacy_results.csv"), df_out)

@info "\n✓ Done!"
@info "  MIA AUC                    = $(round(mia_auc, digits=4))"
@info "  median(DCR synth→real)     = $(round(median_dcr_sr, digits=4))"
@info "  median(DCR real→real)      = $(round(median_rr, digits=4))"
@info "  median(DCR holdout→train)  = $(round(median(holdout_to_train), digits=4))"
@info "  Sanity dcr_sr ≳ rr: $sanity_ok   |  MIA_AUC ≈ 0.5 (±0.15): $mia_ok"
@info "  Output: data/privacy_results.csv"

# ──────────────────────────────────────────────────────────────────────────────
# run_pca_ablation.jl
#
# E1: PCA-space ablation suite. SA's generative loop has two conceptually
# separable pieces: (1) reduce the K=23 real patients to an 18-D PCA memory
# space, (2) sample new points in that space via the Hopfield/Langevin
# sampler, then decode. This script isolates piece (2): it fits FOUR simple,
# well-understood baseline samplers directly in the SAME raw PCA coordinate
# space SA uses (Y = X̂ .* pca_norms', d×23), draws N=100 samples from each,
# decodes them through the identical `decode_sample` pipeline, and evaluates
# all four baselines plus the canonical SA cohort through one identical
# harness:
#   - pooled median marginal relative error (MRE) + bootstrap CI
#   - cross-visit (216×216) correlation-matrix Frobenius distance
#   - BZ2012 mechanistic overlap/KS (TF-only), via the Task-1c helper
#   - Hopfield-memory novelty (1 - max cosine similarity to a stored pattern)
#
# If the baselines are just as good as SA on every axis, the Hopfield/Langevin
# sampler is not pulling its weight. Empirically (see data/pca_ablation_results.csv)
# the four baselines do NOT fail uniformly: PCA+kNN collapses to near-copies of
# the 23 real patients (novelty median 0.05, only 5% of samples novel — the
# classic SMOTE-style memorization failure), PCA+KDE has the worst marginal MRE
# and the weakest mechanistic-plausibility overlap/KS, while PCA+GMM and
# PCA+Copula are close competitors to SA on MRE, mechanistic overlap, and
# cross-visit correlation Frobenius distance. Notably, SA does NOT have the
# best cross-visit correlation match (all four PCA baselines that stay closer
# to the 23 training points reproduce the empirical 216×216 correlation matrix
# more exactly than SA does — SA's Frobenius distance of 0.64 matches the
# already-published, frozen `validate_cross_visit_covariance.jl` value exactly)
# — a bias/novelty-vs-fidelity trade-off: methods that memorize stay closer to
# the empirical correlation structure trivially, at the cost of near-zero
# novelty. SA is the only method that is simultaneously non-memorizing
# (novelty median 0.51, 100% of samples novel) AND competitive on MRE and
# mechanistic plausibility.
#
# Reads:  data/full_longitudinal_memory.jld2      (X̂, pca_concat, std_params_concat,
#                                                   pca_norms, complete_subjects, kept_cols,
#                                                   n_assays, df_clean)
#         data/synthetic_full_longitudinal.csv     (canonical N=100 SA cohort — SA row)
#         data/cleaned_full_data.csv (via df_clean already embedded in the JLD2)
# Writes: data/pca_ablation_results.csv
#
# Real/synthetic cohort selection for the mechanistic axis mirrors
# experiments/test_mechanistic_eval.jl exactly: complete-case (all 3 visits)
# real subjects, TF-only (TM=0.0) BZ2012 ratios via the Task-1c helper
# `mechanistic_eval.jl` (CALIBRATED_OVERRIDES, TF_TGA_COLS, bz2012_ratios,
# overlap_ks) — NOT a re-implementation of the harness.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using HockinMannModel, HypothesisTests
include(joinpath(@__DIR__, "mechanistic_eval.jl"))

using Random, LinearAlgebra, Statistics, Distributions

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Four baseline samplers, fit in the raw PCA coordinate space.
# Y :: d×23 raw PCA coords; returns d×N samples. rng passed for reproducibility.
# ══════════════════════════════════════════════════════════════════════════════

# (a) PCA + diagonal-covariance GMM (1–3 comps, variance floor); robust at n<p.
function sample_gmm(Y, N, rng; ncomp=2, floor_var=1e-3, iters=100)
    d, K = size(Y); pts = collect(eachcol(Y))
    # k-means++ init
    μ = [Y[:, rand(rng, 1:K)] for _ in 1:ncomp]
    σ2 = [fill(var(Y) + floor_var, d) for _ in 1:ncomp]; w = fill(1/ncomp, ncomp)
    for _ in 1:iters
        # E-step
        logr = zeros(K, ncomp)
        for k in 1:K, c in 1:ncomp
            logr[k,c] = log(w[c]) - 0.5*sum(log.(2π.*σ2[c])) - 0.5*sum((pts[k].-μ[c]).^2 ./ σ2[c])
        end
        r = exp.(logr .- maximum(logr, dims=2)); r ./= sum(r, dims=2)
        # M-step
        Nc = vec(sum(r, dims=1)); w = Nc ./ K
        for c in 1:ncomp
            Nc[c] < 1e-6 && continue
            μ[c] = sum(r[k,c].*pts[k] for k in 1:K) ./ Nc[c]
            σ2[c] = sum(r[k,c].*(pts[k].-μ[c]).^2 for k in 1:K) ./ Nc[c] .+ floor_var
        end
    end
    out = zeros(d, N)
    for n in 1:N
        c = rand(rng, Categorical(w)); out[:,n] = μ[c] .+ sqrt.(σ2[c]).*randn(rng, d)
    end
    out
end

# (b) PCA + Gaussian KDE sampler: pick a real point, add isotropic Gaussian noise (Scott bandwidth).
function sample_kde(Y, N, rng)
    d, K = size(Y)
    h = K^(-1/(d+4))                      # Scott's rule factor
    S = cov(Y'; corrected=true)           # d×d empirical cov of the 23 points
    L = cholesky(Symmetric(S) + 1e-8I).L
    out = zeros(d, N)
    for n in 1:N
        j = rand(rng, 1:K); out[:,n] = Y[:,j] .+ h .* (L*randn(rng, d))
    end
    out
end

# (c) PCA + k-NN interpolation (SMOTE-style): real point + λ·(neighbor − point), λ∼U(0,1).
function sample_knn(Y, N, rng; k=3)
    d, K = size(Y); out = zeros(d, N)
    D = [norm(Y[:,i]-Y[:,j]) for i in 1:K, j in 1:K]
    for n in 1:N
        i = rand(rng, 1:K)
        nbrs = sortperm(D[i,:])[2:k+1]     # exclude self
        j = rand(rng, nbrs); λ = rand(rng)
        out[:,n] = Y[:,i] .+ λ.*(Y[:,j].-Y[:,i])
    end
    out
end

# (d) PCA + Gaussian copula: empirical marginal CDF per component + Gaussian correlation.
function sample_copula(Y, N, rng)
    d, K = size(Y)
    U = similar(Y)                        # ranks → (0,1)
    for c in 1:d
        U[c,:] = (invperm(sortperm(Y[c,:])) .- 0.5) ./ K
    end
    Z = quantile.(Normal(), clamp.(U, 1/(2K), 1-1/(2K)))
    R = cor(Z'); R = Symmetric(R) + 1e-6I
    L = cholesky(R).L
    out = zeros(d, N)
    for n in 1:N
        z = L*randn(rng, d); u = cdf.(Normal(), z)
        for c in 1:d                       # inverse empirical marginal (nearest-rank quantile)
            out[c,n] = quantile(Y[c,:], clamp(u[c], 0.0, 1.0))
        end
    end
    out
end

# ══════════════════════════════════════════════════════════════════════════════
# Load canonical SA artifacts
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading canonical SA memory + cohorts …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms complete_subjects kept_cols n_assays df_clean

d_pca, K = size(X̂)
n_visits = length(std_params_concat.col_names) ÷ n_assays
@assert length(kept_cols) == n_assays
@info "  X̂: $(size(X̂))  (d=$d_pca, K=$K real patients)   n_assays=$n_assays  n_visits=$n_visits"

Y = X̂ .* pca_norms'   # d×23 raw (un-normalized) PCA coordinates
@info "  Raw PCA coords Y: $(size(Y))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: self-check (shape / finiteness) — the "failing test" first.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nSelf-check: sampler shape / finiteness …"
rng = MersenneTwister(42)
@assert size(sample_gmm(Y, 100, rng)) == (size(Y, 1), 100)
@assert all(isfinite, sample_kde(Y, 100, rng))
@assert all(isfinite, sample_knn(Y, 100, rng))
@assert all(isfinite, sample_copula(Y, 100, rng))
@info "  ✓ all four samplers produce finite d×100 output"

# ══════════════════════════════════════════════════════════════════════════════
# Shared cohort prep — mirrors experiments/test_mechanistic_eval.jl exactly:
# complete-case (all 3 visits present) real subjects for the mechanistic axis
# and the cross-visit correlation reference.
# ══════════════════════════════════════════════════════════════════════════════
@info "\nShared cohort prep …"
df_real_complete = df_clean[in.(df_clean.SubjectID, Ref(complete_subjects)), :]
@info "  Real complete-case cohort: $(length(complete_subjects)) patients, $(nrow(df_real_complete)) records"

@info "  Running BZ2012 (TF-only, calibrated) on the real cohort …"
real_ratios = bz2012_ratios(df_real_complete, TF_TGA_COLS; TM=0.0)
@info "  real_ratios: $(nrow(real_ratios)) rows"

"""
    longdf_to_concat(df, id_col, kept_cols, n_assays, n_visits) -> Matrix (N × n_assays*n_visits)

Reshape a long-format (one row per patient×visit) DataFrame into a
concatenated-visit matrix, columns ordered V1 block, V2 block, V3 block —
same layout as `run_full_longitudinal.jl`'s `X_concat` / `sa_generate_from_matrix`.
"""
function longdf_to_concat(df::DataFrame, id_col::Symbol, kept_cols::Vector{Symbol},
                           n_assays::Int, n_visits::Int)
    ids = sort(unique(df[!, id_col]))
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

"""
    safe_cor(X) -> Matrix

Complete-case correlation matrix with NaN (from zero-variance columns)
mapped to 0.0. Pattern from validate_cross_visit_covariance.jl:93.
"""
function safe_cor(X)
    C = cor(X)
    C[isnan.(C)] .= 0.0
    return C
end

X_real = longdf_to_concat(df_real_complete, :SubjectID, kept_cols, n_assays, n_visits)
C_real = safe_cor(X_real)
@info "  Real concat matrix: $(size(X_real)); correlation matrix $(size(C_real))"

# ══════════════════════════════════════════════════════════════════════════════
# Decoding + evaluation helpers
# ══════════════════════════════════════════════════════════════════════════════

"""
    decode_to_longdf(Ysamp, pca_model, std_params, kept_cols, n_assays, n_visits)
        -> (df_long, concat)

Decode each column of a d×N raw-PCA-space sample matrix via
`Patient.decode_sample`, then reshape into the 3-visit long format (columns =
`kept_cols` + `SyntheticID` + `Visit`) — identical layout to
`synthetic_full_longitudinal.csv`. Also returns the underlying N×d_concat
matrix (reused for the cross-visit correlation axis).
"""
function decode_to_longdf(Ysamp::Matrix{Float64}, pca_model, std_params,
                           kept_cols::Vector{Symbol}, n_assays::Int, n_visits::Int)
    d, N = size(Ysamp)
    concat = Matrix{Float64}(undef, N, n_assays * n_visits)
    for i in 1:N
        concat[i, :] = decode_sample(Ysamp[:, i], pca_model, std_params)
    end
    df_visits = DataFrame[]
    for v in 1:n_visits
        offset = (v - 1) * n_assays
        cols = concat[:, (offset + 1):(offset + n_assays)]
        df_v = DataFrame(cols, kept_cols)
        df_v.SyntheticID = 1:N
        df_v.Visit = fill(v, N)
        push!(df_visits, df_v)
    end
    df = vcat(df_visits...)
    sort!(df, [:SyntheticID, :Visit])
    return df, concat
end

"""
    pca_directions(concat) -> Matrix (d × N)

Round-trip a decoded N×d_concat measurement-space matrix back through the
canonical standardization + PCA transform to recover its raw PCA coordinates
(exact up to floating point, since `decode_sample` = reconstruct∘destandardize
and this is transform∘standardize — their composition is the identity on the
PCA subspace). Used so novelty is computed identically for SA and all four
baselines, whether or not the raw sampled PCA vector was already at hand.
"""
function pca_directions(concat::Matrix{Float64})
    Z = (concat .- std_params_concat.μ') ./ std_params_concat.σ'
    return MultivariateStats.transform(pca_concat, Z')  # d×N
end

"""
    mre_ci(df_synth) -> (median, ci_lo, ci_hi)

Pooled per-(visit,feature) Mean_Rel_Error via `feature_summary_comparison`
looped over the 3 visits, pooled median with a 10,000-resample bootstrap 95%
CI (percentile indices 250/9750) — replicates paper_summary_table.jl:100-112
(seed 42 before each method's bootstrap, for reproducibility).
"""
function mre_ci(df_synth::DataFrame)
    all_mre = Float64[]
    for v in 1:n_visits
        real_v = df_real_complete[df_real_complete.Visit .== v, :]
        synth_v = df_synth[df_synth.Visit .== v, :]
        comp = feature_summary_comparison(real_v, synth_v, kept_cols)
        append!(all_mre, comp.Mean_Rel_Error)
    end
    n = length(all_mre)
    med = median(all_mre)
    Random.seed!(42)
    bs_medians = [median(rand(all_mre, n)) for _ in 1:10_000]
    sort!(bs_medians)
    return med, bs_medians[250], bs_medians[9750]
end

"""
    novelty_stats(Ypca) -> (median_novelty, frac_novelty_gt_0p2)

For each sample's unit-normalized PCA direction, `sample_novelty(ξ̂, X̂)`
(1 - max cosine similarity to a stored memory pattern).
"""
function novelty_stats(Ypca::Matrix{Float64})
    N = size(Ypca, 2)
    nov = Vector{Float64}(undef, N)
    for i in 1:N
        y = Ypca[:, i]
        ξ̂ = y ./ (norm(y) + 1e-12)
        nov[i] = sample_novelty(ξ̂, X̂)
    end
    return median(nov), mean(nov .> 0.2)
end

"""
    evaluate_method(name, df_synth, concat_synth) -> NamedTuple

Runs the full 4-axis harness (MRE+CI, cross-visit Frobenius, BZ2012
mechanistic overlap/KS TF-only, novelty) on one method's decoded N×d_concat
synthetic cohort.
"""
function evaluate_method(name::String, df_synth::DataFrame, concat_synth::Matrix{Float64})
    med, lo, hi = mre_ci(df_synth)

    C_synth = safe_cor(concat_synth)
    frob = norm(C_synth - C_real) / norm(C_real)

    n_synth = length(unique(df_synth.SyntheticID))
    @info "  [$name] BZ2012 mechanistic eval (TF-only) on $n_synth synthetic patients …"
    synth_ratios = bz2012_ratios(df_synth, TF_TGA_COLS; TM=0.0)
    frac, ksD, ksp = overlap_ks(real_ratios.Ratio, synth_ratios.Ratio)

    Ypca = pca_directions(concat_synth)
    med_nov, frac_nov = novelty_stats(Ypca)

    @info "  [$name] MedianMRE=$(round(100*med,digits=2))%  Frob=$(round(frob,digits=4))  " *
          "Overlap=$(round(frac,digits=3))  KS D=$(round(ksD,digits=3)) p=$(round(ksp,digits=3))  " *
          "Novelty(median)=$(round(med_nov,digits=3)) frac>0.2=$(round(frac_nov,digits=3))"

    return (Method=name, Median_MRE=med, MRE_CI_lo=lo, MRE_CI_hi=hi,
            CrossVisit_Frob=frob, Mech_Overlap_TFonly=frac, Mech_KS_D=ksD, Mech_KS_p=ksp,
            Median_Novelty=med_nov, Frac_Novel_gt_0p2=frac_nov)
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 3/4: Run the harness for SA (canonical N=100 cohort) and the four
# PCA-space baselines (each drawing its own N=100 sample from Y).
# ══════════════════════════════════════════════════════════════════════════════
results = NamedTuple[]

@info "\n=== SA (canonical, data/synthetic_full_longitudinal.csv) ==="
df_synth_sa = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
concat_sa = longdf_to_concat(df_synth_sa, :SyntheticID, kept_cols, n_assays, n_visits)
push!(results, evaluate_method("SA", df_synth_sa, concat_sa))

baseline_samplers = [
    ("PCA+GMM",    (rng) -> sample_gmm(Y, 100, rng)),
    ("PCA+KDE",    (rng) -> sample_kde(Y, 100, rng)),
    ("PCA+kNN",    (rng) -> sample_knn(Y, 100, rng)),
    ("PCA+Copula", (rng) -> sample_copula(Y, 100, rng)),
]

for (name, sampler) in baseline_samplers
    @info "\n=== $name ==="
    rng_m = MersenneTwister(42)   # independent, reproducible per-method seed
    Ysamp = sampler(rng_m)
    @assert size(Ysamp) == (d_pca, 100) "$name: sampler returned $(size(Ysamp)), expected ($d_pca, 100)"
    @assert all(isfinite, Ysamp) "$name: sampler produced non-finite values"
    df_synth, concat_synth = decode_to_longdf(Ysamp, pca_concat, std_params_concat, kept_cols, n_assays, n_visits)
    push!(results, evaluate_method(name, df_synth, concat_synth))
end

# ══════════════════════════════════════════════════════════════════════════════
# Assemble and write the results table
# ══════════════════════════════════════════════════════════════════════════════
df_out = DataFrame(
    Method               = [r.Method for r in results],
    Median_MRE           = [r.Median_MRE for r in results],
    MRE_CI_lo            = [r.MRE_CI_lo for r in results],
    MRE_CI_hi            = [r.MRE_CI_hi for r in results],
    CrossVisit_Frob      = [r.CrossVisit_Frob for r in results],
    Mech_Overlap_TFonly  = [r.Mech_Overlap_TFonly for r in results],
    Mech_KS_D            = [r.Mech_KS_D for r in results],
    Mech_KS_p            = [r.Mech_KS_p for r in results],
    Median_Novelty       = [r.Median_Novelty for r in results],
    Frac_Novel_gt_0p2    = [r.Frac_Novel_gt_0p2 for r in results],
)

CSV.write(joinpath(_PATH_TO_DATA, "pca_ablation_results.csv"), df_out)

println("\n" * "="^100)
println("E1: PCA-space ablation suite — SA vs. baselines fit in the same 18-D raw PCA coordinate space")
println("="^100)
pretty_table(df_out)

# ══════════════════════════════════════════════════════════════════════════════
# Invariant checks (Step 4)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nInvariant checks …"
@assert all(df_out.Median_MRE .> 0) "every Median_MRE must be > 0"
@assert all(0 .<= df_out.Mech_Overlap_TFonly .<= 1) "Mech_Overlap_TFonly must be in [0,1]"
@assert all(0 .<= df_out.Median_Novelty .<= 1) "Median_Novelty must be in [0,1]"

sa_mre = df_out.Median_MRE[findfirst(==("SA"), df_out.Method)]
@info "  SA pooled median MRE = $(round(100*sa_mre, digits=2))% (cross-check target ≈ 1.2%)"
@assert isapprox(sa_mre, 0.012; atol=0.01) "SA MRE cross-check failed: got $(100*sa_mre)%, expected ≈1.2%"

@info "✓ Wrote data/pca_ablation_results.csv"

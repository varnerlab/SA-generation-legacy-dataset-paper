# ──────────────────────────────────────────────────────────────────────────────
# run_sim_lowrank.jl
#
# SIM benchmark, arm 1: low-rank Gaussian ground-truth covariance recovery.
#
# Generates data from a KNOWN population covariance Σpop via a latent factor
# model in the n<p regime, then asks whether SA and MVN can recover Σpop from
# a small training cohort. Because Σpop is independent of the training draw,
# recovering it is genuine generalization — unlike the empirical cross-visit
# covariance metric (E1), which rewards memorization of the real training
# points. SA is expected to degrade gracefully where MVN (rank-deficient
# sample covariance when n<=p) fails outright. Report the truth regardless.
#
# Writes: data/sim_lowrank_results.csv
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "sim_common.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: self-check (TDD) — sample covariance of a big low-rank draw ≈ Σpop
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: self-check — cov_frob(cov(sampler(2000)), Σpop) < 0.15 …"
rng_check = MersenneTwister(42)
sampler_check, Σpop_check = make_lowrank_gt(30, 5, rng_check)
check_val = cov_frob(cov(sampler_check(2000)), Σpop_check)
@info "  cov_frob = $check_val"
@assert check_val < 0.15 "Step-1 self-check FAILED: cov_frob=$check_val ≥ 0.15 — GT or metric mis-wired"
@info "  Step-1 self-check PASSED ✓"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: sweep the grid
#   n ∈ {8,15,23,40,80}, p ∈ {30,120,216}, r ∈ {3,5,10}, tail ∈ {:gaussian,:heavy}
#   90 cells × {SA, MVN}. Canonical hyperparameters: α=0.01, T=2000, pratio=0.95,
#   N=100, β=nothing (recomputed per cohort inside sa_generate_from_matrix).
#   Runtime ~10–30 min is expected (90 cells × 100-chain SA + per-cohort β*
#   search each) — do NOT shrink N, T, or the grid to speed this up.
# ══════════════════════════════════════════════════════════════════════════════
n_grid = [8, 15, 23, 40, 80]
p_grid = [30, 120, 216]
r_grid = [3, 5, 10]
tail_grid = [:gaussian, :heavy]

total_cells = length(n_grid) * length(p_grid) * length(r_grid) * length(tail_grid)
@info "Step 2: sweeping $total_cells cells × {SA, MVN} …"

rows = NamedTuple[]
cell_i = 0
t_start = time()

for (i_p, p) in enumerate(p_grid), (i_r, r) in enumerate(r_grid), (i_tail, tail) in enumerate(tail_grid)
    for (i_n, n) in enumerate(n_grid)
        global cell_i += 1
        seed = cell_seed(i_n, i_p, i_r, i_tail)
        rng = MersenneTwister(seed)

        sampler, Σpop = make_lowrank_gt(p, r, rng; tail=tail)
        Xtrain = sampler(n)
        Xtest = sampler(5000)

        elapsed = round(time() - t_start, digits=1)
        @info "[$cell_i/$total_cells] n=$n p=$p r=$r tail=$tail seed=$seed  (elapsed $(elapsed)s)"

        Xsyn_sa = sa_recover(Xtrain; T=2000, seed=seed)
        Xsyn_mvn = mvn_recover(Xtrain; seed=seed)

        for (method, Xsyn) in (("SA", Xsyn_sa), ("MVN", Xsyn_mvn))
            cf = cov_frob(cov(Xsyn), Σpop)
            mre = marginal_mre(Xsyn, Xtest)
            ks = marginal_ks(Xsyn, Xtest)
            push!(rows, (n=n, p=p, r=r, tail=String(tail), Method=method,
                          Cov_Frob=cf, Marginal_MRE=mre, Marginal_KS=ks))
        end
    end
end

@info "Step 2 complete: $(length(rows)) rows in $(round(time()-t_start,digits=1))s"

results = DataFrame(rows)

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: write CSV + verify the intended contrast
# ══════════════════════════════════════════════════════════════════════════════
out_path = joinpath(_PATH_TO_DATA, "sim_lowrank_results.csv")
CSV.write(out_path, results)
@info "Step 3: wrote $out_path ($(nrow(results)) rows)"

# pivot: mean Cov_Frob by (n,p), averaged over the r×tail replicates, SA vs MVN
agg = combine(groupby(results, [:n, :p, :Method]), :Cov_Frob => mean => :Mean_Cov_Frob)
pivot = DataFrame(n=Int[], p=Int[], SA=Float64[], MVN=Float64[], SA_minus_MVN=Float64[])
for p in p_grid, n in n_grid
    sa_v = only(agg[(agg.n .== n) .& (agg.p .== p) .& (agg.Method .== "SA"), :Mean_Cov_Frob])
    mvn_v = only(agg[(agg.n .== n) .& (agg.p .== p) .& (agg.Method .== "MVN"), :Mean_Cov_Frob])
    push!(pivot, (n=n, p=p, SA=sa_v, MVN=mvn_v, SA_minus_MVN=sa_v - mvn_v))
end
sort!(pivot, [:p, :n])

println()
println("═"^78)
println("Pivot: mean Cov_Frob by (n,p) — SA vs MVN (averaged over r ∈ {3,5,10}, tail ∈ {gaussian,heavy})")
println("═"^78)
pretty_table(pivot)

println()
for row in eachrow(pivot)
    regime = row.n <= row.p ? "n≤p (rank-deficient MVN)" : "n>p"
    verdict = row.SA < row.MVN ? "SA beats MVN" : "MVN beats/ties SA"
    rel = row.MVN > 0 ? round(100 * (row.MVN - row.SA) / row.MVN, digits=1) : NaN
    println("  n=$(row.n) p=$(row.p) [$regime]: SA=$(round(row.SA,digits=4)) MVN=$(round(row.MVN,digits=4)) — $verdict (SA vs MVN: $(rel)% lower Cov_Frob)")
end

# marginal MRE / KS summary, SA vs MVN, gaussian vs heavy tail
println()
println("═"^78)
println("Marginal MRE / KS summary (mean over all cells), by tail")
println("═"^78)
marg = combine(groupby(results, [:tail, :Method]),
               :Marginal_MRE => mean => :Mean_Marginal_MRE,
               :Marginal_KS => mean => :Mean_Marginal_KS)
pretty_table(marg)

@info "run_sim_lowrank.jl complete."

# ──────────────────────────────────────────────────────────────────────────────
# run_sim_nonlinear.jl
#
# SIM benchmark, arm 3 (FINAL): nonlinear-manifold ground truth — the
# linear-PCA applicability boundary.
#
# SA builds its memory subspace with linear PCA. On a NONLINEAR manifold (an
# S-curve), linear methods place generated points in the convex region
# "across the bend," off the true curve. This arm quantifies that: as the
# manifold's curvature rises, SA's off-manifold residual grows, locating the
# honest boundary where the linear-PCA assumption starts to matter (evidences
# reviewer point R2.2 / nonlinearity as a GT probe). At curvature 0 the
# manifold is a straight line and the residual is ~0 (the linear-matched
# baseline). Report the truth: where does the residual become material?
#
# Writes: data/sim_nonlinear_results.csv (kind, curvature, p, Rep, Method,
#                                          OffManifold_Err, Cov_Frob)
#
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "sim_common.jl"))
using LinearAlgebra, Statistics, Random

# ══════════════════════════════════════════════════════════════════════════════
# Ground truth: S-curve manifold + off-manifold-plane metric
# (transcribed verbatim from the task brief — do not deviate)
# ══════════════════════════════════════════════════════════════════════════════

# Golden-section minimizer on [lo,hi] for a unimodal-ish f (used to polish the nearest-t search).
function _golden_min(f, lo, hi; tol=1e-9, iters=200)
    φ = (√5 - 1) / 2
    a, b = lo, hi
    c = b - φ*(b-a); d = a + φ*(b-a); fc = f(c); fd = f(d)
    for _ in 1:iters
        (b-a) < tol && break
        if fc < fd; b,d,fd = d,c,fc; c = b - φ*(b-a); fc = f(c)
        else;       a,c,fc = c,d,fd; d = a + φ*(b-a); fd = f(d); end
    end
    return (a+b)/2
end

# S-curve GT: 1D latent t∈[-1,1] → 2D curve [t, curvature·sin(π t)], embedded into p-dim by a
# random ORTHONORMAL p×2 map Q, plus isotropic noise. curvature=0 ⇒ a straight line (linear baseline).
function make_nonlinear_gt(p, rng; kind=:scurve, curvature=1.0, noise=0.05)
    kind == :scurve || throw(ArgumentError("only :scurve implemented"))
    Q = Matrix(qr(randn(rng, p, 2)).Q)[:, 1:2]          # p×2, orthonormal columns
    mfun(t) = [t, curvature * sin(π * t)]               # 2D manifold coordinate
    function sampler(n)
        t = 2 .* rand(rng, n) .- 1                       # U(-1,1)
        M = reduce(vcat, (mfun(ti)' for ti in t))        # n×2 curve coords
        return M * Q' .+ noise .* randn(rng, n, p)       # n×p (on the embedded curve + noise)
    end
    # population covariance of what sampler emits (non-Gaussian ⇒ estimate from a large draw)
    Σpop = cov(sampler(50_000))
    ctx = (Q=Q, mfun=mfun, tgrid=collect(range(-1, 1; length=2000)))
    return sampler, Σpop, ctx
end

# Off-manifold residual, measured WITHIN the 2D embedding plane (so inherent orthogonal noise is NOT
# counted): project each row to plane coords y=Qᵀx, find the nearest point on the curve (coarse grid +
# golden-section polish), return the mean residual. ~0 for on-curve points, grows as SA interpolates
# across the bend.
function offmanifold_error(X, ctx)
    Q, mfun, tgrid = ctx.Q, ctx.mfun, ctx.tgrid
    Mgrid = [mfun(t) for t in tgrid]
    ds = Vector{Float64}(undef, size(X,1))
    for i in 1:size(X,1)
        y = Q' * @view X[i, :]                           # 2D projection into the embedding plane
        j = argmin([sum(abs2, y .- m) for m in Mgrid])   # coarse nearest grid point
        lo = tgrid[max(j-1,1)]; hi = tgrid[min(j+1,length(tgrid))]
        t★ = _golden_min(τ -> sum(abs2, y .- mfun(τ)), lo, hi)
        ds[i] = norm(y .- mfun(t★))
    end
    return mean(ds)
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: self-check (TDD) — the off-manifold metric detects nonlinearity
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: self-check — offmanifold_error on-curve ≈ 0; SA(curved) ≥ SA(straight) …"

rng = MersenneTwister(42)
sampler_c, Σc, ctx_c = make_nonlinear_gt(30, rng; curvature=1.5)
# on-curve points (no noise) → residual ≈ 0
t = 2 .* rand(rng, 200) .- 1
true_pts = reduce(vcat, (ctx_c.mfun(ti)' for ti in t)) * ctx_c.Q'   # 200×30 exactly on the curve
onc = offmanifold_error(true_pts, ctx_c)
@info "  on-curve residual = $onc"
@assert onc < 1e-6 "Step-1 self-check FAILED: on-curve residual=$onc ≥ 1e-6"

rng2 = MersenneTwister(7)
samp_lin, _, ctx_lin = make_nonlinear_gt(30, rng2; curvature=0.0)
Xtr_lin  = samp_lin(23);  err_lin  = offmanifold_error(sa_recover(Xtr_lin;  seed=7), ctx_lin)
Xtr_curv = sampler_c(23); err_curv = offmanifold_error(sa_recover(Xtr_curv; seed=7), ctx_c)
@info "  err_lin (curvature=0) = $err_lin ; err_curv (curvature=1.5) = $err_curv"
@assert err_curv >= err_lin "Step-1 self-check FAILED: err_curv=$err_curv < err_lin=$err_lin"

@info "  Step-1 self-check PASSED ✓"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: sweep curvature
#   curvature ∈ {0.0, 0.25, 0.5, 1.0, 1.5, 2.0} (0.0 = matched linear baseline),
#   p ∈ {30, 120} (both divisible by 3, so sa_recover's p÷3 reshape works; a
#   quick timing probe (p=120, n=23: sa_recover 2.53s, mvn_recover 0.39s) showed
#   p=120 is cheap here — unlike arms 1-2 there is no per-cohort ODE — so it is
#   included per the brief's "MAY add p=120 if runtime allows"), kind=:scurve,
#   Rep ∈ 1:5, n=23 (operating cohort size). Canonical hyperparameters: α=0.01,
#   T=2000, pratio=0.95, N=100, β=nothing (recomputed per cohort).
# ══════════════════════════════════════════════════════════════════════════════
curvature_grid = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0]
p_grid = [30, 120]
n = 23

"""
    nonlinear_seed(i_curv, i_p, rep; base=42) -> Int

Deterministic seed from 1-based (curvature-index, p-index, Rep) grid-cell
indices (not `rand()`), anchored to the canonical seed 42, mirroring
`sim_common.jl`'s `cell_seed` pattern but sized for this task's 3-D grid.
Weights (10, 100, 10_000) are spaced well above each index's range
(6, 2, 5 respectively) so no two cells collide.
"""
nonlinear_seed(i_curv::Int, i_p::Int, rep::Int; base::Int=42) =
    base + 10_000 * rep + 100 * i_p + 10 * i_curv

total_cells = length(curvature_grid) * length(p_grid) * 5
@info "Step 2: sweeping $total_cells cells × {SA, MVN} …"

rows = NamedTuple[]
cell_i = 0
t_start = time()

for (i_p, p) in enumerate(p_grid), (i_curv, curvature) in enumerate(curvature_grid)
    for rep in 1:5
        global cell_i += 1
        seed = nonlinear_seed(i_curv, i_p, rep)
        rng_cell = MersenneTwister(seed)

        sampler, Σpop, ctx = make_nonlinear_gt(p, rng_cell; kind=:scurve, curvature=curvature)
        Xtrain = sampler(n)

        elapsed = round(time() - t_start, digits=1)
        @info "[$cell_i/$total_cells] curvature=$curvature p=$p Rep=$rep seed=$seed  (elapsed $(elapsed)s)"

        Xsyn_sa = sa_recover(Xtrain; T=2000, seed=seed)
        Xsyn_mvn = mvn_recover(Xtrain; seed=seed)

        for (method, Xsyn) in (("SA", Xsyn_sa), ("MVN", Xsyn_mvn))
            off = offmanifold_error(Xsyn, ctx)
            cf = cov_frob(cov(Xsyn), Σpop)
            push!(rows, (kind="scurve", curvature=curvature, p=p, Rep=rep, Method=method,
                          OffManifold_Err=off, Cov_Frob=cf))
        end
    end
end

@info "Step 2 complete: $(length(rows)) rows in $(round(time()-t_start,digits=1))s"

results = DataFrame(rows)

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: write CSV + verify the intended contrast
# ══════════════════════════════════════════════════════════════════════════════
out_path = joinpath(_PATH_TO_DATA, "sim_nonlinear_results.csv")
CSV.write(out_path, results)
@info "Step 3: wrote $out_path ($(nrow(results)) rows)"

@assert nrow(results) == total_cells * 2 "expected $(total_cells*2) rows, got $(nrow(results))"
@assert all(isfinite, results.OffManifold_Err) && all(results.OffManifold_Err .>= 0)
@assert all(isfinite, results.Cov_Frob) && all(results.Cov_Frob .>= 0)
@info "  Invariant checks passed: all metrics finite and non-negative."

agg = combine(groupby(results, [:curvature, :p, :Method]),
              :OffManifold_Err => mean => :Mean_OffManifold_Err,
              :OffManifold_Err => std  => :SD_OffManifold_Err,
              :Cov_Frob => mean => :Mean_Cov_Frob,
              :Cov_Frob => std  => :SD_Cov_Frob)
sort!(agg, [:p, :Method, :curvature])

println()
println("═"^90)
println("Mean ± SD OffManifold_Err and Cov_Frob by (curvature, p, Method)")
println("═"^90)
pretty_table(agg)

# Locate the boundary curvature: first curvature at which SA's mean
# OffManifold_Err exceeds the curvature-0 baseline by ≥ 2× (report honestly —
# if no such crossing exists within the swept grid, say so, do not force one).
println()
println("═"^90)
println("Boundary curvature: first curvature where SA's mean OffManifold_Err ≥ 2× the curvature=0 baseline")
println("═"^90)
factor = 2.0
for p in p_grid
    sa_rows = sort(agg[(agg.p .== p) .& (agg.Method .== "SA"), :], :curvature)
    baseline = only(sa_rows[sa_rows.curvature .== 0.0, :Mean_OffManifold_Err])
    boundary = nothing
    for row in eachrow(sa_rows)
        if row.curvature > 0.0 && row.Mean_OffManifold_Err >= factor * baseline
            boundary = row.curvature
            break
        end
    end
    if boundary === nothing
        println("  p=$p: SA's OffManifold_Err never reached $(factor)× the curvature=0 baseline " *
                "($(round(baseline,digits=4))) within the swept grid $curvature_grid.")
    else
        println("  p=$p: boundary curvature = $boundary " *
                "(baseline=$(round(baseline,digits=4)), threshold=$(round(factor*baseline,digits=4)))")
    end

    # SA vs MVN at each curvature (both linear methods; report the comparison, not a verdict)
    mvn_rows = sort(agg[(agg.p .== p) .& (agg.Method .== "MVN"), :], :curvature)
    for (rs, rm) in zip(eachrow(sa_rows), eachrow(mvn_rows))
        verdict = rs.Mean_OffManifold_Err < rm.Mean_OffManifold_Err ? "SA lower" : "MVN lower/tied"
        println("    curvature=$(rs.curvature): SA=$(round(rs.Mean_OffManifold_Err,digits=4)) " *
                 "MVN=$(round(rm.Mean_OffManifold_Err,digits=4)) — $verdict")
    end
end

# Cov_Frob trend vs curvature (both methods are linear; both should degrade as the
# GT departs from Gaussian/linear structure — report as-is)
println()
println("═"^90)
println("Cov_Frob trend vs curvature, by p and Method")
println("═"^90)
for p in p_grid, method in ("SA", "MVN")
    sub = sort(agg[(agg.p .== p) .& (agg.Method .== method), :], :curvature)
    vals = round.(sub.Mean_Cov_Frob, digits=4)
    println("  p=$p $method: $(collect(zip(sub.curvature, vals)))")
end

@info "run_sim_nonlinear.jl complete."

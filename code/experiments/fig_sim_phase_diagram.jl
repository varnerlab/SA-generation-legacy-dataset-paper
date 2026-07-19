# ──────────────────────────────────────────────────────────────────────────────
# fig_sim_phase_diagram.jl
#
# Main-text figure for the SIM generalizability story (npj major revision).
# Plots (from already-computed CSVs — no new sampling is run here):
#   Panel A — (n,p) recovery phase diagram (arm 1, gaussian ground-truth cells,
#             averaged over r ∈ {3,5,10}): heatmap of the ratio
#             MVN_Cov_Frob / SA_Cov_Frob over the n × p grid. Ratio > 1 ⇒ SA's
#             covariance-recovery error is lower (SA better); ratio < 1 ⇒ MVN
#             better. The n = p boundary is overlaid — SA's advantage is
#             concentrated in the n < p regime.
#   Panel B — real-23 subsampling degradation curve (arm 2, SA-only):
#             mean ± SD of Cov_Frob_vs_full23 vs subsample size n, against the
#             full-cohort (n=23) published cross-visit covariance-error
#             reference (≈0.64, from validate_cross_visit_covariance.jl / E1).
#             Shows graceful degradation, not collapse, as n shrinks below 23.
#
# Reads:  data/sim_lowrank_results.csv (arm 1), data/sim_subsample_results.csv
#         (arm 2 — SA only; there is no arm-2 MVN).
# Writes: figs/sim_phase_diagram.pdf, figs/sim_phase_diagram.png
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
gr()

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load data
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Loading CSVs …"

df_lowrank = CSV.read(joinpath(_PATH_TO_DATA, "sim_lowrank_results.csv"), DataFrame)
df_sub     = CSV.read(joinpath(_PATH_TO_DATA, "sim_subsample_results.csv"), DataFrame)

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Panel A data — gaussian arm-1 rows, averaged over r, ratio grid
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 2: Building Panel A (n,p) ratio grid (gaussian tail, avg over r) …"

df_g = df_lowrank[df_lowrank.tail .== "gaussian", :]

n_vals = sort(unique(df_g.n))   # 8, 15, 23, 40, 80
p_vals = sort(unique(df_g.p))   # 30, 120, 216
@assert n_vals == [8, 15, 23, 40, 80]
@assert p_vals == [30, 120, 216]

sa_grid  = Matrix{Float64}(undef, length(p_vals), length(n_vals))
mvn_grid = Matrix{Float64}(undef, length(p_vals), length(n_vals))

for (pi, pv) in enumerate(p_vals), (ni, nv) in enumerate(n_vals)
    rows_sa  = df_g[(df_g.n .== nv) .& (df_g.p .== pv) .& (df_g.Method .== "SA"),  :Cov_Frob]
    rows_mvn = df_g[(df_g.n .== nv) .& (df_g.p .== pv) .& (df_g.Method .== "MVN"), :Cov_Frob]
    @assert length(rows_sa) == 3 && length(rows_mvn) == 3 "expected 3 r-replicates (r=3,5,10) per (n,p) cell"
    sa_grid[pi, ni]  = mean(rows_sa)
    mvn_grid[pi, ni] = mean(rows_mvn)
end

ratio_grid = mvn_grid ./ sa_grid   # > 1 ⇒ SA lower Cov_Frob ⇒ SA better

idx_n(v)  = findfirst(==(v), n_vals)
idx_p(v)  = findfirst(==(v), p_vals)

@info "  ratio[n=8,p=216]  = $(round(ratio_grid[idx_p(216), idx_n(8)],  digits=3))  (anchor ≈ 1.54)"
@info "  ratio[n=23,p=216] = $(round(ratio_grid[idx_p(216), idx_n(23)], digits=3))  (anchor ≈ 1.47)"
@info "  ratio[n=40,p=216] = $(round(ratio_grid[idx_p(216), idx_n(40)], digits=3))  (anchor ≈ 1.04)"
@info "  ratio[n=80,p=30]  = $(round(ratio_grid[idx_p(30),  idx_n(80)], digits=3))  (SA < MVN expected ⇒ ratio > 1)"

@assert isapprox(ratio_grid[idx_p(216), idx_n(8)],  1.54; atol=0.05)
@assert isapprox(ratio_grid[idx_p(216), idx_n(23)], 1.47; atol=0.05)
@assert isapprox(ratio_grid[idx_p(216), idx_n(40)], 1.04; atol=0.05)
@assert ratio_grid[idx_p(30), idx_n(80)] > 1.0 "n=80,p=30: SA expected to still beat MVN (ratio > 1)"

@info "  Full ratio grid (rows = p, cols = n):"
for pi in 1:length(p_vals)
    @info "    p=$(p_vals[pi]):  " * join(["n=$(n_vals[ni]): $(round(ratio_grid[pi,ni], digits=2))" for ni in 1:length(n_vals)], "   ")
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Panel B data — subsample degradation curve (arm 2, SA-only)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 3: Building Panel B (subsample degradation curve, SA-only) …"

n_sub = sort(unique(df_sub.n))   # 8, 10, 15, 20 (ascending)
@assert n_sub == [8, 10, 15, 20]

cov_mean = Float64[]
cov_sd   = Float64[]
mre_mean = Float64[]
mre_sd   = Float64[]
for nv in n_sub
    rows = df_sub[df_sub.n .== nv, :]
    push!(cov_mean, mean(rows.Cov_Frob_vs_full23))
    push!(cov_sd,   std(rows.Cov_Frob_vs_full23))
    push!(mre_mean, mean(rows.Marginal_MRE))
    push!(mre_sd,   std(rows.Marginal_MRE))
end

@info "  Cov_Frob_vs_full23 mean±SD by n: " *
      join(["n=$(n_sub[i]): $(round(cov_mean[i],digits=3))±$(round(cov_sd[i],digits=3))" for i in eachindex(n_sub)], "   ")
@info "  (anchors: 20→0.633, 15→0.734, 10→0.976, 8→1.108; SD ≈ 0.03–0.06)"

@assert isapprox(cov_mean[findfirst(==(20), n_sub)], 0.633; atol=0.02)
@assert isapprox(cov_mean[findfirst(==(15), n_sub)], 0.734; atol=0.02)
@assert isapprox(cov_mean[findfirst(==(10), n_sub)], 0.976; atol=0.02)
@assert isapprox(cov_mean[findfirst(==(8),  n_sub)], 1.108; atol=0.02)
@assert all(0.02 .< cov_sd .< 0.07) "SDs expected in the 0.03-0.06 ballpark per the brief"

# Published full-cohort (n=23) reference: SA's cross-visit covariance error on
# the real dataset at full sample size, from validate_cross_visit_covariance.jl
# / E1 (task-2), CrossVisit_Frob ≈ 0.6407-0.6410 (rounded to 0.64 in the
# brief). This is an EXTERNAL published number, not a row of sim_subsample_
# results.csv (that CSV measures distance-from-full-cohort, which is trivially
# 0 at n=23) — it anchors where the subsampling curve is heading as n → 23.
full_cohort_ref = 0.64

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Shared style (author convention — gray panel bg, no grid, boxed
# frame, legends on every panel, bare "A"/"B" labels, no overall title)
# ══════════════════════════════════════════════════════════════════════════════
bg_color  = RGB(0.96, 0.96, 0.96)
sa_color  = RGB(0.20, 0.50, 0.72)   # blue
mvn_color = RGB(0.85, 0.45, 0.25)   # orange

common_attrs = (background_color=bg_color, background_color_subplot=bg_color,
                foreground_color=:black,
                grid=false, framestyle=:box,
                guidefontsize=10, tickfontsize=9, legendfontsize=7,
                margin=7Plots.mm)

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Panel A — (n,p) recovery phase diagram
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 5: Building Panel A ((n,p) phase-diagram heatmap) …"

nx, ny = length(n_vals), length(p_vals)

# diverging colour scale centred exactly at ratio = 1.0
lo, hi = extrema(ratio_grid)
dev = max(1.0 - lo, hi - 1.0) * 1.05
clim_lo, clim_hi = 1.0 - dev, 1.0 + dev

pA = heatmap(1:nx, 1:ny, ratio_grid;
    xticks=(1:nx, string.(n_vals)), yticks=(1:ny, string.(p_vals)),
    xlabel="n (sample size)", ylabel="p (dimensionality)",
    clims=(clim_lo, clim_hi), color=cgrad(:RdBu, rev=false),
    colorbar_title="MVN / SA covariance error (Cov_Frob ratio)",
    title="A", titlelocation=:left, titlefontsize=12,
    common_attrs...)

# cell text labels (ratio value), white text on the darkest (most saturated) cells
for ni in 1:nx, pi in 1:ny
    val = ratio_grid[pi, ni]
    dist = abs(val - 1.0) / dev
    txt_color = dist > 0.55 ? :white : :black
    annotate!(pA, [(ni, pi, text(string(round(val, digits=2)), 8, txt_color, :center))])
end

# n = p boundary (stepped line): within the p=30 row it falls between n=23
# (idx 3) and n=40 (idx 4); for p=120 and p=216 the boundary (n=120, n=216) is
# beyond the largest simulated n=80, so the entire row is n<p and the
# staircase steps out to the right edge of the grid.
x_p30_boundary = (idx_n(23) + idx_n(40)) / 2   # 3.5
x_edge = nx + 0.5                               # 5.5
y_row1_top = idx_p(30) + 0.5                    # 1.5

plot!(pA, [x_p30_boundary, x_p30_boundary, x_edge, x_edge],
          [0.5, y_row1_top, y_row1_top, ny + 0.5];
      color=:black, lw=2.5, ls=:solid, label="n = p boundary", legend=:bottomleft)

annotate!(pA, [(2.0, 2.5, text("n < p\n(SA better)", 9, :black, :center))])

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Panel B — real-23 subsampling degradation curve (SA only)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 6: Building Panel B (subsampling degradation curve) …"

pB = plot(n_sub, cov_mean;
    xlabel="n (subsample size, real cohort)", ylabel="Cov_Frob vs. full n=23 cohort",
    xlims=(6, 25), xticks=(vcat(n_sub, 23), string.(vcat(n_sub, 23))),
    label="SA (mean ± SD)", color=sa_color, lw=2.5,
    marker=:circle, markersize=7, markerstrokecolor=:white, markerstrokewidth=0.8,
    markercolor=sa_color,
    linecolor=sa_color,
    legend=:topright,
    title="B", titlelocation=:left, titlefontsize=12,
    common_attrs...)

# Black error bars drawn as explicit segments (GR ignores `ecolor`, so build them by hand
# to honor the house convention of visible black error bars).
let cap = 0.45
    for (xi, yi, si) in zip(n_sub, cov_mean, cov_sd)
        plot!(pB, [xi, xi], [yi - si, yi + si]; color=:black, lw=1.5, label=nothing)
        plot!(pB, [xi - cap, xi + cap], [yi - si, yi - si]; color=:black, lw=1.5, label=nothing)
        plot!(pB, [xi - cap, xi + cap], [yi + si, yi + si]; color=:black, lw=1.5, label=nothing)
    end
end
# Redraw the SA markers on top of the bars so their centers stay clean.
scatter!(pB, n_sub, cov_mean; marker=:circle, markersize=7, markercolor=sa_color,
    markerstrokecolor=:white, markerstrokewidth=0.8, label=nothing)

hline!(pB, [full_cohort_ref]; color=:black, ls=:dash, lw=1.5, label=nothing)
scatter!(pB, [23], [full_cohort_ref];
    marker=:star5, markersize=12, color=:black, markerstrokecolor=:black,
    label="n=23 full cohort (published ≈0.64)")

annotate!(pB, [(23.3, full_cohort_ref - 0.06, text("full cohort", 7, :black, :left))])
annotate!(pB, [(6.3, 0.615, text("graceful degradation\nas n ↓ below 23", 8, :black, :left))])

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Compose and save (no overall title — author convention)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 7: Composing and saving figure …"

p_all = plot(pA, pB; layout=(1, 2), size=(1200, 500), dpi=200,
    background_color=:white, left_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=4Plots.mm)

savefig(p_all, joinpath(_PATH_TO_FIG, "sim_phase_diagram.pdf"))
savefig(p_all, joinpath(_PATH_TO_FIG, "sim_phase_diagram.png"))
@info "  Saved figs/sim_phase_diagram.pdf and figs/sim_phase_diagram.png"

# ══════════════════════════════════════════════════════════════════════════════
# Final report
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^78)
println("SIM Phase-Diagram Figure — Summary")
println("="^78)
println("  Panel A (gaussian, avg over r=3,5,10), ratio = MVN_Cov_Frob / SA_Cov_Frob:")
for pi in 1:ny
    println("    p=$(p_vals[pi]):  " * join(["n=$(n_vals[ni])→$(round(ratio_grid[pi,ni], digits=2))" for ni in 1:nx], "  "))
end
println("  Panel B: Cov_Frob_vs_full23 mean±SD by n = " *
        join(["$(n_sub[i])→$(round(cov_mean[i],digits=3))±$(round(cov_sd[i],digits=3))" for i in eachindex(n_sub)], ", "))
println("  Full-cohort (n=23) published reference: Cov_Frob ≈ $(full_cohort_ref)")
println("="^78)

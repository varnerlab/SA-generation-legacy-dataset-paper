# ──────────────────────────────────────────────────────────────────────────────
# verify_mechanistic_table.jl
#
# Regenerates every number in Supplementary Table S11 (mechanistic validation
# diagnostics) straight from data/validate_mechanistic_results.csv, using the
# definition the table's caption states and that regen_mechanistic_figures.jl
# implements (lines 104-105):
#
#   ratio   = Predicted / Measured, clamped to [0, 5]
#   overlap = fraction of SYNTHETIC ratios inside the [5th, 95th] percentile
#             band of the REAL ratios, per condition x feature
#   KS      = ApproximateTwoSampleKSTest(real_ratios, synth_ratios)
#
# Why this exists: the published cloud-overlap column drifted out of sync with
# the data while the KS D column did not (KS D reproduced exactly for all ten
# condition x feature pairs, so the ratios themselves never changed). The
# overlap values had been produced by a different code path and then carried
# forward by hand. Run this before submission and paste the printed block, or
# diff it against the table, rather than editing cells individually.
#
# Reads:  data/validate_mechanistic_results.csv
# Writes: nothing (prints a LaTeX-ready block)
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))
using Printf, HypothesisTests, Statistics

df = CSV.read(joinpath(_PATH_TO_DATA, "validate_mechanistic_results.csv"), DataFrame)
df.Ratio = clamp.(df.Predicted ./ df.Measured, 0.0, 5.0)

conditions = unique(df.Condition)
features   = unique(df.Feature)
short = Dict("Lagtime (min)" => "Lagtime", "Peak (nM)" => "Peak",
             "T.Peak (min)" => "T.Peak", "Max Rate (nM/min)" => "Max Rate",
             "ETP (nM·min)" => "ETP")

@info "Supplementary Table S11, regenerated from the canonical CSV:\n"
println("Condition & TGA Feature & Cloud overlap & KS \$D\$ & KS \$p\$ \\\\")
println("\\midrule")

all_real = Float64[]; all_syn = Float64[]
overlaps = Dict{String,Vector{Float64}}()
Ds = Float64[]; ps = Float64[]

for (ci, cond) in enumerate(conditions)
    ov_c = Float64[]
    for feat in features
        sub = df[(df.Condition .== cond) .& (df.Feature .== feat), :]
        r = sub[sub.Source .== "Real", :Ratio]
        s = sub[sub.Source .!= "Real", :Ratio]
        @assert length(r) > 2 && length(s) > 2 "too few records for $cond / $feat"
        lo, hi = quantile(r, 0.05), quantile(r, 0.95)
        ov = mean((s .>= lo) .& (s .<= hi))
        t  = ApproximateTwoSampleKSTest(r, s)
        push!(ov_c, ov); push!(Ds, t.δ); push!(ps, pvalue(t))
        append!(all_real, r); append!(all_syn, s)
        @printf("%-7s & %-16s & \\rone{%.2f} & %.3f & \\rone{%.2f} \\\\\n",
                cond, get(short, feat, feat), ov, t.δ, pvalue(t))
    end
    overlaps[cond] = ov_c
    ci < length(conditions) && println("\\midrule")
end

lo, hi = quantile(all_real, 0.05), quantile(all_real, 0.95)
pooled = mean((all_syn .>= lo) .& (all_syn .<= hi))

println("\\bottomrule")
@info "\nFootnote values:"
for cond in conditions
    @printf("  %s cloud overlap range: %.2f--%.2f\n", cond,
            minimum(overlaps[cond]), maximum(overlaps[cond]))
end
@printf("  Pooled overlap (all conditions and features): %.3f\n", pooled)
@printf("  KS D range: %.3f--%.3f, all p > %.2f\n", minimum(Ds), maximum(Ds), floor(minimum(ps)*100)/100)

# Guard: the pooled overlap is quoted in the main text and in Table S8's SA row.
@assert isapprox(pooled, 0.902; atol = 0.002) "pooled overlap moved to $pooled; update results.tex, Table S8 and Table S11 together"
@info "\n✓ pooled overlap 0.902 still matches the value quoted in the main text and Table S8"

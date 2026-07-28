# Standalone script: condition-specific feature comparison — grouped bar chart
# with a descriptive bootstrap Mann-Whitney non-detection diagnostic
# (subsample synth to match n_real). This is not an equivalence test.
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots, StatsPlots, Printf, Random, StatsBase
using HypothesisTests

const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ── Bootstrap non-detection diagnostic ───────────────────────────────────────
# Subsample synth to n_real, do Mann-Whitney, repeat B times.
# Returns the fraction of replicates where p > α. A high fraction means this
# test did not detect a difference at the available sample size; it does not
# establish equivalence between the real and synthetic distributions.
function bootstrap_mw(x_real, x_synth; B=1000, α=0.05, seed=42)
    rng = MersenneTwister(seed)
    n_real = length(x_real)
    n_synth = length(x_synth)
    n_sub = min(n_real, n_synth)  # subsample to the smaller group size

    n_nonsig = 0
    for _ in 1:B
        idx = sample(rng, 1:n_synth, n_sub; replace=false)
        x_sub = x_synth[idx]
        try
            mw = MannWhitneyUTest(x_real, x_sub)
            p = pvalue(mw)
            if p > α
                n_nonsig += 1
            end
        catch
            # test fails (e.g., all identical) → count as non-significant
            n_nonsig += 1
        end
    end
    return n_nonsig / B
end

# ── Load data ────────────────────────────────────────────────────────────────
# Consume the canonical conditional cohorts produced by run_conditioned_generation.jl
# (α=0.01, seed 42, N=100). The legacy `conditioned_*.csv` files written by
# validate_conditioned_generation.jl use α=0.1 and should not be used for the
# paper figures.
@info "Loading data …"
df_real = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth_u  = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_uncomplicated_cohort.csv"), DataFrame)
df_synth_p  = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_pcos_cohort.csv"), DataFrame)
df_synth_pe = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_developed_pe_cohort.csv"), DataFrame)

df_real_u  = filter(r -> r.DevelopedPE != "Yes", df_real)
df_real_p  = filter(r -> r.PCOS == "Yes", df_real)
df_real_pe = filter(r -> r.DevelopedPE == "Yes", df_real)

@info "Real: Uncomplicated=$(nrow(df_real_u)), PCOS=$(nrow(df_real_p)), Developed PE=$(nrow(df_real_pe))"

# ── Features ─────────────────────────────────────────────────────────────────
features = [:VIII, :vWF, :IX, :V, Symbol("α2-AP"), :VII, :Fbgn, :Plgn]
feat_short = ["FVIII", "vWF", "FIX", "FV", "α₂-AP", "FVII", "Fbgn", "Plgn"]
nf = length(features)

cond_labels = ["Uncomplicated (n=18)", "PCOS (n=3)", "Developed PE (n=5)"]
real_dfs = [df_real_u, df_real_p, df_real_pe]
synth_dfs = [df_synth_u, df_synth_p, df_synth_pe]
cond_real_colors = [RGB(0.20, 0.40, 0.65), RGB(0.85, 0.55, 0.10), RGB(0.75, 0.20, 0.20)]
cond_synth_colors = [RGB(0.50, 0.70, 0.90), RGB(1.0, 0.80, 0.50), RGB(1.0, 0.55, 0.55)]

bg_color = RGB(0.97, 0.97, 0.98)

function draw_errbar!(p, x, y, sd; cap_w=0.10)
    plot!(p, [x, x], [y - sd, y + sd]; lc=:black, lw=1.5, label="")
    plot!(p, [x - cap_w, x + cap_w], [y + sd, y + sd]; lc=:black, lw=1.5, label="")
    plot!(p, [x - cap_w, x + cap_w], [y - sd, y - sd]; lc=:black, lw=1.5, label="")
end

# ── Build figure ─────────────────────────────────────────────────────────────
panels = []
all_results = DataFrame(Condition=String[], Feature=String[], MRE=Float64[],
                        Frac_NoDifferenceDetected=Float64[])

@info "Running bootstrap MW tests (1000 replicates each) …"
for (ci, (cond_label, df_r, df_s)) in enumerate(zip(cond_labels, real_dfs, synth_dfs))
    real_vals = Float64[]
    synth_vals = Float64[]
    real_sd = Float64[]
    synth_sd = Float64[]
    mre_vals = Float64[]
    frac_nonsig = Float64[]

    for (fi, feat) in enumerate(features)
        x_r = Float64.(collect(skipmissing(df_r[!, feat])))
        x_s = Float64.(collect(skipmissing(df_s[!, feat])))

        mr, ms = mean(x_r), mean(x_s)
        sr, ss = std(x_r), std(x_s)
        mre = abs(mr - ms) / max(abs(mr), 1e-12)

        frac = bootstrap_mw(x_r, x_s; B=1000)

        push!(real_vals, mr)
        push!(synth_vals, ms)
        push!(real_sd, sr)
        push!(synth_sd, ss)
        push!(mre_vals, mre)
        push!(frac_nonsig, frac)

        push!(all_results, (cond_label, feat_short[fi], round(mre, digits=3),
                            round(frac, digits=3)))

        @info "  $cond_label $(feat_short[fi]): MRE=$(round(mre*100, digits=1))%, p>0.05 in $(round(frac*100, digits=1))% of replicates"
    end

    x = collect(1:nf)
    w = 0.35

    all_hi = maximum(vcat(real_vals .+ real_sd, synth_vals .+ synth_sd))
    all_lo = minimum(vcat(real_vals .- real_sd, synth_vals .- synth_sd))
    span = all_hi - all_lo
    yhi = all_hi + 0.4 * span
    ylo = max(0, all_lo - 0.15 * span)

    p = bar(x .- w/2, real_vals;
        bar_width=w, color=cond_real_colors[ci], label="Real",
        linecolor=cond_real_colors[ci], linewidth=0.5,
        title=cond_label, titlefontsize=16,
        xticks=(x, feat_short), tickfontsize=12,
        ylabel=ci == 1 ? "Mean (pooled)" : "",
        guidefontsize=15,
        background_color_inside=bg_color,
        grid=false, framestyle=:box,
        legend=:topright,
        legendfontsize=11, margin=5Plots.mm,
        ylim=(ylo, yhi))

    bar!(p, x .+ w/2, synth_vals;
        bar_width=w, color=cond_synth_colors[ci], label="SA synth",
        linecolor=cond_synth_colors[ci], linewidth=0.5)

    for k in 1:nf
        draw_errbar!(p, x[k] - w/2, real_vals[k], real_sd[k]; cap_w=0.10)
        draw_errbar!(p, x[k] + w/2, synth_vals[k], synth_sd[k]; cap_w=0.10)
    end

    # Annotations: MRE + bootstrap result
    for k in 1:nf
        top_val = max(real_vals[k] + real_sd[k], synth_vals[k] + synth_sd[k])

        # MRE
        mre_str = @sprintf("%.1f%%", mre_vals[k] * 100)
        annotate!(p, x[k], top_val + 0.06 * span,
                  text(mre_str, 7, :center, :gray30))

        # Bootstrap MW result
        frac_pct = round(frac_nonsig[k] * 100, digits=0)
        frac_str = @sprintf("p>0.05: %.0f%%", frac_pct)
        txt_color = frac_nonsig[k] >= 0.90 ? :forestgreen : :red
        annotate!(p, x[k], top_val + 0.16 * span,
                  text(frac_str, 6, :center, txt_color))
    end

    push!(panels, p)
end

p_cond = plot(panels...; layout=(1, 3), size=(1800, 600),
    left_margin=10Plots.mm, margin=6Plots.mm)
savefig(p_cond, joinpath(_PATH_TO_FIG, "validate_conditioned_features_v2.pdf"))
savefig(p_cond, joinpath(_PATH_TO_FIG, "validate_conditioned_features_v2.png"))
@info "  Saved validate_conditioned_features_v2.pdf"

# Summary
n_pass = sum(all_results.Frac_NoDifferenceDetected .>= 0.90)
n_total = nrow(all_results)
@info "\n═══ Summary ═══"
@info "  $n_pass / $n_total ($(round(100*n_pass/n_total, digits=1))%) features: p>0.05 in ≥90% of bootstrap replicates"
@info "  Median no-difference-detected fraction: $(round(median(all_results.Frac_NoDifferenceDetected)*100, digits=1))%"

CSV.write(joinpath(_PATH_TO_DATA, "bootstrap_mw_nondetection.csv"), all_results)
@info "  Saved bootstrap_mw_nondetection.csv"

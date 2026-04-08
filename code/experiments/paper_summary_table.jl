# ──────────────────────────────────────────────────────────────────────────────
# paper_summary_table.jl
#
# Per-feature, per-visit summary statistics: real vs synthetic.
# This script does NOT generate synthetic data. It consumes the canonical
# artifacts produced by run_full_longitudinal.jl:
#   - data/synthetic_full_longitudinal.csv  (100 synthetic patients × 3 visits)
#   - data/full_longitudinal_memory.jld2    (kept_cols, complete_subjects, df_clean)
#
# Real cohort = the K=23 complete-case patients used to train SA (filtered by
# complete_subjects). Synthetic cohort = the 100 patients that
# run_full_longitudinal.jl actually wrote out.
#
# Outputs:
#   - data/paper_summary_statistics.csv (per-visit, per-feature MRE)
#   - prints pooled and per-visit aggregates to stdout
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load canonical artifacts
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading canonical SA artifacts …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") complete_subjects kept_cols n_assays df_clean

@info "  Complete-case subjects: $(length(complete_subjects))"
@info "  Assay columns: $(length(kept_cols))"

# Filter df_clean to only the K=23 complete-case patients used to train SA
df_real = df_clean[in.(df_clean.SubjectID, Ref(complete_subjects)), :]
@info "  Real (complete-case) records: $(nrow(df_real))"

# Load canonical synthetic patients
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
@info "  Synthetic records: $(nrow(df_synth))"

# ══════════════════════════════════════════════════════════════════════════════
# Build per-(visit, feature) MRE table
# ══════════════════════════════════════════════════════════════════════════════
@info "\nBuilding per-feature MRE table …"

rows = []
for v in 1:3
    real_v = df_real[df_real.Visit .== v, :]
    synth_v = df_synth[df_synth.Visit .== v, :]
    for col in kept_cols
        real_vals = collect(skipmissing(real_v[!, col]))
        synth_vals = collect(skipmissing(synth_v[!, col]))

        r_mean = mean(real_vals)
        r_std  = std(real_vals)
        s_mean = mean(synth_vals)
        s_std  = std(synth_vals)

        mean_rel_err = abs(r_mean) > 1e-12 ? abs(s_mean - r_mean) / abs(r_mean) : NaN
        std_rel_err  = abs(r_std)  > 1e-12 ? abs(s_std  - r_std)  / abs(r_std)  : NaN

        push!(rows, (
            Visit = v,
            Feature = string(col),
            Real_Mean = round(r_mean, digits=2),
            Real_Std  = round(r_std,  digits=2),
            Synth_Mean = round(s_mean, digits=2),
            Synth_Std  = round(s_std,  digits=2),
            Mean_Rel_Err = round(mean_rel_err, digits=4),
            Std_Rel_Err  = round(std_rel_err,  digits=4),
        ))
    end
end

df_table = DataFrame(rows)
CSV.write(joinpath(_PATH_TO_DATA, "paper_summary_statistics.csv"), df_table)

# ══════════════════════════════════════════════════════════════════════════════
# Per-visit aggregates
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^70)
println("Per-Visit Aggregates")
println("="^70)
for v in 1:3
    vt = df_table[df_table.Visit .== v, :]
    avg_mean_err = mean(filter(!isnan, vt.Mean_Rel_Err))
    avg_std_err  = mean(filter(!isnan, vt.Std_Rel_Err))
    med_mean_err = median(filter(!isnan, vt.Mean_Rel_Err))
    med_std_err  = median(filter(!isnan, vt.Std_Rel_Err))
    println("  Visit $v:  Mean RE = $(round(avg_mean_err*100, digits=2))% (median $(round(med_mean_err*100, digits=2))%)  |  Std RE = $(round(avg_std_err*100, digits=2))% (median $(round(med_std_err*100, digits=2))%)")
end

# ══════════════════════════════════════════════════════════════════════════════
# Pooled (all 3 visits) aggregates with bootstrap CI
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^70)
println("Pooled (All 3 Visits, $(nrow(df_table)) entries)")
println("="^70)

all_mre = filter(!isnan, df_table.Mean_Rel_Err)
n = length(all_mre)

med_pooled  = median(all_mre)
mean_pooled = mean(all_mre)
q25 = quantile(all_mre, 0.25)
q75 = quantile(all_mre, 0.75)
mx  = maximum(all_mre)
below_5 = count(x -> x < 0.05, all_mre)

# Bootstrap 95% CI for pooled median (10,000 reps, fixed seed for reproducibility)
Random.seed!(42)
bs_medians = [median(rand(all_mre, n)) for _ in 1:10_000]
sort!(bs_medians)
ci_lo = bs_medians[250]
ci_hi = bs_medians[9750]

println("  Pooled median MRE  : $(round(med_pooled*100, digits=2))%  (95% bootstrap CI: $(round(ci_lo*100, digits=2))%--$(round(ci_hi*100, digits=2))%)")
println("  Pooled mean MRE    : $(round(mean_pooled*100, digits=2))%")
println("  Pooled 25th pctl   : $(round(q25*100, digits=2))%")
println("  Pooled 75th pctl   : $(round(q75*100, digits=2))%")
println("  Pooled max         : $(round(mx*100, digits=2))%")
println("  Features below 5%  : $below_5 / $n  ($(round(100*below_5/n, digits=1))%)")

@info "\n✓ Saved to data/paper_summary_statistics.csv"

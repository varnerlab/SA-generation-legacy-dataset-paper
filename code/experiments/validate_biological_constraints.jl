# ======================================================================================== #
# validate_biological_constraints.jl — Level 1: Marginal Plausibility
#
# Check whether known physiological relationships hold in synthetic patients:
#   1. AT inversely correlates with thrombin generation (peak, ETP)
#   2. Factor VIII positively correlates with peak/ETP
#   3. Pregnancy progression: fibrinogen, FVIII, vWF increase V1→V2→V3
#   4. Per-patient longitudinal monotonicity matches real data patterns
#   5. No biologically impossible values (negatives, extreme outliers)
# ======================================================================================== #

include(joinpath(@__DIR__, "..", "Include.jl"))

# ── Load data ──────────────────────────────────────────────────────────────────
@info "Loading data …"
df_real  = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

# tag sources
df_real[!, :Source] .= "Real"
df_synth[!, :Source] .= "Synthetic"

# align column names: synthetic has SyntheticID, real has SubjectID
df_synth[!, :SubjectID] = df_synth[!, :SyntheticID]

# ══════════════════════════════════════════════════════════════════════════════
# Check 1: Known correlations (AT vs peak, VIII vs peak)
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Check 1: Known Physiological Correlations ═══"

# Define expected correlations: (factor, TGA feature, expected sign, rationale)
expected_correlations = [
    (:AT, Symbol("TF Initiator Peak (nM)"),     -1, "AT inhibits thrombin → higher AT = lower peak"),
    (:AT, Symbol("TF Initiator ETP (nM·min)"),   -1, "AT inhibits thrombin → higher AT = lower ETP"),
    (:VIII, Symbol("TF Initiator Peak (nM)"),    +1, "FVIII amplifies intrinsic tenase → higher VIII = higher peak"),
    (:VIII, Symbol("TF Initiator ETP (nM·min)"), +1, "FVIII amplifies intrinsic tenase → higher VIII = higher ETP"),
    (:II, Symbol("TF Initiator Peak (nM)"),      +1, "More prothrombin substrate → higher peak"),
    (:II, Symbol("TF Initiator ETP (nM·min)"),   +1, "More prothrombin substrate → higher ETP"),
    (:X, Symbol("TF Initiator Peak (nM)"),       +1, "More FX substrate → more Xa → higher peak"),
]

corr_results = DataFrame(
    Factor       = String[],
    TGA_Feature  = String[],
    Expected     = String[],
    Real_ρ       = Float64[],
    Synth_ρ      = Float64[],
    Real_Correct = Bool[],
    Synth_Correct = Bool[],
)

for (factor, tga, expected_sign, rationale) in expected_correlations
    # compute per-visit correlations and average
    real_ρs = Float64[]
    synth_ρs = Float64[]

    for v in 1:3
        rv = df_real[df_real.Visit .== v, :]
        sv = df_synth[df_synth.Visit .== v, :]

        # filter missing/NaN
        r_x = collect(skipmissing(rv[:, factor]))
        r_y = collect(skipmissing(rv[:, tga]))
        s_x = collect(skipmissing(sv[:, factor]))
        s_y = collect(skipmissing(sv[:, tga]))

        n_r = min(length(r_x), length(r_y))
        n_s = min(length(s_x), length(s_y))

        if n_r >= 5
            push!(real_ρs, corspearman(r_x[1:n_r], r_y[1:n_r]))
        end
        if n_s >= 5
            push!(synth_ρs, corspearman(s_x[1:n_s], s_y[1:n_s]))
        end
    end

    rρ = isempty(real_ρs) ? NaN : mean(real_ρs)
    sρ = isempty(synth_ρs) ? NaN : mean(synth_ρs)
    sign_str = expected_sign > 0 ? "+" : "-"
    r_ok = !isnan(rρ) && sign(rρ) == expected_sign
    s_ok = !isnan(sρ) && sign(sρ) == expected_sign

    push!(corr_results, (string(factor), string(tga), sign_str, rρ, sρ, r_ok, s_ok))

    status_r = r_ok ? "✓" : "✗"
    status_s = s_ok ? "✓" : "✗"
    @info "  $factor vs $(split(string(tga), " ")[end-1]): expected=$sign_str | real=$(round(rρ, digits=3)) $status_r | synth=$(round(sρ, digits=3)) $status_s"
end

n_real_pass = sum(corr_results.Real_Correct)
n_synth_pass = sum(corr_results.Synth_Correct)
n_total = nrow(corr_results)
@info "  Score: Real $(n_real_pass)/$(n_total), Synthetic $(n_synth_pass)/$(n_total)"

# ══════════════════════════════════════════════════════════════════════════════
# Check 2: Pregnancy progression — known increasing factors
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Check 2: Pregnancy Progression Patterns ═══"

# Factors known to increase during pregnancy
progression_factors = [
    (:Fbgn,  "Fibrinogen increases through pregnancy"),
    (:VIII,  "FVIII increases through pregnancy"),
    (:vWF,   "vWF increases through pregnancy"),
    (:VII,   "FVII increases through pregnancy"),
    (Symbol("TF Initiator Peak (nM)"), "Thrombin peak increases (hypercoagulable state)"),
]

progression_results = DataFrame(
    Feature   = String[],
    Source    = String[],
    V1_mean  = Float64[],
    V2_mean  = Float64[],
    V3_mean  = Float64[],
    V1_V2_up = Bool[],
    V2_V3_up = Bool[],
    Monotone = Bool[],
)

for (col, description) in progression_factors
    for (src, df) in [("Real", df_real), ("Synthetic", df_synth)]
        means = Float64[]
        for v in 1:3
            vals = collect(skipmissing(df[df.Visit .== v, col]))
            push!(means, isempty(vals) ? NaN : mean(vals))
        end

        v12_up = means[2] > means[1]
        v23_up = means[3] > means[2]
        mono = v12_up && v23_up

        push!(progression_results, (string(col), src, means..., v12_up, v23_up, mono))
    end

    r_row = filter(r -> r.Source == "Real" && r.Feature == string(col), progression_results)
    s_row = filter(r -> r.Source == "Synthetic" && r.Feature == string(col), progression_results)

    r_trend = r_row.Monotone[1] ? "↑↑" : (r_row.V1_V2_up[1] ? "↑↓" : "↓?")
    s_trend = s_row.Monotone[1] ? "↑↑" : (s_row.V1_V2_up[1] ? "↑↓" : "↓?")

    @info "  $(rpad(string(col), 30)) Real: $(round(r_row.V1_mean[1], digits=1)) → $(round(r_row.V2_mean[1], digits=1)) → $(round(r_row.V3_mean[1], digits=1)) ($r_trend) | Synth: $(round(s_row.V1_mean[1], digits=1)) → $(round(s_row.V2_mean[1], digits=1)) → $(round(s_row.V3_mean[1], digits=1)) ($s_trend)"
end

# ══════════════════════════════════════════════════════════════════════════════
# Check 3: No biologically impossible values
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Check 3: Biologically Impossible Values ═══"

# Concentration-type features should be non-negative
positive_cols = [:II, :V, :VII, :VIII, :IX, :X, :AT, :PC, :Fbgn, :Estradiol,
                 Symbol("TF Initiator Peak (nM)"), Symbol("TF Initiator ETP (nM·min)"),
                 Symbol("TF Initiator Lagtime (min)")]

neg_counts = let
    ncr = 0; ncs = 0; details = []
    for col in positive_cols
        if col in propertynames(df_synth)
            s_vals = collect(skipmissing(df_synth[:, col]))
            n_neg = sum(s_vals .< 0)
            ncs += n_neg
            n_neg > 0 && push!(details, (col, n_neg, minimum(s_vals)))
        end
        if col in propertynames(df_real)
            r_vals = collect(skipmissing(df_real[:, col]))
            ncr += sum(r_vals .< 0)
        end
    end
    (real=ncr, synth=ncs, details=details)
end

@info "  Negative values in positive-definite features:"
@info "    Real: $(neg_counts.real)"
@info "    Synthetic: $(neg_counts.synth)"
if !isempty(neg_counts.details)
    for (col, n, minval) in neg_counts.details
        @info "      $col: $n negative values (min = $(round(minval, digits=2)))"
    end
end

# Check for extreme outliers (>5 SD from real mean)
@info "\n  Extreme outliers (>5 SD from real mean):"
outlier_features = [:II, :V, :VII, :VIII, :IX, :X, :AT, :PC, :Fbgn,
                    Symbol("TF Initiator Peak (nM)"), Symbol("TF Initiator ETP (nM·min)")]

n_outliers_total = let
    total = 0
    for col in outlier_features
        if col in propertynames(df_synth) && col in propertynames(df_real)
            r_vals = collect(skipmissing(df_real[:, col]))
            s_vals = collect(skipmissing(df_synth[:, col]))
            μ_r = mean(r_vals)
            σ_r = std(r_vals)
            if σ_r > 1e-12
                n_out = sum(abs.(s_vals .- μ_r) .> 5 * σ_r)
                total += n_out
                n_out > 0 && @info "    $col: $n_out / $(length(s_vals)) synthetic values beyond 5σ"
            end
        end
    end
    total
end
n_outliers_total == 0 && @info "    None found"

# ══════════════════════════════════════════════════════════════════════════════
# Check 4: Per-patient longitudinal coherence
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Check 4: Longitudinal Coherence ═══"

# For each patient, compute the V1→V3 change direction for key features.
# Compare the fraction of patients with increasing trend in real vs synthetic.

coherence_features = [:VIII, :Fbgn, :vWF, :II, :AT, :X]

@info "  Fraction of patients with V1→V3 increase:"
for col in coherence_features
    if !(col in propertynames(df_synth)) || !(col in propertynames(df_real))
        continue
    end

    # Real patients
    real_subjects = unique(df_real.SubjectID)
    r_inc = 0; r_total = 0
    for s in real_subjects
        v1_rows = df_real[(df_real.SubjectID .== s) .& (df_real.Visit .== 1), col]
        v3_rows = df_real[(df_real.SubjectID .== s) .& (df_real.Visit .== 3), col]
        if length(v1_rows) == 1 && length(v3_rows) == 1
            v1 = v1_rows[1]; v3 = v3_rows[1]
            if !ismissing(v1) && !ismissing(v3)
                r_total += 1
                v3 > v1 && (r_inc += 1)
            end
        end
    end

    # Synthetic patients
    synth_subjects = unique(df_synth.SubjectID)
    s_inc = 0; s_total = 0
    for s in synth_subjects
        v1_rows = df_synth[(df_synth.SubjectID .== s) .& (df_synth.Visit .== 1), col]
        v3_rows = df_synth[(df_synth.SubjectID .== s) .& (df_synth.Visit .== 3), col]
        if length(v1_rows) == 1 && length(v3_rows) == 1
            s_total += 1
            v3_rows[1] > v1_rows[1] && (s_inc += 1)
        end
    end

    r_frac = r_total > 0 ? round(r_inc / r_total, digits=2) : NaN
    s_frac = s_total > 0 ? round(s_inc / s_total, digits=2) : NaN
    @info "    $(rpad(string(col), 10)) Real: $r_frac ($r_inc/$r_total)  Synth: $s_frac ($s_inc/$s_total)"
end

# ══════════════════════════════════════════════════════════════════════════════
# Visualization: correlation scatters for key relationships
# ══════════════════════════════════════════════════════════════════════════════
@info "\nGenerating plots …"

plot_pairs = [
    (:AT, Symbol("TF Initiator Peak (nM)"), "AT (% nominal)", "TF Peak (nM)"),
    (:VIII, Symbol("TF Initiator Peak (nM)"), "FVIII (% nominal)", "TF Peak (nM)"),
    (:II, Symbol("TF Initiator ETP (nM·min)"), "FII (% nominal)", "TF ETP (nM·min)"),
    (:AT, Symbol("TF Initiator ETP (nM·min)"), "AT (% nominal)", "TF ETP (nM·min)"),
]

panels = []
for (xcol, ycol, xlab, ylab) in plot_pairs
    p = plot(; xlabel=xlab, ylabel=ylab, legend=:topright, size=(400, 400))

    for v in 1:3
        rv = df_real[df_real.Visit .== v, :]
        r_x = collect(skipmissing(rv[:, xcol]))
        r_y = collect(skipmissing(rv[:, ycol]))
        n = min(length(r_x), length(r_y))
        scatter!(p, r_x[1:n], r_y[1:n]; label=(v==1 ? "Real" : ""), mc=:steelblue,
                 ms=4, ma=0.6, shape=[:circle, :diamond, :utriangle][v])
    end

    for v in 1:3
        sv = df_synth[df_synth.Visit .== v, :]
        scatter!(p, sv[:, xcol], sv[:, ycol]; label=(v==1 ? "Synth" : ""), mc=:coral,
                 ms=3, ma=0.3, shape=[:circle, :diamond, :utriangle][v])
    end
    push!(panels, p)
end

p_corr = plot(panels...; layout=(2, 2), size=(900, 800),
              plot_title="Level 1: Physiological Correlations", margin=5Plots.mm)
savefig(p_corr, joinpath(_PATH_TO_FIG, "validate_bio_correlations.pdf"))
savefig(p_corr, joinpath(_PATH_TO_FIG, "validate_bio_correlations.png"))
@info "  Saved validate_bio_correlations.pdf"

# ── Pregnancy progression plot ─────────────────────────────────────────────────
prog_cols = [:Fbgn, :VIII, :vWF, :VII, :II, :AT]
prog_panels = []
for col in prog_cols
    if !(col in propertynames(df_synth))
        continue
    end
    r_means = [mean(collect(skipmissing(df_real[df_real.Visit .== v, col]))) for v in 1:3]
    r_sds   = [std(collect(skipmissing(df_real[df_real.Visit .== v, col]))) for v in 1:3]
    s_means = [mean(collect(skipmissing(df_synth[df_synth.Visit .== v, col]))) for v in 1:3]
    s_sds   = [std(collect(skipmissing(df_synth[df_synth.Visit .== v, col]))) for v in 1:3]

    p = plot([1,2,3], r_means; ribbon=r_sds, fillalpha=0.2, label="Real",
             lw=2, mc=:steelblue, marker=:circle, color=:steelblue,
             xlabel="Visit", ylabel="% nominal",
             title=string(col), xticks=([1,2,3], ["BL","1st","3rd"]))
    plot!(p, [1,2,3], s_means; ribbon=s_sds, fillalpha=0.2, label="Synth",
          lw=2, mc=:coral, marker=:diamond, color=:coral)
    push!(prog_panels, p)
end

p_prog = plot(prog_panels...; layout=(2, 3), size=(1000, 600),
              plot_title="Level 1: Pregnancy Progression", margin=5Plots.mm)
savefig(p_prog, joinpath(_PATH_TO_FIG, "validate_pregnancy_progression.pdf"))
savefig(p_prog, joinpath(_PATH_TO_FIG, "validate_pregnancy_progression.png"))
@info "  Saved validate_pregnancy_progression.pdf"

# ── Save results ───────────────────────────────────────────────────────────────
CSV.write(joinpath(_PATH_TO_DATA, "validate_bio_correlations.csv"), corr_results)
CSV.write(joinpath(_PATH_TO_DATA, "validate_pregnancy_progression.csv"), progression_results)
@info "\nDone! Level 1 validation complete."

# ======================================================================================== #
# validate_cross_visit_covariance.jl — Level 2: Joint Structure
#
# Compare the full cross-visit correlation structure in real vs SA vs MVN.
# The concatenated [V1|V2|V3] vector has d_concat = n_assays × 3 dimensions.
# Key test: does the off-diagonal block structure (V1↔V2, V1↔V3, V2↔V3)
# match between real and synthetic data?
#
# MVN estimated a rank-22 covariance from n=23 patients with regularization.
# SA preserves the data manifold geometry directly.
# ======================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames, Statistics, Plots, StatsBase, LinearAlgebra, Random
using Distributions

const _ROOT = joinpath(@__DIR__, "..")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG  = joinpath(_ROOT, "figs")

# ── Load data from CSVs ──────────────────────────────────────────────────────
@info "Loading data from CSVs …"
df_real  = CSV.read(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), DataFrame)
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

# Identify complete real patients (all 3 visits)
all_subjects = unique(df_real.SubjectID)
complete_subjects = [s for s in all_subjects
    if sum(df_real.SubjectID .== s) == 3 &&
       all(v -> any((df_real.SubjectID .== s) .& (df_real.Visit .== v)), 1:3)]

# Identify assay columns (shared between real and synth CSVs)
meta_cols = [:SubjectID, :Visit, :DrawTag, :Condition, :PCOS, :DevelopedPE, :DevelopedHTN]
synth_meta = [:SyntheticID, :Visit]
kept_cols = [c for c in Symbol.(names(df_synth)) if c ∉ synth_meta]
n_assays = length(kept_cols)
n_visits = 3
d_concat = n_assays * n_visits
K = length(complete_subjects)

@info "  $K complete subjects, $n_assays features/visit, $d_concat total dimensions"

# ── Build concatenated matrices ────────────────────────────────────────────────
@info "\nBuilding concatenated data matrices …"

# Real data
X_real = Matrix{Float64}(undef, K, d_concat)
for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_real.SubjectID .== s) .& (df_real.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_real[i, offset + j] = Float64(df_real[row, col])
        end
    end
end

# SA synthetic data
n_synth = length(unique(df_synth.SyntheticID))
X_sa = Matrix{Float64}(undef, n_synth, d_concat)
for i in 1:n_synth
    for v in 1:n_visits
        row_idx = findfirst((df_synth.SyntheticID .== i) .& (df_synth.Visit .== v))
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_sa[i, offset + j] = Float64(df_synth[row_idx, col])
        end
    end
end
@info "  SA: $n_synth patients"

# MVN synthetic data
@info "  Generating MVN synthetic data …"
μ_mvn = vec(mean(X_real, dims=1))
Σ_mvn = cov(X_real)
Σ_mvn = (Σ_mvn + Σ_mvn') / 2
Σ_mvn += 1e-6 * I

Random.seed!(42)
dist_mvn = MvNormal(μ_mvn, Σ_mvn)
X_mvn = rand(dist_mvn, n_synth)'   # n_synth × d_concat
@info "  MVN: $n_synth patients"

# ══════════════════════════════════════════════════════════════════════════════
# Compute correlation matrices
# ══════════════════════════════════════════════════════════════════════════════
@info "\nComputing correlation matrices …"

# Safe correlation: replace NaN (from constant columns) with 0
function safe_cor(X)
    C = cor(X)
    C[isnan.(C)] .= 0.0
    return C
end

C_real = safe_cor(X_real)
C_sa   = safe_cor(X_sa)
C_mvn  = safe_cor(X_mvn)

# ══════════════════════════════════════════════════════════════════════════════
# Compare: full matrix Frobenius distance
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Full Correlation Matrix Distance ═══"

frob_sa  = norm(C_sa - C_real) / norm(C_real)
frob_mvn = norm(C_mvn - C_real) / norm(C_real)
@info "  Relative Frobenius distance to real:"
@info "    SA:  $(round(frob_sa, digits=4))"
@info "    MVN: $(round(frob_mvn, digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Compare: per-block analysis (within-visit vs cross-visit)
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Block-wise Correlation Distance ═══"

function block_distance(C_test, C_ref, v1_range, v2_range)
    block_test = C_test[v1_range, v2_range]
    block_ref  = C_ref[v1_range, v2_range]
    return norm(block_test - block_ref) / max(norm(block_ref), 1e-12)
end

visit_ranges = [(1:n_assays), (n_assays+1:2*n_assays), (2*n_assays+1:3*n_assays)]
block_labels = ["V1-V1", "V2-V2", "V3-V3", "V1-V2", "V1-V3", "V2-V3"]
block_pairs = [(1,1), (2,2), (3,3), (1,2), (1,3), (2,3)]

@info "  $(rpad("Block", 8))  $(rpad("SA", 10))  MVN"
for (label, (i, j)) in zip(block_labels, block_pairs)
    d_sa  = block_distance(C_sa, C_real, visit_ranges[i], visit_ranges[j])
    d_mvn = block_distance(C_mvn, C_real, visit_ranges[i], visit_ranges[j])
    marker_sa  = d_sa < d_mvn ? " ✓" : ""
    marker_mvn = d_mvn < d_sa ? " ✓" : ""
    @info "  $(rpad(label, 8))  $(rpad(round(d_sa, digits=4), 10))$marker_sa  $(round(d_mvn, digits=4))$marker_mvn"
end

# ══════════════════════════════════════════════════════════════════════════════
# Cross-visit feature correlations: does V1 peak predict V3 peak?
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Cross-Visit Feature Correlations ═══"
@info "  (V1 feature vs V3 feature — Spearman ρ)"

cross_visit_features = [:II, :VIII, :AT, :Fbgn, :X,
    Symbol("TF Initiator Peak (nM)"),
    Symbol("TF Initiator ETP (nM·min)")]

@info "  $(rpad("Feature", 30))  $(rpad("Real", 8))  $(rpad("SA", 8))  MVN"
for col in cross_visit_features
    col_idx = findfirst(==(col), kept_cols)
    if isnothing(col_idx)
        continue
    end

    v1_idx = col_idx
    v3_idx = 2 * n_assays + col_idx

    ρ_real = corspearman(X_real[:, v1_idx], X_real[:, v3_idx])
    ρ_sa   = corspearman(X_sa[:, v1_idx], X_sa[:, v3_idx])
    ρ_mvn  = corspearman(X_mvn[:, v1_idx], X_mvn[:, v3_idx])

    # closer to real is better
    err_sa  = abs(ρ_sa - ρ_real)
    err_mvn = abs(ρ_mvn - ρ_real)
    marker = err_sa < err_mvn ? "SA ✓" : "MVN ✓"

    @info "  $(rpad(string(col), 30))  $(rpad(round(ρ_real, digits=3), 8))  $(rpad(round(ρ_sa, digits=3), 8))  $(round(ρ_mvn, digits=3))  $marker"
end

# ══════════════════════════════════════════════════════════════════════════════
# Eigenvalue spectrum comparison
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Covariance Eigenvalue Spectrum ═══"

# standardize before eigendecomposition (same scale)
function standardized_cov(X)
    μ = mean(X, dims=1)
    σ = std(X, dims=1)
    σ[σ .< 1e-12] .= 1.0
    Z = (X .- μ) ./ σ
    return cov(Z)
end

Σ_r = standardized_cov(X_real)
Σ_s = standardized_cov(X_sa)
Σ_m = standardized_cov(X_mvn)

eig_real = sort(eigvals(Σ_r), rev=true)
eig_sa   = sort(eigvals(Σ_s), rev=true)
eig_mvn  = sort(eigvals(Σ_m), rev=true)

# effective rank (number of eigenvalues > 1% of max)
thresh = 0.01 * eig_real[1]
rank_real = sum(eig_real .> thresh)
rank_sa   = sum(eig_sa .> thresh)
rank_mvn  = sum(eig_mvn .> thresh)

@info "  Effective rank (eigenvalues > 1% of max):"
@info "    Real: $rank_real"
@info "    SA:   $rank_sa"
@info "    MVN:  $rank_mvn"

# ══════════════════════════════════════════════════════════════════════════════
# Visualization
# ══════════════════════════════════════════════════════════════════════════════
@info "\nGenerating plots …"

# Correlation heatmaps
clim = (-1, 1)
hm_kw = (c=:RdBu, aspect_ratio=1, xticks=false, yticks=false,
          axis=false, border=:none, colorbar_titlefontsize=8, yflip=true)
# First panel (Real) gets visit labels on axes
mid1 = div(n_assays, 2)
mid2 = n_assays + mid1
mid3 = 2*n_assays + mid1
visit_ticks = ([mid1, mid2, mid3], ["V1", "V2", "V3"])

p1 = heatmap(C_real; title="Real", clim=clim, hm_kw...)
# Add V1/V2/V3 labels as annotations on the Real panel
d = size(C_real, 1)
for (label, mid) in [("V1", mid1), ("V2", mid2), ("V3", mid3)]
    # bottom edge (x-axis labels)
    annotate!(p1, mid, d + 8, text(label, 9, :center, :black))
    # left edge (y-axis labels)
    annotate!(p1, -8, mid, text(label, 9, :center, :black))
end
# Other panels without axis labels
p2 = heatmap(C_sa; title="SA", clim=clim, hm_kw...)
p3 = heatmap(C_mvn; title="MVN", clim=clim, hm_kw...)

# Difference heatmaps (all three pairwise) — same scale as top row
p4 = heatmap(C_sa - C_real; title="SA − Real", clim=clim, hm_kw...)
p5 = heatmap(C_mvn - C_real; title="MVN − Real", clim=clim, hm_kw...)
p6 = heatmap(C_sa - C_mvn; title="SA − MVN", clim=clim, hm_kw...)

# add block dividers
for p in [p1, p2, p3, p4, p5, p6]
    vline!(p, [n_assays + 0.5, 2*n_assays + 0.5]; lc=:black, lw=1, label="")
    hline!(p, [n_assays + 0.5, 2*n_assays + 0.5]; lc=:black, lw=1, label="")
end

p_heat = plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(1400, 900),
              plot_title="Level 2: Cross-Visit Correlation Structure",
              margin=3Plots.mm)
savefig(p_heat, joinpath(_PATH_TO_FIG, "validate_cross_visit_corr.pdf"))
savefig(p_heat, joinpath(_PATH_TO_FIG, "validate_cross_visit_corr.png"))
@info "  Saved validate_cross_visit_corr.pdf"

# Supplementary version: ±0.5 scale on residuals to reveal fine structure
dlim = (-0.5, 0.5)
frob_sa_real = round(norm(C_sa - C_real), digits=2)
frob_mvn_real = round(norm(C_mvn - C_real), digits=2)
frob_sa_mvn = round(norm(C_sa - C_mvn), digits=2)
p4s = heatmap(C_sa - C_real; title="SA − Real (‖Δ‖=$frob_sa_real)", clim=dlim, hm_kw...)
p5s = heatmap(C_mvn - C_real; title="MVN − Real (‖Δ‖=$frob_mvn_real)", clim=dlim, hm_kw...)
p6s = heatmap(C_sa - C_mvn; title="SA − MVN (‖Δ‖=$frob_sa_mvn)", clim=dlim, hm_kw...)
for p in [p4s, p5s, p6s]
    vline!(p, [n_assays + 0.5, 2*n_assays + 0.5]; lc=:black, lw=1, label="")
    hline!(p, [n_assays + 0.5, 2*n_assays + 0.5]; lc=:black, lw=1, label="")
end
p_heat_supp = plot(p1, p2, p3, p4s, p5s, p6s; layout=(2, 3), size=(1400, 900),
                   plot_title="Level 2: Cross-Visit Correlation Structure (±0.5 residual scale)",
                   margin=3Plots.mm)
savefig(p_heat_supp, joinpath(_PATH_TO_FIG, "validate_cross_visit_corr_supp.pdf"))
savefig(p_heat_supp, joinpath(_PATH_TO_FIG, "validate_cross_visit_corr_supp.png"))
@info "  Saved validate_cross_visit_corr_supp.pdf"

# Eigenvalue spectrum — full range to show the MVN hallucination gap
n_show = min(40, length(eig_real))
n_mvn_show = min(n_show, length(eig_mvn))

# Clamp only true zeros to a small floor for log display (keep real values)
eig_floor = 1e-14
eig_r_plot = max.(eig_real[1:n_show], eig_floor)
eig_s_plot = max.(eig_sa[1:n_show], eig_floor)
eig_m_plot = max.(eig_mvn[1:n_mvn_show], eig_floor)

p_eig = plot(1:n_show, eig_r_plot; label="Real (rank 22)", lw=2.5, mc=:steelblue,
             marker=:circle, ms=4, ma=0.8,
             xlabel="Component", ylabel="Eigenvalue",
             title="Covariance Eigenvalue Spectrum",
             yscale=:log10, legend=:right, size=(650, 450),
             ylim=(1e-14, maximum(eig_real) * 3),
             yticks=[1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0, 1e1, 1e2],
             guidefontsize=11, titlefontsize=12, tickfontsize=8, legendfontsize=9)
plot!(p_eig, 1:n_show, eig_s_plot; label="SA (rank 18)", lw=2.5, mc=:coral,
      marker=:diamond, ms=4, ma=0.8)
plot!(p_eig, 1:n_mvn_show, eig_m_plot; label="MVN (rank 23+)", lw=2.5, mc=:forestgreen,
      marker=:utriangle, ms=4, ma=0.8)

# Annotate the rank drop points
vline!(p_eig, [18.5]; lc=:coral, ls=:dash, lw=1.5, label="")
vline!(p_eig, [22.5]; lc=:steelblue, ls=:dash, lw=1.5, label="")
annotate!(p_eig, 18.5, maximum(eig_real) * 1.5, text("PCA\ntrunc.", 7, :center, :coral))
annotate!(p_eig, 22.5, maximum(eig_real) * 1.5, text("K−1\nrank", 7, :center, :steelblue))

# Shade the "spurious MVN variance" region
plot!(p_eig, Shape([23, n_show, n_show, 23],
      [1e-14, 1e-14, maximum(eig_real) * 3, maximum(eig_real) * 3]);
      fillalpha=0.08, fc=:forestgreen, lc=:transparent, label="")

# Label the gap
annotate!(p_eig, 32, 1e-1, text("MVN: spurious\nvariance", 8, :center, :forestgreen))

savefig(p_eig, joinpath(_PATH_TO_FIG, "validate_eigenvalue_spectrum.pdf"))
savefig(p_eig, joinpath(_PATH_TO_FIG, "validate_eigenvalue_spectrum.png"))
@info "  Saved validate_eigenvalue_spectrum.pdf"

@info "\nDone! Level 2 validation complete."

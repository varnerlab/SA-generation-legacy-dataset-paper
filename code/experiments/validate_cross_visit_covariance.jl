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

include(joinpath(@__DIR__, "..", "Include.jl"))

# ── Load pipeline state ───────────────────────────────────────────────────────
@info "Loading pipeline state …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms complete_subjects subject_conditions kept_cols n_assays df_clean

d_pca, K = size(X̂)
n_visits = 3
d_concat = n_assays * n_visits

@info "  $K complete subjects, $n_assays features/visit, $d_concat total dimensions"

# ── Build concatenated matrices ────────────────────────────────────────────────
@info "\nBuilding concatenated data matrices …"

# Real data
X_real = Matrix{Float64}(undef, K, d_concat)
for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_clean.SubjectID .== s) .& (df_clean.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_real[i, offset + j] = df_clean[row, col]
        end
    end
end

# SA synthetic data
df_synth = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)
n_synth = length(unique(df_synth.SyntheticID))
X_sa = Matrix{Float64}(undef, n_synth, d_concat)
for i in 1:n_synth
    for v in 1:n_visits
        row_idx = findfirst((df_synth.SyntheticID .== i) .& (df_synth.Visit .== v))
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_sa[i, offset + j] = df_synth[row_idx, col]
        end
    end
end
@info "  SA: $n_synth patients"

# MVN synthetic data (generate inline — same method as paper_sa_vs_mvn.jl)
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
p1 = heatmap(C_real; title="Real", c=:RdBu, clim=clim, aspect_ratio=1,
             xticks=false, yticks=false, size=(400, 400))
p2 = heatmap(C_sa; title="SA", c=:RdBu, clim=clim, aspect_ratio=1,
             xticks=false, yticks=false, size=(400, 400))
p3 = heatmap(C_mvn; title="MVN", c=:RdBu, clim=clim, aspect_ratio=1,
             xticks=false, yticks=false, size=(400, 400))

# Difference heatmaps
dlim = (-0.5, 0.5)
p4 = heatmap(C_sa - C_real; title="SA - Real", c=:RdBu, clim=dlim, aspect_ratio=1,
             xticks=false, yticks=false, size=(400, 400))
p5 = heatmap(C_mvn - C_real; title="MVN - Real", c=:RdBu, clim=dlim, aspect_ratio=1,
             xticks=false, yticks=false, size=(400, 400))

# add block dividers
for p in [p1, p2, p3, p4, p5]
    vline!(p, [n_assays + 0.5, 2*n_assays + 0.5]; lc=:black, lw=1, label="")
    hline!(p, [n_assays + 0.5, 2*n_assays + 0.5]; lc=:black, lw=1, label="")
end

p_heat = plot(p1, p2, p3, p4, p5; layout=(2, 3), size=(1200, 800),
              plot_title="Level 2: Cross-Visit Correlation Structure")
savefig(p_heat, joinpath(_PATH_TO_FIG, "validate_cross_visit_corr.pdf"))
savefig(p_heat, joinpath(_PATH_TO_FIG, "validate_cross_visit_corr.png"))
@info "  Saved validate_cross_visit_corr.pdf"

# Eigenvalue spectrum
n_show = min(40, length(eig_real))
p_eig = plot(1:n_show, eig_real[1:n_show]; label="Real", lw=2, mc=:steelblue, marker=:circle,
             xlabel="Component", ylabel="Eigenvalue", title="Covariance Eigenvalue Spectrum",
             yscale=:log10, legend=:topright)
plot!(p_eig, 1:n_show, eig_sa[1:n_show]; label="SA", lw=2, mc=:coral, marker=:diamond)
plot!(p_eig, 1:min(n_show, length(eig_mvn)), eig_mvn[1:min(n_show, length(eig_mvn))];
      label="MVN", lw=2, mc=:green, marker=:utriangle)
savefig(p_eig, joinpath(_PATH_TO_FIG, "validate_eigenvalue_spectrum.pdf"))
savefig(p_eig, joinpath(_PATH_TO_FIG, "validate_eigenvalue_spectrum.png"))
@info "  Saved validate_eigenvalue_spectrum.pdf"

@info "\nDone! Level 2 validation complete."

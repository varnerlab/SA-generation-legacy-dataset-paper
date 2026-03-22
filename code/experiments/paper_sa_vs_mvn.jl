# ──────────────────────────────────────────────────────────────────────────────
# paper_sa_vs_mvn.jl
#
# Head-to-head comparison of Stochastic Attention (SA) synthetic patient
# generation versus Multivariate Normal (MVN) generation on the longitudinal
# pregnancy coagulation dataset. MVN mirrors the JuliaCon baseline approach:
# independent per-visit multivariate normals with no cross-visit coupling.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Load saved pipeline artefacts
# ══════════════════════════════════════════════════════════════════════════════
@info "Loading data …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms n_assays df_clean kept_cols

d_pca, K_mem = size(X̂)
n_synth = 500
T_steps = 2000

# ══════════════════════════════════════════════════════════════════════════════
# Helper: two-sample Kolmogorov-Smirnov statistic (no external dependency)
# ══════════════════════════════════════════════════════════════════════════════
function ks_statistic(x, y)
    nx, ny = length(x), length(y)
    xs, ys = sort(x), sort(y)
    all_vals = sort(vcat(xs, ys))
    d_max = 0.0
    for v in all_vals
        f1 = searchsortedlast(xs, v) / nx
        f2 = searchsortedlast(ys, v) / ny
        d_max = max(d_max, abs(f1 - f2))
    end
    return d_max
end

# ══════════════════════════════════════════════════════════════════════════════
# Helper: pairwise-complete covariance (handles missing values)
# ══════════════════════════════════════════════════════════════════════════════
function pairwise_cov(df, cols)
    n = length(cols)
    C = Matrix{Float64}(undef, n, n)
    for i in 1:n
        for j in 1:n
            xi = Float64[]
            xj = Float64[]
            for row in 1:nrow(df)
                a = df[row, cols[i]]
                b = df[row, cols[j]]
                if !ismissing(a) && !ismissing(b)
                    push!(xi, a)
                    push!(xj, b)
                end
            end
            C[i,j] = cov(xi, xj)
        end
    end
    # ensure symmetric positive definite
    C = (C + C') / 2
    C += 1e-6 * I
    return C
end

# ══════════════════════════════════════════════════════════════════════════════
# Helper: pairwise-complete correlation matrix (for evaluation)
# ══════════════════════════════════════════════════════════════════════════════
function pairwise_cor(df, cols)
    n = length(cols)
    C = Matrix{Float64}(undef, n, n)
    for i in 1:n
        for j in 1:n
            if i == j
                C[i,j] = 1.0
            else
                xi = Float64[]
                xj = Float64[]
                for row in 1:nrow(df)
                    a = df[row, cols[i]]
                    b = df[row, cols[j]]
                    if !ismissing(a) && !ismissing(b)
                        push!(xi, a)
                        push!(xj, b)
                    end
                end
                C[i,j] = length(xi) > 2 ? cor(xi, xj) : NaN
            end
        end
    end
    return C
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Generate 500 SA synthetic patients
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Generating $n_synth SA synthetic patients …"

phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star

Random.seed!(123)
sa_concat = Matrix{Float64}(undef, n_synth, length(std_params_concat.col_names))
for i in 1:n_synth
    k0 = rand(1:K_mem)
    ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d_pca)
    ξ₀ ./= norm(ξ₀)
    result = sample(X̂, ξ₀, T_steps; β=β_star, α=0.05)
    ξ_dir = result.Ξ[end, :]
    ξ_dir ./= (norm(ξ_dir) + 1e-12)
    r = rand(pca_norms)
    ξ_rescaled = r .* ξ_dir
    z_std = vec(MultivariateStats.reconstruct(pca_concat, ξ_rescaled))
    sa_concat[i, :] = z_std .* std_params_concat.σ .+ std_params_concat.μ
end

# split SA into per-visit DataFrames
df_sa_visits = DataFrame[]
for v in 1:3
    offset = (v - 1) * n_assays
    cols_data = sa_concat[:, (offset+1):(offset+n_assays)]
    df_v = DataFrame(cols_data, kept_cols)
    df_v.Visit = fill(v, n_synth)
    push!(df_sa_visits, df_v)
end
df_sa = vcat(df_sa_visits...)
@info "  SA generation complete: $(nrow(df_sa)) records"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Generate 500 MVN synthetic patients (concatenated longitudinal)
#
# Fit a single MvNormal on the same K=23 concatenated [V1|V2|V3] patient
# vectors used by SA. This produces coherent longitudinal trajectories
# from the MVN, making it a fair apples-to-apples comparison.
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 2: Generating $n_synth MVN synthetic patients (concatenated) …"

# rebuild concatenated matrix from complete subjects (same as run_full_longitudinal.jl)
subjects = unique(df_clean.SubjectID)
complete_subjects = Int[]
for s in subjects
    visits = sort(unique(df_clean[df_clean.SubjectID .== s, :Visit]))
    if visits == [1, 2, 3]
        push!(complete_subjects, s)
    end
end
K_subj = length(complete_subjects)
n_visits = 3
d_concat = n_assays * n_visits

X_concat = Matrix{Float64}(undef, K_subj, d_concat)
for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_clean.SubjectID .== s) .& (df_clean.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_concat[i, offset + j] = df_clean[row, col]
        end
    end
end
@info "  Concatenated matrix: $K_subj subjects × $d_concat features"

# fit MVN on the concatenated data
μ_mvn = vec(mean(X_concat, dims=1))
Σ_mvn = cov(X_concat)
Σ_mvn = (Σ_mvn + Σ_mvn') / 2  # ensure symmetric
Σ_mvn += 1e-6 * I              # ensure positive definite

Random.seed!(123)
dist_mvn = MvNormal(μ_mvn, Σ_mvn)
mvn_concat = rand(dist_mvn, n_synth)'  # n_synth × d_concat

# split into per-visit DataFrames
df_mvn_visits = DataFrame[]
for v in 1:3
    offset = (v - 1) * n_assays
    cols_data = mvn_concat[:, (offset+1):(offset+n_assays)]
    df_v = DataFrame(cols_data, kept_cols)
    df_v.Visit = fill(v, n_synth)
    push!(df_mvn_visits, df_v)
end
df_mvn = vcat(df_mvn_visits...)
@info "  MVN generation complete: $(nrow(df_mvn)) records"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Evaluate — per-visit metrics
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 3: Computing comparison metrics …"

comparison_rows = []

for v in 1:3
    real_v = df_clean[df_clean.Visit .== v, :]
    sa_v = df_sa[df_sa.Visit .== v, :]
    mvn_v = df_mvn[df_mvn.Visit .== v, :]

    # ── Per-feature means and stds ────────────────────────────────────────
    sa_mean_re = Float64[]
    sa_std_re = Float64[]
    mvn_mean_re = Float64[]
    mvn_std_re = Float64[]
    sa_ks_vals = Float64[]
    mvn_ks_vals = Float64[]

    for col in kept_cols
        real_vals = collect(skipmissing(real_v[!, col]))
        sa_vals = sa_v[!, col]
        mvn_vals = mvn_v[!, col]

        r_mean = mean(real_vals)
        r_std = std(real_vals)
        sa_mean_val = mean(sa_vals)
        sa_std_val = std(sa_vals)
        mvn_mean_val = mean(mvn_vals)
        mvn_std_val = std(mvn_vals)

        # mean relative error
        if abs(r_mean) > 1e-12
            push!(sa_mean_re, abs(sa_mean_val - r_mean) / abs(r_mean))
            push!(mvn_mean_re, abs(mvn_mean_val - r_mean) / abs(r_mean))
        end

        # std relative error
        if abs(r_std) > 1e-12
            push!(sa_std_re, abs(sa_std_val - r_std) / abs(r_std))
            push!(mvn_std_re, abs(mvn_std_val - r_std) / abs(r_std))
        end

        # KS statistic
        push!(sa_ks_vals, ks_statistic(Float64.(real_vals), Float64.(sa_vals)))
        push!(mvn_ks_vals, ks_statistic(Float64.(real_vals), Float64.(mvn_vals)))
    end

    # ── Correlation MAE ───────────────────────────────────────────────────
    C_real = pairwise_cor(real_v, kept_cols)
    C_sa = pairwise_cor(sa_v, kept_cols)
    C_mvn = pairwise_cor(mvn_v, kept_cols)

    # replace NaN with 0 for comparison
    C_real[isnan.(C_real)] .= 0
    C_sa[isnan.(C_sa)] .= 0
    C_mvn[isnan.(C_mvn)] .= 0

    sa_corr_mae = mean(abs.(C_real - C_sa))
    mvn_corr_mae = mean(abs.(C_real - C_mvn))

    push!(comparison_rows, (
        Visit = v,
        SA_Mean_RE = round(mean(sa_mean_re), digits=4),
        MVN_Mean_RE = round(mean(mvn_mean_re), digits=4),
        SA_Std_RE = round(mean(sa_std_re), digits=4),
        MVN_Std_RE = round(mean(mvn_std_re), digits=4),
        SA_Corr_MAE = round(sa_corr_mae, digits=4),
        MVN_Corr_MAE = round(mvn_corr_mae, digits=4),
        SA_KS_Mean = round(mean(sa_ks_vals), digits=4),
        MVN_KS_Mean = round(mean(mvn_ks_vals), digits=4),
        SA_KS_Median = round(median(sa_ks_vals), digits=4),
        MVN_KS_Median = round(median(mvn_ks_vals), digits=4),
    ))
end

df_comparison = DataFrame(comparison_rows)

# save CSV
CSV.write(joinpath(_PATH_TO_DATA, "sa_vs_mvn_comparison.csv"), df_comparison)
@info "  Saved comparison to data/sa_vs_mvn_comparison.csv"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Per-feature KS statistics table (all visits)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 4: Per-feature KS statistics …"

ks_rows = []
for v in 1:3
    real_v = df_clean[df_clean.Visit .== v, :]
    sa_v = df_sa[df_sa.Visit .== v, :]
    mvn_v = df_mvn[df_mvn.Visit .== v, :]
    for col in kept_cols
        real_vals = Float64.(collect(skipmissing(real_v[!, col])))
        sa_vals = Float64.(sa_v[!, col])
        mvn_vals = Float64.(mvn_v[!, col])
        push!(ks_rows, (
            Visit = v,
            Feature = string(col),
            SA_KS = round(ks_statistic(real_vals, sa_vals), digits=4),
            MVN_KS = round(ks_statistic(real_vals, mvn_vals), digits=4),
        ))
    end
end
df_ks = DataFrame(ks_rows)
CSV.write(joinpath(_PATH_TO_DATA, "sa_vs_mvn_ks_by_feature.csv"), df_ks)
@info "  Saved per-feature KS to data/sa_vs_mvn_ks_by_feature.csv"

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: PCA figure (2 x 3 grid — top: SA, bottom: MVN)
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 5: PCA comparison figure …"

# ── Fit per-record PCA on real data (complete rows only) ──────────────────
real_rows = []
real_visits = Int[]
for row in 1:nrow(df_clean)
    vals = Float64[]
    complete = true
    for col in kept_cols
        v = df_clean[row, col]
        if ismissing(v)
            complete = false
            break
        end
        push!(vals, v)
    end
    if complete
        push!(real_rows, vals)
        push!(real_visits, df_clean[row, :Visit])
    end
end
X_real = reduce(hcat, real_rows)  # n_assays × n_real

# standardize (fit on real)
μ_r = vec(mean(X_real, dims=2))
σ_r = vec(std(X_real, dims=2))
for i in eachindex(σ_r)
    σ_r[i] < 1e-12 && (σ_r[i] = 1.0)
end
Z_real = (X_real .- μ_r) ./ σ_r

# fit PCA on real
pca_rec = MultivariateStats.fit(PCA, Z_real; maxoutdim=2)
P_real = MultivariateStats.transform(pca_rec, Z_real)  # 2 × n_real

var1 = round(100 * MultivariateStats.principalvars(pca_rec)[1] / MultivariateStats.tvar(pca_rec), digits=1)
var2 = round(100 * MultivariateStats.principalvars(pca_rec)[2] / MultivariateStats.tvar(pca_rec), digits=1)

# transform SA synthetic
X_sa = Matrix{Float64}(undef, n_assays, nrow(df_sa))
for i in 1:nrow(df_sa)
    for (j, col) in enumerate(kept_cols)
        X_sa[j, i] = df_sa[i, col]
    end
end
Z_sa = (X_sa .- μ_r) ./ σ_r
P_sa = MultivariateStats.transform(pca_rec, Z_sa)

# transform MVN synthetic
X_mvn = Matrix{Float64}(undef, n_assays, nrow(df_mvn))
for i in 1:nrow(df_mvn)
    for (j, col) in enumerate(kept_cols)
        X_mvn[j, i] = df_mvn[i, col]
    end
end
Z_mvn = (X_mvn .- μ_r) ./ σ_r
P_mvn = MultivariateStats.transform(pca_rec, Z_mvn)

sa_visits = df_sa.Visit
mvn_visits = df_mvn.Visit

# ── Colors (matching paper style) ────────────────────────────────────────
synth_colors = [
    RGB(0.55, 0.82, 0.82),  # V1: teal
    RGB(0.55, 0.78, 0.55),  # V2: green
    RGB(0.90, 0.65, 0.50),  # V3: salmon
]
real_colors = [
    RGB(0.05, 0.55, 0.55),  # V1: dark teal
    RGB(0.10, 0.52, 0.10),  # V2: dark green
    RGB(0.80, 0.30, 0.10),  # V3: dark coral
]
visit_labels = ["V1: not pregnant", "V2: first trimester", "V3: third trimester"]

# ── Build 2 × 3 panel figure ─────────────────────────────────────────────
panels = []

# Top row: SA by visit
for v in 1:3
    real_mask = real_visits .== v
    sa_mask = sa_visits .== v

    p = scatter(P_real[1, real_mask], P_real[2, real_mask],
        label="Real",
        color=real_colors[v],
        alpha=0.8, markersize=5,
        markerstrokecolor=:white, markerstrokewidth=0.5,
        markershape=:circle,
        xlabel="PC1 ($(var1)%)", ylabel="PC2 ($(var2)%)",
        title="SA — $(visit_labels[v])",
        titlefontsize=9,
        background_color=RGB(0.96, 0.96, 0.96),
        foreground_color=:black,
        grid=false, framestyle=:box,
        legendfontsize=7,
        margin=3Plots.mm,
    )
    scatter!(p, P_sa[1, sa_mask], P_sa[2, sa_mask],
        label="SA synth",
        color=synth_colors[v],
        alpha=0.4, markersize=3, markerstrokewidth=0,
        markershape=:circle)
    push!(panels, p)
end

# Bottom row: MVN by visit
for v in 1:3
    real_mask = real_visits .== v
    mvn_mask = mvn_visits .== v

    p = scatter(P_real[1, real_mask], P_real[2, real_mask],
        label="Real",
        color=real_colors[v],
        alpha=0.8, markersize=5,
        markerstrokecolor=:white, markerstrokewidth=0.5,
        markershape=:circle,
        xlabel="PC1 ($(var1)%)", ylabel="PC2 ($(var2)%)",
        title="MVN — $(visit_labels[v])",
        titlefontsize=9,
        background_color=RGB(0.96, 0.96, 0.96),
        foreground_color=:black,
        grid=false, framestyle=:box,
        legendfontsize=7,
        margin=3Plots.mm,
    )
    scatter!(p, P_mvn[1, mvn_mask], P_mvn[2, mvn_mask],
        label="MVN synth",
        color=synth_colors[v],
        alpha=0.4, markersize=3, markerstrokewidth=0,
        markershape=:diamond)
    push!(panels, p)
end

p_pca = plot(panels..., layout=(2, 3), size=(1500, 900),
    plot_title="PCA Projection: SA vs MVN by Visit",
    plot_titlefontsize=14,
    background_color=:white,
    margin=5Plots.mm)

savefig(p_pca, joinpath(_PATH_TO_FIG, "sa_vs_mvn_pca_by_visit.pdf"))
savefig(p_pca, joinpath(_PATH_TO_FIG, "sa_vs_mvn_pca_by_visit.png"))
@info "  Saved PCA figure to figs/sa_vs_mvn_pca_by_visit.pdf"

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Bar chart — SA vs MVN on key metrics per visit
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 6: Metrics bar chart …"

sa_color = RGB(0.20, 0.50, 0.72)   # blue
mvn_color = RGB(0.85, 0.45, 0.25)  # orange

metric_names = ["Mean RE", "Std RE", "Corr MAE", "KS (mean)"]
n_metrics = length(metric_names)

bar_panels = []
for v in 1:3
    row = df_comparison[v, :]
    sa_vals = [row.SA_Mean_RE, row.SA_Std_RE, row.SA_Corr_MAE, row.SA_KS_Mean]
    mvn_vals = [row.MVN_Mean_RE, row.MVN_Std_RE, row.MVN_Corr_MAE, row.MVN_KS_Mean]

    x_pos = 1:n_metrics
    p = groupedbar(
        [sa_vals mvn_vals],
        bar_position=:dodge, bar_width=0.7,
        xticks=(1:n_metrics, metric_names),
        xrotation=15,
        ylabel="Error",
        title="Visit $v",
        titlefontsize=11,
        label=["SA" "MVN"],
        color=[sa_color mvn_color],
        alpha=0.85,
        background_color=RGB(0.96, 0.96, 0.96),
        foreground_color=:black,
        grid=false, framestyle=:box,
        legendfontsize=8,
        margin=3Plots.mm,
    )
    push!(bar_panels, p)
end

p_bars = plot(bar_panels..., layout=(1, 3), size=(1400, 450),
    plot_title="SA vs MVN: Metric Comparison by Visit",
    plot_titlefontsize=14,
    background_color=:white,
    margin=5Plots.mm)

savefig(p_bars, joinpath(_PATH_TO_FIG, "sa_vs_mvn_metrics_bar.pdf"))
savefig(p_bars, joinpath(_PATH_TO_FIG, "sa_vs_mvn_metrics_bar.png"))
@info "  Saved bar chart to figs/sa_vs_mvn_metrics_bar.pdf"

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Print summary table
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^90)
println("SA vs MVN Head-to-Head Comparison")
println("="^90)
println()

for v in 1:3
    row = df_comparison[v, :]
    println("  Visit $v:")
    println("    Mean RE:    SA = $(row.SA_Mean_RE)    MVN = $(row.MVN_Mean_RE)    " *
            (row.SA_Mean_RE < row.MVN_Mean_RE ? "← SA wins" : "← MVN wins"))
    println("    Std RE:     SA = $(row.SA_Std_RE)    MVN = $(row.MVN_Std_RE)    " *
            (row.SA_Std_RE < row.MVN_Std_RE ? "← SA wins" : "← MVN wins"))
    println("    Corr MAE:   SA = $(row.SA_Corr_MAE)    MVN = $(row.MVN_Corr_MAE)    " *
            (row.SA_Corr_MAE < row.MVN_Corr_MAE ? "← SA wins" : "← MVN wins"))
    println("    KS (mean):  SA = $(row.SA_KS_Mean)    MVN = $(row.MVN_KS_Mean)    " *
            (row.SA_KS_Mean < row.MVN_KS_Mean ? "← SA wins" : "← MVN wins"))
    println("    KS (med):   SA = $(row.SA_KS_Median)    MVN = $(row.MVN_KS_Median)    " *
            (row.SA_KS_Median < row.MVN_KS_Median ? "← SA wins" : "← MVN wins"))
    println()
end

# overall wins
global sa_wins = 0
global mvn_wins = 0
for v in 1:3
    global sa_wins, mvn_wins
    row = df_comparison[v, :]
    row.SA_Mean_RE < row.MVN_Mean_RE ? (sa_wins += 1) : (mvn_wins += 1)
    row.SA_Std_RE < row.MVN_Std_RE ? (sa_wins += 1) : (mvn_wins += 1)
    row.SA_Corr_MAE < row.MVN_Corr_MAE ? (sa_wins += 1) : (mvn_wins += 1)
    row.SA_KS_Mean < row.MVN_KS_Mean ? (sa_wins += 1) : (mvn_wins += 1)
end
println("  Overall: SA wins $sa_wins / $(sa_wins + mvn_wins) metric-visit matchups")
println("="^90)

@info "\n✓ Done — SA vs MVN comparison complete."
@info "  Figures: figs/sa_vs_mvn_pca_by_visit.pdf, figs/sa_vs_mvn_metrics_bar.pdf"
@info "  Data:    data/sa_vs_mvn_comparison.csv, data/sa_vs_mvn_ks_by_feature.csv"

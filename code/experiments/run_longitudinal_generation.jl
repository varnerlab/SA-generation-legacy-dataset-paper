# ──────────────────────────────────────────────────────────────────────────────
# LEGACY / EXPLORATORY — DO NOT USE FOR PAPER ARTIFACTS.
# This script regenerates synthetic data with non-canonical hyperparameters
# (different α, seed, or N from run_full_longitudinal.jl). Its CSV outputs are
# only consumed by itself and are not used by the paper or supplement.
# Canonical generation: code/experiments/run_full_longitudinal.jl
# ──────────────────────────────────────────────────────────────────────────────
# run_longitudinal_generation.jl
#
# Generate synthetic patient *trajectories* (visit 1 → 2 → 3) using SA.
#
# Key idea: concatenate each patient's 3-visit measurement vectors into a
# single pattern [V1 | V2 | V3]. SA then generates full longitudinal
# trajectories as atomic units, preserving within-patient correlations
# across pregnancy visits.
#
# This is analogous to how protein SA treats a full sequence (L positions ×
# 20 amino acids) as one pattern — here each "position" is a visit and
# each "channel" is a biochemical assay.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Build concatenated visit vectors
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Loading and restructuring data …"
JLD2.@load joinpath(_PATH_TO_DATA, "patient_memory.jld2") df_clean kept_cols

# find subjects with all 3 visits
subjects = unique(df_clean.SubjectID)
complete_subjects = Int[]
for s in subjects
    visits = sort(unique(df_clean[df_clean.SubjectID .== s, :Visit]))
    if visits == [1, 2, 3]
        push!(complete_subjects, s)
    end
end
K = length(complete_subjects)
@info "  Complete subjects (all 3 visits): $K"

n_assays = length(kept_cols)
n_visits = 3
d_concat = n_assays * n_visits
@info "  Concatenated dimension: $n_assays assays × $n_visits visits = $d_concat"

# build concatenated matrix: each row is one subject, columns = [V1_assays | V2_assays | V3_assays]
X_concat = Matrix{Float64}(undef, K, d_concat)
subject_conditions = String[]
subject_pe = String[]
subject_htn = String[]

for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_clean.SubjectID .== s) .& (df_clean.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_concat[i, offset + j] = df_clean[row, col]
        end
    end
    # metadata from first visit
    first_row = findfirst((df_clean.SubjectID .== s) .& (df_clean.Visit .== 1))
    if first_row === nothing
        first_row = findfirst(df_clean.SubjectID .== s)
    end
    push!(subject_conditions, df_clean[first_row, :Condition])
    push!(subject_pe, df_clean[first_row, :DevelopedPE])
    push!(subject_htn, df_clean[first_row, :DevelopedHTN])
end

@info "  Conditions: $(sort(collect(StatsBase.countmap(subject_conditions))))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Standardize + PCA on concatenated vectors
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 2: Building longitudinal memory matrix …"

# create column names for concatenated vector
concat_col_names = Symbol[]
for v in 1:n_visits
    for col in kept_cols
        push!(concat_col_names, Symbol("V$(v)_$(col)"))
    end
end

# standardize
μ_concat = vec(mean(X_concat, dims=1))
σ_concat = vec(std(X_concat, dims=1))
for i in eachindex(σ_concat)
    σ_concat[i] < 1e-12 && (σ_concat[i] = 1.0)
end
Z_concat = (X_concat .- μ_concat') ./ σ_concat'
std_params_concat = StandardizationParams(μ_concat, σ_concat, concat_col_names)
@info "  Standardized: $K subjects × $d_concat features"

# PCA
Z_t = Z_concat'  # d_concat × K
pca_concat = MultivariateStats.fit(PCA, Z_t; pratio=0.95)
d_pca = MultivariateStats.outdim(pca_concat)
X_pca = MultivariateStats.transform(pca_concat, Z_t)  # d_pca × K
var_retained = round(100 * sum(MultivariateStats.principalvars(pca_concat)) /
                     MultivariateStats.tvar(pca_concat), digits=1)
@info "  PCA: $d_concat → $d_pca dimensions ($var_retained% variance)"

# compute norms before unit-normalizing
pca_norms = [norm(X_pca[:, k]) for k in 1:K]
@info "  PCA norms: range=[$(round(minimum(pca_norms),digits=2)), $(round(maximum(pca_norms),digits=2))], mean=$(round(mean(pca_norms),digits=2))"

# unit-norm memory
X̂ = copy(X_pca)
for k in 1:K
    X̂[:, k] ./= (norm(X̂[:, k]) + 1e-12)
end
@info "  Memory: $d_pca × $K (unit-norm)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Phase transition
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 3: Finding β* …"
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80,
                                 β_range=(0.1, 1000.0))
β_star = phase.β_star
@info "  β* = $(round(β_star, digits=2)), √d = $(round(sqrt(d_pca), digits=2))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Generate synthetic trajectories (rescaled)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Generating synthetic patient trajectories …"
n_synth = 100
T_steps = 2000

Random.seed!(42)

synth_concat = Matrix{Float64}(undef, n_synth, d_concat)
synth_pca = Matrix{Float64}(undef, d_pca, n_synth)

for i in 1:n_synth
    # direction from SA
    k0 = rand(1:K)
    ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d_pca)
    ξ₀ ./= norm(ξ₀)

    result = sample(X̂, ξ₀, T_steps; β=β_star, α=0.05)
    ξ_dir = result.Ξ[end, :]
    ξ_dir ./= (norm(ξ_dir) + 1e-12)

    # magnitude from empirical
    r = rand(pca_norms)
    ξ_rescaled = r .* ξ_dir

    synth_pca[:, i] = ξ_rescaled

    # decode: inverse PCA → destandardize
    z_std = vec(MultivariateStats.reconstruct(pca_concat, ξ_rescaled))
    x_orig = z_std .* σ_concat .+ μ_concat
    synth_concat[i, :] = x_orig
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Split back into per-visit DataFrames
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 5: Splitting trajectories into per-visit records …"

df_synth_visits = DataFrame[]
for v in 1:n_visits
    offset = (v - 1) * n_assays
    cols = synth_concat[:, (offset+1):(offset+n_assays)]
    df_v = DataFrame(cols, kept_cols)
    df_v.SyntheticID = 1:n_synth
    df_v.Visit = fill(v, n_synth)
    push!(df_synth_visits, df_v)
end

df_synth_long = vcat(df_synth_visits...)
sort!(df_synth_long, [:SyntheticID, :Visit])
CSV.write(joinpath(_PATH_TO_DATA, "synthetic_longitudinal.csv"), df_synth_long)
@info "  Saved $(nrow(df_synth_long)) records ($(n_synth) patients × 3 visits)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Evaluate
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 6: Evaluation …"

# per-visit mean comparison
@info "  Per-visit feature means (select features):"
show_cols = kept_cols[1:min(8, length(kept_cols))]
for col in show_cols
    print("  $(rpad(string(col), 14))")
    for v in 1:3
        real_vals = collect(skipmissing(df_clean[df_clean.Visit .== v, col]))
        synth_vals = df_synth_long[df_synth_long.Visit .== v, col]
        μ_r = round(mean(real_vals), digits=1)
        μ_s = round(mean(synth_vals), digits=1)
        print("  V$v: $μ_r→$μ_s")
    end
    println()
end

# per-visit feature comparison
for v in 1:3
    real_v = df_clean[df_clean.Visit .== v, :]
    synth_v = df_synth_long[df_synth_long.Visit .== v, :]
    comp = feature_summary_comparison(real_v, synth_v, kept_cols)
    @info "  Visit $v mean relative error: $(round(mean(comp.Mean_Rel_Error), digits=4))"
end

# within-patient visit-to-visit correlation
# for real data: correlation between V1 and V3 of same subject
@info "\n  Within-patient V1→V3 correlation (real vs synthetic):"
for col in show_cols
    # real
    real_v1 = Float64[]
    real_v3 = Float64[]
    for s in complete_subjects
        v1_val = df_clean[(df_clean.SubjectID .== s) .& (df_clean.Visit .== 1), col]
        v3_val = df_clean[(df_clean.SubjectID .== s) .& (df_clean.Visit .== 3), col]
        if !isempty(v1_val) && !isempty(v3_val) && !ismissing(v1_val[1]) && !ismissing(v3_val[1])
            push!(real_v1, v1_val[1])
            push!(real_v3, v3_val[1])
        end
    end
    # synthetic
    synth_v1 = df_synth_long[df_synth_long.Visit .== 1, col]
    synth_v3 = df_synth_long[df_synth_long.Visit .== 3, col]

    r_real = length(real_v1) > 2 ? round(cor(real_v1, real_v3), digits=3) : NaN
    r_synth = round(cor(synth_v1, synth_v3), digits=3)
    println("  $(rpad(string(col), 14))  real=$r_real  synth=$r_synth")
end

# PCA dispersion
X_real_concat_pca = X_pca  # d_pca × K
real_centroid = vec(mean(X_real_concat_pca, dims=2))
real_dists = [norm(X_real_concat_pca[:, k] - real_centroid) for k in 1:K]
synth_centroid = vec(mean(synth_pca, dims=2))
synth_dists = [norm(synth_pca[:, j] - synth_centroid) for j in 1:n_synth]
@info "\n  Dispersion ratio: $(round(mean(synth_dists)/mean(real_dists), digits=3))"
@info "  PC1 σ ratio: $(round(std(synth_pca[1,:])/std(X_real_concat_pca[1,:]), digits=3))"
@info "  PC2 σ ratio: $(round(std(synth_pca[2,:])/std(X_real_concat_pca[2,:]), digits=3))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Figures
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 7: Generating figures …"

# --- Figure 1: PCA space (longitudinal) ---
p1 = scatter(X_real_concat_pca[1, :], X_real_concat_pca[2, :],
    label="Real", color=:steelblue, alpha=0.7, markersize=6,
    xlabel="PC1", ylabel="PC2",
    title="Longitudinal PCA: Real vs Synthetic (K=$K, d=$d_pca)")
scatter!(p1, synth_pca[1, :], synth_pca[2, :],
    label="Synthetic", color=:coral, alpha=0.5, markersize=5, markershape=:diamond)
savefig(p1, joinpath(_PATH_TO_FIG, "longitudinal_pca.pdf"))
savefig(p1, joinpath(_PATH_TO_FIG, "longitudinal_pca.png"))

# --- Figure 2: PCA colored by condition ---
cond_colors = Dict("Healthy" => :steelblue, "Prior PE" => :coral,
                    "PCOS" => :seagreen, "Multi" => :purple)
p2 = plot(xlabel="PC1", ylabel="PC2",
          title="Longitudinal PCA by Condition")
for cond in unique(subject_conditions)
    mask = subject_conditions .== cond
    color = get(cond_colors, cond, :gray)
    scatter!(p2, X_real_concat_pca[1, mask], X_real_concat_pca[2, mask],
        label=cond, color=color, alpha=0.7, markersize=6)
end
savefig(p2, joinpath(_PATH_TO_FIG, "longitudinal_pca_by_condition.pdf"))
savefig(p2, joinpath(_PATH_TO_FIG, "longitudinal_pca_by_condition.png"))

# --- Figure 3: Visit-specific feature trajectories ---
# show how select features evolve across visits (real vs synthetic)
trajectory_cols = [:Estradiol, :Progesterone, :VIII, :Fbgn, Symbol("PAI1"), Symbol("D-dimer")]
trajectory_cols = [c for c in trajectory_cols if c in kept_cols]
n_traj = length(trajectory_cols)

p_trajs = []
for col in trajectory_cols
    real_means = Float64[]
    real_sds = Float64[]
    synth_means = Float64[]
    synth_sds = Float64[]
    for v in 1:3
        rv = collect(skipmissing(df_clean[df_clean.Visit .== v, col]))
        sv = df_synth_long[df_synth_long.Visit .== v, col]
        push!(real_means, mean(rv))
        push!(real_sds, std(rv))
        push!(synth_means, mean(sv))
        push!(synth_sds, std(sv))
    end
    p = plot([1,2,3], real_means, ribbon=real_sds, fillalpha=0.2,
        label="Real", color=:steelblue, lw=2, marker=:circle,
        xlabel="Visit", ylabel=string(col),
        title=string(col), titlefontsize=9,
        xticks=([1,2,3], ["Baseline","1st Tri","3rd Tri"]),
        xrotation=15)
    plot!(p, [1,2,3], synth_means, ribbon=synth_sds, fillalpha=0.2,
        label="Synthetic", color=:coral, lw=2, marker=:diamond)
    push!(p_trajs, p)
end
n_rows = ceil(Int, n_traj / 3)
p3 = plot(p_trajs..., layout=(n_rows, min(3, n_traj)),
          size=(900, 300*n_rows),
          plot_title="Longitudinal Trajectories: Real vs Synthetic")
savefig(p3, joinpath(_PATH_TO_FIG, "longitudinal_trajectories.pdf"))
savefig(p3, joinpath(_PATH_TO_FIG, "longitudinal_trajectories.png"))

# --- Figure 4: Per-visit feature distributions ---
p_visit_dists = []
for v in 1:3
    col = :Estradiol
    rv = collect(skipmissing(df_clean[df_clean.Visit .== v, col]))
    sv = df_synth_long[df_synth_long.Visit .== v, col]
    p = histogram(rv, normalize=:pdf, alpha=0.5, label="Real",
        color=:steelblue, bins=15,
        title="Estradiol — Visit $v", titlefontsize=9)
    histogram!(p, sv, normalize=:pdf, alpha=0.5, label="Synthetic",
        color=:coral, bins=15)
    push!(p_visit_dists, p)
end
p4 = plot(p_visit_dists..., layout=(1, 3), size=(900, 300),
          plot_title="Estradiol Distribution by Visit")
savefig(p4, joinpath(_PATH_TO_FIG, "longitudinal_estradiol_by_visit.pdf"))
savefig(p4, joinpath(_PATH_TO_FIG, "longitudinal_estradiol_by_visit.png"))

# --- Figure 5: Entropy curve ---
p5 = plot(phase.βs, phase.Hs_norm,
    xlabel="β", ylabel="H / log(K)", xscale=:log10,
    title="Entropy Curve (Longitudinal, d=$d_pca, K=$K)",
    lw=2, label="H(β)", color=:steelblue)
vline!(p5, [β_star], ls=:dash, color=:red, lw=2,
       label="β* = $(round(β_star, digits=1))")
vline!(p5, [sqrt(d_pca)], ls=:dot, color=:blue, lw=2,
       label="√d = $(round(sqrt(d_pca), digits=1))")
savefig(p5, joinpath(_PATH_TO_FIG, "longitudinal_entropy.pdf"))
savefig(p5, joinpath(_PATH_TO_FIG, "longitudinal_entropy.png"))

# save metadata
JLD2.@save joinpath(_PATH_TO_DATA, "longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms complete_subjects subject_conditions kept_cols n_assays

@info "\n✓ Done!"
@info "  Output: data/synthetic_longitudinal.csv"
@info "  $(n_synth) synthetic patients × 3 visits = $(3*n_synth) records"
@info "  Figures: figs/longitudinal_*.pdf"

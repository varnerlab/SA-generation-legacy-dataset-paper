# ──────────────────────────────────────────────────────────────────────────────
# run_full_longitudinal.jl
#
# Full pipeline with ALL measurement columns (biochemical + TGA + ROTEM),
# longitudinal generation (concatenated visits), rescaled SA.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load ALL columns
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Loading full dataset (all measurement columns) …"
xlsx_path = joinpath(_PATH_TO_ORIGINAL_DATA, "20200604_R61_Legacy_Biochemical_Measurements.xlsx")
df_raw = load_legacy_data(xlsx_path)
@info "  Loaded $(nrow(df_raw)) records, $(ncol(df_raw)) columns"

acols = assay_columns(df_raw)
@info "  Total assay columns: $(length(acols))"

# clean
df_clean, kept_cols = clean_assay_data(df_raw; max_missing_frac=0.3)
@info "  Clean: $(nrow(df_clean)) records × $(length(kept_cols)) assays"
@info "  Kept columns:"
for c in kept_cols
    println("    ", c)
end

# save
CSV.write(joinpath(_PATH_TO_DATA, "cleaned_full_data.csv"), df_clean)

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Build longitudinal concatenated patterns
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 2: Building longitudinal patterns …"

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
n_assays = length(kept_cols)
n_visits = 3
d_concat = n_assays * n_visits
@info "  Complete subjects: $K"
@info "  Concatenated dim: $n_assays assays × $n_visits visits = $d_concat"

# build concatenated matrix
X_concat = Matrix{Float64}(undef, K, d_concat)
subject_conditions = String[]

for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_clean.SubjectID .== s) .& (df_clean.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_concat[i, offset + j] = df_clean[row, col]
        end
    end
    first_row = findfirst((df_clean.SubjectID .== s) .& (df_clean.Visit .== 1))
    if first_row === nothing
        first_row = findfirst(df_clean.SubjectID .== s)
    end
    push!(subject_conditions, df_clean[first_row, :Condition])
end
@info "  Conditions: $(sort(collect(StatsBase.countmap(subject_conditions))))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Standardize + PCA
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 3: PCA …"

concat_col_names = Symbol[]
for v in 1:n_visits
    for col in kept_cols
        push!(concat_col_names, Symbol("V$(v)_$(col)"))
    end
end

μ_concat = vec(mean(X_concat, dims=1))
σ_concat = vec(std(X_concat, dims=1))
for i in eachindex(σ_concat)
    σ_concat[i] < 1e-12 && (σ_concat[i] = 1.0)
end
Z_concat = (X_concat .- μ_concat') ./ σ_concat'
std_params_concat = StandardizationParams(μ_concat, σ_concat, concat_col_names)

Z_t = Z_concat'
pca_concat = MultivariateStats.fit(PCA, Z_t; pratio=0.95)
d_pca = MultivariateStats.outdim(pca_concat)
X_pca = MultivariateStats.transform(pca_concat, Z_t)
var_retained = round(100 * sum(MultivariateStats.principalvars(pca_concat)) /
                     MultivariateStats.tvar(pca_concat), digits=1)
@info "  PCA: $d_concat → $d_pca dimensions ($var_retained% variance)"

pca_norms = [norm(X_pca[:, k]) for k in 1:K]
X̂ = copy(X_pca)
for k in 1:K
    X̂[:, k] ./= (norm(X̂[:, k]) + 1e-12)
end
@info "  Memory: $d_pca × $K"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: β* and generate
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Phase transition + generation …"
phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
β_star = phase.β_star

n_synth = 100
T_steps = 2000
Random.seed!(42)

synth_concat = Matrix{Float64}(undef, n_synth, d_concat)
synth_pca = Matrix{Float64}(undef, d_pca, n_synth)

for i in 1:n_synth
    k0 = rand(1:K)
    ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d_pca)
    ξ₀ ./= norm(ξ₀)
    result = sample(X̂, ξ₀, T_steps; β=β_star, α=0.05)
    ξ_dir = result.Ξ[end, :]
    ξ_dir ./= (norm(ξ_dir) + 1e-12)
    r = rand(pca_norms)
    ξ_rescaled = r .* ξ_dir
    synth_pca[:, i] = ξ_rescaled
    z_std = vec(MultivariateStats.reconstruct(pca_concat, ξ_rescaled))
    synth_concat[i, :] = z_std .* σ_concat .+ μ_concat
end

# split into per-visit DataFrames
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
CSV.write(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), df_synth_long)
@info "  Generated $(n_synth) patients × 3 visits = $(nrow(df_synth_long)) records"

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Evaluation
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 5: Evaluation …"
for v in 1:3
    real_v = df_clean[df_clean.Visit .== v, :]
    synth_v = df_synth_long[df_synth_long.Visit .== v, :]
    comp = feature_summary_comparison(real_v, synth_v, kept_cols)
    @info "  Visit $v mean relative error: $(round(mean(comp.Mean_Rel_Error), digits=4))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: X vs MCF scatter (1 × 3 grid by visit)
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 6: X vs MCF plot …"

# find the MCF column name (TF Only MCF)
global mcf_col = nothing
for c in kept_cols
    global mcf_col
    if occursin("MCF", string(c)) && occursin("TF Only", string(c))
        mcf_col = c
        @info "  Found MCF (TF Only) column: $c"
        break
    end
end
# fallback: any MCF column
if mcf_col === nothing
    for c in kept_cols
        global mcf_col
        if occursin("MCF", string(c))
            mcf_col = c
            @info "  Found MCF column (fallback): $c"
            break
        end
    end
end
x_col = Symbol("X")
@info "  Using X column: $x_col"

if mcf_col !== nothing && x_col in kept_cols
    panels = []
    visit_labels = ["Visit 1 (Baseline)", "Visit 2 (1st Trimester)", "Visit 3 (3rd Trimester)"]
    for v in 1:3
        # real data
        real_v = df_clean[df_clean.Visit .== v, :]
        real_x = collect(skipmissing(real_v[!, x_col]))
        real_mcf = collect(skipmissing(real_v[!, mcf_col]))
        # need paired data (both non-missing in same row)
        real_pairs_x = Float64[]
        real_pairs_mcf = Float64[]
        for row in 1:nrow(real_v)
            xv = real_v[row, x_col]
            mv = real_v[row, mcf_col]
            if !ismissing(xv) && !ismissing(mv)
                push!(real_pairs_x, xv)
                push!(real_pairs_mcf, mv)
            end
        end

        # synthetic data
        synth_v = df_synth_long[df_synth_long.Visit .== v, :]
        synth_x = synth_v[!, x_col]
        synth_mcf = synth_v[!, mcf_col]

        p = scatter(real_pairs_x, real_pairs_mcf,
            label="Real", color=:steelblue, alpha=0.7, markersize=6,
            xlabel="Factor X (%)", ylabel="MCF (mm)",
            title=visit_labels[v], titlefontsize=10)
        scatter!(p, synth_x, synth_mcf,
            label="Synthetic", color=:coral, alpha=0.5, markersize=5,
            markershape=:diamond)
        push!(panels, p)
    end

    p_xmcf = plot(panels..., layout=(1, 3), size=(1200, 400),
                  plot_title="Factor X vs MCF by Visit")
    savefig(p_xmcf, joinpath(_PATH_TO_FIG, "X_vs_MCF_by_visit.pdf"))
    savefig(p_xmcf, joinpath(_PATH_TO_FIG, "X_vs_MCF_by_visit.png"))
    @info "  Saved X_vs_MCF_by_visit.pdf"
else
    @warn "  Could not find MCF or X column!"
    @info "  Available columns with MCF: $(filter(c -> occursin("MCF", string(c)), string.(kept_cols)))"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Trajectory plots for key features including ROTEM
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 7: Trajectory plots …"

# pick features spanning biochemical, TGA, and ROTEM
trajectory_candidates = [:Estradiol, :VIII, :Fbgn, Symbol("D-dimer"), x_col]
if mcf_col !== nothing
    push!(trajectory_candidates, mcf_col)
end
# find TGA feature
for c in kept_cols
    if occursin("ETP", string(c))
        push!(trajectory_candidates, c)
        break
    end
end
trajectory_cols = [c for c in trajectory_candidates if c in kept_cols]
trajectory_cols = unique(trajectory_cols)

p_trajs = []
for col in trajectory_cols
    real_means = Float64[]
    real_sds = Float64[]
    synth_means = Float64[]
    synth_sds = Float64[]
    for v in 1:3
        rv = collect(skipmissing(df_clean[df_clean.Visit .== v, col]))
        sv = df_synth_long[df_synth_long.Visit .== v, col]
        push!(real_means, isempty(rv) ? NaN : mean(rv))
        push!(real_sds, isempty(rv) ? 0.0 : std(rv))
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
n_rows = ceil(Int, length(p_trajs) / 3)
p_all_traj = plot(p_trajs..., layout=(n_rows, min(3, length(p_trajs))),
                  size=(900, 300*n_rows),
                  plot_title="Longitudinal Trajectories (Full Feature Set)")
savefig(p_all_traj, joinpath(_PATH_TO_FIG, "full_longitudinal_trajectories.pdf"))
savefig(p_all_traj, joinpath(_PATH_TO_FIG, "full_longitudinal_trajectories.png"))

# save pipeline
JLD2.@save joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms complete_subjects subject_conditions kept_cols n_assays df_clean

@info "\n✓ Done!"
@info "  $(length(kept_cols)) assay columns used"
@info "  Output: data/synthetic_full_longitudinal.csv"

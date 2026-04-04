# ──────────────────────────────────────────────────────────────────────────────
# run_conditioned_generation.jl
#
# Generate three conditioned synthetic cohorts using multiplicity-weighted SA:
#   1. Uncomplicated  (18 real patients → 100 synthetic)  [no PE outcome]
#   2. PCOS           (3 real patients  → 100 synthetic)  [cross-cutting comorbidity]
#   3. Developed PE   (5 real patients  → 100 synthetic)  [developed PE during study]
#
# Each synthetic patient is a full longitudinal trajectory (V1|V2|V3).
# MVN cannot do this — you can't fit a covariance matrix from 3 patients.
#
# The multiplicity ratio ρ controls conditioning strength. For each cohort,
# we upweight the designated patients and find β*(ρ) via entropy inflection.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load the saved SA pipeline artefacts
# ══════════════════════════════════════════════════════════════════════════════
@info "Step 1: Loading pipeline artefacts …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms complete_subjects subject_conditions kept_cols n_assays df_clean

d_pca, K = size(X̂)
n_visits = 3
d_concat = n_assays * n_visits
μ_concat = std_params_concat.μ
σ_concat = std_params_concat.σ

@info "  Memory: $d_pca × $K ($(length(complete_subjects)) complete patients)"
@info "  Conditions: $(sort(collect(StatsBase.countmap(subject_conditions))))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Identify subpopulation indices in the memory matrix
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 2: Identifying subpopulations …"

# PE cohort: patients who developed preeclampsia during study pregnancy
pe_indices = Int[]
for (i, s) in enumerate(complete_subjects)
    row = findfirst((df_clean.SubjectID .== s))
    if df_clean[row, :DevelopedPE] == "Yes"
        push!(pe_indices, i)
    end
end

# Uncomplicated cohort: all patients who did NOT develop PE (n=18)
uncomplicated_indices = findall(i -> !(i in pe_indices), 1:K)

# PCOS cohort: cross-cutting comorbidity (n=3, 1 overlaps with PE)
pcos_indices = findall(c -> c == "PCOS", subject_conditions)

@info "  Uncomplicated: $(length(uncomplicated_indices)) patients (indices: $uncomplicated_indices)"
@info "  PCOS:          $(length(pcos_indices)) patients (indices: $pcos_indices)"
@info "  Developed PE:  $(length(pe_indices)) patients (indices: $pe_indices)"

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Fisher separation index — predict conditioning success
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 3: Fisher separation analysis …"

S_uncomplicated = fisher_separation_index(X̂, uncomplicated_indices)
S_pcos = fisher_separation_index(X̂, pcos_indices)
S_pe = fisher_separation_index(X̂, pe_indices)

@info "  Fisher separation index:"
@info "    Uncomplicated vs rest:  S = $(round(S_uncomplicated, digits=4))"
@info "    PCOS vs rest:           S = $(round(S_pcos, digits=4))"
@info "    Developed PE vs rest:   S = $(round(S_pe, digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Configure multiplicity weights and find β*(ρ) for each cohort
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Multiplicity weighting …"

# target: 80% of softmax attention on the designated subpopulation
f_target = 0.80

cohorts = [
    (name="Uncomplicated", indices=uncomplicated_indices),
    (name="PCOS",          indices=pcos_indices),
    (name="Developed PE",  indices=pe_indices),
]

cohort_configs = []
for cohort in cohorts
    K_sub = length(cohort.indices)
    K_bg = K - K_sub

    ρ = ρ_for_target_fraction(K_sub, K_bg, f_target)
    r = multiplicity_vector(K, cohort.indices; ρ=ρ)
    f_eff = effective_subset_fraction(r, cohort.indices)
    K_eff = effective_num_patterns(r)

    # find β* for this weighting
    pt = find_weighted_entropy_inflection(X̂, r; α=0.01, n_betas=80, β_range=(0.1, 1000.0))

    @info "  $(cohort.name):"
    @info "    K_sub=$K_sub, K_bg=$K_bg"
    @info "    ρ = $(round(ρ, digits=2)), f_eff = $(round(f_eff, digits=3)), K_eff = $(round(K_eff, digits=1))"
    @info "    β*(ρ) = $(round(pt.β_star, digits=2))"

    push!(cohort_configs, (
        name    = cohort.name,
        indices = cohort.indices,
        ρ       = ρ,
        r       = r,
        f_eff   = f_eff,
        K_eff   = K_eff,
        β_star  = pt.β_star,
    ))
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Generate conditioned synthetic cohorts
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 5: Generating conditioned cohorts …"

n_synth = 100
T_steps = 2000

all_cohort_results = Dict{String, DataFrame}()

for config in cohort_configs

    @info "  Generating $(config.name) cohort (n=$n_synth, ρ=$(round(config.ρ, digits=1)), β*=$(round(config.β_star, digits=2))) …"

    Random.seed!(42)

    synth_concat = Matrix{Float64}(undef, n_synth, d_concat)
    synth_pca_vecs = Matrix{Float64}(undef, d_pca, n_synth)

    for i in 1:n_synth
        # seed from the designated subpopulation
        k0 = rand(config.indices)
        ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d_pca)
        ξ₀ ./= norm(ξ₀)

        # weighted sample with multiplicity bias
        result = weighted_sample(X̂, ξ₀, T_steps, config.r;
                                  β=config.β_star, α=0.01)

        # extract final direction + rescale by empirical norm
        ξ_dir = result.Ξ[end, :]
        ξ_dir ./= (norm(ξ_dir) + 1e-12)
        r_norm = rand(pca_norms)
        ξ_rescaled = r_norm .* ξ_dir

        synth_pca_vecs[:, i] = ξ_rescaled

        # decode to physical space
        z_std = vec(MultivariateStats.reconstruct(pca_concat, ξ_rescaled))
        synth_concat[i, :] = z_std .* σ_concat .+ μ_concat
    end

    # split into per-visit DataFrames
    df_visits = DataFrame[]
    for v in 1:n_visits
        offset = (v - 1) * n_assays
        cols_data = synth_concat[:, (offset+1):(offset+n_assays)]
        df_v = DataFrame(cols_data, kept_cols)
        df_v.SyntheticID = 1:n_synth
        df_v.Visit = fill(v, n_synth)
        df_v.Cohort = fill(config.name, n_synth)
        push!(df_visits, df_v)
    end
    df_cohort = vcat(df_visits...)
    sort!(df_cohort, [:SyntheticID, :Visit])

    all_cohort_results[config.name] = df_cohort
    @info "    Generated $(nrow(df_cohort)) records ($(n_synth) patients × 3 visits)"
end

# combine all cohorts into one DataFrame
df_all_conditioned = vcat(values(all_cohort_results)...)
sort!(df_all_conditioned, [:Cohort, :SyntheticID, :Visit])
CSV.write(joinpath(_PATH_TO_DATA, "conditioned_synthetic_cohorts.csv"), df_all_conditioned)
@info "  Saved combined conditioned data: $(nrow(df_all_conditioned)) records"

# also save individual cohort files
for (name, df) in all_cohort_results
    fname = "synthetic_$(replace(lowercase(name), " " => "_"))_cohort.csv"
    CSV.write(joinpath(_PATH_TO_DATA, fname), df)
    @info "  Saved $fname"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Evaluation — per-visit feature comparison by cohort
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 6: Evaluation …"

eval_rows = []
for config in cohort_configs
    df_synth = all_cohort_results[config.name]
    sub_ids = complete_subjects[config.indices]

    for v in 1:3
        # real patients from this subpopulation at this visit
        real_mask = (df_clean.Visit .== v) .& [s in sub_ids for s in df_clean.SubjectID]
        real_v = df_clean[real_mask, :]
        synth_v = df_synth[(df_synth.Visit .== v), :]

        if nrow(real_v) == 0
            continue
        end

        # per-feature KS statistics
        ks_vals = Float64[]
        mean_re_vals = Float64[]
        for col in kept_cols
            real_vals = Float64.(collect(skipmissing(real_v[!, col])))
            synth_vals = Float64.(synth_v[!, col])

            if length(real_vals) >= 2
                # KS statistic
                xs, ys = sort(real_vals), sort(synth_vals)
                all_vals = sort(vcat(xs, ys))
                d_max = 0.0
                for val in all_vals
                    f1 = searchsortedlast(xs, val) / length(xs)
                    f2 = searchsortedlast(ys, val) / length(ys)
                    d_max = max(d_max, abs(f1 - f2))
                end
                push!(ks_vals, d_max)

                # mean relative error
                r_mean = mean(real_vals)
                if abs(r_mean) > 1e-12
                    push!(mean_re_vals, abs(mean(synth_vals) - r_mean) / abs(r_mean))
                end
            end
        end

        push!(eval_rows, (
            Cohort = config.name,
            Visit = v,
            N_Real = nrow(real_v),
            N_Synth = nrow(synth_v),
            KS_Mean = round(mean(ks_vals), digits=4),
            KS_Median = round(median(ks_vals), digits=4),
            Mean_RE = round(mean(mean_re_vals), digits=4),
            ρ = round(config.ρ, digits=2),
        ))
    end
end

df_eval = DataFrame(eval_rows)
CSV.write(joinpath(_PATH_TO_DATA, "conditioned_cohort_evaluation.csv"), df_eval)
@info "  Evaluation summary:"
for row in eachrow(df_eval)
    @info "    $(row.Cohort) V$(row.Visit): KS=$(row.KS_Mean), MeanRE=$(row.Mean_RE) ($(row.N_Real) real, $(row.N_Synth) synth)"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: PCA visualization — all three cohorts overlaid with real data
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 7: PCA visualization …"

# fit per-record PCA on real data (all visits, complete rows)
real_rows_mat = []
real_visits_vec = Int[]
real_subject_ids = Int[]
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
        push!(real_rows_mat, vals)
        push!(real_visits_vec, df_clean[row, :Visit])
        push!(real_subject_ids, df_clean[row, :SubjectID])
    end
end
X_real_rec = reduce(hcat, real_rows_mat)

# standardize using real stats
μ_r = vec(mean(X_real_rec, dims=2))
σ_r = vec(std(X_real_rec, dims=2))
for i in eachindex(σ_r)
    σ_r[i] < 1e-12 && (σ_r[i] = 1.0)
end
Z_real_rec = (X_real_rec .- μ_r) ./ σ_r
pca_rec = MultivariateStats.fit(PCA, Z_real_rec; maxoutdim=2)
P_real = MultivariateStats.transform(pca_rec, Z_real_rec)

var1 = round(100 * MultivariateStats.principalvars(pca_rec)[1] / MultivariateStats.tvar(pca_rec), digits=1)
var2 = round(100 * MultivariateStats.principalvars(pca_rec)[2] / MultivariateStats.tvar(pca_rec), digits=1)

# colors
cohort_colors = Dict(
    "Uncomplicated" => RGB(0.20, 0.60, 0.40),
    "PCOS"          => RGB(0.80, 0.30, 0.30),
    "Developed PE"  => RGB(0.30, 0.30, 0.80),
)

# one panel per cohort
panels = []
for config in cohort_configs
    sub_ids = Set(complete_subjects[config.indices])

    p = plot(
        xlabel="PC1 ($(var1)%)", ylabel="PC2 ($(var2)%)",
        title="$(config.name) Cohort (ρ=$(round(config.ρ, digits=1)))",
        titlefontsize=10,
        background_color=:white, grid=false, framestyle=:box,
        legendfontsize=7, margin=3Plots.mm,
    )

    # plot ALL real patients in gray
    scatter!(p, P_real[1, :], P_real[2, :],
        label="Real (all)", color=:gray80, alpha=0.4,
        markersize=4, markerstrokewidth=0)

    # highlight real patients from this subpopulation
    real_sub_indices = findall(sid -> sid in sub_ids, real_subject_ids)

    if !isempty(real_sub_indices)
        scatter!(p, P_real[1, real_sub_indices], P_real[2, real_sub_indices],
            label="Real $(config.name)", color=cohort_colors[config.name],
            alpha=0.9, markersize=7, markerstrokecolor=:white, markerstrokewidth=1)
    end

    # transform synthetic cohort
    df_synth = all_cohort_results[config.name]
    X_synth_rec = Matrix{Float64}(undef, n_assays, nrow(df_synth))
    for i in 1:nrow(df_synth)
        for (j, col) in enumerate(kept_cols)
            X_synth_rec[j, i] = df_synth[i, col]
        end
    end
    Z_synth_rec = (X_synth_rec .- μ_r) ./ σ_r
    P_synth = MultivariateStats.transform(pca_rec, Z_synth_rec)

    scatter!(p, P_synth[1, :], P_synth[2, :],
        label="Synth $(config.name)", color=cohort_colors[config.name],
        alpha=0.3, markersize=3, markerstrokewidth=0, markershape=:diamond)

    push!(panels, p)
end

p_pca = plot(panels..., layout=(1, 3), size=(1500, 450),
    plot_title="Conditioned SA Generation: PCA by Cohort",
    plot_titlefontsize=14,
    background_color=:white, margin=5Plots.mm)

savefig(p_pca, joinpath(_PATH_TO_FIG, "conditioned_cohorts_pca.pdf"))
savefig(p_pca, joinpath(_PATH_TO_FIG, "conditioned_cohorts_pca.png"))
@info "  Saved PCA figure"

# ══════════════════════════════════════════════════════════════════════════════
# Step 8: Key feature trajectory comparison across cohorts
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 8: Trajectory comparison …"

# select key features
trajectory_features = Symbol[]
for c in kept_cols
    s = string(c)
    if s in ["Estradiol", "VIII", "Fbgn", "AT", "PC"] ||
       occursin("TF Initiator Peak", s) ||
       occursin("TF Initiator ETP", s) ||
       occursin("TF Only MCF", s)
        push!(trajectory_features, c)
    end
end
trajectory_features = unique(trajectory_features)
@info "  Trajectory features: $(length(trajectory_features))"

traj_panels = []
for col in trajectory_features
    p = plot(
        xlabel="Visit", ylabel=string(col),
        title=string(col), titlefontsize=8,
        xticks=([1,2,3], ["Baseline","1st Tri","3rd Tri"]),
        xrotation=15, legend=:topright, legendfontsize=6,
        background_color=:white, grid=false, framestyle=:box,
    )

    for config in cohort_configs
        df_synth = all_cohort_results[config.name]
        means = Float64[]
        sds = Float64[]
        for v in 1:3
            vals = df_synth[df_synth.Visit .== v, col]
            push!(means, mean(vals))
            push!(sds, std(vals))
        end
        plot!(p, [1,2,3], means, ribbon=sds, fillalpha=0.15,
            label=config.name, color=cohort_colors[config.name],
            lw=2, marker=:circle, markersize=4)
    end

    push!(traj_panels, p)
end

n_rows = ceil(Int, length(traj_panels) / 3)
p_traj = plot(traj_panels..., layout=(n_rows, min(3, length(traj_panels))),
    size=(1200, 300*n_rows),
    plot_title="Longitudinal Trajectories by Conditioned Cohort",
    plot_titlefontsize=12,
    background_color=:white, margin=5Plots.mm)

savefig(p_traj, joinpath(_PATH_TO_FIG, "conditioned_cohorts_trajectories.pdf"))
savefig(p_traj, joinpath(_PATH_TO_FIG, "conditioned_cohorts_trajectories.png"))
@info "  Saved trajectory figure"

# ══════════════════════════════════════════════════════════════════════════════
# Step 9: Save pipeline state
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 9: Saving …"
JLD2.@save joinpath(_PATH_TO_DATA, "conditioned_generation_state.jld2") cohort_configs all_cohort_results df_eval

@info "\n✓ Conditioned generation complete!"
@info "  Cohorts: Uncomplicated ($(length(uncomplicated_indices)) → $n_synth), PCOS ($(length(pcos_indices)) → $n_synth), Developed PE ($(length(pe_indices)) → $n_synth)"
@info "  Output: data/conditioned_synthetic_cohorts.csv"
@info "  Figures: figs/conditioned_cohorts_pca.pdf, figs/conditioned_cohorts_trajectories.pdf"

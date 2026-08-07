# ======================================================================================== #
# validate_conditioned_generation.jl — Level 3: Conditional Structure
#
# Generate condition-specific synthetic cohorts using multiplicity-weighted SA:
#   Uncomplicated (18 real → 100 synthetic)  [no PE outcome]
#   PCOS          (3 real  → 100 synthetic)  [cross-cutting comorbidity]
#   Developed PE  (5 real  → 100 synthetic)  [developed PE during study]
#
# Validation checks:
#   1. Condition-specific means match real condition-specific means
#   2. PCA separation: do generated cohorts separate as expected?
#   3. Cross-cohort feature differences preserve real between-group differences
#   4. MVN fundamentally cannot do this from n=3 PCOS patients
# ======================================================================================== #

include(joinpath(@__DIR__, "..", "Include.jl"))

# ── Load pipeline state ───────────────────────────────────────────────────────
@info "Step 1: Loading pipeline state …"
JLD2.@load joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2") X̂ pca_concat std_params_concat pca_norms complete_subjects subject_conditions kept_cols n_assays df_clean

d_pca, K = size(X̂)
n_visits = 3
d_concat = n_assays * n_visits
μ_concat = std_params_concat.μ
σ_concat = std_params_concat.σ

@info "  $d_pca PCA dims × $K patients, $n_assays features/visit"
@info "  Conditions: $(sort(collect(StatsBase.countmap(subject_conditions))))"

# ── Identify subpopulations ───────────────────────────────────────────────────
@info "\nStep 2: Identifying subpopulations …"

pe_indices = Int[]
for (i, s) in enumerate(complete_subjects)
    row = findfirst((df_clean.SubjectID .== s))
    if df_clean[row, :DevelopedPE] == "Yes"
        push!(pe_indices, i)
    end
end

# Uncomplicated: all patients who did NOT develop PE (n=18)
uncomplicated_indices = findall(i -> !(i in pe_indices), 1:K)

# PCOS: cross-cutting comorbidity (n=3, 1 overlaps with PE)
pcos_indices = findall(c -> c == "PCOS", subject_conditions)

@info "  Uncomplicated: $(length(uncomplicated_indices)), PCOS: $(length(pcos_indices)), Developed PE: $(length(pe_indices))"

# ── Fisher separation index ──────────────────────────────────────────────────
@info "\nStep 3: Fisher separation analysis …"
S_uncomplicated = fisher_separation_index(X̂, uncomplicated_indices)
S_pcos = fisher_separation_index(X̂, pcos_indices)
S_pe = fisher_separation_index(X̂, pe_indices)

@info "  Fisher separation: Uncomplicated=$(round(S_uncomplicated, digits=4)), PCOS=$(round(S_pcos, digits=4)), Developed PE=$(round(S_pe, digits=4))"

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Generate conditioned cohorts
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 4: Generating conditioned cohorts …"

n_synth_per_cohort = 100
f_target = 0.80
T_langevin = 5000

cohorts = [
    (name="Uncomplicated", indices=uncomplicated_indices),
    (name="PCOS",          indices=pcos_indices),
    (name="Developed PE",  indices=pe_indices),
]

# store results
cohort_data = Dict{String, Matrix{Float64}}()   # name → n_synth × d_concat
cohort_pca  = Dict{String, Matrix{Float64}}()    # name → n_synth × d_pca

for cohort in cohorts
    name = cohort.name
    indices = cohort.indices
    K_sub = length(indices)
    K_bg = K - K_sub

    @info "\n  ── $name ($(K_sub) real → $n_synth_per_cohort synthetic) ──"

    # compute multiplicity ratio ρ
    ρ = ρ_for_target_fraction(K_sub, K_bg, f_target)
    r = multiplicity_vector(K, indices; ρ=ρ)
    eff_frac = effective_subset_fraction(r, indices)
    @info "    ρ = $(round(ρ, digits=2)), component-probability fraction = $(round(eff_frac, digits=3))"

    # select β*(ρ) at the weighted attention-entropy transition
    β_star, _ = find_weighted_entropy_inflection(X̂, r; β_range=(0.1, 50.0), n_betas=200)
    @info "    β*(ρ) = $(round(β_star, digits=2))"

    # generate samples
    pca_samples = Matrix{Float64}(undef, n_synth_per_cohort, d_pca)
    for i in 1:n_synth_per_cohort
        ξ₀ = randn(d_pca)
        ξ₀ ./= norm(ξ₀)
        result = weighted_sample(X̂, ξ₀, T_langevin, r; β=β_star, α=0.01, seed=1000*i)
        pca_samples[i, :] = result.Ξ[end, :]
    end
    cohort_pca[name] = pca_samples

    # rescale: direction × empirical norm
    synth_concat = Matrix{Float64}(undef, n_synth_per_cohort, d_concat)
    for i in 1:n_synth_per_cohort
        direction = pca_samples[i, :]
        direction ./= max(norm(direction), 1e-12)

        # draw norm from empirical distribution
        r_idx = rand(1:length(pca_norms))
        magnitude = pca_norms[r_idx]

        z_pca = magnitude .* direction
        z_concat = MultivariateStats.reconstruct(pca_concat, z_pca)
        x_concat = z_concat .* σ_concat .+ μ_concat
        synth_concat[i, :] = x_concat
    end
    cohort_data[name] = synth_concat

    @info "    Generated $n_synth_per_cohort synthetic patients"
end

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Build DataFrames for analysis
# ══════════════════════════════════════════════════════════════════════════════
@info "\nStep 5: Building analysis DataFrames …"

function concat_to_visits(synth_concat, kept_cols, n_assays, cohort_name)
    n = size(synth_concat, 1)
    dfs = DataFrame[]
    for v in 1:3
        offset = (v - 1) * n_assays
        cols_data = synth_concat[:, (offset+1):(offset+n_assays)]
        df_v = DataFrame(cols_data, kept_cols)
        df_v.SyntheticID = 1:n
        df_v.Visit = fill(v, n)
        df_v.Cohort = fill(cohort_name, n)
        push!(dfs, df_v)
    end
    return vcat(dfs...)
end

df_cohorts = vcat([concat_to_visits(cohort_data[c.name], kept_cols, n_assays, c.name) for c in cohorts]...)
@info "  Combined cohort DataFrame: $(nrow(df_cohorts)) records"

# Build real per-condition DataFrames for comparison
function get_real_condition(df_clean, complete_subjects, subject_conditions, cond_name, pe_check=false)
    rows = DataFrame[]
    for (i, s) in enumerate(complete_subjects)
        if pe_check
            row_idx = findfirst(df_clean.SubjectID .== s)
            if df_clean[row_idx, :DevelopedPE] != "Yes"
                continue
            end
        elseif subject_conditions[i] != cond_name
            continue
        end
        push!(rows, df_clean[df_clean.SubjectID .== s, :])
    end
    return isempty(rows) ? DataFrame() : vcat(rows...)
end

# Uncomplicated: all patients who did NOT develop PE
function get_real_uncomplicated(df_clean, complete_subjects)
    rows = DataFrame[]
    for s in complete_subjects
        row_idx = findfirst(df_clean.SubjectID .== s)
        if df_clean[row_idx, :DevelopedPE] != "Yes"
            push!(rows, df_clean[df_clean.SubjectID .== s, :])
        end
    end
    return isempty(rows) ? DataFrame() : vcat(rows...)
end

df_real_uncomplicated = get_real_uncomplicated(df_clean, complete_subjects)
df_real_pcos    = get_real_condition(df_clean, complete_subjects, subject_conditions, "PCOS")
df_real_pe      = get_real_condition(df_clean, complete_subjects, subject_conditions, "", true)

# ══════════════════════════════════════════════════════════════════════════════
# Validation 1: Condition-specific means
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Validation 1: Condition-Specific Means ═══"

key_features = [:II, :VIII, :AT, :Fbgn, :X, :PC,
    Symbol("TF Initiator Peak (nM)"),
    Symbol("TF Initiator ETP (nM·min)"),
    Symbol("TF Initiator Lagtime (min)")]

for (cond_name, df_real_cond, synth_name) in [
    ("Uncomplicated", df_real_uncomplicated, "Uncomplicated"),
    ("PCOS", df_real_pcos, "PCOS"),
    ("Developed PE", df_real_pe, "Developed PE"),
]
    @info "\n  ── $cond_name ──"
    df_synth_cond = filter(r -> r.Cohort == synth_name, df_cohorts)

    for col in key_features
        if !(col in propertynames(df_synth_cond)) || !(col in propertynames(df_real_cond))
            continue
        end

        r_vals = collect(skipmissing(df_real_cond[:, col]))
        s_vals = collect(skipmissing(df_synth_cond[:, col]))

        if isempty(r_vals) || isempty(s_vals)
            continue
        end

        r_mean = mean(r_vals)
        s_mean = mean(s_vals)
        re = abs(r_mean) > 1e-12 ? abs(s_mean - r_mean) / abs(r_mean) : NaN

        @info "    $(rpad(string(col), 35)) Real=$(round(r_mean, digits=1))  Synth=$(round(s_mean, digits=1))  RE=$(round(re, digits=3))"
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Validation 2: Between-group differences preserved
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Validation 2: Between-Group Differences ═══"
@info "  (Difference in means: Condition - Healthy)"

compare_features = [:VIII, :AT, :II, :Fbgn,
    Symbol("TF Initiator Peak (nM)"),
    Symbol("TF Initiator ETP (nM·min)")]

for (cond_name, df_real_cond, synth_name) in [
    ("PCOS", df_real_pcos, "PCOS"),
    ("Developed PE", df_real_pe, "Developed PE"),
]
    @info "\n  ── $cond_name vs Uncomplicated ──"
    df_synth_cond = filter(r -> r.Cohort == synth_name, df_cohorts)
    df_synth_healthy = filter(r -> r.Cohort == "Uncomplicated", df_cohorts)

    for col in compare_features
        if !(col in propertynames(df_real_cond)) || !(col in propertynames(df_real_uncomplicated))
            continue
        end

        r_cond = collect(skipmissing(df_real_cond[:, col]))
        r_healthy = collect(skipmissing(df_real_uncomplicated[:, col]))
        s_cond = collect(skipmissing(df_synth_cond[:, col]))
        s_healthy = collect(skipmissing(df_synth_healthy[:, col]))

        if isempty(r_cond) || isempty(r_healthy) || isempty(s_cond) || isempty(s_healthy)
            continue
        end

        real_diff = mean(r_cond) - mean(r_healthy)
        synth_diff = mean(s_cond) - mean(s_healthy)

        # same sign = direction preserved
        sign_match = sign(real_diff) == sign(synth_diff)
        marker = sign_match ? "✓" : "✗"

        @info "    $(rpad(string(col), 35)) Real Δ=$(round(real_diff, digits=1))  Synth Δ=$(round(synth_diff, digits=1))  $marker"
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Validation 3: PCA visualization — do cohorts separate?
# ══════════════════════════════════════════════════════════════════════════════
@info "\n═══ Validation 3: PCA Separation ═══"

# Project everything into PCA space (first 2 components for visualization)
# Real patients
Z_real = Matrix{Float64}(undef, K, d_concat)
for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_clean.SubjectID .== s) .& (df_clean.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            Z_real[i, offset + j] = df_clean[row, col]
        end
    end
end

# standardize and project
Z_real_std = (Z_real .- μ_concat') ./ σ_concat'
pca_real = MultivariateStats.transform(pca_concat, Z_real_std')  # d_pca × K

p_pca = plot(; xlabel="PC1", ylabel="PC2",
             title="Level 3: Conditioned Generation — PCA",
             legend=:topright, size=(700, 600))

# Real patients by condition
colors_real = Dict("Uncomplicated" => :steelblue, "PCOS" => :orange, "Developed PE" => :red)
for (cond, color) in colors_real
    if cond == "Developed PE"
        idxs = pe_indices
    elseif cond == "Uncomplicated"
        idxs = uncomplicated_indices
    else
        idxs = findall(c -> c == cond, subject_conditions)
    end
    if !isempty(idxs)
        scatter!(p_pca, pca_real[1, idxs], pca_real[2, idxs];
                 label="Real $cond", mc=color, ms=8, ma=0.9, shape=:circle,
                 markerstrokewidth=2)
    end
end

# Synthetic cohorts
colors_synth = Dict("Uncomplicated" => :lightblue, "PCOS" => :lightsalmon, "Developed PE" => :pink)
shapes_synth = Dict("Uncomplicated" => :diamond, "PCOS" => :utriangle, "Developed PE" => :star5)
for cohort in cohorts
    name = cohort.name
    pca_s = cohort_pca[name]
    # rescale for plotting (same as real)
    pca_rescaled = Matrix{Float64}(undef, n_synth_per_cohort, d_pca)
    for i in 1:n_synth_per_cohort
        direction = pca_s[i, :] ./ max(norm(pca_s[i, :]), 1e-12)
        r_idx = rand(1:length(pca_norms))
        pca_rescaled[i, :] = pca_norms[r_idx] .* direction
    end

    scatter!(p_pca, pca_rescaled[:, 1], pca_rescaled[:, 2];
             label="Synth $name", mc=colors_synth[name], ms=4, ma=0.4,
             shape=shapes_synth[name], markerstrokewidth=0)
end

savefig(p_pca, joinpath(_PATH_TO_FIG, "validate_conditioned_pca.pdf"))
savefig(p_pca, joinpath(_PATH_TO_FIG, "validate_conditioned_pca.png"))
@info "  Saved validate_conditioned_pca.pdf"

# ── Per-visit PCA ──────────────────────────────────────────────────────────────
pca_visit_panels = []
for v in 1:3
    p = plot(; xlabel="PC1", ylabel="PC2", title="Visit $v",
             legend=(v==1 ? :topright : false), size=(400, 400))

    # real by condition
    for (cond, color) in colors_real
        if cond == "Developed PE"
            idxs = pe_indices
        elseif cond == "Uncomplicated"
            idxs = uncomplicated_indices
        else
            idxs = findall(c -> c == cond, subject_conditions)
        end
        if !isempty(idxs)
            scatter!(p, pca_real[1, idxs], pca_real[2, idxs];
                     label=(v==1 ? "Real $cond" : ""), mc=color, ms=7, ma=0.9)
        end
    end

    # synthetic by cohort (project per-visit data)
    for cohort in cohorts
        name = cohort.name
        synth = cohort_data[name]
        offset = (v - 1) * n_assays
        visit_data = synth[:, (offset+1):(offset+n_assays)]

        # need full concat for PCA projection — use the visit data padded with zeros
        # Actually just project the full concat vector
        synth_std = (synth .- μ_concat') ./ σ_concat'
        pca_proj = MultivariateStats.transform(pca_concat, synth_std')

        scatter!(p, pca_proj[1, :], pca_proj[2, :];
                 label=(v==1 ? "Synth $name" : ""), mc=colors_synth[name],
                 ms=3, ma=0.3, shape=shapes_synth[name], markerstrokewidth=0)
    end

    push!(pca_visit_panels, p)
end

p_pca_visits = plot(pca_visit_panels...; layout=(1, 3), size=(1200, 400),
                    plot_title="Level 3: Conditioned Cohorts by Visit")
savefig(p_pca_visits, joinpath(_PATH_TO_FIG, "validate_conditioned_pca_by_visit.pdf"))
savefig(p_pca_visits, joinpath(_PATH_TO_FIG, "validate_conditioned_pca_by_visit.png"))
@info "  Saved validate_conditioned_pca_by_visit.pdf"

# ── Save cohort data ──────────────────────────────────────────────────────────
CSV.write(joinpath(_PATH_TO_DATA, "conditioned_cohorts_all.csv"), df_cohorts)
for cohort in cohorts
    df_c = filter(r -> r.Cohort == cohort.name, df_cohorts)
    CSV.write(joinpath(_PATH_TO_DATA, "conditioned_$(replace(lowercase(cohort.name), " " => "_")).csv"), df_c)
end
@info "\nSaved cohort data to CSV"

@info "\nDone! Level 3 validation complete."

# ──────────────────────────────────────────────────────────────────────────────
# Patient.jl — Data pipeline for patient biochemical records
#
# Provides: Excel loading, cleaning, standardization, PCA pipeline,
# memory matrix construction, inverse decoding, and evaluation metrics.
# ──────────────────────────────────────────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════
# Data I/O: load and clean the legacy biochemical dataset
# ══════════════════════════════════════════════════════════════════════════════

# Column indices in the Excel file (1-indexed, row 3 has headers)
const _HEADER_ROW = 3
const _DATA_START_ROW = 4
const _COL_SUBJECT_ID = 2
const _COL_VISIT = 3
const _COL_DRAW_TAG = 4
const _COL_CONDITION = 5
const _COL_PCOS = 6
const _COL_PE = 7
const _COL_HTN = 8
const _COL_ASSAY_START = 10  # first biochemical assay column (Estradiol)
const _COL_ASSAY_END = 104   # last column (Viscoelastic AUC TF+tPA)
const _SKIP_COLS = Set([56,57,58,59,60,61,62,63,  # Y/N flags (TGA and ROTEM)
                        89, 97])                    # Tag strings (Viscoelastic)

"""
    load_legacy_data(filepath) -> DataFrame

Load the legacy biochemical dataset from the Excel file.
Returns a DataFrame with subject metadata and continuous assay values.
Non-numeric assay values (e.g., '#N/A', empty) are replaced with `missing`.
"""
function load_legacy_data(filepath::String)
    xf = XLSX.readxlsx(filepath)
    ws = xf["Assay Values"]

    # ── Build unique header names ──────────────────────────────────────────
    # Row 1 has group type (Hormones, Coagulation, TGA Parameters, …)
    # Row 2 has sub-group / initiator (TF Initiator, TF Only, TF + tPA, …)
    # Row 3 has the measurement name (Lagtime, Peak, MCF, …)
    # Many measurement names repeat across groups, so we combine them.
    headers = String[]
    assay_cols = Int[]  # actual Excel column indices we keep
    current_type = ""
    current_sub = ""
    for col in _COL_ASSAY_START:_COL_ASSAY_END
        # always update group/sub-group tracking, even for skipped columns
        r1 = ws[1, col]
        r2 = ws[2, col]
        if r1 !== nothing && !ismissing(r1)
            current_type = strip(string(r1))
        end
        if r2 !== nothing && !ismissing(r2)
            current_sub = strip(string(r2))
        end

        col in _SKIP_COLS && continue  # skip Y/N flags and Tag columns

        r3 = ws[3, col]
        meas_name = (r3 === nothing || ismissing(r3)) ? "Col$col" : strip(string(r3))

        # build unique name: for the first 46 assay columns (10-55) the
        # measurement name is already unique; for TGA/Viscoelastic we
        # prefix with sub-group to disambiguate repeated names
        if col <= 55
            header = meas_name
        else
            header = "$(current_sub) $(meas_name)"
        end

        # ensure uniqueness
        if header in headers
            header = "$(header)_$(col)"
        end

        push!(headers, header)
        push!(assay_cols, col)
    end

    n_assays = length(headers)
    records = Dict{String, Vector}()
    records["SubjectID"] = Int[]
    records["Visit"] = Int[]
    records["DrawTag"] = String[]
    records["Condition"] = String[]
    records["PCOS"] = String[]
    records["DevelopedPE"] = String[]
    records["DevelopedHTN"] = String[]
    for h in headers
        records[h] = Union{Float64, Missing}[]
    end

    # track subject ID (only filled in first visit row)
    current_subject = 0
    for row in _DATA_START_ROW:ws.dimension.stop.row_number
        sid_raw = ws[row, _COL_SUBJECT_ID]
        if sid_raw !== nothing && sid_raw !== ""
            current_subject = Int(sid_raw)
        end

        visit_raw = ws[row, _COL_VISIT]
        visit_raw === nothing && continue

        push!(records["SubjectID"], current_subject)
        push!(records["Visit"], Int(visit_raw))

        draw_raw = ws[row, _COL_DRAW_TAG]
        push!(records["DrawTag"], draw_raw === nothing ? "" : string(draw_raw))

        cond_raw = ws[row, _COL_CONDITION]
        push!(records["Condition"], cond_raw === nothing ? "" : string(cond_raw))

        pcos_raw = ws[row, _COL_PCOS]
        push!(records["PCOS"], pcos_raw === nothing ? "" : string(pcos_raw))

        pe_raw = ws[row, _COL_PE]
        push!(records["DevelopedPE"], pe_raw === nothing ? "" : string(pe_raw))

        htn_raw = ws[row, _COL_HTN]
        push!(records["DevelopedHTN"], htn_raw === nothing ? "" : string(htn_raw))

        for (j, col) in enumerate(assay_cols)
            val = ws[row, col]
            if val === nothing || val === "" || ismissing(val) ||
               (val isa String && (val == "#N/A" || val == "---" || val == "No Activity"))
                push!(records[headers[j]], missing)
            elseif val isa Number
                push!(records[headers[j]], Float64(val))
            else
                parsed = tryparse(Float64, string(val))
                push!(records[headers[j]], parsed === nothing ? missing : parsed)
            end
        end
    end

    # build DataFrame with column order preserved
    df = DataFrame()
    df.SubjectID = records["SubjectID"]
    df.Visit = records["Visit"]
    df.DrawTag = records["DrawTag"]
    df.Condition = records["Condition"]
    df.PCOS = records["PCOS"]
    df.DevelopedPE = records["DevelopedPE"]
    df.DevelopedHTN = records["DevelopedHTN"]
    for h in headers
        df[!, Symbol(h)] = records[h]
    end

    # propagate Condition down within each subject
    _propagate_metadata!(df)

    return df
end

"""
    _propagate_metadata!(df)

Fill in empty Condition/PCOS/PE/HTN values by propagating from earlier visits
of the same subject.
"""
function _propagate_metadata!(df::DataFrame)
    for col in [:Condition, :PCOS, :DevelopedPE, :DevelopedHTN]
        for subj in unique(df.SubjectID)
            mask = df.SubjectID .== subj
            vals = df[mask, col]
            fill_val = ""
            for v in vals
                if v != ""
                    fill_val = v
                    break
                end
            end
            if fill_val != ""
                for i in findall(mask)
                    if df[i, col] == ""
                        df[i, col] = fill_val
                    end
                end
            end
        end
    end
end


# ══════════════════════════════════════════════════════════════════════════════
# Data cleaning and standardization
# ══════════════════════════════════════════════════════════════════════════════

"""
    assay_columns(df) -> Vector{Symbol}

Return the names of columns that contain continuous assay measurements.
"""
function assay_columns(df::DataFrame)
    meta_cols = Set(["SubjectID", "Visit", "DrawTag", "Condition", "PCOS", "DevelopedPE", "DevelopedHTN"])
    return [Symbol(c) for c in names(df) if c ∉ meta_cols]
end

"""
    clean_assay_data(df; max_missing_frac=0.3) -> (df_clean, kept_cols)

Remove assay columns with too many missing values, then remove rows with
any remaining missing values. Returns cleaned DataFrame and the list of
retained assay column names.
"""
function clean_assay_data(df::DataFrame; max_missing_frac::Float64=0.3)
    acols = assay_columns(df)
    n_rows = nrow(df)

    # filter columns by missing fraction
    kept_cols = Symbol[]
    for c in acols
        n_miss = count(ismissing, df[!, c])
        if n_miss / n_rows <= max_missing_frac
            push!(kept_cols, c)
        end
    end
    @info "  Columns: $(length(acols)) → $(length(kept_cols)) (removed $(length(acols) - length(kept_cols)) with >$(round(100*max_missing_frac))% missing)"

    # select metadata + kept assay columns
    meta = [:SubjectID, :Visit, :DrawTag, :Condition, :PCOS, :DevelopedPE, :DevelopedHTN]
    df_sub = df[!, vcat(meta, kept_cols)]

    # drop rows with any remaining missing in assay columns
    complete_mask = completecases(df_sub[!, kept_cols])
    df_clean = df_sub[complete_mask, :]
    @info "  Rows: $(nrow(df)) → $(nrow(df_clean)) (removed $(nrow(df) - nrow(df_clean)) incomplete rows)"

    return df_clean, kept_cols
end

"""
    StandardizationParams

Stores the mean and standard deviation for each feature, used to standardize
and de-standardize data.
"""
struct StandardizationParams
    μ::Vector{Float64}
    σ::Vector{Float64}
    col_names::Vector{Symbol}
end

"""
    standardize(df, cols) -> (Z, params)

Z-score standardize the assay columns. Returns the standardized matrix (n × p)
and the standardization parameters for inverse transform.
"""
function standardize(df::DataFrame, cols::Vector{Symbol})
    X = Matrix{Float64}(df[!, cols])
    μ = vec(mean(X, dims=1))
    σ = vec(std(X, dims=1))

    # guard against zero-variance columns
    for i in eachindex(σ)
        σ[i] < 1e-12 && (σ[i] = 1.0)
    end

    Z = (X .- μ') ./ σ'
    params = StandardizationParams(μ, σ, cols)
    return Z, params
end

"""
    destandardize(Z, params) -> Matrix{Float64}

Inverse z-score transform: X = Z * σ + μ
"""
function destandardize(Z::AbstractMatrix{Float64}, params::StandardizationParams)
    return Z .* params.σ' .+ params.μ'
end

function destandardize(z::Vector{Float64}, params::StandardizationParams)
    return z .* params.σ .+ params.μ
end


# ══════════════════════════════════════════════════════════════════════════════
# PCA pipeline: standardize → reduce → normalize → memory matrix
# ══════════════════════════════════════════════════════════════════════════════

"""
    build_memory_matrix(Z; pratio=0.95) -> (X̂, pca_model, d_orig)

Full pipeline: standardized data → PCA → unit-norm columns.

# Arguments
- `Z::Matrix{Float64}`: Standardized data matrix (n × p), where n = number of
  patient records and p = number of assays.

# Keyword Arguments
- `pratio::Float64=0.95`: Fraction of variance to retain in PCA.

# Returns
- `X̂`: Memory matrix (d_pca × K), unit-norm columns in PCA space
- `pca_model`: Fitted PCA model (for inverse transform)
- `d_orig`: Original number of features (p)
"""
function build_memory_matrix(Z::Matrix{Float64}; pratio::Float64=0.95)
    n, p = size(Z)

    # MultivariateStats expects features × samples, so transpose
    Z_t = Z'  # p × n

    # PCA
    pca_model = MultivariateStats.fit(PCA, Z_t; pratio=pratio)
    d_pca = MultivariateStats.outdim(pca_model)
    X_pca = MultivariateStats.transform(pca_model, Z_t)  # d_pca × n
    var_retained = round(100 * sum(MultivariateStats.principalvars(pca_model)) /
                         MultivariateStats.tvar(pca_model), digits=1)
    @info "  PCA: $p → $d_pca dimensions ($var_retained% variance retained)"

    # normalize columns to unit norm
    X̂ = copy(X_pca)
    K = size(X̂, 2)
    for k in 1:K
        nk = norm(X̂[:, k])
        X̂[:, k] ./= (nk + 1e-12)
    end
    @info "  Memory matrix X̂: $d_pca × $K (unit-norm columns)"

    return X̂, pca_model, p
end


"""
    build_memory_matrix_raw(Z; pratio=0.95) -> (X̂, pca_model, d_orig, norms)

Alternative pipeline: standardized data → PCA → **no unit-norm projection**.

For continuous patient data, the PCA projection norms carry real information
(how far a patient is from the population mean). Unit-norming erases this,
producing isotropic samples that lack the anisotropic spread of the real data.

Returns the raw PCA projections plus their norms (for diagnostics).
"""
function build_memory_matrix_raw(Z::Matrix{Float64}; pratio::Float64=0.95)
    n, p = size(Z)
    Z_t = Z'  # p × n

    pca_model = MultivariateStats.fit(PCA, Z_t; pratio=pratio)
    d_pca = MultivariateStats.outdim(pca_model)
    X_pca = MultivariateStats.transform(pca_model, Z_t)  # d_pca × n
    var_retained = round(100 * sum(MultivariateStats.principalvars(pca_model)) /
                         MultivariateStats.tvar(pca_model), digits=1)
    @info "  PCA: $p → $d_pca dimensions ($var_retained% variance retained)"

    K = size(X_pca, 2)
    norms = [norm(X_pca[:, k]) for k in 1:K]
    @info "  Memory matrix (raw): $d_pca × $K — norm range: [$(round(minimum(norms), digits=2)), $(round(maximum(norms), digits=2))], mean=$(round(mean(norms), digits=2))"

    return X_pca, pca_model, p, norms
end


"""
    decode_sample(ξ_pca, pca_model, std_params) -> Vector{Float64}

Decode a PCA-space vector back to the original biochemical measurement space.
Maps through inverse PCA, then de-standardization.
"""
function decode_sample(ξ_pca::Vector{Float64}, pca_model, std_params::StandardizationParams)
    z_standardized = vec(MultivariateStats.reconstruct(pca_model, ξ_pca))
    return destandardize(z_standardized, std_params)
end


"""
    generate_patients(X̂, pca_model, std_params, n_samples, T;
                      β, α=0.1, seed=nothing) -> DataFrame

Generate `n_samples` synthetic patient records using the SA sampler.

Each sample starts from a random initial state drawn uniformly from the memory
patterns (with small noise), runs `T` Langevin steps, and the final state is
decoded back to biochemical measurement space.

Returns a DataFrame with the original assay column names.
"""
function generate_patients(X̂::Matrix{Float64}, pca_model, std_params::StandardizationParams,
                           n_samples::Int, T::Int;
                           β::Float64, α::Float64=0.1,
                           seed::Union{Int, Nothing}=nothing)

    if seed !== nothing
        Random.seed!(seed)
    end

    d, K = size(X̂)
    results = Matrix{Float64}(undef, n_samples, length(std_params.col_names))

    for i in 1:n_samples
        # initialize from a random stored pattern + small perturbation
        k0 = rand(1:K)
        ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d)
        ξ₀ ./= norm(ξ₀)  # re-normalize

        # run SA sampler
        result = sample(X̂, ξ₀, T; β=β, α=α)

        # take the final state
        ξ_final = result.Ξ[end, :]

        # decode back to measurement space
        measurements = decode_sample(ξ_final, pca_model, std_params)
        results[i, :] = measurements
    end

    # build DataFrame
    df_synth = DataFrame(results, std_params.col_names)
    return df_synth
end


"""
    generate_conditioned_patients(X̂, pca_model, std_params, condition_mask,
                                  n_samples, T; β, α=0.1, ρ_ratio=10.0,
                                  seed=nothing) -> DataFrame

Generate synthetic patients conditioned on a clinical subgroup using
multiplicity-weighted SA. Patterns matching `condition_mask` receive
multiplicity weight `ρ_ratio`, others get weight 1.0.
"""
function generate_conditioned_patients(X̂::Matrix{Float64}, pca_model,
                                       std_params::StandardizationParams,
                                       condition_mask::BitVector,
                                       n_samples::Int, T::Int;
                                       β::Float64, α::Float64=0.1,
                                       ρ_ratio::Float64=10.0,
                                       seed::Union{Int, Nothing}=nothing)

    if seed !== nothing
        Random.seed!(seed)
    end

    d, K = size(X̂)
    length(condition_mask) == K || throw(DimensionMismatch(
        "condition_mask length $(length(condition_mask)) ≠ K=$K"))

    # build multiplicity vector
    ρ = ones(K)
    for k in 1:K
        if condition_mask[k]
            ρ[k] = ρ_ratio
        end
    end

    results = Matrix{Float64}(undef, n_samples, length(std_params.col_names))

    for i in 1:n_samples
        # initialize from a random conditioned pattern
        cond_indices = findall(condition_mask)
        k0 = isempty(cond_indices) ? rand(1:K) : rand(cond_indices)
        ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d)
        ξ₀ ./= norm(ξ₀)

        result = weighted_sample(X̂, ξ₀, T, ρ; β=β, α=α)
        ξ_final = result.Ξ[end, :]
        measurements = decode_sample(ξ_final, pca_model, std_params)
        results[i, :] = measurements
    end

    df_synth = DataFrame(results, std_params.col_names)
    return df_synth
end


"""
    generate_patients_rescaled(X̂_unit, pca_model, std_params, pca_norms,
                               n_samples, T; β, α=0.1, seed=nothing) -> DataFrame

Generate synthetic patients using the **direction–magnitude decomposition**:

1. SA samples **directions** on the unit sphere (unit-norm memory, well-behaved
   Hopfield energy, good theory for β*).
2. Each sample is then **rescaled** by a norm drawn from the empirical
   distribution of PCA projection norms from the real data.

This preserves the anisotropic variance structure that unit-norm SA alone erases:
patients far from the population mean (large PCA norm) remain far out, while
patients near the mean stay near the center.

# Arguments
- `X̂_unit`: Unit-norm memory matrix (d × K).
- `pca_model`: Fitted PCA model.
- `std_params`: Standardization parameters.
- `pca_norms`: Vector of PCA projection norms from the real data (length K).
- `n_samples`, `T`: Number of samples and Langevin steps.
"""
function generate_patients_rescaled(X̂_unit::Matrix{Float64}, pca_model,
                                     std_params::StandardizationParams,
                                     pca_norms::Vector{Float64},
                                     n_samples::Int, T::Int;
                                     β::Float64, α::Float64=0.1,
                                     seed::Union{Int, Nothing}=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    d, K = size(X̂_unit)
    results = Matrix{Float64}(undef, n_samples, length(std_params.col_names))

    for i in 1:n_samples
        # Step 1: SA samples a direction on the unit sphere
        k0 = rand(1:K)
        ξ₀ = X̂_unit[:, k0] .+ 0.01 .* randn(d)
        ξ₀ ./= norm(ξ₀)

        result = sample(X̂_unit, ξ₀, T; β=β, α=α)
        ξ_dir = result.Ξ[end, :]
        ξ_dir ./= (norm(ξ_dir) + 1e-12)  # ensure unit norm

        # Step 2: Draw a magnitude from the empirical norm distribution
        r = rand(pca_norms)

        # Step 3: Combine direction × magnitude
        ξ_rescaled = r .* ξ_dir

        # decode
        measurements = decode_sample(ξ_rescaled, pca_model, std_params)
        results[i, :] = measurements
    end

    df_synth = DataFrame(results, std_params.col_names)
    return df_synth
end


"""
    build_concat_matrix(df, id_col, id_values, kept_cols, n_assays;
                        visit_col=:Visit) -> Matrix{Float64}

Rebuilds the canonical `length(id_values) × (n_assays·3)` concatenated-visit
matrix used throughout this repo: row `i` is the `[V1|V2|V3]` concatenation of
`id_values[i]`'s three visits, in `kept_cols` order within each block
(`offset = (v-1)*n_assays`). This mirrors, exactly, the inline construction in
`run_full_longitudinal.jl` (lines ~57–65); that script is left untouched as
the frozen canonical reference (see `experiments/test_build_concat.jl`).

General over `id_col`/`visit_col` so it works both for real long-format
cohorts (`id_col=:SubjectID`) and synthetic ones (`id_col=:SyntheticID`).

# Arguments
- `df::DataFrame`: long-format data with one row per (id, visit).
- `id_col::Symbol`: column identifying each patient (e.g. `:SubjectID`,
  `:SyntheticID`).
- `id_values::Vector`: the ids to build rows for, in order; row `i` of the
  output corresponds to `id_values[i]`.
- `kept_cols::Vector{Symbol}`: assay column names, length `n_assays`; sets
  the within-block column order.
- `n_assays::Int`

# Keyword Arguments
- `visit_col::Symbol=:Visit`: column identifying the visit number (assumed
  `1:3`).

# Returns
- `Matrix{Float64}` of size `length(id_values) × (n_assays*3)`.
"""
function build_concat_matrix(df::DataFrame, id_col::Symbol, id_values::Vector,
                              kept_cols::Vector{Symbol}, n_assays::Int;
                              visit_col::Symbol=:Visit)
    n_visits = 3
    K = length(id_values)
    X = Matrix{Float64}(undef, K, n_assays * n_visits)
    for (i, s) in enumerate(id_values)
        for v in 1:n_visits
            mask = (df[!, id_col] .== s) .& (df[!, visit_col] .== v)
            row = findfirst(mask)
            offset = (v - 1) * n_assays
            for (j, col) in enumerate(kept_cols)
                X[i, offset + j] = df[row, col]
            end
        end
    end
    return X
end


"""
    sa_generate_from_matrix(X_concat, kept_cols, n_assays, n_samples, T;
                            β=nothing, α=0.01, pratio=0.95, seed=42)
        -> (synth_df::DataFrame, mem::NamedTuple)

Faithful extraction of the concatenated-visit SA-generation pipeline inlined
in `run_full_longitudinal.jl` (lines ~53–135), parameterized over an arbitrary
`K × (n_assays·n_visits)` real-data matrix so downstream experiments
(holdout retraining, SIM subsets, …) can reuse it without duplicating the
inline block. `run_full_longitudinal.jl` itself is left untouched as the
frozen canonical reference; this function reproduces its seed-42 output
exactly (see `experiments/test_sa_helper.jl`).

Mirrors, in the same order (so the RNG call sequence matches):
standardize → PCA (`pratio`) → unit-norm columns (+ retain `pca_norms`)
→ β* via `find_entropy_inflection` (only if `β === nothing`)
→ `Random.seed!(seed)` → per-sample ULA `sample(X̂, ξ₀, T; β=β_star, α=α)`
→ rescale the unit direction by a magnitude drawn from `pca_norms`
→ `reconstruct` → destandardize → long-format DataFrame.

# Arguments
- `X_concat::Matrix{Float64}`: `K × (n_assays·n_visits)` concatenated real
  rows; visits are laid out contiguously (columns `1:n_assays` = visit 1,
  `n_assays+1:2n_assays` = visit 2, …).
- `kept_cols::Vector{Symbol}`: assay column names, length `n_assays`.
- `n_assays::Int`
- `n_samples::Int`
- `T::Int`: ULA steps per sample.

# Keyword Arguments
- `β::Union{Float64,Nothing}=nothing`: fixed inverse temperature. If
  `nothing`, computed once via
  `find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1,1000.0))`
  (this inner call always uses `α=0.01`, independent of the `α` below).
- `α::Float64=0.01`: ULA step size used for the per-sample `sample(...)` calls.
- `pratio::Float64=0.95`: PCA variance-retained ratio.
- `seed::Union{Int,Nothing}=42`: seeds the RNG once before the generation
  loop (matching the canonical script's `Random.seed!(42)`); pass `nothing`
  to leave the RNG state as-is.

# Returns
- `synth_df::DataFrame`: long format, columns = `kept_cols` plus
  `SyntheticID` and `Visit`, sorted by `(SyntheticID, Visit)`,
  `n_samples * n_visits` rows.
- `mem::NamedTuple`: `(X̂, pca_model, std_params, pca_norms, β_star, d_pca)`.
"""
function sa_generate_from_matrix(X_concat::Matrix{Float64}, kept_cols::Vector{Symbol},
                                  n_assays::Int, n_samples::Int, T::Int;
                                  β::Union{Float64,Nothing}=nothing,
                                  α::Float64=0.01, pratio::Float64=0.95,
                                  seed::Union{Int,Nothing}=42)

    K, d_concat = size(X_concat)
    n_visits = d_concat ÷ n_assays
    n_visits * n_assays == d_concat || throw(DimensionMismatch(
        "X_concat has $d_concat columns, not divisible by n_assays=$n_assays"))

    # ── Step 3 (mirrors run_full_longitudinal.jl:79-107): standardize + PCA ──
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
    pca_concat = MultivariateStats.fit(PCA, Z_t; pratio=pratio)
    d_pca = MultivariateStats.outdim(pca_concat)
    X_pca = MultivariateStats.transform(pca_concat, Z_t)

    pca_norms = [norm(X_pca[:, k]) for k in 1:K]
    X̂ = copy(X_pca)
    for k in 1:K
        X̂[:, k] ./= (norm(X̂[:, k]) + 1e-12)
    end

    # ── Step 4 (mirrors run_full_longitudinal.jl:112-135): β* and generate ──
    β_star = β
    if β_star === nothing
        phase = find_entropy_inflection(X̂; α=0.01, n_betas=80, β_range=(0.1, 1000.0))
        β_star = phase.β_star
    end

    if seed !== nothing
        Random.seed!(seed)
    end

    synth_concat = Matrix{Float64}(undef, n_samples, d_concat)

    for i in 1:n_samples
        k0 = rand(1:K)
        ξ₀ = X̂[:, k0] .+ 0.01 .* randn(d_pca)
        ξ₀ ./= norm(ξ₀)
        result = sample(X̂, ξ₀, T; β=β_star, α=α)
        ξ_dir = result.Ξ[end, :]
        ξ_dir ./= (norm(ξ_dir) + 1e-12)
        r = rand(pca_norms)
        ξ_rescaled = r .* ξ_dir
        z_std = vec(MultivariateStats.reconstruct(pca_concat, ξ_rescaled))
        synth_concat[i, :] = z_std .* σ_concat .+ μ_concat
    end

    # split into per-visit blocks, stack long-format
    df_synth_visits = DataFrame[]
    for v in 1:n_visits
        offset = (v - 1) * n_assays
        cols = synth_concat[:, (offset+1):(offset+n_assays)]
        df_v = DataFrame(cols, kept_cols)
        df_v.SyntheticID = 1:n_samples
        df_v.Visit = fill(v, n_samples)
        push!(df_synth_visits, df_v)
    end
    synth_df = vcat(df_synth_visits...)
    sort!(synth_df, [:SyntheticID, :Visit])

    mem = (X̂=X̂, pca_model=pca_concat, std_params=std_params_concat,
           pca_norms=pca_norms, β_star=β_star, d_pca=d_pca)

    return synth_df, mem
end


"""
    generate_conditioned_patients_rescaled(X̂_unit, pca_model, std_params,
        pca_norms, condition_mask, n_samples, T; β, α=0.1, ρ_ratio=10.0,
        seed=nothing) -> DataFrame

Conditioned version of rescaled generation. Uses multiplicity-weighted SA
for direction sampling and draws norms from the condition-specific subset.
"""
function generate_conditioned_patients_rescaled(X̂_unit::Matrix{Float64}, pca_model,
                                                 std_params::StandardizationParams,
                                                 pca_norms::Vector{Float64},
                                                 condition_mask::BitVector,
                                                 n_samples::Int, T::Int;
                                                 β::Float64, α::Float64=0.1,
                                                 ρ_ratio::Float64=10.0,
                                                 seed::Union{Int, Nothing}=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    d, K = size(X̂_unit)
    length(condition_mask) == K || throw(DimensionMismatch(
        "condition_mask length $(length(condition_mask)) ≠ K=$K"))

    # multiplicity vector
    ρ = ones(K)
    for k in 1:K
        if condition_mask[k]
            ρ[k] = ρ_ratio
        end
    end

    # condition-specific norms for magnitude sampling
    cond_norms = pca_norms[condition_mask]
    isempty(cond_norms) && (cond_norms = pca_norms)  # fallback to all

    results = Matrix{Float64}(undef, n_samples, length(std_params.col_names))

    for i in 1:n_samples
        cond_indices = findall(condition_mask)
        k0 = isempty(cond_indices) ? rand(1:K) : rand(cond_indices)
        ξ₀ = X̂_unit[:, k0] .+ 0.01 .* randn(d)
        ξ₀ ./= norm(ξ₀)

        result = weighted_sample(X̂_unit, ξ₀, T, ρ; β=β, α=α)
        ξ_dir = result.Ξ[end, :]
        ξ_dir ./= (norm(ξ_dir) + 1e-12)

        r = rand(cond_norms)
        ξ_rescaled = r .* ξ_dir

        measurements = decode_sample(ξ_rescaled, pca_model, std_params)
        results[i, :] = measurements
    end

    df_synth = DataFrame(results, std_params.col_names)
    return df_synth
end


# ══════════════════════════════════════════════════════════════════════════════
# Evaluation metrics for synthetic patient quality
# ══════════════════════════════════════════════════════════════════════════════

"""
    feature_summary_comparison(df_real, df_synth, cols) -> DataFrame

Compare mean and std of each feature between real and synthetic data.
"""
function feature_summary_comparison(df_real::DataFrame, df_synth::DataFrame, cols::Vector{Symbol})
    comparison = DataFrame(
        Feature = String[],
        Real_Mean = Float64[],
        Real_Std = Float64[],
        Synth_Mean = Float64[],
        Synth_Std = Float64[],
        Mean_Rel_Error = Float64[]
    )

    for c in cols
        real_vals = collect(skipmissing(df_real[!, c]))
        synth_vals = df_synth[!, c]

        μ_r = mean(real_vals)
        σ_r = std(real_vals)
        μ_s = mean(synth_vals)
        σ_s = std(synth_vals)
        rel_err = abs(μ_r) > 1e-12 ? abs(μ_s - μ_r) / abs(μ_r) : abs(μ_s - μ_r)

        push!(comparison, (string(c), μ_r, σ_r, μ_s, σ_s, rel_err))
    end

    return comparison
end


"""
    column_correlation_matrix(df, cols) -> Matrix{Float64}

Compute the correlation matrix for the given columns.
"""
function column_correlation_matrix(df::DataFrame, cols::Vector{Symbol})
    X = Matrix{Float64}(coalesce.(df[!, cols], NaN))
    # remove NaN rows
    good = .!any(isnan.(X), dims=2)[:]
    X = X[good, :]
    return cor(X)
end

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
const _COL_ASSAY_END = 55    # last protein assay column (C5b-9), before TGA Y/N flags

"""
    load_legacy_data(filepath) -> DataFrame

Load the legacy biochemical dataset from the Excel file.
Returns a DataFrame with subject metadata and continuous assay values.
Non-numeric assay values (e.g., '#N/A', empty) are replaced with `missing`.
"""
function load_legacy_data(filepath::String)
    xf = XLSX.readxlsx(filepath)
    ws = xf["Assay Values"]

    # read header names from row 3
    headers = String[]
    for col in _COL_ASSAY_START:_COL_ASSAY_END
        val = ws[_HEADER_ROW, col]
        push!(headers, val === nothing ? "Assay_$col" : string(val))
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
        # subject ID: propagate from last non-empty value
        sid_raw = ws[row, _COL_SUBJECT_ID]
        if sid_raw !== nothing && sid_raw !== ""
            current_subject = Int(sid_raw)
        end

        visit_raw = ws[row, _COL_VISIT]
        visit_raw === nothing && continue  # skip empty rows

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

        for (j, col) in enumerate(_COL_ASSAY_START:_COL_ASSAY_END)
            val = ws[row, col]
            if val === nothing || val === "" || (val isa String && (val == "#N/A" || val == "---"))
                push!(records[headers[j]], missing)
            elseif val isa Number
                push!(records[headers[j]], Float64(val))
            else
                # try to parse as number
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

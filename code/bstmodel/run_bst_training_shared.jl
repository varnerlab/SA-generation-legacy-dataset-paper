# ──────────────────────────────────────────────────────────────────────────────
# run_bst_training_shared.jl — BST Training with Shared Parameters Across Visits
#
# Fits ONE set of 12 BST parameters per patient, optimized jointly across all
# 3 visits. The model must explain the patient's changing TGA output purely
# from their changing coagulation factor levels (initial conditions).
#
# This tests whether the BST kinetics are patient-specific but time-invariant:
# the same biochemistry, different inputs at each visit.
#
# Real data: 23 patients with complete longitudinal data (Visits 1, 2, 3)
# Synthetic data: 100 SA-generated patients (Visits 1, 2, 3)
# Total training jobs: 23 + 100 = 123 patient fits (each across 3 visits)
#
# Parallelized via Distributed.jl (pmap).
# ──────────────────────────────────────────────────────────────────────────────

using Distributed

# add workers (leave 1 core for the main process)
const NWORKERS = max(1, Sys.CPU_THREADS - 1)
if nworkers() < NWORKERS
    addprocs(NWORKERS - nworkers() + 1)
    @info "Launched $(nworkers()) workers"
end

# load dependencies on ALL workers
@everywhere begin
    include(joinpath(@__DIR__, "Include.jl"))

    # ──────────────────────────────────────────────────────────────────────────
    # Constants
    # ──────────────────────────────────────────────────────────────────────────
    const PLASMA_CONCENTRATION_NM = Dict(
        :II   => 1400.0,  :V    => 20.0,   :VII  => 10.0,
        :VIII => 0.7,     :IX   => 90.0,   :X    => 170.0,
        :XI   => 30.0,    :XII  => 375.0,  :AT   => 3400.0,
        :PC   => 65.0,    :TFPI => 2.5,
    )

    const SF = 1e9
    const TFPI_COL_SOURCE = "TFPI Free"

    # Training configuration
    const N_RESTARTS = 10
    const ITERATIONS_PER_RESTART = 500
    const TIME_LIMIT_PER_RESTART = 300.0

    # ──────────────────────────────────────────────────────────────────────────
    # Multi-visit loss: sum of normalized relative squared error across visits
    #
    # For each visit, re-setup the model ICs from that visit's data, apply the
    # candidate parameters κ, evaluate, and accumulate the error.
    # ──────────────────────────────────────────────────────────────────────────
    function _multi_visit_loss(κ::Vector{Float64}, model::BSTModel,
                               visit_indices::Vector{Int},
                               visit_dfs::Vector{DataFrame},
                               visit_targets::Vector{Vector{Float64}})

        total_ϵ = 0.0
        n_visits = length(visit_indices)

        for v in 1:n_visits
            idx = visit_indices[v]
            df  = visit_dfs[v]
            Y   = visit_targets[v]

            # set ICs and static factors for this visit
            _setup_model!(model, idx, df)

            # apply candidate parameters
            _apply_parameters!(model, κ)

            # evaluate
            (T, Xₘ) = evaluate(model, tspan=(0.0, 40.0))
            Yₘ = model_output_vector(T, Xₘ[:, 9])

            # normalized relative squared error
            for i in 1:length(Yₘ)
                if abs(Y[i]) > 1e-12
                    total_ϵ += ((Y[i] - Yₘ[i]) / Y[i])^2
                else
                    total_ϵ += (Y[i] - Yₘ[i])^2
                end
            end
        end

        # average over visits so error scale is comparable to per-visit fits
        return total_ϵ / n_visits
    end

    # ──────────────────────────────────────────────────────────────────────────
    # learn_bst_shared_parameters — multi-start optimization across all visits
    # ──────────────────────────────────────────────────────────────────────────
    function learn_bst_shared_parameters(model::BSTModel,
                                          visit_indices::Vector{Int},
                                          visit_dfs::Vector{DataFrame};
                                          pₒ::Union{Nothing,Vector{Float64}} = nothing,
                                          n_restarts::Int = N_RESTARTS,
                                          iterations_per_restart::Int = ITERATIONS_PER_RESTART,
                                          time_limit_per_restart::Float64 = TIME_LIMIT_PER_RESTART,
                                          show_trace::Bool = false)

        # extract target vectors for each visit
        visit_targets = Vector{Vector{Float64}}()
        for v in 1:length(visit_indices)
            idx = visit_indices[v]
            df  = visit_dfs[v]
            Y = Float64[
                df[idx, :Lagtime],
                df[idx, :Peak],
                df[idx, :TPeak],
                df[idx, :Max],
                df[idx, :EPT],
            ]
            push!(visit_targets, Y)
        end

        # parameter bounds
        κ = copy(_DEFAULT_PARAMETER_BOUNDS)
        P = size(κ, 1)
        lb = κ[:, 2]
        ub = κ[:, 3]

        # objective: sum loss across all visits
        obj_function(p) = _multi_visit_loss(p, model, visit_indices, visit_dfs, visit_targets)

        # multi-start optimization
        best_p = nothing
        best_ϵ = Inf
        inner_optimizer = NelderMead()

        for restart in 1:n_restarts

            if restart == 1 && !isnothing(pₒ)
                p_start = clamp.(copy(pₒ), lb, ub)
            elseif restart == 1
                p_start = clamp.(κ[:, 1] .* (1.0 .+ 0.1 .* randn(P)), lb, ub)
            else
                p_start = [10.0^(log10(lb[i]) + rand() * (log10(ub[i]) - log10(lb[i]))) for i in 1:P]
            end

            try
                results = optimize(obj_function, lb, ub, p_start, Fminbox(inner_optimizer),
                    Optim.Options(time_limit=time_limit_per_restart, show_trace=show_trace,
                                  show_every=100, iterations=iterations_per_restart))

                ϵ = Optim.minimum(results)
                if ϵ < best_ϵ
                    best_ϵ = ϵ
                    best_p = Optim.minimizer(results)
                end
            catch e
                continue
            end
        end

        if best_p === nothing
            error("All $n_restarts optimization restarts failed")
        end

        # run final simulations with best parameters to get per-visit results
        per_visit_results = []
        for v in 1:length(visit_indices)
            idx = visit_indices[v]
            df  = visit_dfs[v]

            _setup_model!(model, idx, df)
            _apply_parameters!(model, best_p)
            (T, Xₘ) = evaluate(model, tspan=(0.0, 40.0))
            Yₘ = model_output_vector(T, Xₘ[:, 9])

            # per-visit error
            Y = visit_targets[v]
            ϵ_v = sum(((Y[i] - Yₘ[i]) / max(abs(Y[i]), 1e-12))^2 for i in 1:length(Y)) / length(Y)

            push!(per_visit_results, (T=T, Xm=Xₘ, Ym=Yₘ, Y=Y, error=ϵ_v))
        end

        return (p_best=best_p, total_error=best_ϵ, visits=per_visit_results)
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Worker function: train a single patient across all visits
    # ──────────────────────────────────────────────────────────────────────────
    function train_single_patient_shared(patient_row_indices::Vector{Int},
                                          model_path::String,
                                          visit_dfs::Vector{DataFrame};
                                          pₒ::Union{Nothing,Vector{Float64}} = nothing)

        model = build(model_path)

        try
            result = learn_bst_shared_parameters(model, patient_row_indices, visit_dfs;
                                                   pₒ=pₒ)
            return result
        catch e
            @warn "Training FAILED: $e"
            return nothing
        end
    end
end  # @everywhere

# ──────────────────────────────────────────────────────────────────────────────
# Helper: convert raw data into the BST training format
# ──────────────────────────────────────────────────────────────────────────────
function build_training_dataframe(df_source::DataFrame)::DataFrame

    training_df = DataFrame(
        II      = Float64[], V       = Float64[], VII     = Float64[],
        VIII    = Float64[], IX      = Float64[], X       = Float64[],
        XI      = Float64[], XII     = Float64[], AT      = Float64[],
        PC      = Float64[], TFPI    = Float64[],
        Lagtime = Float64[], Peak    = Float64[], TPeak   = Float64[],
        Max     = Float64[], EPT     = Float64[],
    )

    for i in 1:nrow(df_source)
        row = df_source[i, :]
        push!(training_df, (
            II      = (row["II"]   / 100.0) * PLASMA_CONCENTRATION_NM[:II],
            V       = (row["V"]    / 100.0) * PLASMA_CONCENTRATION_NM[:V],
            VII     = (row["VII"]  / 100.0) * PLASMA_CONCENTRATION_NM[:VII],
            VIII    = (row["VIII"] / 100.0) * PLASMA_CONCENTRATION_NM[:VIII],
            IX      = (row["IX"]   / 100.0) * PLASMA_CONCENTRATION_NM[:IX],
            X       = (row["X"]    / 100.0) * PLASMA_CONCENTRATION_NM[:X],
            XI      = (row["XI"]   / 100.0) * PLASMA_CONCENTRATION_NM[:XI],
            XII     = (row["XII"]  / 100.0) * PLASMA_CONCENTRATION_NM[:XII],
            AT      = (row["AT"]   / 100.0) * PLASMA_CONCENTRATION_NM[:AT],
            PC      = (row["PC"]   / 100.0) * PLASMA_CONCENTRATION_NM[:PC],
            TFPI    = (row[TFPI_COL_SOURCE] / 100.0) * PLASMA_CONCENTRATION_NM[:TFPI],
            Lagtime = row["TF Initiator Lagtime (min)"],
            Peak    = row["TF Initiator Peak (nM)"],
            TPeak   = row["TF Initiator T.Peak (min)"],
            Max     = row["TF Initiator Max Rate (nM/min)"],
            EPT     = row["TF Initiator ETP (nM·min)"],
        ))
    end

    return training_df
end

# ──────────────────────────────────────────────────────────────────────────────
# Helper: parallel training with warm-start
# ──────────────────────────────────────────────────────────────────────────────
function train_patients_shared_parallel(model_path::String,
                                         patient_visit_map::Vector{Vector{Int}},
                                         visit_dfs::Vector{DataFrame};
                                         label::String = "")

    t_start = time()
    n = length(patient_visit_map)
    results = Dict{Int, NamedTuple}()

    # Phase 1: train first patient serially for warm-start
    @info "[$label] Phase 1: training patient 1/$n (serial, for warm-start)..."
    res1 = train_single_patient_shared(patient_visit_map[1], model_path, visit_dfs)

    warm_start = nothing
    if res1 !== nothing
        results[1] = res1
        warm_start = copy(res1.p_best)
        @info "[$label] Patient 1 done — total_err=$(round(res1.total_error; digits=4)). Warm-start acquired."
    else
        @warn "[$label] First patient failed — proceeding without warm-start."
    end

    # Phase 2: train remaining patients in parallel
    remaining_indices = collect(2:n)
    if !isempty(remaining_indices)
        @info "[$label] Phase 2: training $(length(remaining_indices)) patients in parallel ($(nworkers()) workers)..."

        parallel_results = pmap(remaining_indices) do k
            train_single_patient_shared(patient_visit_map[k], model_path, visit_dfs; pₒ=warm_start)
        end

        for (i, res) in enumerate(parallel_results)
            k = remaining_indices[i]
            if res !== nothing
                results[k] = res
            end
        end
    end

    total_time = round((time() - t_start) / 60.0; digits=1)
    @info "[$label] Complete: $(length(results))/$n converged in $(total_time) min"

    return results
end

# ══════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

# --- STEP 1: Load data and filter to complete patients ---
@info "=" ^ 72
@info "STEP 1: Loading data..."

path_to_csv = joinpath(_PATH_TO_DATA, "cleaned_full_data.csv")
df_real = CSV.read(path_to_csv, DataFrame)

# identify complete patients (present at all 3 visits)
gd_real = groupby(df_real, :SubjectID)
complete_ids = sort([g.SubjectID[1] for g in gd_real if length(unique(g.Visit)) == 3])
df_real_complete = filter(row -> row.SubjectID in complete_ids, df_real)
@info "Complete real patients: $(length(complete_ids))"

path_to_synthetic = joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv")
df_synthetic = CSV.read(path_to_synthetic, DataFrame)
syn_ids = sort(unique(df_synthetic.SyntheticID))
@info "Synthetic patients: $(length(syn_ids))"

# --- STEP 2: Build per-visit training DataFrames ---
@info "=" ^ 72
@info "STEP 2: Building per-visit training DataFrames..."

bst_model_path = joinpath(_PATH_TO_MODEL, "Coagulation.bst")

# Real: build training DFs per visit, and map each patient to its row index in each visit's DF
real_visit_dfs = DataFrame[]
real_patient_visit_map = Vector{Vector{Int}}()  # patient_visit_map[k] = [row_in_v1, row_in_v2, row_in_v3]

for visit in [1, 2, 3]
    df_v = filter(row -> row.SubjectID in complete_ids && row.Visit == visit, df_real)
    # sort by SubjectID so patients align across visits
    sort!(df_v, :SubjectID)
    push!(real_visit_dfs, build_training_dataframe(df_v))
end

# each complete patient maps to the same row index across all 3 visit DFs (since sorted by SubjectID)
for k in 1:length(complete_ids)
    push!(real_patient_visit_map, [k, k, k])
end

@info "Real: $(length(real_patient_visit_map)) patients × 3 visits"

# Synthetic: same approach
syn_visit_dfs = DataFrame[]
syn_patient_visit_map = Vector{Vector{Int}}()

for visit in [1, 2, 3]
    df_v = filter(row -> row.Visit == visit, df_synthetic)
    sort!(df_v, :SyntheticID)
    push!(syn_visit_dfs, build_training_dataframe(df_v))
end

for k in 1:length(syn_ids)
    push!(syn_patient_visit_map, [k, k, k])
end

@info "Synthetic: $(length(syn_patient_visit_map)) patients × 3 visits"

# --- STEP 3: Train real patients (shared parameters across visits) ---
@info "=" ^ 72
@info "STEP 3: Training $(length(complete_ids)) real patients (shared params, $(nworkers()) workers)..."

results_real = train_patients_shared_parallel(bst_model_path, real_patient_visit_map,
                                               real_visit_dfs; label="Real-Shared")

# --- STEP 4: Train synthetic patients (shared parameters across visits) ---
@info "=" ^ 72
@info "STEP 4: Training $(length(syn_ids)) synthetic patients (shared params, $(nworkers()) workers)..."

results_synthetic = train_patients_shared_parallel(bst_model_path, syn_patient_visit_map,
                                                    syn_visit_dfs; label="Synth-Shared")

# --- STEP 5: Save results ---
@info "=" ^ 72
@info "STEP 5: Saving results..."

path_results = joinpath(_PATH_TO_DATA, "bst_training_results_shared.jld2")
jldsave(path_results;
    results_real          = results_real,
    results_synthetic     = results_synthetic,
    complete_patient_ids  = complete_ids,
    synthetic_patient_ids = syn_ids,
    n_real_patients       = length(complete_ids),
    n_synthetic_patients  = length(syn_ids),
    n_restarts            = N_RESTARTS,
    iterations_per_restart = ITERATIONS_PER_RESTART,
    n_workers             = nworkers(),
)
@info "Results saved: $path_results"

# --- STEP 6: Summary ---
@info "=" ^ 72
@info "BST Shared-Parameter Training Complete"
@info "  Workers used:       $(nworkers())"
@info "  Real patients:      $(length(results_real)) / $(length(complete_ids)) converged"
@info "  Synthetic patients: $(length(results_synthetic)) / $(length(syn_ids)) converged"

if !isempty(results_real)
    errs = [results_real[k].total_error for k in keys(results_real)]
    @info "  Real total errors:      mean=$(round(mean(errs); digits=4)), " *
          "min=$(round(minimum(errs); digits=4)), max=$(round(maximum(errs); digits=4))"

    # per-visit breakdown
    for v in 1:3
        v_errs = [results_real[k].visits[v].error for k in keys(results_real)]
        @info "    Visit $v: mean=$(round(mean(v_errs); digits=4)), max=$(round(maximum(v_errs); digits=4))"
    end
end

if !isempty(results_synthetic)
    errs = [results_synthetic[k].total_error for k in keys(results_synthetic)]
    @info "  Synthetic total errors: mean=$(round(mean(errs); digits=4)), " *
          "min=$(round(minimum(errs); digits=4)), max=$(round(maximum(errs); digits=4))"

    for v in 1:3
        v_errs = [results_synthetic[k].visits[v].error for k in keys(results_synthetic)]
        @info "    Visit $v: mean=$(round(mean(v_errs); digits=4)), max=$(round(maximum(v_errs); digits=4))"
    end
end

@info "Done."

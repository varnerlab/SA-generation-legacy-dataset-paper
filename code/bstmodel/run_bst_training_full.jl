# ──────────────────────────────────────────────────────────────────────────────
# run_bst_training_full.jl — BST Coagulation Model Training (Full Cohort)
#
# Trains a 9-reaction BST coagulation model (12 trainable parameters:
# 9 alpha values + 3 G-matrix entries) on complete patients across all 3 visits.
#
# Real data: 23 patients with complete longitudinal data (Visits 1, 2, 3)
# Synthetic data: 100 SA-generated patients (Visits 1, 2, 3)
# Total training jobs: (23 + 100) × 3 visits = 369 patient-visit fits
#
# Parallelized via Distributed.jl (pmap) — each patient-visit trains independently.
# Warm-starting: first patient uses defaults, subsequent patients use the
# fitted parameters from the first as a starting point.
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
    # Worker function: train a single patient, return result or nothing
    # ──────────────────────────────────────────────────────────────────────────
    function train_single_patient(idx::Int, model_path::String,
                                  training_df::DataFrame;
                                  pₒ::Union{Nothing,Vector{Float64}} = nothing,
                                  n_restarts::Int = N_RESTARTS)

        model = build(model_path)

        try
            (p_best, T, Xm, Ym, Y, err) = learn_bst_model_parameters(
                idx, model, training_df;
                pₒ = pₒ,
                show_trace = false,
                n_restarts = n_restarts,
                iterations_per_restart = ITERATIONS_PER_RESTART,
                time_limit_per_restart = TIME_LIMIT_PER_RESTART,
            )

            return (
                idx    = idx,
                p_best = p_best,
                T      = T,
                Xm     = Xm,
                Ym     = Ym,
                Y      = Y,
                error  = err,
            )
        catch e
            @warn "Training FAILED for patient $idx: $e"
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
# Helper: parallel training with warm-start strategy
#
# Phase 1 (serial): Train first patient with defaults to get an initial p_best.
# Phase 2 (parallel): Train remaining patients via pmap, warm-started from the
#          Phase 1 result.
# ──────────────────────────────────────────────────────────────────────────────
function train_patients_parallel(model_path::String, training_df::DataFrame,
                                 patient_indices::AbstractVector{Int};
                                 label::String = "")

    t_start = time()
    results = Dict{Int, NamedTuple}()
    n = length(patient_indices)

    # Phase 1: train first patient serially to get warm-start parameters
    @info "[$label] Phase 1: training patient 1/$n (serial, for warm-start)..."
    first_idx = patient_indices[1]
    res1 = train_single_patient(first_idx, model_path, training_df)

    warm_start = nothing
    if res1 !== nothing
        results[res1.idx] = (p_best=res1.p_best, T=res1.T, Xm=res1.Xm,
                              Ym=res1.Ym, Y=res1.Y, error=res1.error)
        warm_start = copy(res1.p_best)
        @info "[$label] Patient $(res1.idx) done — err=$(round(res1.error; digits=4)). Warm-start acquired."
    else
        @warn "[$label] First patient failed — proceeding without warm-start."
    end

    # Phase 2: train remaining patients in parallel
    remaining = patient_indices[2:end]
    if !isempty(remaining)
        @info "[$label] Phase 2: training $(length(remaining)) patients in parallel ($(nworkers()) workers)..."

        parallel_results = pmap(remaining) do idx
            train_single_patient(idx, model_path, training_df; pₒ=warm_start)
        end

        for res in parallel_results
            if res !== nothing
                results[res.idx] = (p_best=res.p_best, T=res.T, Xm=res.Xm,
                                     Ym=res.Ym, Y=res.Y, error=res.error)
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
complete_ids = [g.SubjectID[1] for g in gd_real if length(unique(g.Visit)) == 3]
sort!(complete_ids)
df_real_complete = filter(row -> row.SubjectID in complete_ids, df_real)
@info "Complete real patients: $(length(complete_ids)) ($(nrow(df_real_complete)) total rows across 3 visits)"

path_to_synthetic = joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv")
df_synthetic = CSV.read(path_to_synthetic, DataFrame)
syn_ids = sort(unique(df_synthetic.SyntheticID))
@info "Synthetic patients: $(length(syn_ids)) ($(nrow(df_synthetic)) total rows across 3 visits)"

# --- STEP 2: Build training DataFrames per visit ---
@info "=" ^ 72
@info "STEP 2: Building training DataFrames per visit..."

bst_model_path = joinpath(_PATH_TO_MODEL, "Coagulation.bst")

# store all results keyed by (source, visit)
all_results = Dict{String, Dict{Int, NamedTuple}}()

for visit in [1, 2, 3]

    # --- Real patients for this visit ---
    @info "=" ^ 72
    @info "VISIT $visit — Real patients"

    df_real_visit = filter(row -> row.SubjectID in complete_ids && row.Visit == visit, df_real)
    @info "  Real patients at Visit $visit: $(nrow(df_real_visit))"
    training_df_real = build_training_dataframe(df_real_visit)

    results_real = train_patients_parallel(bst_model_path, training_df_real,
                                            collect(1:nrow(training_df_real));
                                            label="Real-V$visit")
    all_results["real_v$visit"] = results_real

    # --- Synthetic patients for this visit ---
    @info "=" ^ 72
    @info "VISIT $visit — Synthetic patients"

    df_syn_visit = filter(row -> row.Visit == visit, df_synthetic)
    @info "  Synthetic patients at Visit $visit: $(nrow(df_syn_visit))"
    training_df_syn = build_training_dataframe(df_syn_visit)

    results_syn = train_patients_parallel(bst_model_path, training_df_syn,
                                           collect(1:nrow(training_df_syn));
                                           label="Synth-V$visit")
    all_results["synthetic_v$visit"] = results_syn
end

# --- STEP 3: Save results ---
@info "=" ^ 72
@info "STEP 3: Saving results..."

path_results = joinpath(_PATH_TO_DATA, "bst_training_results_full.jld2")
jldsave(path_results;
    all_results         = all_results,
    complete_patient_ids = complete_ids,
    n_real_patients     = length(complete_ids),
    n_synthetic_patients = length(syn_ids),
    n_restarts          = N_RESTARTS,
    iterations_per_restart = ITERATIONS_PER_RESTART,
    n_workers           = nworkers(),
)
@info "Results saved: $path_results"

# --- STEP 4: Summary ---
@info "=" ^ 72
@info "BST Full Training Complete"
@info "  Workers used:       $(nworkers())"
@info "  Real patients:      $(length(complete_ids)) × 3 visits"
@info "  Synthetic patients: $(length(syn_ids)) × 3 visits"

for visit in [1, 2, 3]
    for (source, key) in [("Real", "real_v$visit"), ("Synthetic", "synthetic_v$visit")]
        res = all_results[key]
        if !isempty(res)
            errs = [res[k].error for k in keys(res)]
            @info "  $source V$visit: $(length(res)) converged, " *
                  "err mean=$(round(mean(errs); digits=4)), " *
                  "min=$(round(minimum(errs); digits=4)), " *
                  "max=$(round(maximum(errs); digits=4))"
        else
            @info "  $source V$visit: 0 converged"
        end
    end
end

@info "Done."

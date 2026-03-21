# ──────────────────────────────────────────────────────────────────────────────
# run_bst_training.jl — BST Coagulation Model Training (Proof of Concept)
#
# Trains a 9-reaction BST coagulation model (12 trainable parameters:
# 9 alpha values + 3 G-matrix entries) on patient TGA data. The model takes
# coagulation factor levels as initial conditions and predicts the thrombin
# (FIIa) generation curve. Five TGA summary features are extracted:
# Lagtime, Peak, T.Peak, MaxRate, and ETP.
#
# This script trains on a small subset of real patients and synthetic patients
# as a proof of concept.
# ──────────────────────────────────────────────────────────────────────────────

# load environment and dependencies
include(joinpath(@__DIR__, "Include.jl"))

# ──────────────────────────────────────────────────────────────────────────────
# Constants: physiological plasma concentrations (nM) for % nominal conversion
# Reference: standard plasma concentrations for coagulation factors
# ──────────────────────────────────────────────────────────────────────────────
const PLASMA_CONCENTRATION_NM = Dict(
    :II   => 1400.0,    # FII (prothrombin)
    :V    => 20.0,      # FV
    :VII  => 10.0,      # FVII
    :VIII => 0.7,       # FVIII
    :IX   => 90.0,      # FIX
    :X    => 170.0,     # FX
    :XI   => 30.0,      # FXI
    :XII  => 375.0,     # FXII
    :AT   => 3400.0,    # Antithrombin
    :PC   => 65.0,      # Protein C
    :TFPI => 2.5,       # TFPI
)

# scale factor used throughout the BST model (converts M to nM-like units)
const SF = 1e9

# coagulation factor column names in the source data (% nominal)
const COAG_FACTOR_COLS_SOURCE = [
    "II", "V", "VII", "VIII", "IX", "X", "XI", "XII", "AT", "PC"
]

# TFPI uses a different column name in the source data
const TFPI_COL_SOURCE = "TFPI Free"

# TGA feature column names in the source data
const TGA_COL_MAP = Dict(
    "TF Initiator Lagtime (min)"      => :Lagtime,
    "TF Initiator Peak (nM)"          => :Peak,
    "TF Initiator T.Peak (min)"       => :TPeak,
    "TF Initiator Max Rate (nM/min)"  => :Max,
    "TF Initiator ETP (nM·min)"       => :EPT,
)

# ──────────────────────────────────────────────────────────────────────────────
# Helper: convert a raw data row into the BST training format
# Coag factor levels are converted from % nominal to nM concentrations.
# ──────────────────────────────────────────────────────────────────────────────
"""
    build_training_dataframe(df_source::DataFrame) -> DataFrame

Convert a source DataFrame (with % nominal coagulation factors and TGA features)
into the training DataFrame expected by `learn_bst_model_parameters`. Coagulation
factor columns are converted from % nominal to nM using standard plasma
concentrations. The resulting DataFrame has columns:
    II, V, VII, VIII, IX, X, XI, XII, AT, PC, TFPI, Lagtime, Peak, TPeak, Max, EPT
"""
function build_training_dataframe(df_source::DataFrame)::DataFrame

    @info "Building training DataFrame from source data ($(nrow(df_source)) rows)..."

    # initialize the output DataFrame with the column schema expected by Learn.jl
    training_df = DataFrame(
        II      = Float64[],
        V       = Float64[],
        VII     = Float64[],
        VIII    = Float64[],
        IX      = Float64[],
        X       = Float64[],
        XI      = Float64[],
        XII     = Float64[],
        AT      = Float64[],
        PC      = Float64[],
        TFPI    = Float64[],
        Lagtime = Float64[],
        Peak    = Float64[],
        TPeak   = Float64[],
        Max     = Float64[],
        EPT     = Float64[],
    )

    for i in 1:nrow(df_source)

        row = df_source[i, :]

        # convert coagulation factors from % nominal to nM
        # formula: concentration_nM = (percent_nominal / 100) * plasma_concentration_nM
        II_nM   = (row["II"]   / 100.0) * PLASMA_CONCENTRATION_NM[:II]
        V_nM    = (row["V"]    / 100.0) * PLASMA_CONCENTRATION_NM[:V]
        VII_nM  = (row["VII"]  / 100.0) * PLASMA_CONCENTRATION_NM[:VII]
        VIII_nM = (row["VIII"] / 100.0) * PLASMA_CONCENTRATION_NM[:VIII]
        IX_nM   = (row["IX"]   / 100.0) * PLASMA_CONCENTRATION_NM[:IX]
        X_nM    = (row["X"]    / 100.0) * PLASMA_CONCENTRATION_NM[:X]
        XI_nM   = (row["XI"]   / 100.0) * PLASMA_CONCENTRATION_NM[:XI]
        XII_nM  = (row["XII"]  / 100.0) * PLASMA_CONCENTRATION_NM[:XII]
        AT_nM   = (row["AT"]   / 100.0) * PLASMA_CONCENTRATION_NM[:AT]
        PC_nM   = (row["PC"]   / 100.0) * PLASMA_CONCENTRATION_NM[:PC]
        TFPI_nM = (row[TFPI_COL_SOURCE] / 100.0) * PLASMA_CONCENTRATION_NM[:TFPI]

        # extract TGA features (already in physical units from the assay)
        lagtime = row["TF Initiator Lagtime (min)"]
        peak    = row["TF Initiator Peak (nM)"]
        tpeak   = row["TF Initiator T.Peak (min)"]
        maxrate = row["TF Initiator Max Rate (nM/min)"]
        etp     = row["TF Initiator ETP (nM·min)"]

        # append to training DataFrame
        push!(training_df, (
            II      = II_nM,
            V       = V_nM,
            VII     = VII_nM,
            VIII    = VIII_nM,
            IX      = IX_nM,
            X       = X_nM,
            XI      = XI_nM,
            XII     = XII_nM,
            AT      = AT_nM,
            PC      = PC_nM,
            TFPI    = TFPI_nM,
            Lagtime = lagtime,
            Peak    = peak,
            TPeak   = tpeak,
            Max     = maxrate,
            EPT     = etp,
        ))
    end

    @info "Training DataFrame built: $(nrow(training_df)) rows, $(ncol(training_df)) columns."
    return training_df
end

# ──────────────────────────────────────────────────────────────────────────────
# Helper: train BST model on a set of patients and collect results
# ──────────────────────────────────────────────────────────────────────────────
"""
    train_patients(model_path, training_df, patient_indices; label="")

Train the BST model on each patient index in `patient_indices`. Returns a
dictionary mapping patient index to a named tuple of results.
"""
function train_patients(model_path::String, training_df::DataFrame,
                        patient_indices::AbstractVector{Int}; label::String = "",
                        show_trace::Bool = false, n_restarts::Int = 10)

    results = Dict{Int, NamedTuple}()

    for (k, idx) in enumerate(patient_indices)

        @info "[$label] Training patient $k/$(length(patient_indices)) (row index $idx)..."

        # rebuild the model fresh for each patient (avoids state leakage)
        model = build(model_path)

        try
            (p_best, T, Xm, Ym, Y, err) = learn_bst_model_parameters(idx, model, training_df;
                show_trace=show_trace, n_restarts=n_restarts)

            results[idx] = (
                p_best = p_best,
                T      = T,
                Xm     = Xm,
                Ym     = Ym,       # model-predicted TGA features
                Y      = Y,        # measured TGA features
                error  = err,
            )

            @info "[$label] Patient $idx trained — error = $(round(err; digits=4))"
            @info "[$label]   Measured:  Lagtime=$(Y[1]), Peak=$(Y[2]), TPeak=$(Y[3]), Max=$(Y[4]), EPT=$(Y[5])"
            @info "[$label]   Predicted: Lagtime=$(round(Ym[1];digits=2)), Peak=$(round(Ym[2];digits=2)), TPeak=$(round(Ym[3];digits=2)), Max=$(round(Ym[4];digits=2)), EPT=$(round(Ym[5];digits=2))"

        catch e
            @warn "[$label] Training failed for patient $idx: $e"
        end
    end

    return results
end

# ──────────────────────────────────────────────────────────────────────────────
# STEP 1: Load real patient data
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "STEP 1: Loading real patient data..."

# load the cleaned full longitudinal dataset
path_to_memory = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
@info "Loading JLD2 memory file: $path_to_memory"
memory = load(path_to_memory)
df_clean = memory["df_clean"]
kept_cols = memory["kept_cols"]
@info "Loaded df_clean: $(nrow(df_clean)) rows, $(ncol(df_clean)) columns."

# also load the CSV version for column name inspection
path_to_csv = joinpath(_PATH_TO_DATA, "cleaned_full_data.csv")
@info "Loading cleaned CSV: $path_to_csv"
df_real = CSV.read(path_to_csv, DataFrame)
@info "Loaded real patient data: $(nrow(df_real)) rows, $(ncol(df_real)) columns."

# filter to Visit 1 for the proof-of-concept training
df_real_v1 = filter(row -> row.Visit == 1, df_real)
@info "Filtered to Visit 1: $(nrow(df_real_v1)) patients."

# ──────────────────────────────────────────────────────────────────────────────
# STEP 2: Load synthetic patient data
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "STEP 2: Loading synthetic patient data..."

# load the synthetic longitudinal data generated by the SA pipeline
path_to_synthetic = joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv")
@info "Loading synthetic data: $path_to_synthetic"
df_synthetic = CSV.read(path_to_synthetic, DataFrame)
@info "Loaded synthetic patient data: $(nrow(df_synthetic)) rows, $(ncol(df_synthetic)) columns."

# filter to Visit 1 for proof-of-concept
df_synthetic_v1 = filter(row -> row.Visit == 1, df_synthetic)
@info "Filtered synthetic data to Visit 1: $(nrow(df_synthetic_v1)) patients."

# ──────────────────────────────────────────────────────────────────────────────
# STEP 3: Build training DataFrames (convert % nominal -> nM)
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "STEP 3: Building training DataFrames..."

# real patients
@info "Converting real patient coagulation factors from % nominal to nM..."
training_df_real = build_training_dataframe(df_real_v1)

# synthetic patients
@info "Converting synthetic patient coagulation factors from % nominal to nM..."
training_df_synthetic = build_training_dataframe(df_synthetic_v1)

# ──────────────────────────────────────────────────────────────────────────────
# STEP 4: Build BST model and verify
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "STEP 4: Building BST coagulation model..."

bst_model_path = joinpath(_PATH_TO_MODEL, "Coagulation.bst")
@info "BST model file: $bst_model_path"
test_model = build(bst_model_path)
@info "BST model built successfully."
@info "  Dynamic states:  $(test_model.number_of_dynamic_states)"
@info "  Static factors:  $(length(test_model.static_factors_array))"
@info "  Reactions:       $(length(test_model.α))"
@info "  Species list:    $(test_model.total_species_list)"

# ──────────────────────────────────────────────────────────────────────────────
# STEP 5: Train on real patients (proof of concept — first 5 at Visit 1)
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "STEP 5: Training BST model on real patients (proof of concept)..."

n_real_train = min(2, nrow(training_df_real))
real_patient_indices = collect(1:n_real_train)
@info "Training on $n_real_train real patients (rows 1:$n_real_train of Visit 1 data)."

results_real = train_patients(bst_model_path, training_df_real, real_patient_indices;
                              label = "Real", show_trace = false, n_restarts = 10)

@info "Real patient training complete. $(length(results_real))/$n_real_train patients converged."

# ──────────────────────────────────────────────────────────────────────────────
# STEP 6: Train on synthetic patients (proof of concept — first 5 at Visit 1)
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "STEP 6: Training BST model on synthetic patients (proof of concept)..."

n_syn_train = min(2, nrow(training_df_synthetic))
syn_patient_indices = collect(1:n_syn_train)
@info "Training on $n_syn_train synthetic patients (rows 1:$n_syn_train of Visit 1 data)."

results_synthetic = train_patients(bst_model_path, training_df_synthetic, syn_patient_indices;
                                   label = "Synthetic", show_trace = false, n_restarts = 10)

@info "Synthetic patient training complete. $(length(results_synthetic))/$n_syn_train patients converged."

# ──────────────────────────────────────────────────────────────────────────────
# STEP 7: Plot — overlay real vs model-predicted FIIa curves
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "STEP 7: Generating plots..."

# --- Plot A: Real patient FIIa curves ---
@info "Plotting real patient FIIa curves..."
p_real = plot(
    title  = "BST Model Fit — Real Patients (Visit 1)",
    xlabel = "Time (min)",
    ylabel = "FIIa (nM)",
    legend = :topright,
    size   = (800, 500),
    dpi    = 300,
)

color_palette = palette(:tab10)
for (k, idx) in enumerate(sort(collect(keys(results_real))))
    res = results_real[idx]
    T_sim = res.T
    FIIa_sim = res.Xm[:, 9]  # FIIa is the 9th dynamic state
    plot!(p_real, T_sim, FIIa_sim,
          label = "Patient $idx (err=$(round(res.error; digits=2)))",
          lw = 2, color = color_palette[k])
end

path_real_plot = joinpath(_PATH_TO_FIG, "bst_fit_real_patients.pdf")
savefig(p_real, path_real_plot)
@info "Saved real patient plot: $path_real_plot"

# --- Plot B: Synthetic patient FIIa curves ---
@info "Plotting synthetic patient FIIa curves..."
p_syn = plot(
    title  = "BST Model Fit — Synthetic Patients (Visit 1)",
    xlabel = "Time (min)",
    ylabel = "FIIa (nM)",
    legend = :topright,
    size   = (800, 500),
    dpi    = 300,
)

for (k, idx) in enumerate(sort(collect(keys(results_synthetic))))
    res = results_synthetic[idx]
    T_sim = res.T
    FIIa_sim = res.Xm[:, 9]
    plot!(p_syn, T_sim, FIIa_sim,
          label = "Synth $idx (err=$(round(res.error; digits=2)))",
          lw = 2, color = color_palette[k])
end

path_syn_plot = joinpath(_PATH_TO_FIG, "bst_fit_synthetic_patients.pdf")
savefig(p_syn, path_syn_plot)
@info "Saved synthetic patient plot: $path_syn_plot"

# --- Plot C: Combined comparison (real vs synthetic side by side) ---
@info "Plotting combined comparison..."
p_combined = plot(p_real, p_syn, layout = (1, 2), size = (1600, 500), dpi = 300)
path_combined_plot = joinpath(_PATH_TO_FIG, "bst_fit_comparison.pdf")
savefig(p_combined, path_combined_plot)
@info "Saved combined comparison plot: $path_combined_plot"

# ──────────────────────────────────────────────────────────────────────────────
# STEP 8: Save results to JLD2
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "STEP 8: Saving results..."

path_results = joinpath(_PATH_TO_DATA, "bst_training_results.jld2")
@info "Saving training results to: $path_results"

jldsave(path_results;
    results_real       = results_real,
    results_synthetic  = results_synthetic,
    training_df_real   = training_df_real,
    training_df_synthetic = training_df_synthetic,
    n_real_trained     = length(results_real),
    n_synthetic_trained = length(results_synthetic),
)

@info "Results saved successfully."

# ──────────────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────────────
@info "=" ^ 72
@info "BST Training Pipeline Complete"
@info "  Real patients trained:      $(length(results_real)) / $n_real_train"
@info "  Synthetic patients trained:  $(length(results_synthetic)) / $n_syn_train"

# print error summary for real patients
if !isempty(results_real)
    real_errors = [results_real[k].error for k in keys(results_real)]
    @info "  Real patient errors:  mean=$(round(mean(real_errors); digits=4)), " *
          "min=$(round(minimum(real_errors); digits=4)), max=$(round(maximum(real_errors); digits=4))"
end

# print error summary for synthetic patients
if !isempty(results_synthetic)
    syn_errors = [results_synthetic[k].error for k in keys(results_synthetic)]
    @info "  Synthetic patient errors:  mean=$(round(mean(syn_errors); digits=4)), " *
          "min=$(round(minimum(syn_errors); digits=4)), max=$(round(maximum(syn_errors); digits=4))"
end

@info "Plots saved to: $_PATH_TO_FIG"
@info "Results saved to: $path_results"
@info "Done."

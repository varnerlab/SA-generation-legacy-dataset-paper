# ──────────────────────────────────────────────────────────────────────────────
# test_sa_helper.jl
#
# Equivalence check (not a package test — this repo has no test suite):
# confirms that `sa_generate_from_matrix` (code/src/Patient.jl) reproduces
# the canonical seed-42 cohort produced inline by run_full_longitudinal.jl
# (data/synthetic_full_longitudinal.csv), byte-for-byte up to floating point
# tolerance. run_full_longitudinal.jl itself is NOT modified or re-run here;
# we rebuild its X_concat input from the saved pipeline memory
# (full_longitudinal_memory.jld2) and feed it to the new helper.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

@info "Loading canonical pipeline memory and cleaned data …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
d = JLD2.load(mem_path)
df_clean = d["df_clean"]
kept_cols = d["kept_cols"]
n_assays = d["n_assays"]
complete_subjects = d["complete_subjects"]

# ── Rebuild X_concat exactly as run_full_longitudinal.jl does (lines 53-71) ──
K = length(complete_subjects)
n_visits = 3
d_concat = n_assays * n_visits
X_concat = Matrix{Float64}(undef, K, d_concat)

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
@info "  Rebuilt X_concat: $(size(X_concat))"

# ── Load the canonical output ────────────────────────────────────────────────
canon = CSV.read(joinpath(_PATH_TO_DATA, "synthetic_full_longitudinal.csv"), DataFrame)

# ── Call the new helper ──────────────────────────────────────────────────────
@info "Calling sa_generate_from_matrix(...; seed=42) …"
synth_df, mem = sa_generate_from_matrix(X_concat, kept_cols, n_assays, 100, 2000; seed=42)

# align column order + row order (SyntheticID, Visit) then compare numeric block
sort!(synth_df, [:SyntheticID, :Visit])
sort!(canon, [:SyntheticID, :Visit])

@assert names(synth_df) == vcat(string.(kept_cols), ["SyntheticID", "Visit"]) ||
        Set(names(synth_df)) == Set(names(canon)) "column sets differ"

A = Matrix(synth_df[!, kept_cols])
B = Matrix(canon[!, kept_cols])
max_rel_err = maximum(abs.(A .- B) ./ (abs.(B) .+ 1e-8))
@info "  Max relative error vs canonical cohort: $max_rel_err"

@assert isapprox(A, B; rtol=1e-6, atol=1e-8) "helper output diverges from canonical cohort (max rel err = $max_rel_err)"
@assert isapprox(mem.β_star, 2.94; atol=0.05) "β* = $(mem.β_star) does not match canonical β* ≈ 2.94"

println("helper reproduces canonical cohort ✓  (β* = ", round(mem.β_star, digits=3), ")")

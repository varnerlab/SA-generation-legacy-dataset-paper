# ──────────────────────────────────────────────────────────────────────────────
# test_build_concat.jl
#
# Equivalence check (not a package test — this repo has no test suite):
# confirms that `build_concat_matrix` (code/src/Patient.jl) reproduces the
# canonical [V1|V2|V3] concatenated matrix built inline by
# run_full_longitudinal.jl:57-65. run_full_longitudinal.jl itself is NOT
# modified or re-run here; we rebuild its inline X_concat construction
# locally from the saved pipeline memory (full_longitudinal_memory.jld2)
# and compare it against the new helper's output.
# ──────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "Include.jl"))

@info "Loading canonical pipeline memory and cleaned data …"
mem_path = joinpath(_PATH_TO_DATA, "full_longitudinal_memory.jld2")
d = JLD2.load(mem_path)
df_clean = d["df_clean"]
kept_cols = d["kept_cols"]
n_assays = d["n_assays"]
complete_subjects = d["complete_subjects"]

# ── Rebuild X_concat exactly as run_full_longitudinal.jl does (lines 57-65) ──
K = length(complete_subjects)
n_visits = 3
d_concat = n_assays * n_visits
X_ref = Matrix{Float64}(undef, K, d_concat)

for (i, s) in enumerate(complete_subjects)
    for v in 1:n_visits
        mask = (df_clean.SubjectID .== s) .& (df_clean.Visit .== v)
        row = findfirst(mask)
        offset = (v - 1) * n_assays
        for (j, col) in enumerate(kept_cols)
            X_ref[i, offset + j] = df_clean[row, col]
        end
    end
end
@info "  Rebuilt X_ref (inline): $(size(X_ref))"

# ── Call the new helper ──────────────────────────────────────────────────────
X_new = build_concat_matrix(df_clean, :SubjectID, complete_subjects, kept_cols, n_assays)
@info "  build_concat_matrix output: $(size(X_new))"

@assert size(X_new) == (length(complete_subjects), n_assays * 3)
@assert isapprox(X_new, X_ref; atol=1e-12) "build_concat_matrix diverges from the canonical inline layout"

println("build_concat_matrix reproduces canonical concat ✓")

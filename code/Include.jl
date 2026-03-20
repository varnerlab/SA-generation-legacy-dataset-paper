# ──────────────────────────────────────────────────────────────────────────────
# Include.jl — SA Synthetic Patient Generation Study
# Sets up the environment and loads all source modules
# ──────────────────────────────────────────────────────────────────────────────

# setup paths
const _ROOT = @__DIR__
const _PATH_TO_SRC = joinpath(_ROOT, "src")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIG = joinpath(_ROOT, "figs")
const _PATH_TO_ORIGINAL_DATA = joinpath(_ROOT, "original_legacy_data")

# activate environment (Manifest.toml pins exact dependency versions)
using Pkg
Pkg.activate(_ROOT)
if !isfile(joinpath(_ROOT, "Manifest.toml"))
    @warn "Manifest.toml not found — resolving and instantiating dependencies. " *
          "Pin versions by committing the generated Manifest.toml."
    Pkg.resolve()
    Pkg.instantiate()
end

# load required packages
using XLSX
using DataFrames
using CSV
using Statistics
using LinearAlgebra
using MultivariateStats
using NNlib
using Random
using Distributions
using Plots
using StatsPlots
using Colors
using PrettyTables
using StatsBase
using JLD2

# NOTE: No global Random.seed!() here. Each experiment script and sampling
# function manages its own RNG seed explicitly via seed= arguments.

# include the source modules
include(joinpath(_PATH_TO_SRC, "Compute.jl"))
include(joinpath(_PATH_TO_SRC, "Utilities.jl"))
include(joinpath(_PATH_TO_SRC, "Patient.jl"))

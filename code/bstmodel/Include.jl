# ──────────────────────────────────────────────────────────────────────────────
# Include.jl — BST Coagulation Model Training Pipeline
# Sets up paths, activates the project environment, and loads all dependencies
# ──────────────────────────────────────────────────────────────────────────────

# setup paths
const _ROOT = @__DIR__
const _PATH_TO_SRC = joinpath(_ROOT, "src")
const _PATH_TO_MODEL = joinpath(_ROOT, "model")
const _PATH_TO_DATA = joinpath(dirname(_ROOT), "data")
const _PATH_TO_FIG = joinpath(dirname(_ROOT), "figs")

# activate the parent project environment (shared Project.toml at code/ level)
using Pkg
Pkg.activate(dirname(_ROOT))

# check for Manifest.toml — resolve and instantiate if missing
if !isfile(joinpath(dirname(_ROOT), "Manifest.toml"))
    @warn "Manifest.toml not found — resolving and instantiating dependencies. " *
          "Pin versions by committing the generated Manifest.toml."
    Pkg.resolve()
    Pkg.instantiate()
end

# load required packages ──────────────────────────────────────────────────────
@info "Loading packages for BST coagulation model training..."

# data handling
using DataFrames
using CSV
using JLD2

# math and statistics
using Statistics
using LinearAlgebra
using Random
using Distributions
using NumericalIntegration

# optimization
using Optim

# plotting
using Plots
using Colors

# BST model toolkit
using BSTModelKit

@info "All packages loaded successfully."

# include bstmodel source modules ─────────────────────────────────────────────
@info "Loading BST model source files..."
include(joinpath(_PATH_TO_SRC, "Utility.jl"))
include(joinpath(_PATH_TO_SRC, "Compute.jl"))
include(joinpath(_PATH_TO_SRC, "Learn.jl"))
include(joinpath(_PATH_TO_SRC, "Files.jl"))
@info "BST model source files loaded."

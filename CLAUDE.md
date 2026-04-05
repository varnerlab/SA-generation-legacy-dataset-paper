# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research paper and codebase for validated synthetic patient generation using Stochastic Attention (SA) applied to a small longitudinal pregnancy coagulation dataset (N=23 patients, 72 assays, 3 visits). The SA algorithm uses modern Hopfield network energy with Unadjusted Langevin sampling to generate biologically plausible synthetic patients.

## Common Commands

### LaTeX Paper
```bash
cd paper && make          # Build main.pdf (pdflatex → bibtex → pdflatex x2)
cd paper && make clean    # Remove build artifacts
```

### Julia Code
All Julia scripts assume the working directory is `code/` and begin with `include("Include.jl")` to load the environment.
```bash
cd code && julia Include.jl                              # Load environment / precompile
cd code && julia experiments/<script_name>.jl             # Run a specific experiment
cd code && julia experiments/run_full_longitudinal.jl     # Generate 100 synthetic patients
cd code && julia experiments/run_conditioned_generation.jl # Generate condition-specific cohorts
```

### Python Baseline
```bash
cd code && python3 experiments/run_ctgan_baseline.py      # CTGAN neural baseline
```

## Architecture

### Core Julia Modules (`code/src/`)
- **Compute.jl** — SA sampling engine: `sample()` (ULA on Hopfield energy), `weighted_sample()` (multiplicity-weighted for conditional generation), `hopfield_energy()`, `attention_entropy()`
- **Patient.jl** — Data I/O: `load_legacy_data()` reads the Excel source, cleans and structures the 72-assay panel across 3 visits
- **Utilities.jl** — Helpers for PCA, similarity metrics

### Entry Point
- **Include.jl** — Loads `Project.toml` environment, imports all packages, includes all `src/` modules. Every experiment script starts with `include("Include.jl")`.

### Experiment Scripts (`code/experiments/`)
Scripts fall into categories:
- **Generation:** `run_full_longitudinal.jl`, `run_conditioned_generation.jl`
- **Validation (4 levels):** `validate_biological_constraints.jl`, `validate_cross_visit_covariance.jl`, `validate_conditioned_generation.jl`, `validate_mechanistic_plausibility.jl`
- **Baselines:** `paper_sa_vs_mvn.jl`, `run_ctgan_baseline.py`
- **Figures:** `paper_pca_by_visit.jl`, `paper_correlation_heatmaps.jl`, `regen_*_figures.jl`
- **Mechanistic model:** `calibrate_hockin_mann.jl`, `validate_hockin_mann*.jl`
- **Sensitivity:** `run_beta_sweep.jl`

### BST Model (`code/bstmodel/`)
Separate mechanistic model package with its own `Include.jl` and training scripts.

### Paper (`paper/`)
LaTeX manuscript with modular sections in `paper/sections/`. ArXiv mirror in `arxiv/`.

## Key Algorithm Parameters
- **beta (inverse temperature):** phase transition at ~2.94; controls interpolation vs memorization
- **alpha (step size):** 0.1, controls ULA dynamics
- **d_PCA:** 18 principal components used for memory representation
- **T:** 1000+ ULA iterations per synthetic sample

## Data
- Source: `code/original_legacy_data/20200604_R61_Legacy_Biochemical_Measurements.xlsx`
- Processed outputs in `code/data/` (CSV for tables, JLD2 for serialized Julia objects)
- 23 patients across 3 visits: V1 (pre-pregnancy), V2 (end 1st trimester), V3 (mid-3rd trimester)
- Patient subgroups by outcome: Uncomplicated (18), PCOS (3), Developed PE (5)

## Conventions
- Julia scripts use `include("Include.jl")` rather than module imports; always run from `code/`
- Generated figures go to `code/figs/` as PDF/PNG
- Manifest.toml is gitignored; `Project.toml` is the dependency spec

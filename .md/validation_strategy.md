# Synthetic Patient Validation Strategy

Multi-level validation framework for SA-generated synthetic patients. Each level tests a progressively stronger property. Separate scripts per level.

## Level 1: Marginal Plausibility (biological constraints)
**Script:** `code/experiments/validate_biological_constraints.jl`

Do known physiological relationships hold in synthetic patients?
- AT inversely correlates with thrombin generation (peak, ETP)
- Factor VIII deficiency → reduced peak/ETP
- Pregnancy progression: fibrinogen, FVIII, vWF increase V1→V2→V3
- Per-patient longitudinal monotonicity patterns match real data
- Per-feature distributions match real data (already done in existing scripts)

**Why it matters:** Basic sanity check. If synthetic patients violate known biology, nothing else matters.

## Level 2: Joint Structure (cross-visit covariance)
**Script:** `code/experiments/validate_cross_visit_covariance.jl`

Does the V1↔V2↔V3 correlation structure hold in synthetic vs real data?
- Compare full cross-visit correlation matrices (real vs SA vs MVN)
- Cross-visit feature correlations: does V1 peak predict V3 peak?
- Off-diagonal block structure of the [V1|V2|V3] covariance
- MVN estimated a rank-22 covariance from n=23 patients — should fail here

**Why it matters:** SA preserves the data manifold geometry directly. MVN's regularized singular covariance cannot faithfully represent the joint cross-visit structure. This is the key structural argument for SA over MVN.

## Level 3: Conditional Structure (conditioned generation)
**Script:** `code/experiments/validate_conditioned_generation.jl`

Can SA amplify rare subpopulations while preserving condition-specific signatures?
- Generate condition-specific cohorts: Healthy (14→100), PCOS (3→100), PE (5→100)
- Check condition-specific means match real condition-specific means
- PCA separation: do generated cohorts separate as expected?
- Cross-cohort feature differences match real between-group differences
- MVN fundamentally cannot condition on n=3 PCOS patients

**Why it matters:** Clinical utility of synthetic data requires generating specific patient populations. SA's multiplicity-weighted sampling enables this; MVN has no mechanism for it.

## Level 4: Mechanistic Consistency (Hockin/Mann BZ2012)
**Script:** `code/experiments/validate_mechanistic_plausibility.jl` (main validation)
**Supporting:** `code/experiments/calibrate_hockin_mann.jl` (calibration), `code/experiments/validate_hockin_mann.jl` (exploratory), `code/experiments/validate_hockin_mann_analysis.jl` (IIa-only investigation)

Do synthetic patients' factor-to-TGA relationships fall within the same mechanistically-predicted range as real patients?
- Hockin/Mann BZ2012 coagulation model (58 species, 64 rate constants) via HockinMannModel.jl package
- Two conditions per patient: TF-only (TM=0) and TF+TM (TM=1nM)
- Global calibration: 5 rate constants fit on V1 real patients (prothrombinase_kcat, intrinsic/extrinsic_xase_kcat, PC_activation_kcat, mIIa_conversion_k), remaining 59 at literature values
- Runs ALL real patients (23 × 3 visits) + ALL synthetic patients, both conditions
- Population-level check: synthetic patients produce the same predicted-vs-measured cloud as real patients
- Cloud overlap metric: % of synthetic pred/meas ratios within real 5th–95th percentile
- Rank correlations (Spearman ρ) by source, visit, condition, and feature

**Calibration results (from calibrate_hockin_mann.jl):**
- prothrombinase_kcat: 63.5 → 21.94 s⁻¹ (0.345x)
- intrinsic_xase_kcat: 8.2 → 0.1693 s⁻¹ (0.021x)
- extrinsic_xase_kcat: 6.0 → 93.21 s⁻¹ (15.5x)
- PC_activation_kcat: 0.41 → 0.8896 s⁻¹ (2.17x)
- mIIa_conversion_k: 2.3e8 → 1.153e9 M⁻¹s⁻¹ (5.01x)

**Current results:**
- TF-only calibration generalizes across visits (pred/meas ≈ 1.0x for V2/V3)
- TF+TM calibration does not generalize (overcorrects PC pathway)
- Rank correlations are weak — model captures population-level but not inter-patient variability
- ETP is the best-predicted feature (Spearman ρ ≈ 0.6-0.8 for TF-only)

**Plots generated:**
- `validate_mechanistic_tf_only.pdf` / `tf_tm.pdf` — predicted vs measured scatter (real + synth, colored by visit)
- `validate_mechanistic_ratios_tf_only.pdf` / `tf_tm.pdf` — pred/meas ratio distributions
- `validate_mechanistic_rankcorr_tf_only.pdf` / `tf_tm.pdf` — rank correlation bar charts
- `validate_mechanistic_generalization_tf_only.pdf` / `tf_tm.pdf` — visit transfer boxplots

**Framing:** Not a patient-level predictor, but an independent mechanistic plausibility check. Synthetic patients land in the same cloud as real patients → the SA-generated factor-to-TGA mapping is biologically reasonable. Note: BST model was considered but not used; HockinMann BZ2012 is the mechanistic validator.

## Paper Narrative Arc

SA-generated patients pass validation at every level:
1. **Marginals** — individual features are biologically plausible
2. **Joint structure** — cross-visit correlations are preserved (MVN fails here)
3. **Conditional structure** — rare subpopulations can be amplified (MVN cannot do this)
4. **Mechanistic consistency** — an independent ODE model confirms factor-to-TGA plausibility

MVN matches level 1 by construction (moment matching) but fails at levels 2-4.

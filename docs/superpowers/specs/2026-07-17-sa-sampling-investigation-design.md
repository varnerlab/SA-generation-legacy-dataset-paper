# Design: SA Sampling-Correctness Investigation (diagnostic + MCMC ladder)

**Date:** 2026-07-17
**Manuscript:** Validated Synthetic Patient Generation for Small Longitudinal Cohorts (npj SBA, major revision)
**Origin:** Reviewer R1.5 (privacy). The E3 experiment found a membership-inference leak at the published operating point (MIA AUC ≈ 0.97 across 40 splits; DCR synth→real 13.78 < real→real 15.86). The generator initializes each ULA chain **at** a stored memory and returns a single endpoint, so the memorization may be a sampling/initialization artifact rather than a property of the target distribution. This investigation disentangles the two.
**Companion docs:** `docs/superpowers/specs/2026-07-16-sa-paper-major-revision-design.md` (parent revision), `docs/superpowers/plans/2026-07-16-sa-paper-major-revision-implementation.md`.

---

## 1. Goal and hypothesis

Determine how much of the E3 memorization is a **sampling artifact** (fixable by sampling the Hopfield target correctly) versus an **intrinsic property of π(ξ) ∝ exp(−βE) at β\***.

The current per-sample scheme (`Patient.jl:394-404`, `sa_generate_from_matrix`, `run_full_longitudinal.jl`) has three separable, suspect features:
1. **Anchor initialization** — `ξ₀ = X̂[:,k0] + 0.01·randn`, i.e. each chain starts *at* a stored pattern.
2. **Single-endpoint extraction** — keep only `Ξ[end,:]`; no burn-in discard, no pooling of post-burn-in draws (100 chains × 1 sample).
3. **ULA discretization bias** — no Metropolis correction.

We attribute the leak across these via a four-rung ladder, and separately verify with a mixing diagnostic whether T=2000 from an anchor even reaches the stationary distribution.

## 2. Conditions held constant (so only the sampler varies)

- **β = β\* = 2.94** (fixed across all rungs — co-varying β would confound the sampling attribution; the β–privacy trade-off is a separate axis, out of scope here).
- **Memory** X̂ (18-dim unit-normalized PCA memory of the 23 reals), from `full_longitudinal_memory.jld2`.
- **Decode** unchanged: each sampled state → unit direction (`ξ/‖ξ‖`) → × a magnitude **drawn from the 23 empirical PCA-norms** → inverse PCA → destandardize (`decode_sample`). *The magnitude draw reuses a real patient's radial coordinate and is a secondary, orthogonal memorization channel; it is held fixed here (1-D vs the 18-D direction) and explicitly flagged as a limitation, not silently ignored.*
- **N = 100** synthetic patients per rung.
- **Evaluation harness** identical across rungs (§5).
- **Reproducible seeds** throughout.

## 3. Deliverable 1 — mixing diagnostic

`code/experiments/run_mcmc_mixing_diagnostic.jl`, building on the existing `run_sampling_diagnostics.jl` (which already computes ULA/MALA energy traces, integrated autocorrelation time τ_int, ESS, MALA acceptance). Additions:

- **Chains from uniform-sphere inits:** `ξ₀ = randn(d); ξ₀ ./= norm(ξ₀)` (over-dispersed relative to the 23 memories — the correct start for a convergence check).
- **Gelman–Rubin R̂** across the sphere-init chains, computed on the energy scalar `E(ξ_t)` (the scalar the existing diagnostic already traces), reported vs. step.
- **Anchor-distance trace:** distance from `ξ_t` to its nearest stored memory (and to its own init anchor) vs. t — does a chain leave its starting basin?
- **Init-dependence comparison:** endpoint distribution from **anchor-init** vs **sphere-init** at the same T — if they differ, T=2000 has not reached the stationary law (init still matters).
- **Outputs:** `code/data/mcmc_mixing_diagnostic.csv` (R̂ vs step, τ_int, ESS, anchor-distance summaries, MALA acceptance) + a diagnostic figure (energy traces with burn-in marker, R̂ curve, anchor-distance curve) → `code/figs/mcmc_mixing_diagnostic.pdf`.
- **Feeds the ladder:** the burn-in length **B** (first step where R̂ < 1.05 and energy plateaus) and thinning lag **τ ≈ ⌈τ_int⌉** used by V2/V3 come from this diagnostic — pooling parameters are principled, not arbitrary.

## 4. Deliverable 2 — the ladder

`code/experiments/run_mcmc_ladder.jl`. Four N=100 cohorts, each decoded and evaluated through the identical harness. Only the **direction sampler** differs:

| Rung | Init | Sampler | Extraction |
|---|---|---|---|
| **V0** (published baseline) | at a stored memory (`X̂[:,k0]+0.01·randn`) | ULA (`sample`) | endpoint; 100 chains × 1 |
| **V1** | uniform sphere | ULA (`sample`) | endpoint; 100 chains × 1 |
| **V2** | uniform sphere | ULA (`sample`) | burn-in **B** discard + thin at **τ**, pool; 20 chains × 5 draws |
| **V3** | uniform sphere | MALA (`mala_sample`) | burn-in **B** discard + thin at **τ**, pool; 20 chains × 5 draws |

- **V0** reuses the frozen published scheme (equivalently, `sa_generate_from_matrix` on all 23) — no new generation logic, it is the reference.
- **V1** changes *only* the initialization (isolates the init effect, V0→V1).
- **V2** adds burn-in + thinned multi-chain pooling (isolates single-endpoint extraction, V1→V2). Default **20 chains × 5 thinned post-burn-in draws = 100**; each chain runs `B + 5τ` steps. *(Adjustable — flagged open default.)*
- **V3** swaps ULA for the Metropolis-corrected `mala_sample` (isolates ULA discretization bias, V2→V3).
- All rungs seeded and reproducible; every rung uses the same β\*, the same magnitude-draw decode, and the same N=100.

## 5. Metrics and evaluation harness (identical per rung)

Reuse existing, already-reviewed harness code:
- **Fidelity:** pooled median MRE + bootstrap CI (`feature_summary_comparison` + `paper_summary_table.jl` bootstrap pattern); cross-visit correlation Frobenius (`safe_cor`, `norm(C−C_real)/norm(C_real)`); mechanistic cloud overlap + KS (`mechanistic_eval.jl` `bz2012_ratios`/`overlap_ks`, TF-only); novelty (`sample_novelty`).
- **Privacy — DCR (all four rungs):** per-synthetic min distance to the 23 reals vs. real→real leave-one-out, Euclidean in the canonical standardized space (the E3 DCR computation). Cheap, split-independent — carries the per-rung comparison.
- **Privacy — repeated-split MIA (V0 and the best-DCR rung only):** the reviewed E3b machinery (`run_privacy_split_sensitivity.jl`), reporting the AUC distribution. Running the full retrain-based MIA on every rung is expensive; DCR ranks the rungs, MIA confirms the endpoints. *(Adjustable — flagged open default.)*

**Output:** `code/data/mcmc_ladder_results.csv` — one row per rung with all fidelity + DCR columns, plus a companion `mcmc_ladder_mia.csv` for the V0/best-rung MIA distributions. A rung × metric summary table for the supplement.

## 6. Decision logic (what each outcome means)

Reading the ladder + diagnostic together:
- **DCR/MIA improves at V0→V1** → memorization is mainly an **initialization** artifact (anchored chains never mixed).
- **improves at V1→V2** → the **single-endpoint** extraction was the driver.
- **improves at V2→V3** → **ULA discretization** bias.
- **stays high through V3 with fidelity intact** → memorization is a **genuine property of π at β\*** (sharp modes at the stored patterns); privacy then requires lowering β (a real trade-off), and the honest R1.5 answer is "fundamental at this operating point."
- The diagnostic's init-dependence result independently corroborates: if anchor-init and sphere-init endpoints differ at T=2000, the published sampler was not sampling π (supports the "artifact" reading).

**Fidelity is tracked at every rung.** A rung that fixes privacy but collapses MRE/mechanistic fidelity is a **trade-off**, not a free lunch — the table quantifies the cost.

## 7. Scope and framing

- The ladder variants V1–V3 are **diagnostic probes to characterize the published method**, NOT proposed replacements. The published N=100 cohort and every existing result stay frozen. This keeps the work within the parent revision's "no new SA variant" scope guard (probes ≠ redesign).
- **If** a rung is strictly better (privacy fixed, fidelity retained), that is surfaced to the author as a finding and a *possible* recommended-sampling note — not silently swapped in.
- **Where it lands:** the R1.5 response (privacy — the leak + whether it is sampler-driven), a new **supplement section** (mixing diagnostic + ladder table + figures), and a short Methods/Discussion note on sampling. Not a new main-text headline method.

## 8. Files and reuse map

**New:**
- `code/experiments/run_mcmc_mixing_diagnostic.jl` → `mcmc_mixing_diagnostic.csv`, `figs/mcmc_mixing_diagnostic.pdf`
- `code/experiments/run_mcmc_ladder.jl` → `mcmc_ladder_results.csv`, `mcmc_ladder_mia.csv`
- (supplement figure) `code/experiments/fig_mcmc_ladder.jl` → `figs/mcmc_ladder.pdf`

**Reuse (do not reinvent):** `Compute.sample` (ULA), `Compute.mala_sample` (MALA), `find_entropy_inflection` (β\*), `decode_sample`, `sample_novelty`, `mechanistic_eval.jl` (`bz2012_ratios`/`overlap_ks`), the E3 DCR computation, the E3b repeated-split MIA (`run_privacy_split_sensitivity.jl`), `run_sampling_diagnostics.jl` (ESS/ACF/MALA-acceptance primitives), and `sa_generate_from_matrix` for the V0 baseline. The recurring `rebuild_concat_matrix`/concatenation helper (currently copy-pasted in E3/E3b) should be **hoisted to `src/Patient.jl`** as part of this work so the ladder reuses it (closes a standing Minor review finding).

## 9. Open defaults (chosen; adjustable at spec review)

1. **Privacy metric coverage:** DCR on all four rungs + repeated-split MIA on V0 and the best-DCR rung only. *(Alternative: full MIA on every rung — more retrains.)*
2. **V2/V3 pooling allocation:** 20 chains × 5 thinned post-burn-in draws = 100. *(Alternative: e.g. 10×10, or more chains for better between-chain independence.)*
3. **R̂ scalar:** computed on the energy `E(ξ_t)`. *(Alternative: per-coordinate R̂ or distance-to-nearest-memory.)*

## 10. Success criteria

- Mixing diagnostic reports R̂, τ_int, ESS, MALA acceptance, and a clear anchor-vs-sphere init-dependence verdict; B and τ derived and passed to the ladder.
- Ladder produces a complete rung × {MRE, cross-visit Frobenius, mechanistic overlap/KS, novelty, DCR} table + MIA on V0 and the best rung.
- The attribution question is answered: the leak is assigned to init / pooling / discretization, or declared intrinsic to π at β\* — with the fidelity cost of any privacy improvement quantified.
- Published cohort untouched; results framed as diagnostic probes; supplement section + R1.5 response text drafted.

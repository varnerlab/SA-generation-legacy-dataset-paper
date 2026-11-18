# Round-2 fact check and reference check (2026-11-18)

Two independent checks ran on each question. Codex (`codex exec`, read-only) checked the
round-2 additions and the letter against `code/data/` and `code/experiments/`. Codex and three
Claude subagents each checked all 40 cited references. The raw reports are in the session
scratchpad (`refcheck/codex-references.md`, `refcheck/claude-refs-{A,B,C}.md`,
`factcheck/codex-factcheck.md`). This file keeps the findings.

## Fixed in this session (uncommitted)

**Round-2 text and letter**

- Language pass on every `\rthree{}` block and the letter. Methods round-2 paragraph split into
  three paragraphs, each opening with its question. "The kernel-density baseline above"
  replaced by naming the baseline. "Constraint term" renamed "positivity penalty" everywhere.
  Letter: LaTeX macro names, "below", and "bounded statement" removed; date set to
  18 November 2026.
- Langevin cost claim corrected. The text said "the cost of Langevin sampling does not depend on
  that pass rate". In `constrained_positivity_demo.csv` the step count is 2,000 at λ=0 and
  959,108 at λ=100 (T = ⌈20/α⌉, α ∝ 1/λ). New text: rejection kept 50–100% of draws, while the
  number of Langevin steps grew in proportion to λ.
- Codex fact check: every number in Supplementary Tables 17 and 18 and in the prose matches the
  CSVs. Nine wording problems, all fixed:
  1. "a non-quadratic term does not [keep a Gaussian mixture]" scoped to this penalty.
  2. "Langevin agreed with exact sampling" changed to "did not differ detectably" (KS tests on six
     hormone marginals cannot show equivalence).
  3. Table 18: "independent draw" changed to "separate cohort"; the caption now says both KS
     samplers used the same magnitudes.
  4. "Removing the negative tail raised the hormone means" removed (signed means were not saved).
  5. "every negative value" changed to "every negative hormone value" (PAI values can be negative).
  6. Step rule stated as ⌈20/α⌉ steps.
  7. 15.9 named as the median nearest-neighbor distance between real patients.
  8. "results in the main text are unchanged" changed to "the reported synthetic cohort and its
     results are unchanged".
  9. Data Availability (not fixed by text; see open item 10).

**Bib metadata (Crossref-confirmed)**

Azziz2004 Timothy J. Key; Hockin2002 Matthew F. Hockin; Ribba2012 Linda Dainese and
Čajavec-Bernard; Bernstein2016 adds Carole A. McBride; Grouzi2022 issue 1, article 60;
Abalos2014 issue S1; McLean2012TGA hyphenated title and Kathleen E.; Ramsauer2021 Günter;
Goodfellow2014 Ian J.; Toy2023WSJ URL.

## Pre-round-2 text (fixed later on 2026-11-18, with the author's approval)

The author approved fixing items 1-9 and pushing the repository (item 10). All fixes are in
`\rthree{}` and are green in the marked manuscript. Decisions: the storage temperature stays at
-80 C (Hale2012 is no longer cited for it); the Kuhnel2024 time-trend clause was dropped rather
than verified; the 2-8% preeclampsia figure became 4.6% of deliveries, citing the new Abalos2013
entry (Crossref and PubMed 23746796 checked). The MVN fix rests on a scratch rerun of
`validate_cross_visit_covariance.jl` and `regen_pca_figure.jl` (outputs redirected to the
scratchpad, committed figures untouched): MVN closer to the real correlations in every block;
MVN eigenvalues beyond component 22 hold 0.3% of standardized variance; in the PC1-PC2 plots
the root-mean-square distance from the real centroid is real 3.3/4.0/4.0, MVN 2.8/3.7/3.9 and
SA 2.2/2.7/3.5 at Visits 1/2/3. The letter's opening note discloses the MVN correction and
"13 citations". The table below is the original record.

| # | Location | Problem | Evidence | Suggested fix |
|---|---|---|---|---|
| 1 | methods.tex:280; results.tex:104, 116; floats.tex:123-124; supplementary.tex:483, 568 | Text says the MVN baseline used Ledoit-Wolf shrinkage. No script does. | Every MVN fit is `cov(X) + 1e-6*I`: `paper_sa_vs_mvn.jl:171-173`, `validate_cross_visit_covariance.jl:78-80`, `regen_pca_figure.jl:89`, `run_convex_hull_analysis.jl:364-366`. | Describe the 1e-6 ridge and drop the Ledoit-Wolf attribution. The comparison claim is reversed, not just mislabeled: the supplementary correlation figure titles (same file as `code/figs/`) show SA-Real Frobenius distance 38.66 and MVN-Real 20.74, so MVN reproduced the real correlation matrix more closely. A Python re-draw of the unshrunk MVN gave off-diagonal-block error 13.2-22.1 over 200 seeds vs SA 31.8, and mean off-diagonal abs. correlation 0.216 (real 0.199, SA 0.185), so MVN did not shrink correlations toward zero. The ridge is at most 0.14% of any non-constant feature's variance, so "spurious variance in 194 dimensions" is not supported. Rewrite results.tex:102-120, floats.tex:122-126, discussion.tex:7-8, supplementary.tex:483-485 and 567-569, methods.tex:279-281; check the PCA dispersion and convex-hull observations before keeping them. |
| 2 | introduction.tex:5-8 | "disorders related to pregnancy ... were the leading underlying cause of death" | Chen2025 abstract (PMID 40202759): "cardiovascular disease was the leading cause of the overall pregnancy-related deaths". The 17.4% (1,094/6,283) is correct. | "...disorders related to pregnancy, including hemorrhage and venous complications, accounted for 17.4% of 6,283 pregnancy-related deaths". |
| 3 | introduction.tex:25 | "complicates 2-8% of pregnancies worldwide" cites Abalos2014 | Abalos2014 reports 2.16% in a 29-country facility survey. The range comes from Abalos et al. 2013, Eur J Obstet Gynecol Reprod Biol 170:1-7, doi:10.1016/j.ejogrb.2013.05.005 (PMID 23746796): 4.6%, uncertainty range 2.7-8.2%. | Cite the 2013 review for the range, or state the 2.16% survey estimate. |
| 4 | introduction.tex:29-33 | Bravo2012 and BZ2012 roles swapped | Bravo2012 is the APC-mediated FVa inactivation model (no thrombomodulin module). BZ2012 adds thrombomodulin-protein C activation and the Bravo module to Hockin-Mann. | "...through a module for APC-mediated Factor Va inactivation (Bravo2012) to the BZ2012 model, which adds thrombomodulin-dependent protein C activation and that module to the Hockin-Mann framework (BrummelZiedins2012)." |
| 5 | introduction.tex:14-16 | Feedback "via Factors V, VIII, and XI" and the thrombomodulin-protein C pathway cite Hockin2002 | Hockin2002 names FV, FVIII, FVII and has no thrombomodulin or protein C. Mann2011 covers FV, FVIII, FXI and the TM-PC pathway. | Cite Mann2011 for this sentence. |
| 6 | discussion.tex:72-74 (and :44) | "the same geometry-based operating-point rule" | Varner2026SAProtein defines beta* as the entropy inflection (d2H/dbeta2 = 0) and generates at 2beta*. This paper picks the most negative second derivative against log beta and samples at beta*. Varner2026Multiplicity v2 matches this paper. | "related entropy-based rules", and cite Multiplicity v2 for the exact rule. |
| 7 | discussion.tex:80-82 | Kuhnel2024: "cross-visit correlation error near 0.70" and "weaker time trends were not reproduced" | 0.70 is the full-correlation-matrix error across all variables and visits (both checks agree). At n = 1,312, 22% of synthetic datasets reproduced a significant time trend; at n = 10,000, 100% did (one checker only; confirm in the paper). | "correlation error near 0.70 across all variables and visits", and "weaker time trends were reproduced in only 22% of synthetic datasets at n = 1,312". |
| 8 | methods.tex:15 | Protocol citations | Bernstein2016 and Psoinos2021 describe the cohort and visit design, not plasma processing or coagulation assays. Hale2012 stored samples at -70 C (text says -80 C). McLean2012TGA supports the collection details. | Cite Bernstein2016/Psoinos2021 for the cohort, McLean2012TGA (and Hale2012) for collection and assays. Confirm -80 C. |
| 9 | introduction.tex:20-28, 46-58 | Claims stronger than the cited sources | Hastie2009 shows cross-validation working at n=20, p=500 (section 7.10.3), so "cross-validation is unreliable, and overfitting is nearly inevitable" is overstated. Goodfellow2014, Kingma2014, Xu2019CTGAN introduce methods and do not show small-n failure; DeMitri2025 covers image GANs only. Luan2007/2010 are not patient-specific. Palomba2015 reports elevated PAI-1 and related markers, not "impaired fibrinolysis". Uncited "antithrombin ... decline" conflicts with Hellgren2003 ("mainly unchanged"). Ledoit2004 does not support "introducing bias that distorts the joint distribution". | Soften each claim to what the source shows, or add a source that states it. |
| 10 | main.tex:103 | Data Availability says all intermediate outputs are public | Five round-2 CSVs are not on public `main`: `kde_bandwidth_check.csv`, `constrained_positivity_demo.csv`, `constrained_positivity_agreement.csv`, `log_transform_check.csv`, `real_strict_check.csv`. | Push the repository before submitting. |

Minor: Ramsauer2021 author list differs between arXiv (16 authors, current bib) and the ICLR
2021 record (13 authors); Xu2019CTGAN pages conflict across records (7335-7345 vs 7333-7343);
Varner2023JuliaCon note says "Conference proceedings talk" but no proceedings paper was found.

## Not verified

- Varner2023JuliaCon claim that earlier work sampled each visit with an MVN (the video could not
  be opened; the talk abstract does not say it).
- The definitive ICLR byline for Ramsauer2021 (OpenReview blocked automated access).
- Full texts for Luan2010, Bremme2003, Azziz2004, Ribba2012 (judged from abstracts).
- Peer-reviewed versions of the four arXiv preprints (negative search only).
- Runtime: the Langevin-versus-rejection cost comparison rests on step and draw counts, not
  wall-clock measurements.

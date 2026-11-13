# Response to Reviewers, Round 2

**Manuscript:** Validated Synthetic Patient Generation for Small Longitudinal Cohorts: Coagulation Dynamics Across Pregnancy\
**Journal:** npj Systems Biology and Applications\
**Decision:** Minor revision\
**Date:** 16 November 2026

---

## Note on this document

Each reviewer comment is quoted below. Each comment is followed by our response and the location of the corresponding manuscript change. Text added in this round appears in green in the marked manuscript (`\rthree{...}`). Text added in the first round now appears in black. A `\reviewmode` setting in the preamble changes all colored text to black for the final version.

All new supplementary material is in one new section of the Supplementary Information, "Supplementary Bandwidth and Plausibility Results" (pages 14--17). The section contains the new results and Supplementary Tables 17 and 18. No other part of the supplement changed. The supplement is supplied without color.

We also corrected the Data Availability statement. The de-identified real patient data are in the public project repository, and the statement now says so.

---

# Reviewer 1

> The reviewer thanks the authors for the effort in making the revisions. The authors have answered all my questions and the manuscript is now suitable for publication.

We thank Reviewer 1 for the careful reading in both rounds. The comments improved the privacy analysis, the robustness tests, and the description of the method.

---

# Reviewer 3

We thank Reviewer 3 for the careful reading and for the three remaining comments. The new analyses for comments 1 and 3 are reported in the new Supplement section, "Supplementary Bandwidth and Plausibility Results" (pages 14--17).

### R3.1: Advantage over a Gaussian mixture, and the relation to kernel-density estimation

> The revised results show that simpler approaches, particularly the Gaussian mixture and copula models, achieve comparable or better performance on several metrics. The authors should clarify the specific methodological advantage of stochastic attention, particularly given its exact representation as a finite Gaussian mixture. The relationship between stochastic attention and conventional Gaussian kernel-density estimation deserves further discussion.

*Stochastic attention did not perform better than the fitted Gaussian mixture or the copula on these metrics, and we do not claim that it did.* Its advantages are in how the distribution is built and how it can be extended. The stored profiles define the distribution directly, with no fitting step: there is one mixture component centered on each stored profile, and every component has the same isotropic variance. The inverse temperature is the only parameter, and we chose it from the attention-entropy transition without held-out data. The fitted two-component mixture instead estimated two component means, two diagonal variance vectors, and the mixing weights from the $23$ profiles. The energy form also accepts added terms, which we tested with a positivity term (below and R3.3).

The relation to kernel-density estimation is exact. With equal multiplicity weights, the stochastic-attention distribution is an isotropic Gaussian kernel-density estimate on the $23$ unit-normalized profiles, with kernel width $h=\beta^{-1/2}$. With unequal weights, it is a weighted kernel-density estimate, which is the form we used for conditional generation. Stochastic attention differs from conventional kernel-density estimation in two ways. First, it places the kernels on unit-normalized profiles and draws each profile's magnitude separately. Second, it sets the kernel width from the entropy transition instead of a bandwidth rule.

To test the second difference, we compared the entropy transition with three standard bandwidth rules on the same profiles (new Supplement section, Supplementary Table 17). The entropy transition gave $h=0.58$. Leave-one-out likelihood gave $0.29$, Scott's rule gave $0.21$, and Silverman's rule gave $0.19$. These narrower kernels corresponded to $\beta$ values $4.1$--$9.0$ times $\beta^*$. Profiles drawn with the narrower kernels stayed closer to the stored profiles. Median novelty fell from $0.50$ to $0.21$--$0.34$, and the median distance to the closest real record fell from $13.9$ to $9.9$--$11.9$. Cross-visit correlation error fell from $0.60$ to $0.43$--$0.50$, and median relative error rose from $1.7\%$ to $2.1$--$2.2\%$. Mechanistic overlap was $0.91$--$0.92$ for every width. Of the four widths, the entropy transition gave the widest kernel, the most novel profiles, and the lowest marginal error, but the highest cross-visit correlation error.

The Gaussian-mixture identity also shows how the method extends. Stochastic attention defines an energy on the stored profiles, and the Langevin sampler acts on that energy. With no other terms in the energy, the target is the Gaussian mixture. We used this case in the first round to check the Langevin sampler against exact draws. A constraint term $\lambda C$ can be added to the energy, and the same sampler applies. A quadratic term keeps the target a Gaussian mixture, but a non-quadratic term does not. We tested a non-quadratic positivity term (R3.3). With this term, Langevin sampling agreed with exact rejection sampling from the same target at every constraint weight. Rejection was cheaper in this example because about half of the draws passed. Rejection slows as fewer draws pass, but the cost of Langevin sampling does not depend on that pass rate.

**Changes:** Methods, "Multiplicity-Weighted Stochastic Attention" (one sentence stating the kernel-density form and kernel width); Methods, "Baseline and Validation Framework" (new paragraph describing the bandwidth comparison, the real-patient check, the positivity term, and the log representations); Results, "Marginal Plausibility" (how the kernel-density baseline differed from the kernel-density form of stochastic attention); Discussion (three sentences added to the paragraph on the PCA representation and the generator); new Supplement section "Supplementary Bandwidth and Plausibility Results," first paragraph and Supplementary Table 17. Green. `\rthree{...}`

---

### R3.2: Membership-inference AUC and sharing synthetic data

> The membership-inference analysis reports an AUC of approximately 0.97, indicating substantial privacy risk. The authors should clarify the implications for sharing synthetic patient data and distinguish exploratory computational applications from privacy-preserving release.

*We agree. An AUC of $0.97$ means that the synthetic cohort does not protect membership, so stochastic attention as used here is not a privacy-preserving release method.* Someone who already holds a candidate de-identified profile can usually tell whether that profile was stored in the memory matrix. The revised Discussion now states this directly. The synthetic cohorts support exploratory computational work, such as hypothesis generation, mechanistic modeling, and workflow testing. A cohort generated from restricted data should be shared under the same terms as those data. Releasing synthetic data from sensitive cohorts would require a formal privacy method, such as differential privacy, which we did not test.

The membership signal came from the target distribution, not from the numerical sampler. Direct samples from the exact mixture gave a mean AUC of $0.98$. In the first round, we also varied the inverse temperature over a $40$-fold range. The median distance from synthetic profiles to the closest real record plateaued near $14.3$, below the median distance of $15.9$ between real patients. No setting of the inverse temperature gave a formal privacy guarantee.

In this study, the de-identified source data are public in the project repository, together with the list of the $23$ patients used to build the memory matrix. The synthetic cohorts therefore reveal nothing about the participants beyond what the repository already shows.

**Changes:** Discussion (four sentences added to the limitations paragraph, after the bounded statement of what the membership test showed); Data Availability (corrected). Green. `\rthree{...}`

---

### R3.3: Biological plausibility and positivity

> Only 54% of the generated profiles satisfy the biological plausibility criteria, and 24% satisfy both biological and strict mechanistic criteria. The authors should emphasize this limitation and discuss whether positivity-preserving transformations or constrained reconstruction are necessary for future applications.

*We agree that this is a limitation, and the revised Discussion now states it with the numbers. Future applications that use hormone values need a step that keeps those values positive, and the log transform was the better of the two remedies we tested.* We tested both remedies the reviewer named. The results are in the new Supplement section, paragraphs 2--4 and Supplementary Table 18.

Every profile that failed the biological check had a negative estradiol or progesterone value. Linear decoding has unbounded support, and the hormones lie closer to zero relative to their variation than the coagulation factors do. The strict mechanistic check is also demanding for real patients. It requires all $15$ predicted-to-recorded TGA ratios (five features at three visits) to fall within the real 5th to 95th percentile bands. With the bands built from all $23$ real profiles, which were the bands applied to the synthetic profiles, $8$ of $23$ real profiles passed. When each real profile was checked against bands built from the other $22$, $7$ of $23$ passed. The $24\%$ of synthetic profiles that passed both checks should be read against this real-patient rate of $30$--$35\%$.

We first tested constrained reconstruction. We added a positivity term $\lambda C$ to the energy. The term $C$ sums the squared negative parts of the six decoded hormone values, each divided by its real-cohort standard deviation. We sampled with Langevin dynamics at $\lambda=0$, $1$, $10$, and $100$, and compared each run with exact rejection sampling from the same target. Langevin and exact samples agreed at every $\lambda$: two-sample Kolmogorov-Smirnov tests on $1{,}000$ draws per sampler gave $p\ge0.002$ for all $24$ comparisons, above the Bonferroni threshold of $0.0004$. As $\lambda$ increased, the number of profiles with a negative hormone value fell from about $50$ at $\lambda=0$ to $4$--$6$ at $\lambda=100$. A hard constraint removed every negative value and kept $51\%$ of exact draws. The constraint also distorted the hormone distributions. Removing the negative tail raised the hormone means, so the median hormone error rose from about $5\%$ at $\lambda=0$ to $22$--$27\%$ at $\lambda=100$ and $30\%$ with the hard constraint. The median relative error across all features rose from $1.2$--$1.3\%$ to $2.4$--$2.5\%$.

We then tested a positivity-preserving transform. We replaced the six hormone values with their logarithms before standardization and PCA, and exponentiated them after decoding. All $100$ profiles passed the biological check, and $45$ also passed the strict mechanistic check. Without the transform, sampled the same way, $48$ profiles passed the biological check and $20$ passed both checks. Median relative error was $1.4\%$, compared with $1.7\%$ without the transform. Median hormone error was $6.2\%$, compared with $6.4\%$. Cross-visit correlation error was $0.58$, compared with $0.60$, and mechanistic overlap was $0.90$, compared with $0.91$. Log-transforming all $72$ assays also removed the negative values and gave a median relative error of $1.3\%$, but mechanistic overlap fell to $0.88$.

The log transform of the hormones therefore removed the negative values with similar fidelity. The constraint term reduced them at a cost in hormone accuracy. The reported cohort and the results in the main text are unchanged. The log representation is reported in the Supplement and recommended in the Discussion for future applications.

**Changes:** Results, "Interpolation and Extrapolation" (replaced the first-round sentence on positive-variable transformations with the real-patient pass rate and the log-transform result); Discussion (new sentences in the paragraph on the linear representation, stating the $54/100$ and $24/100$ rates, the real-patient rate, both remedies, and the recommendation); new Supplement section "Supplementary Bandwidth and Plausibility Results," paragraphs 2--4 and Supplementary Table 18. Green. `\rthree{...}`

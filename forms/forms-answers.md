# Reporting Summary — Fill-in Answers

Open `nr-reporting-summary-flat.pdf` in **Adobe Reader** (not Preview).

---

## Header

- **Corresponding author(s):** Jeffrey D. Varner
- **Last updated by author(s):** April 2026

---

## Statistics (checkboxes)

| Item | Check |
|---|---|
| Exact sample size for each group | **Confirmed** |
| Distinct samples vs repeated measures | **Confirmed** |
| Statistical tests and one/two-sided | **Confirmed** |
| Covariates tested | **n/a** |
| Assumptions or corrections | **Confirmed** |
| Central tendency and variation | **Confirmed** |
| Null hypothesis test statistics with P values | **Confirmed** |
| Bayesian analysis | **n/a** |
| Hierarchical/complex designs | **n/a** |
| Effect sizes | **n/a** |

---

## Software and code

**Data collection:**
Blood samples collected and assayed at University of Vermont as previously described (McLean et al. 2012, Hale et al. 2012). Coagulation factor activity levels by Stago clotting assays, fibrinolysis by ELISA, TGA on Synergy4 plate reader, ROTEM Delta Instrument (Werfen).

**Data analysis:**
Julia v1.12 (custom SA pipeline), Python 3 with sdv library for CTGAN/TVAE baselines. ODE integration via DifferentialEquations.jl. PCA via MultivariateStats.jl. All code available at https://github.com/varnerlab/SA-generation-legacy-dataset-paper.

---

## Data

Synthetic datasets and all validation outputs are available at https://github.com/varnerlab/SA-generation-legacy-dataset-paper. Real patient data were collected under IRB approval at the University of Vermont; de-identified data are available upon reasonable request to the corresponding author.

---

## Human participants

**Sex and gender:**
All participants were women desiring pregnancy. Sex was a criterion for enrollment given the study focus on pregnancy-related coagulation changes. Gender identity was not collected.

**Race, ethnicity:**
Race is reported in Table 1 (22 White, 1 not disclosed). The cohort reflects the demographics of the recruitment site (University of Vermont). Race-stratified analysis was not performed due to insufficient diversity in this small cohort.

**Population characteristics:**
23 women, mean age 30.2 +/- 5.1 years, mean prepregnancy BMI 26.6 +/- 5.3. 16 nulliparous, 7 parous. Enrollment: 14 healthy nulliparous, 6 prior PE, 3 PCOS. Full demographics in Table 1.

**Recruitment:**
Women desiring pregnancy were recruited at the University of Vermont Medical Center under three enrollment conditions (healthy nulliparous, prior preeclampsia, PCOS). Described in prior publications (McLean et al. 2012, Hale et al. 2012, Bernstein et al. 2016).

**Ethics oversight:**
University of Vermont Institutional Review Board. Written informed consent obtained from all participants.

---

## Field-specific reporting

Select: **Life sciences**

---

## Life sciences study design

**Sample size:**
K=23 patients with complete longitudinal data across 3 visits (from ~50 initially enrolled). N=100 synthetic patients generated per cohort. No formal power analysis was performed; the study purpose was to validate a generative method on a small existing cohort.

**Data exclusions:**
Patients with incomplete visit data or >30% missing assay values were excluded, yielding 23 of ~50 subjects. Exclusion criteria are stated in Methods.

**Replication:**
All computational results are reproducible via fixed random seed (42). Code and data are publicly available. ULA vs MALA control experiment (10 chains x 5,000 iterations) confirmed sampler equivalence (Supplementary Table S8).

**Randomization:**
This is a computational/observational study. Synthetic patient generation used random initialization on the unit sphere. No experimental group randomization was performed.

**Blinding:**
Blinding was not relevant to this computational study. All analyses were performed on the complete dataset with no subjective outcome assessment.

---

## Materials & experimental systems

All **n/a**:
- Antibodies: n/a
- Eukaryotic cell lines: n/a
- Palaeontology and archaeology: n/a
- Animals and other organisms: n/a
- Clinical data: n/a (secondary analysis of existing biospecimens, not a clinical trial)
- Dual use research of concern: n/a
- Plants: n/a

---

## Methods

All **n/a**:
- ChIP-seq: n/a
- Flow cytometry: n/a
- MRI-based neuroimaging: n/a

---

## Dual use / Hazards

All **No**:
- Public health: No
- National security: No
- Crops and/or livestock: No
- Ecosystems: No
- Any other significant area: No

## Experiments of concern

All **No**:
- Render a vaccine ineffective: No
- Confer resistance to antibiotics/antivirals: No
- Enhance virulence of a pathogen: No
- Increase transmissibility: No
- Alter host range: No
- Enable evasion of diagnostics: No
- Enable weaponization: No
- Any other potentially harmful combination: No

---

## Editorial Policy Checklist

**No longer required.** The PDF states: "This form is no longer required for Nature Portfolio submissions and has been removed."

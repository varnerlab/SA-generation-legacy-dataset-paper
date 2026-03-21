---
name: Hockin/Mann coagulation model package
description: Plan to build a standalone Julia package implementing the Hockin/Mann (2002) detailed coagulation model as an independent mechanistic validator for synthetic patient data
type: project
---

Build a standalone Julia package `HockinMannModel.jl` in its own repo implementing the Hockin/Mann (2002) detailed coagulation model (~30+ ODEs, literature rate constants).

**Why:** The BST model (9 reactions, 12 fitted parameters) cannot serve as an independent validator for synthetic patients because we fit parameters to each patient — there's no independent prediction. The Hockin/Mann model uses **literature rate constants with no fitting**. You supply only patient factor levels as initial conditions, and the model predicts the thrombin curve. Comparing predicted TGA features to SA-generated TGA features is a genuine independent validation with no circularity.

**How to apply:** Build as a separate callable package in its own repo. This paper's repo adds it as a dependency. The API should be simple: supply patient factor levels (nM), get back a thrombin curve and extracted TGA features (lagtime, peak, tpeak, max_rate, ETP).

**Proposed repo structure:**
```
HockinMannModel.jl/
├── src/
│   ├── HockinMannModel.jl    # module definition
│   ├── Model.jl              # ODEs, species, stoichiometry
│   ├── Parameters.jl         # literature rate constants (frozen)
│   └── Simulate.jl           # solve + extract TGA features
├── test/
│   └── runtests.jl           # validate against Hockin et al. 2002 figures
├── Project.toml
└── README.md
```

**Proposed API:**
```julia
using HockinMannModel
result = simulate(patient_factors; tspan=(0.0, 40.0))
tga = extract_tga_features(result)  # lagtime, peak, tpeak, max_rate, etp
```

**Validation strategy:** The model with frozen literature parameters predicts what a patient's thrombin curve *should* look like given their factor levels. For synthetic patients, compare the Hockin/Mann prediction to the SA-generated TGA features. Divergence indicates the synthetic patient's factor-to-TGA relationship is biologically implausible.

**Reference:** Hockin MF, Jones KC, Everse SJ, Mann KG. A model for the stoichiometric regulation of blood coagulation. J Biol Chem. 2002;277(21):18322-18333.

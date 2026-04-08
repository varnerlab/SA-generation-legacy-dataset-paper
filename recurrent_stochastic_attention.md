# Recurrent Stochastic Attention (SA-RNN)

**Status:** brainstorm — not yet implemented or formalized.
**Date captured:** 2026-04-08
**Source:** conversation with Claude during the SA legacy-dataset paper session.

---

## The core idea

A vanilla Elman RNN has a hidden state vector $h_t \in \mathbb{R}^d$ that's
updated each step:

$$h_t = f(W_x x_t + W_h h_{t-1} + b)$$

**The nutty proposal:** replace $h_t$ with a *position in a modern Hopfield
energy landscape*, and replace the recurrence with one or more
**Langevin steps biased by the new input**. The recurrent state becomes a
particle wandering on a learned semantic landscape, with the input field
steering it.

This is a domain-agnostic architecture. The only requirement is that you
have:
1. A small set $\{m_k\}$ of well-characterized exemplars in some vector space.
2. A vectorization that captures the essential structure.
3. A need to generate sequentially (predict next, refine, walk).

---

## Three flavors

### Flavor A — Hidden state = current SA sample in a fixed landscape

Patterns $\{m_k\}$ are learned (or just data) and frozen. At each timestep:

$$\hat{\xi}_t = \text{LangevinStep}(\hat{\xi}_{t-1}, \{m_k\};\ \text{input bias from } x_t)$$

The input modulates either $\beta$ (retrieval sharpness) or the multiplicity
vector $\mathbf{r}$. The whole sequence is one Langevin trajectory on a fixed
energy surface.

**Closest analogue:** an RNN whose recurrence is a physically meaningful
Markov chain on a learned energy landscape.

### Flavor B — Hidden state = the entire memory bank, which evolves

Patterns themselves are the recurrent state. Each timestep, the input rewrites
or modulates the bank.

$$\{m_k\}^{(t)} = g(\{m_k\}^{(t-1)}, x_t),\quad \hat{\xi}_t \sim \text{SA}(\{m_k\}^{(t)})$$

State is now $K \cdot d$-dimensional. **Closest analogue:** Neural Turing
Machine / Differentiable Neural Computer, but with explicit Hopfield physics
instead of learned read/write heads.

### Flavor C — Multiplicity weights ARE the recurrent state ⭐

Patterns are fixed (typically the data itself). Recurrent state is just the
$K$-dimensional multiplicity vector $\mathbf{r}^{(t)} \in \mathbb{R}_+^K$.
Each timestep updates the multiplicity based on input:

$$\mathbf{r}^{(t)} = h(\mathbf{r}^{(t-1)}, x_t),\quad \hat{\xi}_t \sim \text{SA}(\{m_k\},\ \mathbf{r}^{(t)})$$

This is the cleanest version and the one most aligned with existing SA
infrastructure. The recurrent state is a learned distribution over training
exemplars, evolving over time, with theoretically grounded retrieval physics
underneath.

**Recommendation:** start with Flavor C.

---

## Why this is interesting

1. **The hidden state is interpretable as a distribution over training
   examples.** $r_k^{(t)}$ tells you "how relevant is exemplar $k$ to the
   current step." Built-in attention map over your training data, available
   for free at every timestep, with the same units across the whole training
   run. You can literally watch the model "think" through which exemplars
   it's drawing on.

2. **Built-in stochasticity with controllable temperature.** $\beta$ gives
   one knob from "snap to nearest exemplar" (high $\beta$, memorization) to
   "smooth interpolation between exemplars" (low $\beta$, exploration). Anneal
   $\beta$ within a sequence — confident first, exploratory later, or vice
   versa.

3. **The $n \ll p$ regime is its native environment.** Standard RNNs need
   lots of training sequences. SA-RNN where the memory bank IS the training
   set just needs *enough patterns to define a useful landscape*.

4. **Continual learning is structurally easier.** Add a new exemplar = add a
   new pattern to the bank. Landscape automatically deforms. No catastrophic
   forgetting because the old patterns are still there as fixed attractors.

5. **Counterfactual trajectories via multiplicity injection.** Mid-sequence,
   intervene on $\mathbf{r}^{(t)}$ — "from now on, upweight the PCOS
   patterns." The trajectory carries the intervention forward. Dynamic version
   of the conditional generation lever from the static SA paper.

6. **A clean recurrent analogue of in-context learning.** Transformers do
   in-context learning by attention over the prompt. SA-RNN with the prompt
   as the pattern bank does the same thing, but with theoretical grounding
   (Hopfield retrieval physics) and a recurrent structure that handles long
   sequences without quadratic cost.

---

## The unifying claim across all SA application domains

Existing SA papers from the Varner lab:
- **FBA / CFPS:** https://github.com/varnerlab/SA-CFPS-Paper (cell-free
  protein synthesis design)
- **Protein binding:** https://arxiv.org/abs/2603.20115
- **Protein sequences:** https://arxiv.org/abs/2603.14717
- **Patients (longitudinal coagulation):** the current paper

Stripped of domain-specific framing, all four are doing one thing:

> Given a small set $\{m_k\}$ of characterized exemplars in some vector space,
> generate novel samples that respect the geometry of those exemplars, with
> optional conditional steering via multiplicity weights.

Recurrence turns this from one-shot sampling into stepwise generation
conditioned on a partial trajectory, with the recurrent state being a learned
distribution over which exemplars are currently in scope. Same machinery,
sequential.

**SA-RNN is therefore a single architecture with four ready-made
application demonstrations from day one.**

---

## What recurrence buys each domain

### FBA / CFPS

- **Time-dependent flux trajectories.** Model how the flux distribution
  evolves over the course of a CFPS run (substrate depletion, enzyme
  inhibition, product accumulation), with each step's hidden state being
  "which characterized reactions does my current state most resemble."
- **Active design / sequential experimentation.** At each step, the SA-RNN
  proposes the *next* experiment given the trajectory of previous experiments.
  Multiplicity reweights to favor unexplored regions or high-yield
  neighborhoods. Bayesian-optimization-like, but with a Hopfield prior instead
  of a Gaussian process — works in regimes where GPs choke (high dimensions,
  small libraries, structured designs).
- **Batch-to-batch evolution.** Long fermentation / continuous CFPS runs.
  The recurrent state tracks drift from the original calibrated design.

### Protein binding

- **Directed mutation walks.** Start from a seed, take SA-Langevin steps in
  sequence-feature space, where each step reweights which library binders
  are most relevant to the current candidate. The trajectory IS the design
  path — stop, branch, or backtrack.
- **Multiplicity-constrained optimization.** Hold multiplicity high on
  "specificity" exemplars and "affinity" exemplars simultaneously, sample
  sequences that satisfy both. Constraint propagation through a learned
  geometry rather than through a loss function.
- **Pose / conformation sequences.** If binders have docking-pose
  vectorizations, the recurrent SA could walk a candidate binder through a
  sequence of conformational states biased toward known bound poses.

### Protein sequences

This is the one with the most striking pitch.

- **Residue-by-residue generation from a small family alignment.** Patterns
  are $K = 20$ to $K = 200$ family sequences. Recurrent state at position $t$
  is "which family members does my partial sequence so far most resemble."
  Each new residue is sampled from the SA distribution restricted to the
  current multiplicity weights.
- **Compare to PLMs.** ESM and friends learn from billions of sequences.
  SA-RNN over a family alignment doesn't *learn* anything — it samples
  directly from the geometry of your handful of sequences. For **orphan
  families, ancient lineages with few extant homologs, designed families with
  a handful of seeds, niche viruses with limited sequence depth** — anywhere
  you can't bring billions of sequences to bear — SA-RNN could be the right
  tool. Pitch: "PLM-style sequence generation but you don't need a PLM."
- **Counterfactual residues.** Mid-sequence, inject "now favor the
  thermophile members of this family" via multiplicity weights, and watch
  the rest of the sequence shift toward thermostability.

### Patients (current legacy-dataset paper)

- **Sequential trajectory forecasting from partial histories.** Predict the
  rest of a held-out patient's trajectory given their V1+V2 measurements.
  The recurrent state tells you "this incoming patient is starting to look
  like our PCOS subgroup" by visit 2.
- **Counterfactual treatment simulators.** Run a patient forward from
  baseline under different interventions by modifying $\mathbf{r}^{(t)}$
  mid-sequence.
- **Mechanistic-data hybrid models.** Mix real patients and ODE-simulated
  patients in the memory bank; the multiplicity dynamics let the SA-RNN
  draw on both as needed.

---

## Shape of a research program

This is potentially a **single methodological contribution with four
published application demonstrations**, all pre-existing in the Varner lab.

1. **Theory paper.** *"Stochastic Attention as a Recurrent State for
   Few-Shot Sequential Generation."* Establishes the SA-RNN architecture
   (Flavor C), gives convergence / capacity properties from underlying
   Hopfield theory, validates on one illustrative dataset.

2. **Application 1: Patients (extension of the current paper).** Sequential
   clinical-trajectory forecasting from $K = 23$. Shows SA-RNN works where
   vanilla LSTMs can't (they can't be trained on $N = 23$ sequences).

3. **Application 2: Protein sequences.** Family-conditioned residue-by-residue
   generation from a small Pfam family. Compare to ESM at matched sequence
   count. SA-RNN should win in the small-$N$ regime.

4. **Application 3: CFPS.** Sequential experiment design / batch-evolution
   prediction. Show that SA-RNN proposes useful next-experiments given a
   partial design history.

5. **Application 4: Binder design.** Directed mutation walks from a seed,
   with multiplicity constraints for specificity + affinity. Novel
   high-affinity binders generated by trajectory rather than by static
   optimization.

6. **Infrastructure deliverable.** A single Julia package
   `RecurrentStochasticAttention.jl` that takes a memory bank and an update
   rule and produces an iterator over samples — totally domain-agnostic, with
   adapters for each of the four application domains. Infrastructure as
   research artifact.

---

## Concrete first move

Cleanest first experiment is probably the **protein-sequence one**, because:

- Discrete-output domain — recurrent step (one residue at a time) is natural.
- Family alignments give $K$ in the range where SA shines (20–200).
- Baselines (PLMs at small $N$, profile HMMs, PSSMs) are well established.
- Result is binary: either SA-RNN generates valid family members at small
  $N$ or it doesn't.
- An SA implementation for protein sequences already exists.

If that works, patients / CFPS / binders are straight ports of the same
recurrent update rule with different vector spaces and different decoders.

---

## Sketch of the Flavor-C update rule

Pseudocode for one timestep, given a memory bank
$M = [m_1, \ldots, m_K] \in \mathbb{R}^{d \times K}$ (fixed throughout the
sequence) and an input $x_t$:

```
function sa_rnn_step(r_prev, x_t, M; β, α, T_langevin)
    # 1. Update multiplicity from input
    #    Many possible update rules — start simple:
    similarities = M' * x_t                    # K-vector of dot products
    bias = exp(γ * similarities)                # γ controls input strength
    r_t = normalize(r_prev .* bias)             # element-wise mass update

    # 2. Sample from the multiplicity-weighted Hopfield landscape
    #    (this is exactly the existing weighted_sample from Compute.jl)
    ξ_init = sample_init_from(M, r_t)           # warm-start near a likely pattern
    ξ_t = weighted_sample(M, ξ_init, T_langevin; β=β, α=α, r=r_t)

    # 3. Decode the sample to whatever the output space is
    #    (per-domain: residue logits / patient features / flux vector / ...)
    output = decoder(ξ_t)

    return r_t, ξ_t, output
end
```

Things to vary / explore:

- **Update rule for $\mathbf{r}^{(t)}$.** Multiplicative (above), additive
  in log-space, attention-style softmax, learned gating (LSTM-like). Start
  multiplicative.
- **Input encoder.** How does $x_t$ become a query against $M$? Could be
  the raw vector, a learned linear projection, or a feature extractor.
- **Annealing $\beta$.** Across timesteps within a sequence; or as a
  function of confidence (entropy of $\mathbf{r}^{(t)}$).
- **Number of Langevin steps per timestep.** Probably small (1–10), since
  we're not trying to converge to the stationary distribution — we're using
  it as an inductive bias.
- **Decoder.** Domain-specific. For sequences, a softmax over residue
  alphabet. For patients / CFPS, inverse PCA + destandardization.

---

## Failure modes to flag upfront

- **Backprop is awkward.** Langevin steps include noise; reparameterization
  helps but isn't free. Start non-trainable (patterns are data, no learnable
  weights) and only later add a learned similarity function or input encoder.
- **Cost per step is higher than vanilla RNN.** A standard RNN cell is one
  matrix-vector product; one Langevin step is several. For long sequences
  this matters. Mitigation: only run a few Langevin steps per timestep.
- **Recurrent mode collapse.** If $\mathbf{r}^{(t)}$ gets pulled too sharply
  toward a single pattern, the rest of the trajectory just memorizes that
  one exemplar. Fix: regularize $\mathbf{r}$ toward uniform, or use a
  participation-ratio constraint like the static SA paper does for
  conditional generation.
- **What if exemplars aren't sequential?** The current SA papers concatenate
  all positions/visits into one long vector per pattern. For an SA-RNN you
  probably want patterns that ARE sequences (or per-step slices that the
  recurrence walks through). Slight reframing of data layout, not a
  deal-breaker.
- **Output decoding for discrete domains.** For protein sequences the output
  is a residue, not a continuous vector. Need to decide whether the SA sample
  in PCA space gets decoded to residue logits via a learned head or by some
  template-matching scheme.

---

## Open questions worth thinking about later

- **Theoretical capacity.** Modern Hopfield networks have proven retrieval
  capacity. What's the analogous notion of "expressive capacity" for SA-RNN
  as a sequence model? How does it scale with $K$?
- **Convergence of the multiplicity dynamics.** Does $\mathbf{r}^{(t)}$
  converge to a fixed point under typical update rules? Is the fixed point
  meaningful?
- **Connection to score-based sequence models.** Diffusion language models
  also use noise-based sampling for sequence generation. Where does SA-RNN
  fit in that landscape? Is it a special case of something? Or a genuinely
  different inductive bias?
- **Equivalence to known architectures.** Is SA-RNN with a particular update
  rule equivalent to a known recurrent attention mechanism? (E.g.,
  fast-weight RNNs, retentive networks, state-space models.) If so, the
  contribution is the *small-data Hopfield framing*, not the architecture
  itself.
- **What's the right benchmark for "few-shot sequential generation"?** The
  ML community doesn't have a standard small-sample sequence-generation
  benchmark — the SA-RNN paper might need to define one (probably using one
  of the four application domains as the canonical test).

---

## TL;DR

Take SA, make the multiplicity vector the recurrent state of an
Elman-style RNN, update it each timestep based on input, sample the next
output from the multiplicity-weighted Hopfield landscape. Result: a
domain-agnostic architecture for few-shot sequential generation that should
work in the same small-$N$ regime where current SA already works, with the
added ability to forecast trajectories, propose next experiments, walk
mutation paths, and counterfactually steer mid-sequence.

Cleanest first experiment: protein sequences from a Pfam family, compared
against PLMs at matched sequence count. If that works, the patient / CFPS /
binder ports are straightforward.

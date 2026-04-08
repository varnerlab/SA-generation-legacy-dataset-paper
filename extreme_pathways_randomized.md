# Approximate / Random Extreme Pathways for SA Memory Generation

**Status:** brainstorm — not yet implemented or formalized.
**Date captured:** 2026-04-08
**Source:** conversation with Claude during the SA legacy-dataset paper session.
**Sibling note:** [recurrent_stochastic_attention.md](recurrent_stochastic_attention.md)

---

## The real question

Long-standing interest: implement extreme pathways (ExPa) for stoichiometric
networks. Long-standing block: convex independence step in the double
description method has insanely bad combinatorial scaling. Can we get a
useful approximation via random / sampling-based methods?

But the deeper motivation, made explicit now, is **not** "I want extreme
pathways for their own sake." It's:

> **Where do the memory patterns for SA on flux-balance / metabolic-network
> problems come from?**

Existing SA application: https://github.com/varnerlab/SA-CFPS-Paper —
SA applied to cell-free protein synthesis design, working over the flux cone.
The current memory bank in that work is presumably whatever experimental /
sampled / hand-curated flux vectors were available. **The principled answer
to "what should the memory bank be" is: the extreme pathways of the
network.** They are, by definition, the unique minimal convex basis for the
steady-state flux cone — every achievable flux distribution is a non-negative
combination of them. They are *the* natural memories for an SA system that
samples in the flux cone.

So the bind is:

- SA needs a memory bank to define its energy landscape.
- The right memory bank for flux problems is the set of extreme pathways.
- Extreme pathways are intractable to compute exactly for genome-scale
  networks.
- Therefore SA on metabolic networks is currently working with ad hoc
  memories instead of principled ones.

A randomized / approximate ExPa algorithm would unlock principled SA memory
banks for any metabolic network, including genome-scale ones. **That is the
real motivation for the question.**

---

## Why extreme pathways are the right memories

The steady-state flux cone is

$$C = \{ v \in \mathbb{R}^n_{\geq 0} : Sv = 0 \}$$

This is a polyhedral convex cone in $\mathbb{R}^n$. Its extreme rays are the
unique minimal set $\{p_1, \ldots, p_K\}$ such that every $v \in C$ can be
written as

$$v = \sum_k \alpha_k p_k, \qquad \alpha_k \geq 0$$

In the SA framing:

| SA concept                  | Flux cone equivalent                                |
| --------------------------- | --------------------------------------------------- |
| Memory pattern $m_k$        | Extreme pathway $p_k$                               |
| Stored pattern set          | Set of extreme pathways                             |
| Hopfield energy minimum     | Vertex of the cone (extreme ray)                    |
| Langevin sample             | Novel flux distribution near a combination of rays  |
| Multiplicity-weighted bank  | Reweighted basis favoring particular rays           |
| Direction-magnitude split   | Direction in cone × scale of total flux             |

Under this mapping, **SA on the flux cone with extreme pathways as memories
is exactly "sample novel flux distributions that respect the geometric
structure of the cone, biased toward known basic pathways."** That's a
mechanistically meaningful generative model in a way that no other choice of
memory bank is.

The conditional generation lever maps too: multiplicity-weighting a subset of
extreme pathways means "generate fluxes that primarily use these basic
pathways" — i.e., subgroup-conditioned flux generation, where the subgroups
are *biological* (pathways using glycolysis vs gluconeogenesis vs pentose
phosphate, etc.) rather than statistical.

---

## Why this is hard

Schilling–Letscher–Palsson (2000) defined ExPa via a double description method
on the constraint cone. The bottlenecks:

1. **Pivoting blowup.** Each constraint is processed by combining pairs of
   rays from the previous iteration. The number of candidate pairs is
   $\binom{K^{(t)}}{2}$ at iteration $t$, and $K^{(t)}$ itself can grow
   exponentially.

2. **Convex independence is per-pair.** For each candidate pair, you check
   whether their combination is a true new extreme ray (i.e., not in the
   conic hull of the existing set). The check itself is polynomial, but it
   has to be done for combinatorially many pairs.

3. **The output is exponential.** For genome-scale networks the number of
   extreme pathways is itself exponential in the number of reactions. *Any*
   algorithm that enumerates them all is at least linear in the output size,
   so it's at least exponential. There is no polynomial-time complete
   algorithm in principle.

Counting EFMs (a closely related set) is at least #P-hard via reduction from
monotone-2SAT counting. ExPa inherits the same intractability.

**Conclusion:** any useful algorithm has to give up completeness. The
question is what to give up *for*.

---

## Reframing the goal: not "complete" but "useful"

Three relaxations of "complete extreme pathway set":

1. **Coverage relaxation.** Find a set $\mathcal{E}$ such that every reaction
   is used by at least one pathway in $\mathcal{E}$, and every "biologically
   interesting subspace" is represented. Loses uniqueness but keeps utility
   for FBA, knockout analysis, pathway-existence reasoning.

2. **Importance relaxation.** Find the $K$ "most important" extreme pathways
   under some weighting (shortest, highest yield, lowest cost). De Figueiredo
   et al. (2009) did this for EFMs with k-shortest enumeration.

3. **Probabilistic relaxation.** Find a random sample of $K$ extreme pathways
   from some target distribution. With high probability, the sample covers
   the high-mass regions of the true extreme ray set.

For SA memory generation, **option 1 is what we actually want.** SA doesn't
need uniqueness; it needs a memory bank that spans the geometry of the
biologically relevant flux distributions. Coverage is the right success
criterion.

---

## What's already been tried (so we don't reinvent)

- **DDM with smart pivoting:** efmtool (Marashi/von Kamp/Klamt), polco,
  efmlrs, lrs. Better constants, same asymptotic wall.
- **k-shortest EFM enumeration via MILP:** Pey et al., De Figueiredo et al.
  Gives a structurally meaningful subset.
- **Flux cone sampling:** CHRR (Haraldsdóttir et al.), ACHR (Schellenberger),
  hit-and-run. These sample interior fluxes, not extreme rays directly.
- **Random EFM sampling:** Gerstl et al. (~2019) — unbiased random
  sampling of EFMs. Doesn't give completeness but gives an unbiased subset.
  **Worth a literature scope here — this line might be more developed than I
  remember.**
- **Minimal cut sets / elementary flux patterns:** Klamt et al. — alternative
  representations that scale better in some cases.
- **Network compression / lumping:** preprocess the network before ExPa.
  Helps in practice but doesn't change asymptotic behavior.

What I am **not aware of** (and would be embarrassed to miss in a literature
scope):

- A randomized algorithm for ExPa with provable coverage guarantees.
- A learning-based proposer that uses past extreme rays to generate new
  candidate rays much faster than DDM enumeration of pairs.
- Any work that frames extreme pathway discovery as the **memory-bank
  generation problem for a downstream sampler** (this is the SA-aligned
  framing and it's potentially novel as a research direction).

**First action item: literature scope.** Search Klamt, Marashi, Bockmayr,
Müller, Avis, Fukuda, Gerstl, Pey, Haraldsdóttir, von Kamp, Wagner.
Check workshops on polyhedral computation and computational metabolic
engineering. Last 5 years especially.

---

## Approaches worth trying

Listed in order from "engineering-y" to "research-y."

### Approach 1 — Random-objective LP

For each iteration:

1. Sample a random direction $w \sim \mathcal{N}(0, I_n)$.
2. Solve the LP: $\max\ w^T v$ s.t. $Sv = 0$, $v \geq 0$, $\mathbf{1}^T v = 1$.
3. The optimum is an extreme ray of the flux cone (LP optima are vertices).
4. Add to collection if new.

**Why it works:** every extreme ray is the optimum of *some* LP. With enough
random objectives, you sample the rays whose basin of optimality is large.

**Why it's approximate:** rays with small basins are exponentially less
likely to be hit. You miss the rare ones — usually the boring/degenerate
ones, but no guarantees.

**Cost:** one LP per iteration, polynomial-time per LP. Massively cheaper
than DDM.

**Status of prior art:** random-objective vertex enumeration is a standard
technique in computational geometry. Whether anyone has applied it
specifically to ExPa enumeration as a publication, I'm not certain — worth
searching.

### Approach 2 — Sample-then-cluster

Use CHRR / ACHR to sample $N$ flux distributions from the cone interior.
Project each onto the unit sphere to get directions. Cluster the directions
(spherical k-means or DBSCAN on cosine distance). Cluster centroids are
candidate extreme ray directions. Verify each candidate is extreme via local
LP.

**Why it might work:** dense interior sampling reveals the cone's vertex
structure as density modes near the boundary.

**Why it might not:** uniform polytope sampling overwhelmingly hits the
*interior*, not the boundary. You'd need importance weights that favor
boundary points — possibly via simulated annealing with an "extremity bonus."

### Approach 3 — Active learning over extreme rays

Maintain a discovered set $\mathcal{E}$. Each iteration:

1. Find a direction not well-covered by $\mathcal{E}$ (e.g., a direction
   orthogonal to the conic hull of current rays).
2. Solve a small LP in that direction to find a new extreme ray.
3. Add to $\mathcal{E}$.
4. Stop when discovery rate falls below a threshold.

This is essentially **directional flux variability analysis as a discovery
loop**. Each iteration is one LP. Natural stopping criterion. Might miss
"thin spikes" — extreme rays in unusual directions that random probing
misses.

### Approach 4 — Submodular cover (more theoretical)

Define a coverage function $f(\mathcal{E})$ — e.g., the volume of the cone
covered by $\text{cone}(\mathcal{E})$, or the fraction of biologically
realizable flux distributions within $\varepsilon$ of some combination of
rays in $\mathcal{E}$. If $f$ is monotone and submodular, greedy gives a
$(1 - 1/e)$-approximation.

The bottleneck becomes "find the single best ray to add at each step," which
is itself hard but potentially LP-tractable. **I don't know if this has been
done; if not, it's a clean theoretical contribution waiting.**

### Approach 5 — SA as the proposer (the wheelhouse-aligned approach) ⭐

This is the one that ties everything together.

**The chicken-and-egg.** SA needs a memory bank. The principled memory bank
is the extreme pathway set. Extreme pathways are intractable. So we need to
*bootstrap*: start with a small seed of extreme rays (cheaply obtained), use
them as the SA memory bank, run SA-Langevin to propose new candidate flux
vectors, verify extremity via LP, add successful candidates to the memory
bank, repeat.

The loop:

```
seed = compute_exact_expa(small_subnetwork)        # cheap, ~10-50 rays
memory_bank = seed
loop:
    candidate = sa_langevin_step(memory_bank)      # propose novel flux
    candidate = project_onto_cone(candidate)       # enforce Sv=0, v>=0
    if is_extreme_ray(candidate):                  # one feasibility LP
        memory_bank = memory_bank ∪ {candidate}
    if discovery_rate < threshold or coverage > target:
        break
return memory_bank
```

**Why this is interesting beyond just ExPa:**

1. **It closes a real gap in the SA research program.** Currently SA on flux
   problems uses ad hoc memories. After this, SA on flux problems uses the
   principled ones. Every downstream SA-CFPS analysis becomes more
   defensible.

2. **It's self-referential in a satisfying way.** SA generates its own
   memory bank by running on a seed memory bank and using LP verification as
   a filter. The proposal distribution becomes more accurate as the bank
   grows, which makes proposals more likely to be true extreme rays, which
   grows the bank faster — positive feedback toward convergence.

3. **Multiplicity weighting steers the discovery process.** Want to find
   extreme rays that use the pentose phosphate pathway? Set $\rho > 1$ for
   patterns currently using PPP reactions; SA will preferentially propose
   candidates in that region. This is **active learning with biological
   priors built in**, not just random exploration.

4. **The recurrent SA-RNN idea plugs in directly.** Each iteration of the
   discovery loop is one timestep; the recurrent state is the multiplicity
   over discovered rays; the output is a candidate flux vector to verify.
   Cross-pollinates with the sibling brainstorm.

5. **Generalizes beyond flux cones.** Any time SA is applied to a vector
   space defined by linear constraints (or any constraint geometry with
   extreme points), the same bootstrapping approach applies. Polytope
   vertices, conformer libraries, allowed gene-expression states under
   regulatory constraints — anywhere the constraint geometry has well-defined
   extreme points.

**Why this might not work:**

- **Projection onto the cone is a non-trivial operation.** SA-Langevin
  produces unconstrained samples; you'd need to project them onto
  $\{v : Sv = 0, v \geq 0\}$ before checking extremity. LP projection works
  but adds cost per candidate.
- **Verifying extremity is itself an LP.** Cheap individually but not free
  — you do it once per candidate, and most candidates won't be extreme.
- **Convergence isn't guaranteed.** The discovery rate might plateau before
  full coverage, especially if extreme rays cluster in regions far from the
  seed bank. You'd need a separate "exploration injection" mechanism — e.g.,
  periodic random-objective LPs to find rays in unexplored directions.

---

## The bigger reframing

This brainstorm is no longer really about extreme pathways. It's about a
**research program in SA memory generation**:

> **Memory-generation problem for SA.** SA's behavior is determined by its
> memory bank. For SA to be useful in a domain, the memory bank must
> faithfully represent the geometry of the relevant region of the vector
> space. Where does that bank come from?

For each existing SA application:

- **Patients (legacy dataset):** memory bank = K=23 real patient profiles.
  **Source:** the data itself. No memory generation problem because the data
  is the memory.
- **Protein sequences:** memory bank = family alignment sequences.
  **Source:** Pfam, MSA tools. No generation problem; the alignment is the
  memory.
- **Protein binding:** memory bank = characterized binders.
  **Source:** experimental affinity datasets. Usually small but obtainable.
- **CFPS / FBA:** memory bank = ??? Currently ad hoc (sampled fluxes, FBA
  optima, hand-curated). **The principled answer is extreme pathways, but
  they're intractable.**

The CFPS / FBA case is the only one where the memory generation problem is
genuinely open. Every other domain gets memories "for free" from data. The
flux cone is special because:

- The constraint geometry is non-trivial (polyhedral cone in high dimensions).
- The natural "data" (experimentally measured fluxes) doesn't span the cone
  uniformly — you mostly get fluxes near "typical" growth conditions, not
  the extreme rays.
- The principled basis (extreme pathways) is computable in theory but not in
  practice.

Solving the FBA memory generation problem closes the SA research program for
constraint-geometric domains. And the bootstrapping approach (Approach 5) is
the kind of solution that would actually generalize: any future SA
application to a polyhedrally-constrained vector space inherits the same
machinery.

---

## Concrete first move

Pick a medium-sized model. Suggested: **E. coli core (95 reactions, 72
metabolites)** — small enough that DDM gives ground truth (~20k extreme
pathways), big enough that the question is non-trivial.

Implement:

1. **Approach 1 (random-objective LP)** as a baseline. Measure coverage of
   the true ExPa set vs iteration count.
2. **Approach 5 (SA-bootstrap)** seeded with a tiny initial set (10–20 rays
   from DDM on a subnetwork).
3. Compare:
   - Coverage of true ExPa set (fraction of true rays found)
   - Time to discover the first 10%, 50%, 90% of true rays
   - Approximation quality of $\text{cone}(\mathcal{E})$ vs the true cone
     (Hausdorff distance, volume ratio)
   - Wall-clock cost vs DDM

If SA-bootstrap dominates random-objective LP — or even matches it with much
better biological steering — that's the headline result. If it underperforms,
the failure mode itself tells us something about the geometry.

After E. coli core, scale up to **iJO1366** (E. coli iJO1366, ~2500 reactions,
known ExPa-intractable). At that scale, DDM gives no ground truth, so
evaluation has to be coverage-based: "what fraction of biologically known
pathways did we recover?"

---

## Open questions worth thinking about

- **Connection to tropical geometry.** Tropical metabolic networks
  (Allamigeon, Gaubert) reinterpret the cone in $(\mathbb{R}, \max, +)$. Are
  tropical extreme rays cheaper to compute, and how do they relate to
  classical extreme rays?
- **Importance vs uniformity.** If you weight extreme pathways by some
  biological importance (yield, robustness, knockout sensitivity), can you
  define a sampler that targets the *important* rays rather than uniformly
  sampling all rays?
- **Active learning theory.** If SA-bootstrap is "active learning with an
  LP verifier," there's existing theory on active learning sample
  complexity. Can it give convergence guarantees for the discovery loop?
- **The right benchmark.** No standard "approximate ExPa" benchmark exists.
  E. coli core works for ground truth; iJO1366 works for scale. Anything in
  between?

---

## Status summary

| Question | Answer |
|---|---|
| Is computing complete ExPa intractable? | Yes, provably (#P-hard counting). |
| Is approximate / random ExPa intractable? | Probably not, depending on what "approximate" means. Open as far as I know. |
| Has anyone tried random-objective LP for ExPa? | Almost certainly yes for general vertex enumeration; not sure about ExPa specifically. **Literature scope needed.** |
| Has anyone framed this as a memory-generation problem for a downstream sampler? | Not that I know of. **Potentially novel.** |
| Does the SA-bootstrap loop have prior art? | The "active learning with LP verifier" pattern is general; the SA-as-proposer specialization is novel-ish. |

---

## TL;DR

The real question isn't "compute extreme pathways" — it's "where do SA memory
banks come from for flux problems." Extreme pathways are the principled
answer but they're intractable. Bootstrap by seeding SA with a small ExPa
subset, using SA-Langevin as the proposer for new candidate rays, and
verifying extremity via LP. This makes SA on metabolic networks defensible
in a way it currently isn't, gives the SA research program a unifying
"memory-generation" thread, and is potentially publishable both as a
methods paper (approximate ExPa) and as a research-program paper
(memory generation for constraint-geometric SA).

Cleanest first experiment: run on E. coli core where DDM gives ground truth,
compare random-objective LP baseline against SA-bootstrap, report coverage
vs iteration. If SA-bootstrap wins, scale up to iJO1366 where there is no
ground truth and coverage becomes the only success criterion.

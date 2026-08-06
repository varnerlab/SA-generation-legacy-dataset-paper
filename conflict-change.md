# Conflict: introduction mechanism vs. response-to-reviewers

**Date:** 6 August 2026
**Status:** Resolved on disk. Not committed.

## The conflict

The committed manuscript and the response document gave opposite explanations
for why SA samples are not copies of stored patients.

The introduction at `HEAD` (commit `70edfb2`,
`paper/sections/introduction.tex`) attributed the variation to overlap between
memory basins:

> continuous energy landscape\rone{, and the attention mechanism makes that
> landscape a distribution centered on the stored patterns. Sampling it yields
> new profiles rather than copies, because the memory-centered components
> **overlap in the reduced space**}.

The response document
(`peer-review-feedback/response-to-reviewers.md`, line 125, under
**R1.9 / R2.minor2**) explicitly rejects that explanation:

> The Introduction now explains that the Hopfield energy induces a continuous
> probability distribution—here, exactly a finite Gaussian mixture with one
> nonzero-variance component centered on each stored pattern—and that this
> **within-component variance, rather than overlap between memory basins**,
> produces non-identical samples.

This matters because the response document promises reviewers a specific
mechanism. The committed introduction described a different one.

## Resolution

The working-tree edit to `paper/sections/introduction.tex` replaces the overlap
explanation with the within-component-variance explanation:

> continuous energy landscape\rone{. This energy induces a continuous
> probability distribution that is exactly a finite Gaussian mixture in the
> reduced space, with one nonzero-variance component centered on each stored
> pattern. Sampling from this distribution therefore produces continuous
> variations around the stored patterns rather than exact copies}.

The manuscript now matches what the response document claims. The response
document needed no change. It was already written against the corrected
mechanism.

The same correction is consistent with the R1.5 privacy analysis, which
establishes the target as a finite isotropic Gaussian mixture with covariance
beta-inverse times the identity.

## Outstanding

1. The fix is uncommitted. It sits in the working tree along with the
   proof-of-concept rewording and the removal of "presented in the next
   section."
2. `arxiv/sections/introduction.tex` carries this same correction. Full
   paper-to-arxiv sync is deferred until the main-paper edits are done.
3. Minor and cosmetic: the `**Changes:**` line for R1.9 / R2.minor2 reads
   "Violet." but lists both `\rone{...}` and `\rboth{...}`, and the introduction
   passage above is marked `\rone{}` (blue). Elsewhere the color word and the
   macro agree. Either drop the word "Violet" or re-mark the passage as
   `\rboth{}`, since R1.9 / R2.minor2 is a joint comment.

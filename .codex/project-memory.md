# Project Goal

Develop MiSo into a statistics paper centered on interpretable discovery of
recurring factor-usage motifs in NMF with posterior uncertainty.

Publication name for manuscript authorship: Dat Do.

The intended paper is not just another NMF method. The main contribution
should be a probabilistic model that learns repeated factor-usage motifs,
clusters samples by those motifs, and quantifies motif-level uncertainty through posterior distributions such as `gamma_bar[s, d, k]`.

## Core Framing

MiSo represents recurring factor-usage patterns as low-dimensional
submanifolds spanned by subsets of potentially dense factor rows. It does not
impose sparsity on F, and it provides posterior uncertainty over which factors
define each motif.

## Immediate Roadmap

Validate MiSo against MF-Poisson-SuSiE and simple baselines across
controlled simulations:

- correct `K`, `D`, and known or well-initialized `F`
- learned `F` with correct `K` and `D`
- overspecified `K`
- overspecified `D`
- correlated or redundant factor rows
- hard tree/no-anchor scenarios where no group occupies a single factor

Key diagnostics:

- motif clustering accuracy from `omega`
- factor recovery from cosine matching of fitted `F` to true `F0`
- support/motif recovery from posterior motif scores
- posterior uncertainty in `gamma_bar`
- entropy and top-two gaps of `gamma_bar[s, d, ]`
- stability across random initializations

## Working Hypothesis

The paper becomes compelling if MiSo improves motif recovery and
uncertainty interpretation over per-sample MF-Poisson-SuSiE, especially in
settings with overfitted or correlated factors where hard assignments are
unstable or arbitrary.

## Current Evidence Snapshot

The strongest current positive signal is the tree/no-anchor simulation: in the
small three-seed pilot, MiSo has much better support-pattern recovery
than Poisson NMF + k-means or MF-Poisson-SuSiE + k-means, while hard clustering
is comparable.

The block simulation supports feasibility, but it is not yet a clean
support-recovery superiority result across seeds. MiSo has strong hard
clustering and factor recovery, while support accuracy is close to the baseline
pipelines in the current three-seed pilot.

Clean correctly specified settings tend to produce nearly point-mass
`gamma_bar`, which is expected under strong motif identification. The
posterior-uncertainty claim should be tested mainly in ambiguous settings:
overspecified `K`, duplicated or correlated factors, overspecified `D`, weak
motifs, and unstable hard assignments.

The first multi-seed uncertainty pilot is
`analysis/miso-uncertainty-multiseed.Rmd`. With fixed ambiguous dictionaries,
overspecified `K` shows duplicate-pair posterior splitting in selected motif
dimensions, with max split about `0.46` but mean split only about `0.08` across
all dimensions. Correlated factors show stronger targeted pair uncertainty, with mean
targeted score about `0.14` and max about `0.48`, but hard clustering is weak
at about `0.61`. These are useful uncertainty diagnostics, not final recovery
claims.

For overspecified `D`, use EB/ARD dimension activity as a primary diagnostic:
the manuscript notation is `alpha^0_{sd} / beta^0_{sd}` for the prior mean
loading size of dimension `d` in motif `s`. The code still stores these as
`alpha0[s,d] / lambda0[s,d]`. The code now has `motif_slot_activity()` and
`select_motif_dimensions()` in `code/miso.R`. In the current single-seed
over-`D` notebook, repeated-surplus initialization keeps all five dimensions active,
whereas diffuse-surplus initialization gives motif-specific selected dimensions
closer to the matched true dimensions. Dimension-scale shrinkage can reduce
posterior loading mass but leave EB prior fractions above a simple threshold, so
scale and ARD should be coupled carefully before making this a final
model-selection claim.

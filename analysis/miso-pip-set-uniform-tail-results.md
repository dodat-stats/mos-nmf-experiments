# Known-K PIP-set initialization experiment

## Setup

The experiment uses the anchor and tree simulations from
`analysis/6.miso-anchor-tree.Rmd`, with the fitted factor count fixed at its
true value:

- anchor: `K = 5`, `D = 3`, true `S = 8`;
- tree: `K = 7`, `D = 3`, true `S = 6`.

Vanilla Poisson NMF is followed by row-wise Poisson SuSiE with the NMF factor
matrix frozen. The number of fitted motifs is not supplied. It is the number
of distinct truncated PIP sets. A supported dimension starts near a point mass
on its ordered factor, while a dimension beyond the truncated support starts
uniformly over all `K` factors. The aligned row-wise Poisson-SuSiE posterior
moments initialize the MiSo Gamma priors. Each MiSo fit uses 30 iterations and
updates both the factor matrix and Gamma priors.

These are single-seed diagnostic experiments, not Monte Carlo summaries.

## Main results

| Scenario | PIP level | True S | Inferred S | Clusters with a uniform tail | Soft motif accuracy | Exact support accuracy | Support Jaccard | Mean factor cosine | Deviance per entry | Final ELBO |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Anchor | 0.95 | 8 | 10 | 0 | 0.953 | 0.992 | 0.992 | 0.990 | 0.1673 | -31,263 |
| Anchor | 0.50 | 8 | 11 | 10 | 0.883 | 0.985 | 0.988 | 0.967 | 0.1686 | -31,439 |
| Tree | 0.95 | 6 | 35 | 1 | 0.675 | 0.846 | 0.916 | 0.793 | 0.2752 | -51,682 |
| Tree | 0.50 | 6 | 39 | 22 | 0.603 | 0.808 | 0.901 | 0.772 | 0.2749 | -51,552 |

At the proposed 0.95 truncation level, the uniform-tail rule is essentially
not used: every anchor cluster has support size three, and 34 of the 35 tree
clusters have support size three. The 0.50 runs deliberately exercise the
uniform-tail rule; they are an ablation rather than a proposed default.

## What happens to uniform tail dimensions?

| Scenario | PIP level | Number of tail slots | Initial median tail mean | Initial median supported mean | Final median tail mean | Final median supported mean | Final mean tail gamma entropy |
|---|---:|---:|---:|---:|---:|---:|---:|
| Anchor | 0.50 | 10 | 1.869 | 6.578 | 0.005 | 9.827 | 0.600 |
| Tree | 0.50 | 23 | 4.047 | 16.979 | 0.014 | 11.322 | 1.000 |

This is encouraging. The empirical-Bayes updates make the unsupported uniform
slots nearly inactive, while supported dimensions retain substantial loading
mass. In the tree fit, the tail gamma entropy remains one and its mean maximum
factor probability is 0.145, essentially the uniform value `1 / 7`. MiSo can
therefore represent an inactive dimension as both small in loading magnitude
and uncertain in factor identity.

## Comparison with the existing initializations

| Scenario | Initialization | Fitted S | Soft motif accuracy | Exact support accuracy | Support Jaccard | Mean factor cosine | Deviance per entry | Final ELBO |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Anchor | Distinct | 8 | 0.684 | 0.994 | 0.994 | 0.989 | 0.1672 | -31,541 |
| Anchor | Threshold/recycle | 8 | 0.992 | 0.992 | 0.992 | 0.990 | 0.1684 | -31,344 |
| Anchor | Aligned MF | 8 | 0.992 | 0.992 | 0.992 | 0.990 | 0.1674 | -31,286 |
| Anchor | PIP sets, 0.95 | 10 | 0.953 | 0.992 | 0.992 | 0.990 | 0.1673 | -31,263 |
| Tree | Distinct | 6 | 0.660 | 0.988 | 0.992 | 0.819 | 0.2771 | -51,868 |
| Tree | Threshold/recycle | 6 | 0.968 | 0.642 | 0.793 | 0.755 | 0.2780 | -51,767 |
| Tree | Aligned MF | 6 | 0.987 | 0.990 | 0.993 | 0.773 | 0.2771 | -51,600 |
| Tree | PIP sets, 0.95 | 35 | 0.675 | 0.846 | 0.916 | 0.793 | 0.2752 | -51,682 |

The PIP-set method remains competitive for factor and observation-support
recovery in the anchor scenario. It is substantially worse for tree motif
recovery because the initial partition is fragmented. Its reconstruction
deviance remains good, showing that predictive fit alone does not diagnose
motif recovery. The 0.50 tree fit even has a relatively high ELBO despite poor
motif recovery, so the current empirical-Bayes ELBO should not be used by
itself to select the number of unregularized mixture components.

## Main diagnosis

Ordinary SuSiE PIPs ignore loading magnitude. Every Poisson-SuSiE slot must
select some factor, even when its Gamma loading has been shrunk close to zero.
Consequently, an inactive or diffuse slot can still contribute substantial PIP
to several factors.

The exact-set clustering exposes this problem combinatorially:

- at 0.95, the anchor initialization produces all 10 possible three-factor
  subsets of five factors;
- at 0.95, the tree initialization produces 34 different three-factor subsets
  of seven factors, nearly all of the 35 possible triples, plus one pair;
- at 0.50, the anchor initialization produces all 10 possible pairs plus one
  triple;
- at 0.50, the tree initialization produces all 21 possible pairs, 17 triples,
  and one singleton.

Thus uniform tail initialization is behaving as intended, but ordinary PIP
truncation followed by exact equality is not a stable way to determine `S` in
the tree scenario. A natural next experiment is to form sets from normalized
posterior expected loading

\[
R_{ik}
=
\sum_{d=1}^D
\bar\gamma^{\mathrm{PS}}_{idk}
\frac{\alpha^{\mathrm{PS}}_{idk}}{\beta^{\mathrm{PS}}_{idk}},
\]

rather than from unweighted PIPs. This retains the useful SuSiE uncertainty
but prevents a dimension with negligible loading from creating spurious
support membership. Approximate matching or merging of closely related sets
would address the separate instability caused by requiring exact equality.

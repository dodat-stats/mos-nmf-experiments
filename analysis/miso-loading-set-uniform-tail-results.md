# Expected-loading score initialization for MiSo

## Change being tested

The support score for observation $i$ and NMF factor $k$ is

\[
R_{ik}
=
\sum_{d=1}^{D}
\bar\gamma^{\mathrm{PS}}_{idk}
\frac{\alpha^{\mathrm{PS}}_{idk}}
     {\beta^{\mathrm{PS}}_{idk}}.
\]

Thus $R_{ik}$ is the posterior expected loading assigned to factor $k$,
not merely the posterior probability that factor $k$ is selected.  For each
observation, factors are ordered by $R_{ik}$, and the smallest leading set
whose cumulative share is at least \(\tau\) is retained.  Everything after
this substitution is identical to the PIP experiment: exact-set clustering,
near-point-mass supported slots, uniform unsupported tail slots, transferred
Gamma moments, simulation seeds, and MiSo fitting settings.

## Pre-MiSo support sets

The exact and Jaccard columns compare each observation's truncated initial set
with its true factor support after matching the NMF factors to the true factors.

| Scenario | Score | \(\tau\) | \(S^*\) | Exact | Precision | Recall | Jaccard |
|---|---:|---:|---:|---:|---:|---:|---:|
| Anchor | PIP | 0.95 | 10 | 0.373 | 0.583 | 0.999 | 0.582 |
| Anchor | \(R\) | 0.95 | 21 | 0.883 | 0.940 | 0.999 | 0.940 |
| Anchor | PIP | 0.50 | 11 | 0.000 | 0.686 | 0.873 | 0.561 |
| Anchor | \(R\) | 0.50 | 12 | 0.615 | 0.991 | 0.839 | 0.838 |
| Tree | PIP | 0.95 | 35 | 0.244 | 0.673 | 0.755 | 0.591 |
| Tree | \(R\) | 0.95 | 53 | 0.169 | 0.772 | 0.650 | 0.556 |
| Tree | PIP | 0.50 | 39 | 0.081 | 0.753 | 0.584 | 0.521 |
| Tree | \(R\) | 0.50 | 18 | 0.021 | 0.904 | 0.381 | 0.375 |

The expected-loading score is clearly better for the anchor supports.  For the
tree supports it favors the largest contributors and omits weaker true factors,
especially at \(\tau=0.50\).  This gives high precision but inadequate recall.

## MiSo results after convergence or 100 iterations

| Scenario | Score | \(\tau\) | \(S^*\) | Soft motif accuracy | Exact support accuracy | Support Jaccard | Mean factor cosine | ELBO | Deviance/entry |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Anchor | PIP | 0.95 | 10 | 0.955 | 0.992 | 0.992 | 0.990 | -31259.5 | 0.16727 |
| Anchor | \(R\) | 0.95 | 21 | 0.979 | 0.996 | 0.996 | 0.990 | -31259.3 | 0.16721 |
| Anchor | PIP | 0.50 | 11 | 0.884 | 0.988 | 0.989 | 0.967 | -31434.9 | 0.16856 |
| Anchor | \(R\) | 0.50 | 12 | 0.889 | 0.992 | 0.992 | 0.989 | -31264.0 | 0.16716 |
| Tree | PIP | 0.95 | 35 | 0.607 | 0.852 | 0.926 | 0.770 | -51428.2 | 0.27418 |
| Tree | \(R\) | 0.95 | 53 | 0.659 | 0.983 | 0.992 | 0.778 | -51410.5 | 0.27367 |
| Tree | PIP | 0.50 | 39 | 0.638 | 0.788 | 0.892 | 0.767 | -51387.5 | 0.27395 |
| Tree | \(R\) | 0.50 | 18 | 0.723 | 0.473 | 0.718 | 0.760 | -51413.0 | 0.27463 |

Larger ELBO and smaller deviance are better. These are single-seed diagnostic
runs. Both anchor fits at \(\tau=0.95\) converged in 36--39 iterations and the
two fits at \(\tau=0.50\) converged in 39--46 iterations. For tree, PIP at
\(\tau=0.50\), \(R\) at \(\tau=0.50\), and \(R\) at \(\tau=0.95\) converged
in 94, 94, and 89 iterations, respectively. PIP at \(\tau=0.95\) reached the
100-iteration limit without satisfying the stopping rule.

## Does the largest ELBO select the best fit?

For anchor, essentially yes. The largest ELBO is attained by \(R\) at
\(\tau=0.95\), although its advantage over PIP at \(\tau=0.95\) is only 0.17
ELBO units. It has the best soft motif accuracy, exact support accuracy,
support Jaccard, and mean factor cosine. PIP at \(\tau=0.95\) has slightly
better hard motif accuracy, and \(R\) at \(\tau=0.50\) has negligibly smaller
deviance per entry.

For tree, no. PIP at \(\tau=0.50\) has the largest ELBO, but \(R\) at
\(\tau=0.50\) has the best soft and hard motif accuracies, while \(R\) at
\(\tau=0.95\) has by far the best support recovery, factor cosine, and
deviance. The ELBO winner is close to the best predictive deviance but is not
close to the best structural recovery.

## Conclusion

Replacing PIP by $R_{ik}$ is a useful improvement for the anchor simulation:
the initial support sets are much more faithful and the final fit is at least as
good on every factor/support metric. The result is not uniformly better for the
tree simulation. At low \(\tau\), $R_{ik}$ substantially reduces the number of
candidate submanifolds but obtains that reduction by dropping weaker members
of the true support. At high \(\tau\), MiSo eventually recovers the support very
well, but the number of distinct exact sets grows to 53.

Therefore $R_{ik}$ is informative, but exact row-wise score-set clustering
still needs a stabilization or merging step for unequal multi-factor loadings.
The most promising next comparison is a hybrid rule: use PIP to decide whether
a factor belongs to the support and use $R_{ik}$ to order or weight the
supported factors.

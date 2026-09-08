# Simple-tree simulation with Gamma(2, 1/9) loadings

## Design and realized geometry

The simple-tree design is unchanged:

- `N = 480`, `M = 500`;
- true and fitted `K = 3`, `D = 2`;
- true supports `{1, 2}` and `{1, 3}`;
- active loadings are independent `Gamma(shape = 2, rate = 1/9)` draws;
- the true mean and variance of each active loading are 18 and 162.

The realized active loading mean is 17.764 and the mean total loading is
35.528. The empirical 5%, 50%, and 95% quantiles of the within-motif
composition are `(0.158, 0.504, 0.870)`, compared with `(0.368, 0.501,
0.634)` in the previous `Gamma(18, 1)` simulation.

The uncentered second-to-first singular-value ratios of the two true loading
matrices are 0.440 and 0.442. They were 0.169 and 0.154 under the previous
simulation. The corresponding ratios for the noiseless Poisson means are
0.415 and 0.400. Thus the new simulation is substantially more
two-dimensional in the nonnegative-cone geometry.

## Main recovery results

All 24 combinations of mixture-weight treatment, initialization score, and
threshold converged.

| Fit selected for | Pi treatment | Score and level | Initial `S*` | Active `S` | Two-dimensional active clusters | True-mean deviance / entry | Relative mean error | Soft motif accuracy | Mean factor cosine | Aggregate prior-mean error |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Best reconstruction | Point | PIP, 0.60 | 6 | 5 | 2 | **0.01595** | **0.21732** | 0.844 | 0.984 | 0.109 |
| Best motif recovery | Point | PIP, 0.70 | 3 | 3 | 2 | 0.01623 | 0.21744 | **0.929** | 0.987 | 0.082 |
| Dirichlet comparison | Dirichlet(1) | PIP, 0.70 | 3 | 3 | 2 | 0.01622 | 0.21738 | 0.927 | 0.987 | 0.084 |
| Best factor and prior recovery | Dirichlet(1) | Expected loading, 0.50 | 3 | 3 | 3 | 0.01638 | 0.21735 | 0.846 | **0.989** | **0.067** |

Here an active cluster has allocation mass above 0.01. A cluster is called
two-dimensional in this diagnostic when both learned Gamma prior means exceed
1.8, which is 10% of the true active-loading mean.

Vanilla rank-3 NMF has true-mean deviance 0.02611, relative mean error 0.23567,
and mean factor cosine 0.981. The best MiSo reconstruction therefore reduces
true-mean Poisson deviance by 38.9% and relative Frobenius error by 7.8%.

## Factor and Gamma-prior recovery

The broader loading distribution largely resolves the earlier factor
identifiability problem. Vanilla NMF's mean factor cosine rises from 0.804 in
the `Gamma(18, 1)` simulation to 0.981 here. The best MiSo factor cosine is
0.989.

For the structurally clearest PIP fit at level 0.70, the two dominant clusters
have masses 0.489 and 0.455. Their learned Gamma parameters are:

| Cluster mass | Dimension | Fitted alpha | Fitted beta | Fitted mean | True alpha | True beta | True mean |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.489 | 1 | 2.365 | 0.130 | 18.210 | 2 | 0.111 | 18 |
| 0.489 | 2 | 1.655 | 0.102 | 16.181 | 2 | 0.111 | 18 |
| 0.455 | 1 | 1.799 | 0.099 | 18.255 | 2 | 0.111 | 18 |
| 0.455 | 2 | 2.122 | 0.113 | 18.728 | 2 | 0.111 | 18 |

This is meaningful recovery of both the loading means and the individual
Gamma shape and rate parameters. The remaining cluster has mass 0.056 and
prior means `(33.367, 0.012)`, so it is effectively one-dimensional.

After aligning fitted factors and averaging over redundant clusters, the same
fit estimates the true motif-specific factor means as

```text
true support {1, 2}: true (18, 18, 0), fitted (16.627, 19.872, 0.072)
true support {1, 3}: true (18, 0, 18), fitted (16.246, 0.467, 17.744)
```

Its relative aggregate prior-mean error is 0.082, versus a best value of 0.717
under the previous concentrated-loading simulation.

## Mixture-component recovery

The new simulation improves motif recovery but still does not recover exactly
two mixture components. For PIP levels 0.70 through 0.95:

```text
point pi allocation masses:       (0.489, 0.455, 0.056)
Dirichlet(1) allocation masses:    (0.489, 0.453, 0.058)
```

The two dominant clusters correspond closely to the two true motifs, while a
small one-dimensional cluster remains. The Dirichlet(1) update does not remove
it and slightly increases its mass in this run.

Lower PIP thresholds and most expected-loading thresholds initialize more
clusters. These fits generally retain four to six active clusters, although
only two of them are genuinely two-dimensional. Thus broad loading directions
solve much of the factor and loading-prior recovery problem, but not automatic
selection of `S`.

The highest-ELBO fits are also not the structurally best fits. The point-pi
ELBO selects expected-loading level 0.70, with four active clusters and soft
motif accuracy 0.793. The Dirichlet-pi ELBO selects expected-loading level
0.80, with five active clusters and soft motif accuracy 0.791.

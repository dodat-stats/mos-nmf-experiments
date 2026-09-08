# Learned Gamma prior means in the simple-tree Dirichlet-pi experiment

## True generating parameters

Every active loading was generated in shape--rate parameterization as

```text
L[i,k] ~ Gamma(alpha = 18, beta = 1),
```

so its true mean is 18 and its true variance is 18. Inactive loadings are
exactly zero. The two true motifs have supports `{1, 2}` and `{1, 3}`.

The realized simulation is close to this population distribution:

| True motif | Active factor | Empirical mean | Empirical variance | Moment alpha | Moment beta |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 17.518 | 17.796 | 17.244 | 0.984 |
| 1 | 2 | 17.931 | 16.408 | 19.596 | 1.093 |
| 2 | 1 | 18.173 | 17.531 | 18.839 | 1.037 |
| 2 | 3 | 17.744 | 16.884 | 18.647 | 1.051 |

Thus failure to recover `(18, 1)` cannot be attributed to an atypical draw of
the true loadings.

For cluster `s` and dimension `d`, the learned loading prior is

```text
lambda[i,s,d] ~ Gamma(alpha0[s,d], beta0[s,d]),
E[lambda[i,s,d]] = alpha0[s,d] / beta0[s,d].
```

Clusters below are ordered by decreasing final allocation mass. Each pair is
the learned prior mean for dimensions 1 and 2. The number before each pair is
the cluster's rank by mass, not its original label.

| Initialization score | Level | Allocation mass and `(dimension 1, dimension 2)` prior means |
|---|---:|---|
| Expected loading | 0.50 | 1: mass 0.476, (35.774, 0.169)<br>2: mass 0.368, (25.776, 11.169)<br>3: mass 0.132, (34.421, 0.351)<br>4: mass 0.024, (19.584, 18.323) |
| Expected loading | 0.60 | 1: mass 0.435, (36.073, 0.017)<br>2: mass 0.396, (26.913, 10.400)<br>3: mass 0.104, (32.472, 0.251)<br>4: mass 0.065, (21.371, 14.244)<br>5: mass 0.000, (13.479, 3.323) |
| Expected loading | 0.70 | 1: mass 0.424, (36.398, 0.009)<br>2: mass 0.379, (27.439, 9.895)<br>3: mass 0.121, (33.092, 0.162)<br>4: mass 0.076, (20.878, 13.034)<br>5: mass 0.000, (13.421, 3.311) |
| Expected loading | 0.80 | 1: mass 0.423, (36.401, 0.008)<br>2: mass 0.365, (27.760, 9.668)<br>3: mass 0.131, (33.406, 0.017)<br>4: mass 0.077, (21.510, 12.397)<br>5: mass 0.004, (33.150, 0.454)<br>6: mass 0.000, (13.888, 2.824) |
| Expected loading | 0.90 | 1: mass 0.423, (36.407, 0.005)<br>2: mass 0.361, (27.972, 9.471)<br>3: mass 0.139, (33.459, 0.006)<br>4: mass 0.077, (21.865, 12.016)<br>5: mass 0.000, (14.962, 1.679) |
| Expected loading | 0.95 | 1: mass 0.425, (36.345, 0.005)<br>2: mass 0.361, (27.983, 9.450)<br>3: mass 0.139, (33.499, 0.005)<br>4: mass 0.075, (21.914, 12.232)<br>5: mass 0.000, (14.820, 1.792) |
| PIP | 0.50 | 1: mass 0.425, (36.336, 0.005)<br>2: mass 0.364, (27.815, 9.599)<br>3: mass 0.136, (33.451, 0.018)<br>4: mass 0.075, (21.961, 12.244)<br>5: mass 0.000, (14.028, 2.444)<br>6: mass 0.000, (13.586, 3.130) |
| PIP | 0.60 | 1: mass 0.425, (36.342, 0.005)<br>2: mass 0.358, (28.025, 9.377)<br>3: mass 0.142, (33.668, 0.004)<br>4: mass 0.075, (21.902, 12.254)<br>5: mass 0.000, (16.215, 1.441) |
| PIP | 0.70 | 1: mass 0.500, (27.527, 8.804)<br>2: mass 0.359, (36.248, 0.409)<br>3: mass 0.141, (27.113, 7.283) |
| PIP | 0.80 | 1: mass 0.500, (27.527, 8.804)<br>2: mass 0.359, (36.248, 0.409)<br>3: mass 0.141, (27.113, 7.283) |
| PIP | 0.90 | 1: mass 0.500, (27.527, 8.804)<br>2: mass 0.359, (36.248, 0.409)<br>3: mass 0.141, (27.113, 7.283) |
| PIP | 0.95 | 1: mass 0.500, (27.527, 8.804)<br>2: mass 0.359, (36.248, 0.409)<br>3: mass 0.141, (27.113, 7.283) |

## Main pattern

For the expected-loading initializations from 0.60 through 0.95, the four
active clusters have a stable pattern:

```text
large cluster:       approximately (36, 0)
large cluster:       approximately (27, 10)
smaller cluster:     approximately (33, 0)
smaller cluster:     approximately (21, 12--14)
```

Thus the Gamma VEB update often switches off the second dimension within a
cluster, even though the Dirichlet update does not remove the entire redundant
cluster. This is dimension pruning inside a motif rather than mixture-component
pruning.

The sum of the two learned prior means is nevertheless usually close to 36.
This agrees with the true expected total loading, `18 + 18 = 36`. For example,
the dominant cluster with prior means approximately `(36, 0)` has concentrated
almost all of that total scale in one fitted dimension rather than losing the
signal.

The means for a cluster with essentially zero allocation mass are not
identified by the data and should not be interpreted. Once its responsibility
sum becomes numerically negligible, the current VEB update leaves its previous
`alpha0` and `beta0` values unchanged.

Because the fitted factors only have cosine similarities around 0.76--0.86
with the true factors, a fitted dimension with near-zero loading does not by
itself imply that the corresponding true biological loading is absent. A
fitted factor can partially absorb signal from more than one true factor.

## Direct recovery comparison for the best fit

For the expected-loading initialization at level 0.60, the fitted Gamma
parameters are:

| Cluster rank | Mass | `alpha0[,1]` | `beta0[,1]` | Mean 1 | `alpha0[,2]` | `beta0[,2]` | Mean 2 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.435 | 22.887 | 0.634 | 36.073 | 0.737 | 44.026 | 0.017 |
| 2 | 0.396 | 16.641 | 0.618 | 26.913 | 22.579 | 2.171 | 10.400 |
| 3 | 0.104 | 27.293 | 0.841 | 32.472 | 9.867 | 39.296 | 0.251 |
| 4 | 0.065 | 24.655 | 1.154 | 21.371 | 19.440 | 1.365 | 14.244 |
| 5 | 0.000 | 28.740 | 2.132 | 13.479 | 21.705 | 6.532 | 3.323 |

These do not recover the common true values `alpha = 18` and `beta = 1`.
Some individual shapes or rates are close, but no pair of fitted dimensions in
the active clusters jointly reproduces the two true Gamma distributions.

A fairer comparison averages over the redundant fitted clusters and maps the
fitted factors to the true-factor ordering. This recovers total loading scale
but not its allocation between the shared and branch-specific factors:

| True motif | True factor-loading means | Aggregated fitted prior means |
|---:|---|---|
| support `{1, 2}` | (18, 18, 0) | (1.852, 34.173, 0.003) |
| support `{1, 3}` | (18, 0, 18) | (8.245, 0.052, 28.065) |

The fitted totals are 36.03 and 36.36, close to the true total 36, but the
factor-specific prior means are not recovered. Since a mixture of the fitted
cluster-specific Gamma distributions is not itself a single Gamma
distribution, there is no meaningful allocation-weighted `alpha` and `beta`
pair to compare with `(18, 1)`; only moments such as the mean aggregate
directly.

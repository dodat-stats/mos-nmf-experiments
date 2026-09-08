# Doubling both N and M: merge-and-refit result

## Experiment

- `N = 960`, `M = 1000`, `K = 3`, `D = 2`.
- True supports are `{1, 2}` and `{1, 3}`, with equal group probability.
- Each active loading is independently Gamma(shape = 2, rate = 1/9), so its
  true prior mean is 18.
- The three true factor rows are independently generated from Gamma(0.1, 0.01)
  and normalized to sum to one.
- A rank-three Poisson NMF and row-wise Poisson-SuSiE are fitted first.
- MiSo is initialized using PIP truncation level `0.70`.
- The same simulation and fitting seeds (`52` and `3`) are used as in the
  smaller experiment.

PIP truncation produces `S* = 3` clusters, of sizes 468, 460, and 32.

## Was MiSo refitted after merging?

Yes. The posterior responsibilities of the source cluster are added to the
target cluster. The source dimensions are first matched to the target
dimensions. The Gamma priors are then moment matched, and the posterior
factor probabilities are averaged with weights proportional to the expected
loading contributed by each cluster. This only creates an `S = 2`
initialization.

Starting from this initialization, the entire MiSo model is optimized again:
`F`, `gamma_bar`, `alpha0`, `beta0`, `omega`, and `pi` (or the variational
Dirichlet parameter) are all updated.

## Which clusters were merged?

In both versions of the mixing-weight model, internal cluster 3 has the
smallest fitted mass. Its expected factor-loading signature is nearly
orthogonal to cluster 1 but has cosine similarity about 0.73 with cluster 2.
It is therefore merged into cluster 2.

| mixing weights | source mass | target mass | cosine(source, cluster 1) | cosine(source, cluster 2) |
|---|---:|---:|---:|---:|
| point estimate | 0.04895 | 0.45048 | 0.00017 | 0.73226 |
| Dirichlet(1,1,1) | 0.05387 | 0.44555 | 0.00021 | 0.73103 |

The small cluster is a split of the `{F1, F3}` motif: one dimension has prior
mean 35.38 (35.01 under Dirichlet mixing weights) and selects true `F3`; its
other dimension has prior mean 0.013 (0.016) and diffuse posterior factor
probabilities. The merge therefore combines this nearly one-dimensional `F3`
cluster with the main `{F1, F3}` cluster.

## Results

| mixing weights | model | ELBO | change from S=3 | data deviance / entry | true-mean deviance / entry | relative true-mean error | soft motif accuracy | exact support accuracy | mean F cosine |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| point estimate | unmerged S=3 | -99933.41 | 0 | 0.141400 | 0.007990 | 0.213443 | 0.94038 | 0.99271 | 0.99132 |
| point estimate | merged and refit S=2 | -99937.84 | -4.44 | 0.141451 | 0.008341 | 0.212899 | 0.98924 | 0.99063 | 0.99295 |
| Dirichlet(1,1,1) | unmerged S=3 | -99940.28 | 0 | 0.141399 | 0.008000 | 0.213721 | 0.93545 | 0.99271 | 0.99106 |
| Dirichlet(1,1,1) | merged and refit S=2 | -99941.78 | -1.50 | 0.141450 | 0.008342 | 0.212932 | 0.98923 | 0.99063 | 0.99291 |

All four fits converged before the 100-iteration limit, and their recorded
ELBO sequences were nondecreasing.

## Conclusion

The phenomenon persists after doubling both sample size and observed
dimension. The fitted `S = 3` model keeps a non-negligible cluster of about
5% mass, and its ELBO remains higher than that of the merged-and-refitted
`S = 2` model. At the same time, merging produces the more faithful structural
description: soft motif accuracy rises from about 0.94 to 0.989, factor
recovery improves slightly, and the relative error of the group-level prior
means decreases from 0.106 to 0.076 in the point-estimate model.

Thus this larger run strengthens the interpretation from the smaller run:
the third component is a small, specialized split that the variational ELBO
prefers to retain, even though merging it gives a simpler and more accurate
recovery of the two true motifs. This is evidence from one simulated data set,
not yet a statement about the probability of over-splitting over repeated
samples.

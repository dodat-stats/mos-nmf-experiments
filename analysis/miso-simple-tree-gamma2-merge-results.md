# Merging the small cluster in the Gamma(2, 1/9) simple-tree fit

## Which factor does the small cluster represent?

In the point-pi PIP-0.70 fit, the optimal learned-to-true factor permutation is
`(3, 2, 1)`. The cluster masses and factor-aligned prior means are:

```text
mass 0.489: approximately true support {F1, F3}
mass 0.455: approximately true support {F1, F2}
mass 0.056: approximately (0, 33.37, 0)
```

Thus the small cluster represents true `F2`, not the shared true `F1`. It
receives 11.2% of the responsibility within the true `{F1, F2}` motif and
essentially none within the `{F1, F3}` motif.

## Merge initialization

The small cluster was merged separately into each dominant cluster.
Responsibilities were added exactly. Slot labels were aligned by their
expected factor-loading signatures, the categorical factor probabilities were
averaged to preserve expected factor loading, and each merged Gamma prior was
moment matched. MiSo was then refit with `S = 2`, updating `F`, `gamma_bar`, and
the Gamma prior parameters.

More precisely, let `t` be the target cluster and `s` the cluster being
removed. Define

```text
mu[c,d] = alpha0[c,d] / beta0[c,d],
v[c,d,k] = mu[c,d] * gamma_bar[c,d,k],
n[c] = sum_i omega[i,c].
```

First, the source slots are permuted to maximize

```text
sum_d sum_k v[t,d,k] * v[s,permutation(d),k].
```

For the correct merge here, the selected permutation is `(1, 2)`. The merged
responsibility is

```text
omega_new[i,t] = omega[i,t] + omega[i,s],
```

and the source column is deleted. Let

```text
w_t = n[t] / (n[t] + n[s]),
w_s = n[s] / (n[t] + n[s]).
```

For each aligned slot, the Gamma mixture is moment matched:

```text
mu_new = w_t * mu_t + w_s * mu_s,

second_moment_new =
  w_t * (alpha_t / beta_t^2 + mu_t^2)
  + w_s * (alpha_s / beta_s^2 + mu_s^2),

variance_new = second_moment_new - mu_new^2,
alpha_new = mu_new^2 / variance_new,
beta_new = mu_new / variance_new.
```

The factor-assignment probabilities are merged by preserving expected factor
loading:

```text
gamma_bar_new[d,k] =
  (w_t * mu_t * gamma_bar[t,d,k]
   + w_s * mu_s * gamma_bar[s,permutation(d),k]) / mu_new.
```

The other cluster is unchanged, `pi_init = colMeans(omega_new)`, and the
already fitted `F` is used as the starting factor matrix. All MiSo parameter
blocks are then updated normally.

For the correct merge, the source is original cluster 2 with mass 0.055855 and
the target is original cluster 1 with mass 0.454967. Hence
`(w_t, w_s) = (0.890657, 0.109343)`.

## Parameters before and after the correct merge

All `gamma_bar` vectors below are shown in true-factor order `(F1, F2, F3)`.
The original point-pi `S = 3` posterior is:

| Mass | Dimension | Alpha | Beta | Prior mean | Posterior `gamma_bar` |
|---:|---:|---:|---:|---:|---|
| 0.489178 | 1 | 2.365226 | 0.129886 | 18.210074 | (0, 0, 1) |
| 0.489178 | 2 | 1.655031 | 0.102281 | 16.181231 | (1, 0, 0) |
| 0.454967 | 1 | 1.799272 | 0.098564 | 18.254794 | (0, 1, 0) |
| 0.454967 | 2 | 2.122131 | 0.113311 | 18.728374 | (1, 0, 0) |
| 0.055855 | 1 | 2.981983 | 0.089368 | 33.367426 | (0, 1, 0) |
| 0.055855 | 2 | 0.236575 | 20.559163 | 0.011507 | (0.334071, 0.333191, 0.332739) |

Moment matching produces the following `S = 2` initialization:

| Initial mass | Dimension | Alpha | Beta | Prior mean | Initial `gamma_bar` |
|---:|---:|---:|---:|---:|---|
| 0.510822 | 1 | 1.737972 | 0.087303 | 19.907261 | (0, 1, 0) |
| 0.510822 | 2 | 1.534699 | 0.091998 | 16.681808 | (0.999950, 0.000025, 0.000025) |
| 0.489178 | 1 | 2.365226 | 0.129886 | 18.210074 | (0, 0, 1) |
| 0.489178 | 2 | 1.655031 | 0.102281 | 16.181231 | (1, 0, 0) |

After refitting with point pi, the learned posterior is:

| Final mass | Dimension | Alpha | Beta | Prior mean | Posterior `gamma_bar` |
|---:|---:|---:|---:|---:|---|
| 0.509914 | 1 | 1.784741 | 0.093158 | 19.158267 | (0, 1, 0) |
| 0.509914 | 2 | 1.546309 | 0.088619 | 17.448953 | (1, 0, 0) |
| 0.490086 | 1 | 2.344245 | 0.130795 | 17.922992 | (0, 0, 1) |
| 0.490086 | 2 | 1.732589 | 0.105318 | 16.450949 | (1, 0, 0) |

The Dirichlet-pi posterior is almost identical:

| Final mass | Dimension | Alpha | Beta | Prior mean | Posterior `gamma_bar` |
|---:|---:|---:|---:|---:|---|
| 0.509911 | 1 | 1.784687 | 0.093195 | 19.149982 | (0, 1, 0) |
| 0.509911 | 2 | 1.548487 | 0.088702 | 17.457110 | (1, 0, 0) |
| 0.490089 | 1 | 2.344290 | 0.130865 | 17.913839 | (0, 0, 1) |
| 0.490089 | 2 | 1.735110 | 0.105412 | 16.460233 | (1, 0, 0) |

## Results

| Pi treatment | Fit | Final ELBO | Change from S=3 | Data deviance / entry | True-mean deviance / entry | Relative mean error | Soft motif accuracy | Mean factor cosine | Prior-mean error |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Point | Unmerged S=3 | **-37712.44** | -- | **0.202736** | **0.016227** | 0.21744 | 0.92940 | 0.98734 | 0.08218 |
| Point | Merge into wrong `{F1,F3}` motif | -37717.14 | -4.71 | 0.202861 | 0.016290 | **0.21497** | 0.98546 | **0.99041** | 0.05757 |
| Point | Merge into correct `{F1,F2}` motif | -37717.19 | -4.76 | 0.202863 | 0.016267 | 0.21504 | **0.98565** | 0.99039 | **0.05740** |
| Dirichlet(1) | Unmerged S=3 | **-37718.48** | -- | **0.202735** | **0.016219** | 0.21738 | 0.92697 | 0.98715 | 0.08408 |
| Dirichlet(1) | Merge into wrong `{F1,F3}` motif | -37720.14 | -1.66 | 0.202860 | 0.016291 | 0.21498 | 0.98547 | 0.99038 | 0.05790 |
| Dirichlet(1) | Merge into correct `{F1,F2}` motif | -37719.92 | -1.44 | 0.202864 | 0.016280 | **0.21502** | **0.98565** | **0.99044** | **0.05709** |

The point-pi allocation changes from `(0.489, 0.455, 0.056)` to approximately
`(0.510, 0.490)`. The Dirichlet-pi result is nearly identical.

Both merge directions converge to essentially the same two-cluster solution,
because subsequent updates can reorganize the factor probabilities and
responsibilities. The geometrically correct merge is slightly better for the
Dirichlet ELBO and motif recovery, but initialization direction is not decisive
in this example.

## Interpretation

Merging does not increase ELBO. The point-pi ELBO decreases by about 4.7 and
the Dirichlet-pi ELBO decreases by 1.4--1.7. The unmerged model also achieves a
slightly lower deviance on the observed counts, explaining the ELBO preference:
the small cluster buys a modest amount of in-sample fit.

The two-cluster solution is nevertheless substantially cleaner. It nearly
perfectly recovers the motifs, improves factor cosine and factor-specific prior
means, and lowers relative Frobenius error against the true mean. Its
true-mean Poisson deviance is slightly worse, so no single recovery metric
uniformly dominates. The learned Gamma parameters reported above are
reasonably close to the true `alpha = 2`, `beta = 1/9`, and mean 18.

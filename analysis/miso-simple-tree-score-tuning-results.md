# Simple-tree initialization-score experiment

## Design

- `N = 480`, `M = 500`.
- True and fitted dimensions: `K = 3`, `D = 2`.
- Two equally likely submanifolds with supports `{1, 2}` and `{1, 3}`.
- Nonzero loadings are independent `Gamma(shape = 18, rate = 1)` draws.
- Each row of the true factor matrix is generated independently from
  `Gamma(shape = 0.1, rate = 0.01)` and normalized to sum to one.
- A shared rank-3 vanilla Poisson NMF and row-wise fixed-factor
  Poisson-SuSiE fit are used before initializing MiSo.
- The PIP and posterior expected-loading scores are tested at truncation
  levels `0.50, 0.60, 0.70, 0.80, 0.90, 0.95`.
- MiSo learns `F`, `gamma_bar`, and the Gamma prior parameters. All 12 fits
  converged.

## Reconstruction results

The main reconstruction target is the known noiseless Poisson mean
`Lambda_true = L_true %*% F_true`.

| Fit | Score level | True-mean Poisson deviance / entry | Relative Frobenius error | Data deviance / entry | Final ELBO |
|---|---:|---:|---:|---:|---:|
| Vanilla rank-3 NMF | -- | 0.02924 | 0.28577 | 0.22147 | -- |
| Best expected-loading reconstruction | 0.60 | **0.01680** | 0.22551 | 0.22367 | -40761.84 |
| Best PIP Poisson reconstruction | 0.50 | 0.01783 | 0.22712 | 0.22351 | -40750.89 |
| Best PIP Frobenius reconstruction | 0.70 | 0.01807 | **0.22306** | 0.22420 | -40823.74 |
| Highest-ELBO fit (expected loading) | 0.80 | 0.01773 | 0.22713 | 0.22351 | **-40750.87** |

Relative to vanilla NMF, the expected-loading fit at level 0.60 reduces
true-mean Poisson deviance by 42.5%. The PIP fit at level 0.70 reduces relative
Frobenius error by 21.9%.

The raw data deviance should not be used alone here: vanilla NMF attains lower
in-sample deviance than even the true Poisson mean because it also fits some
realized Poisson noise.

## Interpretation

If reconstruction of the Poisson mean is the primary target, the current best
choice is the posterior expected-loading score at level 0.60. If only ELBO is
available for tuning on real data, the expected-loading fit at level 0.80 and
the PIP fit at level 0.50 are essentially tied in this run: their ELBOs differ
by only 0.025, and both have near-best reconstruction.

The experiment does not recover the true number of submanifolds automatically.
The expected-loading fit at level 0.60 initializes five clusters and finishes
with four mixture weights above 0.01. PIP levels from 0.70 through 0.95
initialize three clusters and retain three effective clusters. Thus this run
supports the reconstruction method, but not selection of the true `S = 2`.

## Motif, factor, and mixture-weight recovery

Here `S*` is the number of distinct clusters produced by the initialization.
The three cosine columns compare the fitted rows of `F` with true factors 1,
2, and 3 after finding the permutation that maximizes the sum of their cosine
similarities. The same learned-to-true permutation, `(3, 2, 1)`, is selected in
every fit. The final mixture weights are sorted decreasingly because the
cluster labels are arbitrary.

| Score | Level | Initial `S*` | Soft motif accuracy | Cosine F1 | Cosine F2 | Cosine F3 | Mean cosine | Final sorted `pi` |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Expected loading | 0.50 | 4 | **0.844** | 0.833 | 0.742 | 0.774 | 0.783 | (0.487, 0.357, 0.143, 0.013) |
| Expected loading | 0.60 | 5 | 0.831 | **0.865** | 0.759 | 0.794 | **0.806** | (0.435, 0.395, 0.105, 0.065, 0.000) |
| Expected loading | 0.70 | 5 | 0.806 | 0.847 | 0.760 | 0.779 | 0.795 | (0.424, 0.381, 0.119, 0.076, 0.000) |
| Expected loading | 0.80 | 6 | 0.791 | 0.839 | 0.759 | 0.771 | 0.790 | (0.423, 0.367, 0.133, 0.077, 0.000, 0.000) |
| Expected loading | 0.90 | 5 | 0.778 | 0.821 | 0.758 | 0.763 | 0.781 | (0.418, 0.360, 0.140, 0.082, 0.000) |
| Expected loading | 0.95 | 5 | 0.776 | 0.819 | 0.758 | 0.762 | 0.779 | (0.418, 0.358, 0.142, 0.082, 0.000) |
| PIP | 0.50 | 6 | 0.786 | 0.832 | 0.758 | 0.768 | 0.786 | (0.423, 0.363, 0.137, 0.077, 0.000, 0.000) |
| PIP | 0.60 | 5 | 0.766 | 0.794 | 0.756 | 0.755 | 0.768 | (0.411, 0.355, 0.145, 0.089, 0.000) |
| PIP | 0.70 | 3 | 0.763 | 0.640 | 0.764 | 0.740 | 0.715 | (0.500, 0.263, 0.237) |
| PIP | 0.80 | 3 | 0.763 | 0.640 | 0.764 | 0.740 | 0.715 | (0.500, 0.263, 0.237) |
| PIP | 0.90 | 3 | 0.763 | 0.640 | 0.764 | 0.740 | 0.715 | (0.500, 0.263, 0.237) |
| PIP | 0.95 | 3 | 0.763 | 0.640 | 0.764 | 0.740 | 0.715 | (0.500, 0.263, 0.237) |

For reference, vanilla rank-3 NMF has factor cosines `(0.827, 0.797,
0.789)` in true-factor order, with mean `0.804`. Only the expected-loading fit
at level 0.60 slightly improves this mean, to `0.806`. Therefore, the much
larger improvement in reconstruction of the Poisson mean primarily comes from
the structured loading fit rather than markedly better recovery of `F`.

The soft motif accuracy is

```text
(1 / N) max over one-to-one matches h from true motifs to fitted clusters
        sum_i omega[i, h(true_motif[i])].
```

Consequently, redundant fitted clusters are treated as unmatched. The final
mixture weights explain the sub-perfect soft accuracies. For example, the
expected-loading fit at level 0.60 approximately splits the two equally sized
true motifs into weights `(0.435, 0.065)` and `(0.395, 0.105)`. Its two dominant
matched clusters therefore contain approximately `0.435 + 0.395 = 0.830` of
the posterior mass, agreeing with its soft motif accuracy of `0.831`.

The reported final mixture weights use the current maximum-likelihood update
`pi[s] = mean(omega[, s])`; they do not yet use a Dirichlet prior on `pi`.

## Code names

The public `miso()` arguments are:

- `init_pip_truncate_level` for PIP-based initialization;
- `init_loading_truncate_level` for posterior expected-loading initialization.

The lower-level initializer `init_pip_set_miso()` currently calls the shared
argument `pip_truncate_level`; `support_score = "pip"` selects PIP and
`support_score = "expected_loading"` selects the expected-loading score. In the
experiment script, the common grid column is called `truncate_level`.

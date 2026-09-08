# Simple-tree experiment with a Dirichlet prior on `pi`

## Variational update

The new implementation uses

```text
pi ~ Dirichlet(1, ..., 1),
q(pi) = Dirichlet(delta_tilde),
q(z_i) = Categorical(omega_i).
```

The coordinate updates are

```text
delta_tilde[s] = 1 + sum_i omega[i, s],

omega[i, s] proportional to exp(
  component_ELBO[i, s]
  + digamma(delta_tilde[s])
  - digamma(sum_s delta_tilde[s])
).
```

The returned posterior mean and allocation mass are, respectively,

```text
E_q[pi[s]] = delta_tilde[s] / (N + S*),
allocation_mass[s] = sum_i omega[i, s] / N.
```

The Dirichlet prior and entropy terms are included in the ELBO.

## Does the prior recover the true `S = 2`?

No. It removes components that already have very little support, but every
initialization still has either three or four components with allocation mass
above 0.01.

| Score | Level | Initial `S*` | Active with point `pi` | Active with Dirichlet `pi` | Final allocation masses, decreasing |
|---|---:|---:|---:|---:|---|
| Expected loading | 0.50 | 4 | 4 | 4 | (0.476, 0.368, 0.132, 0.024) |
| Expected loading | 0.60 | 5 | 4 | 4 | (0.435, 0.396, 0.104, 0.065, 0.000) |
| Expected loading | 0.70 | 5 | 4 | 4 | (0.424, 0.379, 0.121, 0.076, 0.000) |
| Expected loading | 0.80 | 6 | 4 | 4 | (0.423, 0.365, 0.131, 0.077, 0.004, 0.000) |
| Expected loading | 0.90 | 5 | 4 | 4 | (0.423, 0.361, 0.139, 0.077, 0.000) |
| Expected loading | 0.95 | 5 | 4 | 4 | (0.425, 0.361, 0.139, 0.075, 0.000) |
| PIP | 0.50 | 6 | 4 | 4 | (0.425, 0.364, 0.136, 0.075, 0.000, 0.000) |
| PIP | 0.60 | 5 | 4 | 4 | (0.425, 0.358, 0.142, 0.075, 0.000) |
| PIP | 0.70--0.95 | 3 | 3 | 3 | (0.500, 0.359, 0.141) |

Here an active component has allocation mass greater than 0.01. The true
mixture weights are `(0.5, 0.5)`.

The posterior mean of `pi[s]` cannot be exactly zero under a proper
`Dirichlet(1, ..., 1)` prior. A component with zero expected observations has
posterior mean `1 / (N + S*)`, which is approximately 0.0021 in this
experiment. The allocation mass can be zero, so it is the more direct measure
of whether a component has been emptied.

## Recovery performance

The best Poisson-mean reconstruction remains the expected-loading
initialization at level 0.60.

| Quantity | Point estimate of `pi` | Dirichlet `q(pi)` |
|---|---:|---:|
| True-mean Poisson deviance / entry | 0.016800 | 0.016778 |
| Relative true-mean Frobenius error | 0.225511 | 0.225564 |
| Soft motif accuracy | 0.830817 | 0.831619 |
| Mean optimally matched factor cosine | 0.806107 | 0.806235 |

For the PIP initialization at levels 0.70 through 0.95, the Dirichlet update
makes the redundant split more unequal but does not eliminate it:

```text
point pi:             (0.500, 0.263, 0.237)
Dirichlet allocation: (0.500, 0.359, 0.141)
```

Its soft motif accuracy improves from 0.763 to 0.859, but the third component
still carries 14.1% of the observations.

## Interpretation

The `Dirichlet(1, ..., 1)` prior is uniform on the simplex; it does not directly
prefer a boundary point. The variational expected-log-weight update can empty a
component that is already weak, but it cannot overcome the likelihood and
component-prior gain from a stable split of a true motif in this run. Thus this
experiment does not provide evidence that this prior alone selects `S`.

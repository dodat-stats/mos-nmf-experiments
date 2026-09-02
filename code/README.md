# MiSo R implementation

The frozen public API has three fitting functions:

1. `poisson_susie(y, F, alpha0, beta0, ...)` fits Poisson SuSiE to one
   observation with a fixed factor matrix.
2. `poisson_susie_nmf(Y, K, D, ...)` fits Poisson-SuSiE-NMF and learns the
   shared factor matrix.
3. `miso(Y, K, S, D, ...)` fits the mixture-of-submanifolds model.

For an interactive session, source `code/miso.R`; it loads all three public
functions. The main implementation files are `poisson-susie.R`,
`joint-learn-susie-poi-F.R`, and `miso.R`.

## MiSo initialization

By default, `miso()` clusters row-normalized preliminary loadings from the
Poisson-SuSiE-NMF fit and initializes each motif with the top `D` distinct
factors in its cluster center. Thus the default requires `D <= K` and avoids
collapsing several dimensions onto one repeated factor. Use
`motif_initialization = "threshold"` only to reproduce the older
threshold-and-recycle initializer. For that legacy initializer,
`motif_min_share` controls the threshold and `surplus_slots` controls whether
recycled dimensions are repeated or initialized uniformly.

## Notation

`lambda` denotes a nonnegative latent loading. `beta`, `beta0`, and `beta_sd`
denote Gamma rates. New fit objects use only the rate names `beta`, `beta0`,
and `beta_sd`. Post-processing helpers can still read pre-freeze saved fits
whose rate fields were named `lambda`, `lambda0`, or `lambda_sd`.

The old names `poi_susie()`, `mf_poi_susie()`, and
`mf_poi_susie_fixed_F()` remain as compatibility wrappers, but they are not
part of the documented public API.

`mf-poi-susie.R`, `mv-poi-susie.R`, and `miso-poi-susie.R` are historical
experiment prototypes. They are retained only so old notebooks remain
reproducible and are excluded from the frozen implementation.

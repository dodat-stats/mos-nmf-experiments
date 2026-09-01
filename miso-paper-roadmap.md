---
editor_options: 
  markdown: 
    wrap: 72
---

# MiSo Paper Roadmap

Publication name for manuscript authorship: Dat Do.

## Paper Thesis

The strongest paper is about interpretable discovery of recurring factor-usage
motifs in NMF, not only about improving reconstruction. The central claim
should be:

> MiSo learns repeated factor-usage motifs, clusters samples by
> those motifs, and quantifies motif-level uncertainty over the factors
> that define each motif.

The model-level distinction is:

-   MF-Poisson-SuSiE gives sample-level sparse loading uncertainty.
-   MiSo pools that structure across samples through motif-level
    latent submanifolds.
-   MiSo represents motif dimensions with posterior motif
    probabilities `gamma_bar[s, d, k]`.

## Current Model Stack

| Layer | Main artifact | Role |
|------------------------|------------------------|------------------------|
| Single Poisson-SuSiE | `analysis/1.poisson-susie.Rmd` | Establishes sparse Poisson single-vector inference with SuSiE-like factor uncertainty. |
| MF-Poisson-SuSiE | `code/joint-learn-susie-poi-F.R`, `analysis/2.joint-learn-susie-poi-F.Rmd` | Learns `F` and per-sample sparse loading supports. Provides initialization for MiSo. |
| MiSo | `tex/miso-main.tex`, `code/miso.R`, `analysis/miso-alternating.Rmd` | Target mixture-of-submanifolds model with global motif posterior `gamma_bar[s, d, k]`. |

## Evidence Map

### Correct Specification / Learned F

Primary artifacts:

-   `analysis/joint-learn-susie-poi-F.Rmd`
-   `analysis/miso-alternating.Rmd`
-   `analysis/miso-baseline-benchmark.Rmd`
-   `analysis/miso-multiseed-benchmark.Rmd`

Current evidence:

-   MF-Poisson-SuSiE can recover factors and loading supports in the
    block scenario after Poisson NMF initialization.
-   MiSo improves the block scenario in the single-seed
    benchmark: clustering, support recovery, and factor recovery are all
    better than Poisson NMF + k-means and MF-Poisson-SuSiE + k-means.
-   The three-seed pilot is more mixed in the block scenario. MiSo
    has the best mean hard clustering and factor cosine, but
    support recovery is close to or slightly below the baseline
    pipelines in that smaller pilot.
-   In clean settings, `gamma_bar` often collapses to point masses. This
    is not a problem by itself; with many samples per motif, motif
    identities can be strongly identified.

Paper status:

-   This supports the basic feasibility claim.
-   It does not yet prove a robust support-recovery advantage across
    seeds or settings.

Next checks:

-   Expand the pilot benchmark beyond three random seeds.
-   Add held-out Poisson log likelihood or reconstruction deviance so
    the paper does not rely only on clustering/support metrics.

### Hard Tree / No-Anchor Scenario

Primary artifacts:

-   `analysis/tree-scenario.Rmd`
-   `analysis/miso-alternating.Rmd`
-   `analysis/miso-baseline-benchmark.Rmd`
-   `analysis/miso-multiseed-benchmark.Rmd`

Current evidence:

-   The no-anchor tree scenario is hard for MF-Poisson-SuSiE factor
    recovery.
-   In the compact benchmark, MiSo improves support accuracy in
    the tree scenario relative to both baselines.
-   However, hard cluster accuracy from `max.col(omega)` is worse than
    k-means on the baseline loading scores in the current seeded
    benchmark.
-   In the three-seed pilot, the support-recovery advantage is more
    stable: MiSo averages about `0.76` support accuracy, while
    both baseline pipelines are around `0.41`. Hard clustering is
    comparable in that smaller pilot.

Paper status:

-   This is currently the strongest positive recovery setting for
    MiSo. It suggests motif sharing helps most when sparsity patterns
    are combinations of factors and no group is explained by a single
    anchor factor.
-   The hard-label responsibility story still needs more diagnostics
    because the single-seed and multi-seed pilots do not give identical
    cluster conclusions.

Next checks:

-   Evaluate soft cluster quality using group responsibility matrices,
    not only hard labels.
-   Add branch-aware diagnostics for shared tree structure.
-   Test initialization variants for `omega` and `gamma_bar`.

### Overspecified K

Primary artifacts:

-   `analysis/overspecify-K-mf-poi-susie.Rmd`
-   `analysis/miso-uncertainty-stress.Rmd`
-   `analysis/miso-uncertainty-multiseed.Rmd`

Current evidence:

-   MF-Poisson-SuSiE with overfitted factors can create
    duplicate/competing rows.
-   MiSo can express uncertainty across duplicate candidate
    factor rows in an overcomplete fixed dictionary.
-   The strongest diagnostic is posterior mass split over known
    duplicate pairs, not raw entropy alone.
-   In a three-seed fixed-dictionary pilot, most motif dimensions remain
    concentrated as expected, but the largest duplicate-pair posterior
    split reaches about `0.46`; the mean split across all dimensions is about
    `0.08` because many dimensions do not involve the duplicated factors.

Paper status:

-   This is one of the most promising uncertainty results.
-   It directly supports the claim that soft motif posteriors can
    represent ambiguity that hard motif assignments cannot.
-   The evidence should be presented using pair-specific diagnostics
    over relevant dimensions, rather than global mean entropy over every
    motif dimension.

Next checks:

-   Repeat over multiple duplicate strengths and larger seed grids.
-   Compare against MF-Poisson-SuSiE summaries and hard support summaries
    using the same duplicate-pair uncertainty metric.

### Overspecified D

Primary artifacts:

-   `analysis/overspecify-D-mf-poi-susie.Rmd`
-   `analysis/miso-overspecify-D.Rmd`

Current evidence:

-   MF-Poisson-SuSiE handles overspecified `D` reasonably when extra
    effects have low loading mass.
-   MiSo with repeated surplus motif dimensions does not
    automatically shrink extra dimensions. It can split one true dimension
    across repeated fitted dimensions.
-   Diffuse surplus initialization greatly reduces extra effective beta
    mass, but can reduce support accuracy.
-   The MiSo code now exposes an ARD-style diagnostic through
    the empirical-Bayes prior mean `alpha^0_{sd} / beta^0_{sd}` in the
    manuscript notation; in code this is stored as `alpha0[s,d] / beta0[s,d]`. In the
    current single-seed over-`D` notebook, this diagnostic correctly
    shows that repeated surplus initialization leaves all five dimensions
    active, while diffuse-surplus initialization selects a smaller
    motif-specific `D` that is close to the matched true dimension for
    several submanifolds.
-   The experimental dimension-scale update suppresses extra dimensions further
    and makes surplus `gamma_bar` diffuse, but in the current seeded run
    it badly degrades hard clustering accuracy. It also shows that
    dimension-scale shrinkage and ARD prior means must be coupled carefully:
    posterior loading mass can be small even when the EB prior fraction remains above
    a simple active threshold.

Paper status:

-   This is currently a model-development gap, not a solved result.
-   Overspecifying `D` should be presented carefully until an explicit
    shrinkage mechanism is stable.
-   The ARD diagnostic is a promising model-selection summary: each
    submanifold can have its own selected `D_hat_s`, but the threshold
    and interaction with initialization need multi-seed validation.

Next checks:

-   Validate the ARD selector over multiple seeds and thresholds.
-   Compare ARD selection with the diffuse-posterior heuristic and with
    true motif dimensions after matching submanifolds to groups.
-   Avoid relying only on initialization for surplus-dimension shrinkage.
-   Track both effective beta mass and support recovery when testing any
    shrinkage mechanism.

### Correlated Factors

Primary artifact:

-   `analysis/miso-uncertainty-stress.Rmd`
-   `analysis/miso-uncertainty-multiseed.Rmd`

Current evidence:

-   Correlated factors are a difficult uncertainty stress test.
-   Raw entropy can be misleading because some dimensions become globally
    diffuse, which is not the same as targeted uncertainty between the
    correlated factors.
-   Pair-specific diagnostics such as `mass_factor_1_or_2` and
    `split_between_1_and_2` are more meaningful.
-   In a three-seed fixed-dictionary pilot, the mean targeted
    pair-uncertainty score is about `0.14` and the maximum reaches about
    `0.48`, but hard clustering accuracy is only about `0.61`. This
    supports the diagnostic value of `gamma_bar`, not a recovery-success
    claim.

Paper status:

-   The current result is not yet a clean positive demonstration.
-   It is still useful for defining the right uncertainty diagnostics.
-   The paper should treat correlated factors as a
    calibration/stress-test setting unless a cleaner simulation sweep
    produces controlled targeted ambiguity with better recovery.

Next checks:

-   Sweep factor correlation strength.
-   Reduce signal or sample size to create controlled ambiguity rather
    than total nonidentification.
-   Compare targeted pair uncertainty to hard assignment instability
    across random initializations.

## Proposed Simulation Table For Paper

| Simulation | Main claim | Current status |
|------------------------|------------------------|------------------------|
| Block, correct `K,D`, learned `F` | Basic recovery and motif clustering | Positive in current seeded benchmark. |
| Tree/no-anchor | Motif sharing helps support recovery without anchor groups | Mixed: support improves, hard clustering weaker. |
| Overspecified `K` with duplicate factors | `gamma_bar` captures duplicate-factor uncertainty | Promising; three-seed fixed-dictionary pilot added. |
| Overspecified `D` | Extra dimensions should be inactive or uncertain | ARD diagnostic added; shrinkage and thresholding still need validation. |
| Correlated factors | `gamma_bar` captures targeted ambiguity | Multi-seed stress pilot added; uncertainty diagnostic works better than recovery. |
| Multiple seeds | Robustness and stability | Pilots added for block/tree and uncertainty stress tests; needs larger seed grid. |

## Methodological Gaps Before Paper

1.  Stabilize or formalize surplus-dimension shrinkage.
2.  Expand multi-seed benchmarks with compact tables.
3.  Add a real-data example where factor-usage motifs have substantive
    meaning.
4.  Decide model-selection guidance for `S`, `D`, and `K`.
5.  Clarify the variational objective when using practical updates such
    as alternating `F` or dimension scales.
6.  Add uncertainty calibration/stability diagnostics, not only point
    recovery.

## Suggested Manuscript Structure

1.  Introduction: NMF lacks motif-level interpretation of recurring
    factor-usage patterns.
2.  Background: Poisson NMF, sparse NMF, SuSiE-style support
    uncertainty.
3.  Model: mixture of low-dimensional submanifolds, spanned by subsets of
    potentially dense factors, with random motif indicators
    `gamma_sd`.
4.  Variational inference: `omega`, `gamma_bar`, Gamma loading
    posteriors, and blockwise `xi` computation.
5.  Simulations:
    -   basic recovery
    -   no-anchor/tree setting
    -   overspecified `K`
    -   overspecified `D`
    -   correlated factors
6.  Real data application.
7.  Discussion: motif uncertainty, model selection, identifiability, and
    limitations.

## Near-Term Priority

The next highest-value technical step is to stabilize dimension
selection for overspecified `D`, using EB/ARD prior means together with
posterior loading mass and diffuse-dimension diagnostics. The next
highest-value empirical step is to expand the multi-seed benchmarks
beyond the current small pilots: larger seed grids, duplicate-strength
sweeps, correlation-strength sweeps, ARD-threshold sweeps, and direct
comparison with MF-Poisson-SuSiE uncertainty summaries and hard support
summaries.

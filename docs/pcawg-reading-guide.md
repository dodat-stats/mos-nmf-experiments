# PCAWG SBS96: data, scientific goal, and NMF analysis

This is a focused reading guide for the MiSo pilot in
`analysis/real-data-exploration/12.pcawg-miso.R`.

## 1. What the data are

The fitted input is an integer count matrix

\[
Y \in \mathbb{N}_0^{N \times 96},
\]

where a row is a tumour and a column is a single-base-substitution (SBS)
category. The 96 categories combine six pyrimidine-oriented substitutions
(`C>A`, `C>G`, `C>T`, `T>A`, `T>C`, and `T>G`) with the 16 possible immediate
5-prime and 3-prime sequence contexts.

The downloaded matrix contains 2,780 PCAWG whole-genome spectra from 37 cancer
types. The current pilot retains 2,438 tumours from 36 types after requiring
500--50,000 SBS mutations and excluding skin melanoma. This exploratory filter
prevents sparse samples and the most extreme UV/hypermutated regime from
dominating the first fit; it is not part of the MiSo model.

An entry `Y[i,j]` is the observed number of mutations of category `j` in tumour
`i`. We fit these raw counts. We do not normalize each tumour to sum to one.

## 2. The scientific question

Each tumour accumulates mutations from several biological processes, such as
age-associated damage, tobacco exposure, APOBEC activity, or defective DNA
repair. Each process produces a characteristic distribution over the 96
categories, called a mutational signature.

Standard signature analysis asks:

1. Which signature profiles are present?
2. How many mutations in each tumour are attributable to each signature?

MiSo adds a third question:

3. Are tumour exposure vectors organized into recurring low-dimensional
   submanifolds (motifs), and what relationships among mutational processes do
   those motifs reveal?

This is a natural MiSo problem because exposure magnitudes should remain
unconstrained: a tumour with twice the mutational burden should not be forced
onto the same simplex as a lower-burden tumour.

## 3. The usual NMF fit

With the orientation used in this repository,

\[
Y \approx L F,
\]

where `F` is a `K x 96` matrix of mutational-signature profiles and `L` is an
`N x K` matrix of nonnegative, absolute signature exposures. Equivalently, the
Poisson model is

\[
Y_{ij} \sim \operatorname{Poisson}\!\left(\sum_{k=1}^K L_{ik}F_{kj}\right).
\]

SigProfiler uses NMF with generalized Kullback--Leibler divergence, many random
restarts, and human assessment of solution stability and reconstruction
accuracy over candidate ranks. SignatureAnalyzer instead uses Bayesian NMF
with automatic relevance determination. PCAWG applied the methods
hierarchically and treated low-burden and hypermutated samples separately to
reduce signature bleeding.

The PCAWG authors emphasize that rank selection is not fully automatic. In
their low-burden analysis, SigProfiler found 31 SBS signatures and
SignatureAnalyzer found 35. This motivates `K = 36` as a transparent pilot
choice here; it is not a claimed estimate of the true rank.

## 4. What this pilot does

- Fit only the observed SBS96 counts. Cancer labels and published PCAWG
  signatures/exposures are not supplied to NMF or MiSo.
- Diagnose `K` with whole-tumour validation: fit factors on training tumours,
  infer validation loadings from half of each validation spectrum, and score
  the other half.
- Fit vanilla Poisson NMF and MiSo to the same 80% Poisson thinning of all
  retained spectra and score the remaining 20%.
- Use an intentionally over-specified `S = 30` and `D = 5` for MiSo.
- Only after fitting, match learned factor profiles to published PCAWG
  signatures by cosine similarity, compare exposures by Spearman correlation,
  and use cancer type as a descriptive label.

The validation curve continues improving through the largest displayed
candidate (`K = 64`). An almost category-specific basis can predict the held-out
half of a tumour once its loading is inferred from the other half, so predictive
likelihood alone does not solve rank selection in this setting.

## 5. Initial result and its limitations

The full pilot (`K = 36`, `S = 30`, `D = 5`) currently shows:

- median cosine similarity 0.881 between learned factors and their closest
  published PCAWG signatures, versus 0.854 for vanilla NMF; 27 of 36 MiSo
  factors have cosine at least 0.85 (19 of 36 for NMF);
- only 12 distinct closest reference matches for MiSo, versus 24 for NMF,
  showing that MiSo currently duplicates a smaller set of broad processes;
- 1.40 active MiSo factors per tumour for 90% loading mass, compared with 18.88
  for vanilla NMF;
- a mixture-weighted post-hoc 90%-mass motif dimension of 1.45;
- held-out deviance per entry 2.015 for MiSo versus 1.290 for vanilla NMF;
- median Spearman correlation 0.270 after aggregating learned factors matched
  to the same published signature; and
- cancer-type NMI 0.436, burden-decile NMI 0.195, and motif-wise log-burden
  R-squared 0.557.

The result is promising for visualizing continuous radial and angular variation
inside recurring two-ray motifs, but it is not yet publication-ready. In
particular, overfitted `K` splits broad signatures (14 learned factors match
SBS40 and 7 match SBS5), some fitted motifs are redundant, held-out prediction
is worse than unconstrained NMF, and the current Gamma-prior rule does not
automatically truncate the fitted `D = 5` even though the post-hoc 90%-mass
summary is usually one or two. Factor truncation/merging and motif merging must
therefore be finalized before drawing biological conclusions.

## 6. Reading order

1. **Alexandrov et al. (2020), “The repertoire of mutational signatures in
   human cancer.”** This is the main paper. Read the Abstract; “Mutational
   signature analysis”; Figures 1--4; Methods sections for SigProfiler and
   SignatureAnalyzer; and Data availability. It defines SBS96, explains the
   biological goal, describes both NMF fits and their rank-selection issues,
   and links the exact observed spectra, signatures, and exposures used here.
   <https://doi.org/10.1038/s41586-020-1943-3>

2. **PCAWG Consortium (2020), “Pan-cancer analysis of whole genomes.”** Read
   the Abstract and “The pan-cancer analysis of whole genomes” for cohort
   construction, uniform variant calling, quality control, and the distinction
   between the 2,658 donor-level flagship cohort and the 2,780 spectra used by
   the mutational-signature analysis.
   <https://doi.org/10.1038/s41586-020-1969-6>

3. **PCAWG7 data package.** This is the convenient open mirror used by the
   script for the observed SBS96 matrix and published validation files.
   <https://github.com/steverozen/PCAWG7>

4. **COSMIC Mutational Signatures.** Use the SBS catalogue to interpret the
   biological process and evidence associated with a matched signature. Treat
   unknown or flat signatures, especially SBS5 and SBS40, cautiously.
   <https://cancer.sanger.ac.uk/signatures/sbs/>

5. **SigProfilerMatrixGenerator.** Read this only if starting from raw VCF/MAF
   calls rather than the published SBS96 matrix. It explains how variants are
   classified into mutational count matrices.
   <https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator>

6. **SignatureAnalyzer implementation.** This is useful for comparing MiSo
   with the PCAWG Bayesian ARD-NMF baseline and for understanding their sparse
   attribution workflow.
   <https://github.com/getzlab/SignatureAnalyzer>

The original open observed-spectra accession cited by the paper is
<https://www.synapse.org/#!Synapse:syn11801889>. The analysis script uses the
PCAWG7 mirror because it gives stable, checksum-verifiable individual files.

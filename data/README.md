# Data

Raw datasets are downloaded into subdirectories here and are deliberately not
tracked by git.

## PCAWG SBS96 pilot

Run from the project root:

```sh
Rscript --vanilla analysis/12.pcawg-miso.R
```

The script creates `data/pcawg/` and downloads checksum-verified copies of:

- the observed 96-category SBS spectra for 2,780 PCAWG whole genomes;
- the PCAWG sample sheet;
- published SigProfiler SBS signatures and per-tumour exposures; and
- the COSMIC v3.2 short aetiology annotations.

The observed spectra are used for fitting. Published signatures, exposures,
cancer labels, and aetiologies are used only after fitting for interpretation
and external validation. The files are mirrored by the
[PCAWG7 data package](https://github.com/steverozen/PCAWG7), whose purpose is to
distribute data from Alexandrov et al. (2020).

Use `--quick` for a short end-to-end check:

```sh
Rscript --vanilla analysis/12.pcawg-miso.R --quick
```

See [`docs/pcawg-reading-guide.md`](../docs/pcawg-reading-guide.md) for the data
definition, NMF formulation, pilot design, and recommended reading.

## version 1.1.0

---

### Disassortative mating

- `am_simulate()` and `am_covariance_structure()` now support negative
  (disassortative) cross-mate correlations. The equilibrium covariance becomes
  diagonal minus low rank in that case, tracked by an `attr(U, "sign")` that
  `rb_dplr()` honors automatically.
- fixed `am_covariance_structure()` returning `NaN` at `r = 0`
- `rb_dplr()` gains a `sign` argument, and its infeasibility error now names
  the offending locus and suggests concrete remedies
- documented that negative `r` leaves the Bahadur feasible region sooner than
  positive `r`, and that the magnitude of `r` that samples reliably depends on
  sample size as well as `min_MAF`; see `?am_covariance_structure` for
  measured envelopes

### Genotype output

- `am_simulate()` gains `path`, `format`, and `batch_size` for streaming
  genotypes to disk in batches, removing the full-matrix allocation
- new `write_genotypes()` and `read_genotypes()` supporting individual-major
  int8, variant-major int8, and PLINK bed
- README documents reading the int8 and bed output into Python

### Command line interface

- new `rbahadur` executable, shipped in `exec/`. `rbahadur simulate` streams a
  simulation to disk without opening R, and `rbahadur info` inspects an
  existing run and verifies the data file against its metadata. Exit status
  separates usage errors (1) from runtime failures (2).
- `rbahadur_cli_path()` returns the location of that script so it can be put
  on the search path; `rbahadur_main()` exposes the same interface from R
- `--csv` writes portable `_pheno.csv` and `_variants.csv` sidecars alongside
  the R-only `.rds`, for pipelines that continue outside R

### Other

- added a testthat suite
- `utils` added to Imports

## version 1.0.0

---

- citation updated after publication
- typos fixed in documentation

## version 0.9.2

---

- am_simulate() can use user-specified allele frequencies
- am_simulate() can now return optionally return haplotypes

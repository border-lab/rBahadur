## version 1.1.0

---

### Disassortative mating

- `am_simulate()` and `am_covariance_structure()` now support negative
  (disassortative) cross-mate correlations. The equilibrium covariance becomes
  diagonal minus low rank in that case, tracked by an `attr(U, "sign")` that
  `rb_dplr()` honors automatically.
- fixed `am_covariance_structure()` returning `NaN` at `r = 0`
- fixed `h2_eq()` returning `NaN` at `r = 0`, where the answer is just `h2_0`,
  since no assortment leaves heritability where it started. The closed form
  divided by `2r` when the bracket it multiplied also vanished there, and the
  same cancellation made small `r` ill conditioned: at `r = 1e-12` it was wrong
  in the fifth digit. It is now defined through `vg_eq()`, as equilibrium
  genetic variance over equilibrium total variance, which divides by nothing
  that vanishes. `rbahadur simulate --r 0` previously reported
  `equilibrium h2 NaN`.
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
- genotype matrices are now checked before integer conversion, so fractional
  or non-finite dosages are rejected rather than silently truncated
- genotype metadata is now parsed strictly: unknown layouts, incompatible
  dtypes, invalid dimensions, duplicate fields, and unsupported versions are
  reported as corruption instead of being guessed; int8 readers and
  `rbahadur info` also detect dosage bytes outside 0, 1, and 2
- PLINK output from `am_mosaic()` now carries the reference panel's chromosome,
  physical position, genetic map, marker ID, and allele metadata when present

### Local LD

- new `am_mosaic()`, which combines the genome-wide linkage disequilibrium
  induced by assortative mating with the local linkage disequilibrium induced
  by limited recombination. Causal variants are drawn with `rb_dplr()` and the
  intervening markers are filled by copying contiguous blocks from a reference
  panel, following Algorithm S4 of the supplementary note. Unlike the published
  vignette, block boundaries are drawn from a genetic map rather than uniformly,
  so breakpoints concentrate where recombination occurs.
- `am_mosaic()` accepts the same `path`, `format`, and `batch_size` arguments as
  `am_simulate()`. The variant-major layouts stream over markers and the
  individual-major layout streams over people, so no full genotype matrix is
  held either way.
- new `kg_reference()` returning a bundled 1000 Genomes panel (520 haplotypes
  across 2500 common SNVs in a 1 Mb window of chromosome 22, GRCh38, with
  genetic map positions), so examples and tests exercise real LD offline
- new `vcf_to_panel()` to build a panel from a phased VCF and a PLINK genetic
  map, and `download_1kg_panel()` to fetch a region of 1000 Genomes directly
- `vcf_to_panel()` also accepts unphased calls, with a prominent warning that
  their supplied allele order is treated as phase and can create artificial LD
- VCF records are now width-checked independently, symbolic alleles are not
  mistaken for SNVs, marker positions must be ordered, and genetic-map
  interpolation is restricted to the retained chromosome
- streamed reference downloads are written to temporary files and published
  to the cache only after the requested interval has demonstrably completed;
  map archives and extracted maps are likewise published atomically
- new vignette walking through the method, verifying that the assortative
  mating structure is preserved exactly and that local LD matches the panel,
  and documenting the tradeoff between block length and the infinitesimal
  limit underlying `vg_eq()`

### Performance

Except for the corrected one- and two-variable `rb_dplr()` edge cases, output
on previously tested simulation dimensions is unchanged: it is bit-identical
under the same seed across all formats and batch sizes and leaves the random
number stream in the same state.

- `am_mosaic()` expands markers 3 to 6 times faster and in roughly a third of
  the memory. Measured on a 2,000 haplotype by 100,000 marker panel with
  n = 1,000 and 400 causal variants: in memory 19.3s to 5.7s, individual-major
  18.7s to 7.3s, variant-major 31.0s to 5.9s, and bed 36.3s to 6.0s, with peak
  resident memory falling from 4.0 GB to 1.3 GB.
- reference panels are now held as one byte per allele rather than being
  expanded to integer on entry, which was four times the memory and, on a
  genome-scale panel, more time than the simulation it was preparing for
- the variant-major and bed layouts no longer rescan every individual at every
  marker to find who has crossed a block boundary; boundaries are compiled once
  into a marker-ordered schedule
- the individual-major layout copies blocks out of a transposed panel, so a
  haplotype is a sequential read rather than a stride over the panel
- `rb_dplr()` draws its uniforms one locus at a time rather than allocating an
  n by m matrix of them up front, the same equivalence `am_simulate()`'s
  streaming path already relied on
- PLINK bed packing recodes by lookup and packs a block of variants per call
  instead of one call per variant, which also speeds up `write_genotypes()`

### Command line interface

- new `rbahadur` executable, shipped in `exec/`. `rbahadur simulate` streams a
  simulation to disk without opening R, and `rbahadur info` inspects an
  existing run and verifies the data file against its metadata. Exit status
  separates usage errors (1) from runtime failures (2).
- `rbahadur_cli_path()` returns the location of that script so it can be put
  on the search path; `rbahadur_main()` exposes the same interface from R
- `--csv` writes portable `_pheno.csv` and `_variants.csv` sidecars alongside
  the R-only `.rds`, for pipelines that continue outside R
- unknown or misspelled command-line options and extra positional arguments
  are now rejected as usage errors instead of being silently ignored
- duplicate options are rejected, and `--seed 0` is accepted as a valid R seed

### Other

- added a testthat suite
- `utils` added to Imports
- simulations with fewer than 50 causal variants remain supported but now warn
  that the equilibrium variance and heritability formulas are large-locus
  targets; the warning can be disabled with
  `options(rBahadur.warn_small_m = FALSE)`
- declared and preflight-checked the external commands used by
  `download_1kg_panel()`
- simulation, Bernoulli, equilibrium, reference-panel, metadata, and PLINK
  inputs now receive explicit domain and integer validation; the broken
  one- and two-variable recursion edges in `rb_dplr()` and `rb_unstr()` are
  handled directly

## version 1.0.0

---

- citation updated after publication
- typos fixed in documentation

## version 0.9.2

---

- am_simulate() can use user-specified allele frequencies
- am_simulate() can now return optionally return haplotypes

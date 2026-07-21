# Design: negative assortative mating and binary int8 genotype output

Date: 2026-07-21
Package: rBahadur (currently 1.0.0, on CRAN)

## Summary

Two independent features:

1. Support negative (disassortative) mating, meaning cross-mate phenotypic
   correlation `r < 0`, which currently produces `NaN` and aborts.
2. Write simulated genotypes to disk as binary int8, with an optional streaming
   path so simulations are no longer bounded by available RAM.

A third, smaller item is folded in because it shares a code path: `r = 0`
(panmixia) is also broken today and returns all `NaN`.

## Part 1: negative assortative mating

### Current behavior

`am_covariance_structure(beta, AF, r)` computes the rank-one loading vector `U`
such that the equilibrium haploid LD matrix is `C = D + U U^T`. The expression
contains `sqrt(r)` in a denominator, so every entry of `U` is `NaN` when
`r < 0`, and `am_simulate()` then fails inside `rb_dplr()` with
"missing value where TRUE/FALSE needed". When `r = 0` the same expression
divides by zero and yields `NaN` as well.

### Why a sign flip is required, not a formula patch

Let `a_k` be the raw (unstandardized) haploid substitution effects. Under the
DPLR structure the genetic variance at equilibrium satisfies, up to an `O(1/m)`
diagonal correction,

    vg_eq - vg_0 = (sum_k a_k U_k)^2

and the Nagylaki equilibrium relation gives `vg_eq - vg_0 = vg_eq * rg_eq`.
Therefore

    (sum_k a_k U_k)^2 = vg_eq * rg_eq.

For `r < 0` we have `rg_eq < 0`, so the right hand side is negative and no real
`U` exists. This is not a numerical artifact. A positive semidefinite rank-one
term `U U^T` can only ever *increase* genetic variance, whereas negative
assortment *reduces* it. The low-rank term must therefore be subtracted:

    C = D - V V^T   for r < 0.

### The correct V is the analytic continuation of the existing formula

Substituting `sqrt(r) = i*sqrt(|r|)` into the existing expression makes `U`
purely imaginary, `U = i*V`, hence `U U^T = -V V^T`, which is exactly the
required diagonal-minus-low-rank structure. Writing it out as a real quantity:

    V_k = sqrt(vtot/2) / (2*b_k*sqrt(|r|))
          * ( (1-rg_eq) - sqrt((1-rg_eq)^2 - 4*b_k^2*|r|/vtot) )
          * sd_k

where `b_k` is the standardized haploid effect, `sd_k = sqrt(AF*(1-AF))`, and
`vtot = vg_eq + (1-h2_0)`.

This was verified numerically against complex arithmetic: the real part of the
continued expression is exactly zero, and its imaginary part equals `V` to
machine precision. The identity `-(sum_k a_k V_k)^2 = vg_eq - vg_0` reproduces
the theoretical value to four decimal places at m = 5000, the residual being
the expected `O(1/m)` diagonal term.

### Changes

**`am_covariance_structure(beta, AF, r)`**, signature unchanged:

- `r > 0`: unchanged.
- `r < 0`: use the `V` expression above.
- `r = 0`: return a vector of zeros, since panmixia induces no LD. This fixes an
  existing bug.
- Attach `attr(U, "sign") <- sign(r)`, using `1` when `r = 0`.
- Guard the discriminant `(1-rg_eq)^2 - 4*b_k^2*|r|/vtot >= 0` and raise an
  informative error naming `r`, `m`, and the largest offending effect size. In
  practice this bound is permissive: it holds down to `r = -0.999` even at
  m = 20, so it is a safety net rather than a routine constraint.

The return value stays a plain numeric vector. An attribute does not change how
the vector behaves in arithmetic, indexing, or printing, so existing user code
that calls `am_covariance_structure()` directly continues to work unchanged.
This matters because the package is on CRAN.

**`rb_dplr(n, mu, U, sign = NULL)`**, one new optional argument:

- `sign = NULL` resolves to `attr(U, "sign")`, falling back to `+1` when the
  attribute is absent. Old code and hand-built `U` vectors keep working.
- The sampler draws from `C = D + sign * U U^T`. The only change to the
  recursion is the two conditional-probability lines, which become
  `p <- mu[m] + sign * x * U[m]`.

The accumulator updates are deliberately left alone. Under `U = i*V` induction
gives `x = i * x_V` at every step, so `p = mu - x_V * V`, and the factor `c`
stays real throughout. Flipping the sign at the two probability lines is thus
exactly equivalent to running the recursion in complex arithmetic.

**`am_simulate()`**: no logic change. Once the two functions above accept
negative `r` it flows through. Documentation updated to state that `r` lies in
`(-1, 1)` and to describe the disassortative case.

**`rb_unstr()`**: no change. It accepts a full correlation matrix `C` directly,
so negative correlations are already expressible.

### Validation

An end-to-end prototype confirmed the approach at m = 2000, n = 5000. Empirical
heritability and genetic variance track the theoretical `h2_eq` and `vg_eq`
across `r` in `[-0.6, 0.6]`, and realized allele frequencies correlate with
target frequencies at 0.9998.

| r | empirical h2 | h2_eq | empirical vg | vg_eq |
|-------|--------|--------|--------|--------|
| -0.60 | 0.4342 | 0.4415 | 0.3914 | 0.3953 |
| -0.40 | 0.4554 | 0.4580 | 0.4337 | 0.4226 |
| -0.20 | 0.4670 | 0.4772 | 0.4507 | 0.4564 |
| +0.20 | 0.5269 | 0.5279 | 0.5665 | 0.5590 |
| +0.60 | 0.6101 | 0.6126 | 0.7830 | 0.7906 |

Note that `h2_eq()`, `rg_eq()`, and `vg_eq()` already return correct values for
negative `r` with no modification, and `h2_eq` equals `vg_eq/vtot` exactly.

## Part 2: binary int8 genotype output

### Motivation

`rb_dplr()` allocates `matrix(NaN, nrow = n, ncol = 2*m)`, a double matrix over
haplotype columns. At n = 10,000 and m = 100,000 that is roughly 16 TB, so the
package cannot currently reach the scale its own design targets. Writing int8
to disk costs one byte per genotype, and streaming removes the full-matrix
allocation entirely.

### API

`am_simulate()` gains three optional arguments:

    am_simulate(h2_0, r, m, n, afs = NULL, min_MAF = .1, haplotypes = FALSE,
                path = NULL,
                format = c("individual", "variant", "bed"),
                batch_size = NULL)

- `path = NULL` keeps today's behavior exactly, returning `X` in memory.
- `path` supplied streams genotypes to disk and returns the same list without
  `X`, adding `path`, `format`, `n`, and `m`.

Two exported helpers cover the in-memory export case and reading results back:

- `write_genotypes(X, path, format = "individual")`
- `read_genotypes(path)`

`read_genotypes()` reads the sidecar metadata and reconstructs an `n` by `m`
integer matrix regardless of on-disk layout, so the three formats are
interchangeable from R.

### On-disk formats

`"individual"` is the default.

- **individual**: flat int8, values 0/1/2, one individual's `m` loci
  contiguous. File size `n*m` bytes.
- **variant**: flat int8, values 0/1/2, one locus' `n` individuals contiguous.
  File size `n*m` bytes. This matches PLINK variant ordering.
- **bed**: PLINK 2-bit variant-major, magic bytes `0x6c 0x1b 0x01`, four
  samples per byte. Roughly four times smaller and directly readable by plink,
  GCTA, and bigsnpr. Dosage of the effect allele maps to PLINK codes as
  2 -> `00`, 1 -> `10`, 0 -> `11`, with the effect allele written as A1 in the
  `.bim`. This convention is documented in the help page because A1/A2 handling
  is a common source of error.

### Batching

The batching dimension follows the layout, because each layout wants sequential
writes and the two layouts disagree about what "sequential" means.

- **individual**: batch over individuals. Rows are independent draws, so a
  batch runs the full locus recursion on a subset of rows, and its output is
  written contiguously. Setting `batch_size = n` draws `runif(m*n)` in exactly
  the order the current code does, so the streaming path reproduces the
  in-memory path bit for bit under the same seed. This becomes a test.
- **variant** and **bed**: batch over loci, carrying the recursion state `x`
  and `c` across blocks. This is the current algorithm with the full matrix
  replaced by a block buffer, so it is a small change.

`batch_size = NULL` picks a batch targeting roughly a 128 MB working buffer,
clamped to at least 1 and at most the relevant dimension.

Because the two layouts consume random numbers in different orders, the
streaming `"variant"` path is not bit-identical to the `"individual"` path
under the same seed. This is documented rather than engineered around.

### Buffers

Streaming buffers use `integer` storage rather than the current
`matrix(NaN, ...)` double, halving working memory before the int8 conversion.

### Sidecar metadata

No new package dependencies. `Depends` stays at `R` and `stats`, using base
`writeBin`, `readBin`, `saveRDS`, and `writeLines`.

- `<prefix>.meta`: plain text `key: value` lines recording `n`, `m`, `format`,
  `dtype`, and byte order.
- `<prefix>.rds`: the simulation objects, meaning `AF`, `beta_std`, `beta_raw`,
  `y`, `g`, and the call parameters.
- `bed` additionally writes `<prefix>.bim` and `<prefix>.fam` as PLINK requires.

### Interaction with `haplotypes = TRUE`

When streaming, `haplotypes = TRUE` is rejected with a clear error. Returning
the full `2*m` haplotype matrix in memory would defeat the purpose of
streaming. Streaming haplotypes to a second file is deliberately out of scope
for this change.

## Part 3: tests

The package currently has no `tests/` directory. Add `tests/testthat` and
`testthat` under `Suggests`, covering:

- `am_covariance_structure()` for `r < 0` matches the imaginary part of the
  complex continuation.
- `r = 0` returns zeros rather than `NaN`.
- The `sign` attribute is set, and is honored by `rb_dplr()`.
- Empirical `h2` from `am_simulate()` tracks `h2_eq()` within Monte Carlo
  tolerance for negative, zero, and positive `r`.
- Realized allele frequencies track target frequencies for negative `r`.
- Round trip `write_genotypes()` then `read_genotypes()` for all three formats.
- Streaming with `batch_size = n` reproduces the in-memory result exactly.
- `.bed` bit packing is correct, checked against a hand-built reference block.
- The discriminant guard raises an informative error rather than producing
  `NaN`.

## Documentation

- Roxygen updates for the changed functions, plus new pages for
  `write_genotypes()` and `read_genotypes()`; regenerate `NAMESPACE` and `man/`.
- README gains a short negative-assortment example and an on-disk output
  example.
- `NEWS.md` gains a 1.1.0 section. Bump `DESCRIPTION` to 1.1.0.
- `DESCRIPTION` `Description:` field updated, since it currently says only
  "assortative mating".

## Out of scope

- Streaming haplotypes.
- Multivariate or bivariate assortative mating.
- Compression of the int8 output.
- Cleaning the 17 stray lines in `.gitignore` that were pasted from
  `rBahadur.Rproj` and match nothing. Noted separately, not part of this work.

## Open question

The names `write_genotypes()` and `read_genotypes()` do not follow either
existing prefix convention in the package, which uses `rb_` for Bahadur
samplers and `am_` for assortative mating helpers. They are genuinely neither.
The alternative is `am_write_genotypes()` and `am_read_genotypes()` for
namespace consistency at the cost of a slightly misleading prefix.

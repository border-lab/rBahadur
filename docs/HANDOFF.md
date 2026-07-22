# rBahadur handoff

Written 2026-07-22. Everything below is committed on `main`. The package went
from 1.0.0, which had no automated tests at all, to 1.1.0 with three new
feature areas and 67 test blocks.

## Current state

| | |
|---|---|
| branch | `main`, clean tree, no untracked files |
| version | 1.1.0 |
| `R CMD check --as-cran --run-donttest` | **Status: OK**, zero errors, warnings, and notes |
| tests | 67 blocks, 0 failed, 0 skipped |
| exported functions | 15 |
| tarball | 122,229 bytes |

Feature branches `negative-am-int8`, `cli`, and `local-ld` are all merged into
`main` and can be deleted whenever you like.

## How to verify

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'

R CMD build .
_R_CHECK_SYSTEM_CLOCK_=0 R CMD check --as-cran --run-donttest rBahadur_1.1.0.tar.gz
```

Two things about that recipe are load-bearing.

**Run the tests with `NOT_CRAN=true`.** Several carry `skip_on_cran()` and are
silently skipped otherwise. Early in this work a test reported as passing was
in fact being skipped and genuinely failing.

**Check with `--as-cran` on a built tarball, not `devtools::check(cran = TRUE)`.**
The latter reports clean while missing real problems. It missed all three of
these, which are now fixed:

- the bundled `.rds` used serialization version 3, which silently forced
  `R >= 3.5.0`, and `R CMD build` was rewriting the declared `R (>= 3.3.0)` in
  the tarball, so source and built package disagreed about the floor
- `URL:` pointed at `github.com/rborder/rBahadur`, which redirects to
  `border-lab/`, and moved URLs are flagged
- the `rb_unstr` example took 10.4s, over CRAN's 5s budget, because it builds a
  dense 400x400 matrix. Pre-existing since 1.0.0. Now 0.6s.

`_R_CHECK_SYSTEM_CLOCK_=0` suppresses an intermittent `future file timestamps`
NOTE, which is the check failing to reach the world clock API rather than
anything about the package.

It is also worth installing the tarball and exercising it, rather than relying
on `load_all()`, which hides installation-time questions such as whether
`exec/rbahadur` actually lands executable. It does: `-rwxrwxr-x`.

## What shipped in 1.1.0

### 1. Disassortative mating

`am_covariance_structure()` previously returned all `NaN` for negative
cross-mate correlation, and divided by zero at `r = 0`.

The fix is structural rather than cosmetic. At equilibrium
`(sum_k a_k U_k)^2 = vg_eq * rg_eq`, and for `r < 0` the right side is
negative, so no real loading vector exists. A rank-one term can only *increase*
genetic variance while disassortment *reduces* it, so the covariance must
become `D - VV^T`, with `V` the analytic continuation of the positive branch.
The sign travels on `attr(U, "sign")`, which `rb_dplr()` honors through a new
`sign` argument.

The formula was later collapsed to one expression using
`sign(r) * sqrt(abs(r))`, verified bit-identical to the branched version across
44 cases. The `r == 0` early return has to stay, because `sign(0)` is 0.

**Feasibility envelope, measured.** Negative `r` leaves the Bahadur order-2
feasible region much sooner than positive `r`, and worsens as `n` grows,
because infeasibility is a tail event across individuals. Positive `r` was
feasible in every configuration tried. At `h2_0 = 0.5`:

| configuration | result |
|---|---|
| `r = -0.6`, m = 800, n = 1500 | 92% of seeds feasible |
| `r = -0.4`, m = 1500, n = 2000 | 6 of 6 seeds |
| `r = -0.4`, m = 1500, n = 4000 | fails for some seeds |
| `r = -0.3`, m = 1500, n = 4000 | 6 of 6 seeds |

Raising `min_MAF` widens it. Documented under "Feasibility under negative
assortment" in `?am_covariance_structure`. `rb_dplr()` aborts the whole call
rather than discarding the offending individual, which is deliberate:
discard-and-resample would silently bias the sample toward the interior of the
feasible region.

**Independently validated.** A from-scratch forward-time simulator (25
generations from panmixia, m = 600, n = 4000, sharing no code with the package)
converges to `vg_eq` from below for negative `r` and from above for positive,
within 1 to 4 percent in every case.

### 2. Binary genotype output

`write_genotypes()` and `read_genotypes()` cover individual-major int8,
variant-major int8, and PLINK `.bed`. `am_simulate()` gained `path`, `format`,
and `batch_size` for streaming.

- `path = NULL` is bit-identical to the 1.0.0 release across all six returned
  components, verified against a git worktree at the pre-change commit.
- Variant-major and bed reproduce the in-memory genotypes at *any* batch size,
  because R fills matrices column-major and drawing `runif(n)` per locus
  consumes the stream in the same order. Individual-major matches only at
  `batch_size >= n`.
- `.bed` output validated against real PLINK v1.9.0-b.7.7: a known matrix round
  trips exactly including a padded final byte, allele frequencies match, and
  missing calls report correctly under `--missing`.
- The README documents reading both int8 layouts into `numpy`, including
  `mmap`. Note the documented gotcha: `pandas_plink` counts the opposite allele
  by default and returns `2 - X` unless you pass `ref="a0"`.

### 3. Command line interface

`exec/rbahadur` with `simulate` and `info`. Logic lives in `R/cli.R` so it is
unit-testable; the shipped script is a shim converting the return value into an
exit code. Exit 0 for success, 1 for a usage error, 2 for a runtime failure, so
a caller can distinguish a typo from an infeasible parameter set.

```bash
ln -s $(Rscript -e "cat(rBahadur::rbahadur_cli_path())") ~/bin/rbahadur
```

### 4. Local LD via `am_mosaic()`

Implements Algorithm S4 of the supplementary note. Causal variants come from
`rb_dplr()`, giving the dense genome-wide assortative mating structure;
intervening markers are filled by copying contiguous blocks from a reference
panel, giving local LD for free because the blocks are real haplotypes. Each
block's donor must carry the allele already drawn at that block's causal locus,
so the AM structure is preserved **exactly**, not approximately.

**What was improved over the published vignette.** The supplement specifies
recombination-weighted breakpoints, but the vignette says outright: *"For
simplicity, we draw haplotype block boundaries uniformly at random between
causal variants, though a genetic map could also be used."* That is now
implemented, drawing each boundary from a categorical weighted by the cM
distance of each interval.

Bundled data: `kg_reference()` returns a 52 KB real 1000 Genomes panel
(chr22:20-21 Mb, GRCh38, 520 haplotypes across 2500 common SNVs, with
interpolated map positions), so examples and tests exercise real LD offline.
`vcf_to_panel()` and `download_1kg_panel()` build panels at realistic scale.

**The tradeoff, documented in the vignette.** Block length falls as causal
variants increase, but `vg_eq()` is an infinitesimal-limit result and needs
many causal variants. Measured shortfall in `vg`: 17.7% at m = 10, 10.3% at
m = 25, and within replicate noise from m = 50 up. In a 1 Mb window these
pressures conflict; at genome scale they do not, which is the regime the method
is for.

## Outstanding

### Release chores

- `cran-comments.md` still reads "This is a new release. 0 errors | 0 warnings
  | 1 note" and needs rewriting for 1.1.0. It is now accurate to claim zero
  notes under `--as-cran`.
- If you ever regenerate `inst/extdata/kg_chr22_panel.rds`, save it with
  `version = 2` or the R dependency floor silently jumps to 3.5.0 again.

### Deferred, none can produce wrong output

Triaged as follow-up during review rather than merge blockers:

- The `.bed` magic bytes and per-variant pack loop are duplicated between
  `R/genotype_io.R` and `R/am_simulate.R` rather than shared.
- A mid-stream failure leaves a partial data file with no `.meta`, so it is
  unreadable rather than dangerous, but nothing unlinks it. At real scale that
  is tens of GB of litter.
- The auto `batch_size` heuristic documents "roughly 128 MB" but real peak is
  about 3x that. The documentation was corrected rather than the heuristic.
- `format` and `batch_size` are silently ignored when `path = NULL`.
- No test covers `am_covariance_structure()`'s infeasibility guard or its `r`
  type/length/NA validation branch.
- `am_simulate(m = 1)` fails with an opaque "missing value where TRUE/FALSE
  needed". The cause is pre-existing and unrelated: `scale(rnorm(1))` is `NaN`
  because the standard deviation of one value is undefined. Relatedly
  `rb_dplr`'s `for (m in 2:(M-1))` iterates `c(2, 1)` when `M == 2`; it returns
  output rather than crashing, and only bites at exactly one causal variant.

### Possible correction to the published supplement

Three defects in the code in Additional file 1, beyond the uniform-breakpoint
simplification already noted:

- maternal and paternal block endpoints are constructed asymmetrically,
  `rep(0, n)` versus `rep(0, 2)`
- the final gather indexes with `as.vector(ii_maternal)[1:n]` against an index
  vector of length `M*n`, so it recycles
- it builds full `M x n` matrices, roughly 14 GB at the 173,545 by 10,000 scale
  the vignette itself describes

Worth deciding whether that warrants an erratum, given the vignette is the
documented path for combining assortative mating with local LD.

## Environment notes

- **There is a stale rBahadur 1.0.0 in your personal library**
  (`~/R/x86_64-pc-linux-gnu-library/4.5`). A bare `library(rBahadur)` outside
  this project loads that, not the new build, and it exports 7 functions rather
  than 15. Reinstall over it when you are ready, or you will get confusing
  "could not find function am_mosaic" errors.
- `bspm` on this machine redirects `install.packages()` to apt, which silently
  fails for anything not in the apt repo. Call `bspm::disable()` first.
- PLINK 1.9 was used for the `.bed` interop check but lived in session scratch
  and will not persist. Re-fetch from
  `https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip`.
- The design spec and implementation plan are in `docs/superpowers/`.
  Per-task execution notes are in `.superpowers/sdd/`, gitignored and local.
- `docs/` is in `.Rbuildignore`, so nothing here ships in the tarball.

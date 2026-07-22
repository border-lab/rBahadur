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

### 5. Expansion performance

`am_mosaic()` was rewritten around two expansion kernels, one per output
orientation, and the reference panel now stays `raw`. Output is unchanged:
bit-identical to the previous implementation under the same seed across all
formats and batch sizes, with the random number stream left in the same state.
That was verified over 45 cases spanning both signs of `r`, explicit and
adjacent causal indices, raw and integer panels, panels with a hotspot map and
with no map at all, every layout at four batch sizes, `am_simulate()`,
`write_genotypes()`, and `rb_dplr()` down to the degenerate `M = 2`.

Measured on a 2,000 haplotype by 100,000 marker panel, `n = 1,000`, 400 causal
variants, timing the call alone against a prebuilt panel:

| layout | before | after | peak RSS |
|---|---|---|---|
| in memory | 19.3s | 5.7s | 4.0 GB to 1.3 GB |
| individual-major | 18.7s | 7.3s | 4.0 GB to 1.3 GB |
| variant-major | 31.0s | 5.9s | 3.6 GB to 1.4 GB |
| bed | 36.3s | 6.0s | 3.6 GB to 1.5 GB |

What changed, in rough order of how much it was worth:

- **The panel is no longer expanded to integer.** `.mosaic_check_panel()` used
  to convert the whole panel on entry, which is four bytes per allele where one
  will do, and at genome scale costs more time than the simulation. It stays
  `raw` and is converted a block of markers at a time where it is used.
- **Marker-ordered sweep.** The old variant-major gather rescanned all `n`
  individuals at every marker, in a `repeat` loop, to find who had crossed a
  block boundary, building two-column index matrices to do it. Boundaries are
  now compiled once into a schedule listing who changes at each marker, so a
  marker costs one gather per parental copy. Both parental copies share one
  pass, so the panel block is converted once rather than twice.
- **Transposed panel for individual-major.** A haplotype is a contiguous read
  rather than a stride over the panel. This is the one path that holds a second
  copy of the panel, and it is still half the old footprint, since the old copy
  was integer.
- **`.gt_pack_bed()` takes a block.** It recodes by lookup instead of four
  masked assignments, and packs a block of variants per call rather than one
  call per variant. `write_genotypes()` and `am_simulate()` benefit too.
- **`rb_dplr()` draws uniforms per locus** instead of allocating an `n` by `M`
  matrix of them, the same column-major equivalence `.rb_dplr_stream()` already
  documented and is tested on. `M == 2` keeps the up-front draw, because the
  recursion counts down there and reads locus 1 twice; that is the degenerate
  case noted under Deferred below, and its output is deliberately unchanged.

Two things are worth knowing before optimising further. The remaining time is
roughly half in the sweep itself, which is near the floor for vectorised R at
one gather per individual per marker, and the schedule build is another 15 to
20 percent. And the naive `n` by `p` donor index matrix is exactly what to
avoid: at the 173,545 by 10,000 scale it is 7 GB, which is the same trap the
published supplement falls into.

## Outstanding

### Release chores

- `cran-comments.md` has been rewritten for 1.1.0: a minor release rather than
  a new one, 0/0/0, the three `--as-cran` fixes listed, and no reverse
  dependencies, which was confirmed against the CRAN package database rather
  than assumed.
- **Still to do before submitting: check on a platform other than this one.**
  `cran-comments.md` lists only local Ubuntu, because that is all that has
  actually been run. Add win-builder and macOS results, or R-hub, before
  submission.
- Locally the check reports one NOTE, "Skipping checking HTML validation: no
  command 'tidy' found". That is this machine lacking `tidy`, not the package.
  With `_R_CHECK_RD_VALIDATE_RD2HTML_=FALSE` the status is OK with no notes at
  all, which is what `cran-comments.md` claims.
- If you ever regenerate `inst/extdata/kg_chr22_panel.rds`, save it with
  `version = 2` or the R dependency floor silently jumps to 3.5.0 again.

### Deferred, none can produce wrong output

Triaged as follow-up during review rather than merge blockers:

- The `.bed` magic bytes are still written from three places rather than
  shared. The per-variant pack loops that sat alongside them are gone:
  `.gt_pack_bed()` now takes a whole block.
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

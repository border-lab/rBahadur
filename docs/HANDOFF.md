# rBahadur handoff

Written 2026-07-22. Everything below is on `main`, committed, with
`R CMD check --cran` reporting 0 errors, 0 warnings, 0 notes at version 1.1.0,
and 67 testthat blocks passing with zero skips.

## Verify the state

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'

R CMD build .
_R_CHECK_SYSTEM_CLOCK_=0 R CMD check --as-cran rBahadur_1.1.0.tar.gz
```

**Use `--as-cran` on a built tarball, not `devtools::check(cran = TRUE)`.** The
latter reports clean while missing real problems. It missed all three of these:

- the bundled `.rds` used serialization version 3, which silently forced
  `R >= 3.5.0`, and `R CMD build` was rewriting the declared `R (>= 3.3.0)` in
  the tarball so source and built package disagreed about the floor
- `URL:` pointed at `github.com/rborder/rBahadur`, which now redirects to
  `border-lab/`, and moved URLs are flagged
- the `rb_unstr` example took 10.4s, over CRAN's 5s example budget

All three are fixed. `--as-cran` on the tarball is now `Status: OK` with zero
notes. Also worth building and installing the tarball rather than relying on
`load_all()`, which hides installation-time problems such as whether
`exec/rbahadur` actually lands executable.

`_R_CHECK_SYSTEM_CLOCK_=0` matters: without it the check intermittently emits a
`future file timestamps` NOTE, which is `R CMD check` failing to reach the world
clock API rather than anything about the package.

Tests must be run with `NOT_CRAN=true`. Several carry `skip_on_cran()` and are
silently skipped otherwise. Early in this work a test reported as passing was in
fact being skipped and genuinely failing.

## What shipped in 1.1.0

The package went from 1.0.0 (no tests at all) to three feature areas.

### 1. Disassortative mating

`am_covariance_structure()` previously returned all `NaN` for negative
cross-mate correlation and divided by zero at `r = 0`.

The fix is structural, not cosmetic. At equilibrium
`(sum_k a_k U_k)^2 = vg_eq * rg_eq`, and for `r < 0` the right side is negative,
so no real loading vector exists. A rank-one term can only *increase* genetic
variance while disassortment *reduces* it, so the covariance must become
`D - VV^T`. `V` is the analytic continuation of the positive branch. The sign
travels on `attr(U, "sign")`, which `rb_dplr()` honors through a `sign`
argument.

The formula was later collapsed to a single expression using
`sign(r) * sqrt(abs(r))`, verified bit-identical to the branched version across
44 cases. The `r == 0` early return has to stay, because `sign(0)` is 0.

**Feasibility envelope, measured.** Negative `r` leaves the Bahadur order-2
feasible region much sooner than positive `r`, and worsens as `n` grows because
infeasibility is a tail event across individuals. Positive `r` was feasible in
every configuration tried. At `h2_0 = 0.5`:

| configuration | result |
|---|---|
| `r = -0.6`, m = 800, n = 1500 | 92% of seeds feasible |
| `r = -0.4`, m = 1500, n = 2000 | 6 of 6 seeds |
| `r = -0.4`, m = 1500, n = 4000 | fails for some seeds |
| `r = -0.3`, m = 1500, n = 4000 | 6 of 6 seeds |

Raising `min_MAF` widens it. This is documented in
`?am_covariance_structure` under "Feasibility under negative assortment".
`rb_dplr()` aborts the whole call rather than discarding the offending
individual, which is deliberate: discard-and-resample would silently bias the
sample toward the interior of the region.

**Independently validated.** A from-scratch forward-time simulator (25
generations from panmixia, m = 600, n = 4000, sharing no code with the package)
converges to `vg_eq` from below for negative `r` and from above for positive,
within 1 to 4 percent in every case.

### 2. Binary genotype output

`write_genotypes()` and `read_genotypes()` cover individual-major int8,
variant-major int8, and PLINK `.bed`. `am_simulate()` gained `path`, `format`,
and `batch_size` to stream to disk.

- `path = NULL` is bit-identical to the 1.0.0 release, verified against a git
  worktree at the pre-change commit across all six returned components.
- Variant-major and bed reproduce the in-memory genotypes at *any* batch size,
  because R fills matrices column-major and drawing `runif(n)` per locus
  consumes the stream in the same order. Individual-major matches only at
  `batch_size >= n`.
- `.bed` output was validated against real PLINK v1.9.0-b.7.7: a known matrix
  round-trips exactly including a padded final byte, plink's A1 frequencies
  match, and missing calls report correctly under `--missing`.

### 3. Command line interface

`exec/rbahadur` with `simulate` and `info`. Logic lives in `R/cli.R` so it is
unit-testable; the script is a shim converting the return value to an exit code.
Exit 0 success, 1 usage error, 2 runtime failure, so a caller can tell a typo
from an infeasible parameter set. Put it on PATH with:

```bash
ln -s $(Rscript -e "cat(rBahadur::rbahadur_cli_path())") ~/bin/rbahadur
```

### 4. Local LD via `am_mosaic()`

Implements Algorithm S4 of the supplementary note. Causal variants come from
`rb_dplr()` (dense genome-wide AM structure); intervening markers are filled by
copying contiguous blocks from a reference panel (local LD, free, because the
blocks are real haplotypes). Each block's donor must carry the allele already
drawn at that block's causal locus, so the AM structure is preserved *exactly*,
not approximately.

**What was improved over the published vignette.** The supplement specifies
recombination-weighted breakpoints, but the vignette says outright: *"For
simplicity, we draw haplotype block boundaries uniformly at random between
causal variants, though a genetic map could also be used."* That is now
implemented, drawing each boundary from a categorical weighted by the cM
distance of each interval. Also: the published code has asymmetric maternal and
paternal endpoint construction (`rep(0, n)` versus `rep(0, 2)`), an index-length
mismatch in the final gather (`as.vector(ii_maternal)[1:n]` against an index
vector of length `M*n`), and builds full `M x n` matrices, roughly 14 GB at the
scale its own example describes.

**The real tradeoff, documented in the vignette.** Block length falls as causal
variants increase, but `vg_eq()` is an infinitesimal-limit result and needs many
causal variants. Measured shortfall in `vg`: 17.7% at m = 10, 10.3% at m = 25,
within replicate noise from m = 50 up. In a 1 Mb window these conflict; at
genome scale they do not, which is the regime the method is for.

Bundled data: `kg_reference()` returns a 52 KB real 1000 Genomes panel
(chr22:20-21 Mb, GRCh38, 520 haplotypes x 2500 common SNVs, with interpolated
map positions) so examples and tests exercise real LD offline.
`vcf_to_panel()` and `download_1kg_panel()` build panels at realistic scale.

## Outstanding, in priority order

### Release chores

- `cran-comments.md` still reads "This is a new release. 0 errors | 0 warnings
  | 1 note" and needs rewriting for a 1.1.0 submission. It is now accurate to
  say 0 errors, 0 warnings, 0 notes under `--as-cran`.
- Decide whether this ships as 1.1.0 or gets split. Everything currently sits
  in one 1.1.0 section of `NEWS.md`.
- If you ever regenerate `inst/extdata/kg_chr22_panel.rds`, save it with
  `version = 2` or the R dependency floor silently jumps to 3.5.0 again.

### Deferred, none can produce wrong output

These were triaged as follow-up during review, not merge blockers:

- The `.bed` magic bytes and the per-variant pack loop are duplicated between
  `R/genotype_io.R` and `R/am_simulate.R` rather than shared.
- A mid-stream failure leaves a partial data file with no `.meta`, so it is
  unreadable rather than dangerous, but nothing unlinks it. At real scale that
  is tens of GB of litter.
- The auto `batch_size` heuristic documents "roughly 128 MB" but real peak is
  about 3x that; the doc was corrected rather than the heuristic.
- `format` and `batch_size` are silently ignored when `path = NULL`.
- No test covers `am_covariance_structure()`'s infeasibility guard or its `r`
  type/length/NA validation branch.
- `am_simulate(m = 1)` fails with an opaque "missing value where TRUE/FALSE
  needed". Root cause is pre-existing and unrelated to this work:
  `scale(rnorm(1))` is `NaN` because the standard deviation of one value is
  undefined. Relatedly `rb_dplr`'s `for (m in 2:(M-1))` iterates `c(2, 1)` when
  `M == 2`; it returns output rather than crashing and only bites at exactly one
  causal variant.

## Two findings that belong to other projects

### xftsim silently ignores negative assortative mating

**This is a bug in your other package and it is silent.** Requested `r = +0.60`
realizes 0.585, 0.602, 0.607 across generations. Requested `r = -0.60` realizes
-0.000, -0.023, +0.026, which is random mating, with no error or warning.

Cause, in `xftsim/mate.py`, `LinearAssortativeMatingRegime.mate()`:

```python
R = np.sum(cross_cov) / np.sqrt(np.sum(within_cov1)*(np.sum(within_cov2)))
mating_score1 = sum_scaled_mate1 * np.sqrt(R) + np.sqrt(np.abs(1-R))*np.random.normal(...)
```

When `r < 0`, `R < 0`, so `np.sqrt(R)` is `NaN`, the whole mating score vector
becomes `NaN`, and `np.argsort` on an all-`NaN` array returns the identity
permutation, so mates pair arbitrarily. The `assert -1 <= r <= 1` in the
constructor accepts the value, so nothing complains.

It is the same class of bug rBahadur had, and wants the same fix:
`np.sign(R)*np.sqrt(np.abs(R))`, applied to one side's score only so high pairs
with low.

This is why the negative-r validation was done with a from-scratch simulator
instead of xftsim.

**Status: fixed and merged.** Sasha Gusev had already reported it as
border-lab/xftsim#18 with essentially this fix. It is now on both `main`
(PR #21) and `dev` (PR #23), both at version 0.3.1, each verified from a fresh
clone to realize both signs correctly. `dev` was cherry-picked rather than
merged from `main`, so that `main`'s deletion of `claude.md` and
`devtools/claude.md` did not ride along with a bug fix.

border-lab/xftsim#22 remains open for the `v0.9alpha` branch. That branch needs
**no code change**: it already takes `abs(self.r)` before the square root and
negates one sex, which is the same fix by another route. What it needs is the
regression tests ported, since it has nothing asserting a negative `r` is
actually realized, and care when the branches converge (156 commits diverged)
so the fix is not clobbered.

### The published supplement's vignette

The three code defects listed under `am_mosaic()` above are in the published
Additional file 1. Worth deciding whether that warrants a correction, given the
vignette is the documented path for combining AM with local LD.

## Environment notes

- A conda env `xftsim-rbahadur` was created at `~/miniforge3/envs` (Python 3.9,
  needed because xftsim pins `numba==0.56.4`). Remove with
  `conda env remove -n xftsim-rbahadur` if you do not want it.
- PLINK 1.9 was downloaded to a session scratch directory and will not survive.
  Re-fetch from `https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip`
  if you want to redo the `.bed` interop check.
- `bspm` redirects `install.packages()` to apt on this machine, which silently
  fails for packages not in the apt repo. Use `bspm::disable()` first.
- Design spec and implementation plan for the 1.1.0 work are in
  `docs/superpowers/`. Per-task execution notes are in `.superpowers/sdd/`,
  which is gitignored and local only.

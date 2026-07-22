## Submission

This is a minor release, updating rBahadur from 1.0.0 to 1.1.0. It adds
disassortative mating, binary genotype output as int8 or PLINK bed, a command
line interface, and simulation of realistic local linkage disequilibrium by
copying haplotype blocks from a reference panel. NEWS.md has the full list.

Three issues raised by `--as-cran` since 1.0.0 have been fixed:

* the bundled `.rds` in `inst/extdata` is now written with `version = 2`. It
  previously used serialization version 3, which silently raised the effective
  R dependency to 3.5.0 while `DESCRIPTION` continued to declare 3.3.0, so the
  source package and the built tarball disagreed about the floor.
* `URL:` now points at <https://github.com/border-lab/rBahadur>. The address
  given in 1.0.0 redirected there.
* the `rb_unstr` example took 10.4 seconds, over the 5 second budget, because
  it built a dense 400 by 400 matrix. The slowest example is now 0.6 seconds.

## R CMD check results

0 errors | 0 warnings | 0 notes

Checked with `R CMD check --as-cran --run-donttest` against a built tarball
rather than the source directory, and with `NOT_CRAN=true` so that no test is
skipped.

## Test environments

* Ubuntu 24.04, R 4.5.2, x86_64-pc-linux-gnu (local)

## Reverse dependencies

There are none on CRAN.

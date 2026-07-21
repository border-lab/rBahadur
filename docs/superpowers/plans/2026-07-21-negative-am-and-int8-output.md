# Negative Assortative Mating and Binary int8 Genotype Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `rBahadur` simulate under negative (disassortative) mating, and write simulated genotypes to disk as binary int8 or PLINK `.bed` with an optional batched streaming path that removes the full-matrix allocation.

**Architecture:** At assortative mating equilibrium the haploid LD matrix is `C = D + s*U U^T` where `s = sign(r)`. Negative `r` needs `s = -1`, because a rank-one term can only increase genetic variance while disassortment reduces it. The correct loading vector is the analytic continuation of the existing formula, and the sampler change is a sign multiplier on two conditional-probability lines. Separately, genotype output gains three on-disk layouts plus a streaming driver whose batching dimension follows the layout.

**Tech Stack:** Base R only (`stats`, plus `writeBin`, `readBin`, `readLines`, `writeLines`, `saveRDS`). Tests use `testthat` 3rd edition under `Suggests`. Docs via roxygen2 7.2.3.

## Global Constraints

- No new hard dependencies. `Depends` stays `R (>= 3.3.0), stats`. `testthat` goes under `Suggests` only.
- The package is on CRAN. `am_covariance_structure()` must keep returning a plain numeric vector, and `rb_dplr()` must keep working when called with three positional arguments.
- `am_simulate()` with `path = NULL` must remain bit-for-bit identical to the current release under the same seed. Do not reorder any `rnorm()` or `runif()` call in that path.
- Genotype values are `0`, `1`, `2`. Only the `bed` format can represent missing.
- `r` lies in the open interval `(-1, 1)`.
- Style: no em dashes and no non-ASCII characters in R source, roxygen, or Markdown.
- Target version is 1.1.0.
- Negative `r` has a tighter Bahadur feasibility envelope than positive `r`,
  measured during execution. Positive `r` was feasible in every configuration
  tried. Negative `r` degrades as `n` grows, because infeasibility is a tail
  event across individuals: it only takes one individual to leave the region.
  Measured points, all at `h2_0 = 0.5`: `r = -0.6` at m = 800, n = 1500 is 92
  percent of seeds; `r = -0.4` at m = 1500, n = 2000 is 6 of 6 seeds but at
  n = 4000 it fails for some seeds; `r = -0.3` at m = 1500, n = 4000 is 6 of 6
  seeds. Raising `min_MAF` widens the envelope. Tests must stay inside the
  reliable range and use fixed seeds. Verify feasibility across every seed a
  test actually uses, not just one.

## File Structure

| File | Status | Responsibility |
|---|---|---|
| `R/am_covariance_structure.R` | modify | Loading vector `U` for positive, zero, and negative `r`; sets the `sign` attribute |
| `R/rb_dplr.R` | modify | DPLR sampler; gains `sign` argument and the `.rb_sign()` resolver |
| `R/am_simulate.R` | modify | Simulation driver; gains `path`, `format`, `batch_size` |
| `R/genotype_io.R` | create | On-disk encoding: int8 layouts, bed packing, sidecar metadata, `write_genotypes()`, `read_genotypes()` |
| `R/am_stream.R` | create | `.rb_dplr_stream()`, the locus-blocked recursion that emits column blocks |
| `tests/testthat.R` | create | testthat entry point |
| `tests/testthat/test-*.R` | create | One file per unit under test |
| `DESCRIPTION`, `NEWS.md`, `README.md`, `NAMESPACE`, `man/` | modify | Metadata and docs |

Tasks 1 to 3 deliver negative assortment and are independently shippable. Tasks 4 to 7 deliver genotype output. Task 8 is the release wrap-up.

---

### Task 1: Negative and zero `r` in `am_covariance_structure()`

**Files:**
- Modify: `R/am_covariance_structure.R`
- Modify: `DESCRIPTION` (add `Suggests: testthat (>= 3.0.0)` and `Config/testthat/edition: 3`)
- Create: `tests/testthat.R`
- Test: `tests/testthat/test-am_covariance_structure.R`

**Interfaces:**
- Consumes: `rg_eq(r, h2_0)`, `vg_eq(r, vg_0, h2_0)` from `R/am_equilibrium_parameters.R`.
- Produces: `am_covariance_structure(beta, AF, r)` returning a numeric vector of length `2*length(beta)` carrying `attr(x, "sign")` equal to `1` or `-1`.

- [ ] **Step 1: Create the testthat entry point**

Create `tests/testthat.R`:

```r
library(testthat)
library(rBahadur)

test_check("rBahadur")
```

- [ ] **Step 2: Declare testthat in DESCRIPTION**

Append these two lines to `DESCRIPTION`, after `RoxygenNote: 7.2.3`:

```
Suggests:
    testthat (>= 3.0.0)
Config/testthat/edition: 3
```

- [ ] **Step 3: Write the failing test**

Create `tests/testthat/test-am_covariance_structure.R`:

```r
make_arch <- function(m = 400, h2_0 = 0.5, seed = 1) {
  set.seed(seed)
  list(beta = as.vector(scale(rnorm(m))) * sqrt(h2_0 / m),
       AF = runif(m, 0.1, 0.9),
       h2_0 = h2_0)
}

test_that("positive r is unchanged and tagged with sign 1", {
  a <- make_arch()
  U <- am_covariance_structure(a$beta, a$AF, 0.5)
  expect_true(all(is.finite(U)))
  expect_identical(attr(U, "sign"), 1)
  expect_length(U, 2 * length(a$beta))
})

test_that("r = 0 returns zeros rather than NaN", {
  a <- make_arch()
  U <- am_covariance_structure(a$beta, a$AF, 0)
  expect_true(all(U == 0))
  expect_identical(attr(U, "sign"), 1)
})

test_that("negative r matches the imaginary part of the complex continuation", {
  a <- make_arch()
  for (r in c(-0.8, -0.5, -0.2)) {
    b <- rep(a$beta, each = 2)
    sdh <- rep(sqrt(a$AF * (1 - a$AF)), each = 2)
    h20 <- sum(a$beta^2)
    rg <- rg_eq(r, h20)
    vtot <- vg_eq(r, h20, h20) + (1 - h20)
    rc <- as.complex(r)
    Uc <- sqrt(vtot / 2) / (2 * b * sqrt(rc)) *
      (sqrt(4 * b^2 * rc / vtot + (1 - rg)^2) - (1 - rg)) * sdh

    U <- am_covariance_structure(a$beta, a$AF, r)
    expect_equal(max(abs(Re(Uc))), 0)
    expect_equal(as.vector(U), Im(Uc), tolerance = 1e-12)
    expect_identical(attr(U, "sign"), -1)
  }
})

test_that("negative r reproduces the equilibrium variance deficit", {
  a <- make_arch(m = 4000)
  r <- -0.5
  h20 <- sum(a$beta^2)
  a_hap <- rep(a$beta / sqrt(2 * a$AF * (1 - a$AF)), each = 2)
  V <- am_covariance_structure(a$beta, a$AF, r)
  expect_equal(-sum(a_hap * V)^2, vg_eq(r, h20, h20) - h20, tolerance = 1e-3)
})

test_that("r outside (-1, 1) is rejected", {
  a <- make_arch()
  expect_error(am_covariance_structure(a$beta, a$AF, 1), "open interval")
  expect_error(am_covariance_structure(a$beta, a$AF, -1), "open interval")
})
```

- [ ] **Step 4: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-am_covariance_structure.R")'`

Expected: the `r = 0`, negative `r`, and rejection tests FAIL. The `r = 0` and negative cases fail because `U` is all `NaN`. The rejection test fails because no validation exists yet.

- [ ] **Step 5: Rewrite the function body**

Replace the body of `am_covariance_structure()` in `R/am_covariance_structure.R` (keep the existing roxygen block above it, which Task 8 updates):

```r
am_covariance_structure <- function(beta, AF, r) {
  if (!is.numeric(r) || length(r) != 1L || is.na(r)) {
    stop("`r` must be a single non-missing numeric value")
  }
  if (r <= -1 || r >= 1) {
    stop("`r` must lie in the open interval (-1, 1)")
  }
  ## obtain haploid substitution effects, variances
  beta_hap <- rep(beta, each = 2)
  sd_hap <- rep(sqrt(AF * (1 - AF)), each = 2)
  h2_0 <- sum(beta**2)

  ## panmixia induces no linkage disequilibrium
  if (r == 0) {
    U <- rep(0, length(beta_hap))
    attr(U, "sign") <- 1
    return(U)
  }

  ## compute equilibrium variance components
  rgeq <- rg_eq(r = r, h2_0)
  vgeq <- vg_eq(r = r, h2_0, h2_0)
  vtot <- vgeq + (1 - h2_0)
  abs_r <- abs(r)
  scale_term <- sqrt(vtot / 2) / (2 * beta_hap * sqrt(abs_r))

  if (r > 0) {
    ## C = D + U U^T
    U <- scale_term *
      (sqrt(4 * beta_hap**2 * abs_r / vtot + (1 - rgeq)^2) - (1 - rgeq)) * sd_hap
    attr(U, "sign") <- 1
  } else {
    ## C = D - U U^T; this branch is the analytic continuation of the branch
    ## above, which is purely imaginary for r < 0
    disc <- (1 - rgeq)^2 - 4 * beta_hap**2 * abs_r / vtot
    if (any(disc < 0)) {
      stop(sprintf(
        paste0("Infeasible negative-assortment structure at r = %g with %d causal ",
               "variants: the largest squared standardized effect (%g) violates ",
               "4*beta^2*|r|/vtot <= (1 - rg_eq)^2. Reduce |r| or increase `m`."),
        r, length(beta), max(beta_hap[disc < 0]**2)))
    }
    U <- scale_term * ((1 - rgeq) - sqrt(disc)) * sd_hap
    attr(U, "sign") <- -1
  }
  return(U)
}
```

Note that `sd_hap` is a plain numeric vector, so `U` stays a plain numeric vector and the attribute rides along without changing arithmetic or indexing behavior.

- [ ] **Step 6: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-am_covariance_structure.R")'`

Expected: PASS, 5 tests, 0 failures.

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION tests/testthat.R tests/testthat/test-am_covariance_structure.R R/am_covariance_structure.R
git commit -m "support negative and zero r in am_covariance_structure"
```

---

### Task 2: Signed low-rank term in `rb_dplr()`

**Files:**
- Modify: `R/rb_dplr.R`
- Test: `tests/testthat/test-rb_dplr.R`

**Interfaces:**
- Consumes: `attr(U, "sign")` set by Task 1.
- Produces: `rb_dplr(n, mu, U, sign = NULL)` returning an `n` by `length(mu)` numeric 0/1 matrix, and the internal helper `.rb_sign(U)` returning `1` or `-1`, used by Task 7.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-rb_dplr.R`:

```r
test_that(".rb_sign resolves the attribute with a positive fallback", {
  expect_identical(rBahadur:::.rb_sign(structure(c(1, 2), sign = -1)), -1)
  expect_identical(rBahadur:::.rb_sign(structure(c(1, 2), sign = 1)), 1)
  expect_identical(rBahadur:::.rb_sign(c(1, 2)), 1)
})

test_that("an explicit sign argument overrides the attribute", {
  set.seed(4)
  m <- 40
  beta <- as.vector(scale(rnorm(m))) * sqrt(0.5 / m)
  AF <- runif(m, 0.2, 0.8)
  U <- am_covariance_structure(beta, AF, 0.4)
  mu <- rep(AF, each = 2)

  set.seed(5); a <- rb_dplr(50, mu, U, sign = -1)
  set.seed(5); b <- rb_dplr(50, mu, structure(as.vector(U), sign = -1))
  expect_identical(a, b)

  set.seed(5); pos <- rb_dplr(50, mu, U)
  expect_false(identical(a, pos))
})

test_that("sign = 1 is bit-identical to passing no sign at all", {
  set.seed(6)
  m <- 40
  beta <- as.vector(scale(rnorm(m))) * sqrt(0.5 / m)
  AF <- runif(m, 0.2, 0.8)
  U <- as.vector(am_covariance_structure(beta, AF, 0.4))
  mu <- rep(AF, each = 2)

  set.seed(8); a <- rb_dplr(30, mu, U)
  set.seed(8); b <- rb_dplr(30, mu, U, sign = 1)
  expect_identical(a, b)
})

test_that("an invalid sign is rejected", {
  expect_error(rb_dplr(5, rep(0.5, 6), rep(0.01, 6), sign = 0), "either 1 or -1")
})

test_that("negative-r draws induce negative linkage disequilibrium", {
  skip_on_cran()
  set.seed(9)
  m <- 800
  beta <- as.vector(scale(rnorm(m))) * sqrt(0.5 / m)
  AF <- runif(m, 0.2, 0.8)
  mu <- rep(AF, each = 2)

  ## r is held at 0.4 in magnitude: negative assortment leaves the Bahadur
  ## feasible region well before positive assortment does, and r = -0.6 with
  ## this n fails for a meaningful fraction of seeds
  Uneg <- am_covariance_structure(beta, AF, -0.4)
  Upos <- am_covariance_structure(beta, AF, 0.4)
  bu <- beta / sqrt(2 * AF * (1 - AF))

  gv <- function(U) {
    H <- rb_dplr(1500, mu, U)
    X <- H[, seq(1, 2 * m, 2)] + H[, seq(2, 2 * m, 2)]
    var(as.vector(X %*% bu))
  }
  ## disassortment strips genetic variance, assortment inflates it
  expect_lt(gv(Uneg), 0.5)
  expect_gt(gv(Upos), 0.5)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-rb_dplr.R")'`

Expected: FAIL with `unused argument (sign = -1)` and `could not find function ".rb_sign"`.

- [ ] **Step 3: Add the sign resolver**

Add to the top of `R/rb_dplr.R`, above the roxygen block:

```r
## Resolve the low-rank sign for a loading vector, defaulting to +1 so that
## hand-built `U` vectors and pre-1.1.0 code keep working.
.rb_sign <- function(U) {
  s <- attr(U, "sign")
  if (is.null(s)) 1 else s
}
```

- [ ] **Step 4: Add the argument and apply it**

Change the signature of `rb_dplr()` from `function(n, mu, U)` to:

```r
rb_dplr <- function(n, mu, U, sign = NULL) {

  if (is.null(sign)) sign <- .rb_sign(U)
  if (length(sign) != 1L || !sign %in% c(-1, 1)) {
    stop("`sign` must be either 1 or -1")
  }
  ## bind to a local name so the argument does not shadow base::sign()
  s <- sign

  M <- length(mu)
```

Then change exactly two lines in the existing body. Inside the recursion loop:

```r
    p <- mu[m] + s * x * U[m]
```

and the final step after the loop:

```r
  p <- mu[M] + s * x * U[M]
```

Leave the `x` and `c` accumulator updates untouched. Under `U = i*V` induction gives `x = i * x_V` at every step while `c` stays real, so negating at these two lines is exactly equivalent to running the recursion in complex arithmetic.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-rb_dplr.R")'`

Expected: PASS, 5 tests, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add R/rb_dplr.R tests/testthat/test-rb_dplr.R
git commit -m "add signed low-rank term to rb_dplr"
```

---

### Task 3: Negative assortment end to end through `am_simulate()`

**Files:**
- Test: `tests/testthat/test-am_simulate-negative.R`

**Interfaces:**
- Consumes: Tasks 1 and 2. No source change is required, since `r` already flows from `am_simulate()` into `am_covariance_structure()` and the resulting `sign` attribute already flows into `rb_dplr()`.
- Produces: nothing new. This task exists to prove the integration and to gate the feature before genotype output work begins.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-am_simulate-negative.R`:

```r
test_that("am_simulate runs for negative, zero, and positive r", {
  set.seed(21)
  for (r in c(-0.5, 0, 0.5)) {
    d <- am_simulate(h2_0 = 0.5, r = r, m = 200, n = 300)
    expect_true(all(is.finite(d$X)))
    expect_true(all(d$X %in% c(0, 1, 2)))
    expect_equal(dim(d$X), c(300L, 200L))
  }
})

test_that("empirical heritability tracks h2_eq across the sign range", {
  skip_on_cran()
  set.seed(22)
  h2_0 <- 0.5
  ## negative r is held within the reliable Bahadur feasible range; see the
  ## Global Constraints note on the negative-assortment envelope
  for (r in c(-0.3, -0.2, 0.3, 0.6)) {
    d <- am_simulate(h2_0 = h2_0, r = r, m = 1500, n = 4000)
    emp <- var(as.vector(d$g)) / var(as.vector(d$y))
    expect_equal(emp, h2_eq(r, h2_0), tolerance = 0.05)
  }
})

test_that("negative r reduces genetic variance and positive r inflates it", {
  skip_on_cran()
  set.seed(23)
  h2_0 <- 0.5
  vneg <- var(as.vector(am_simulate(h2_0, -0.3, 1500, 4000)$g))
  vpos <- var(as.vector(am_simulate(h2_0, 0.6, 1500, 4000)$g))
  expect_lt(vneg, h2_0)
  expect_gt(vpos, h2_0)
  expect_equal(vneg, vg_eq(-0.3, h2_0, h2_0), tolerance = 0.05)
  expect_equal(vpos, vg_eq(0.6, h2_0, h2_0), tolerance = 0.05)
})

test_that("allele frequencies are preserved under negative r", {
  skip_on_cran()
  set.seed(24)
  d <- am_simulate(0.5, -0.3, 1500, 4000)
  expect_gt(cor(d$AF, colMeans(d$X) / 2), 0.99)
})
```

- [ ] **Step 2: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-am_simulate-negative.R")'`

Expected: PASS, 4 tests, 0 failures. If any fail, the defect is in Task 1 or Task 2, not here. Do not add compensating logic to `am_simulate()`.

- [ ] **Step 3: Run the whole suite**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`

Expected: PASS, 14 tests, 0 failures.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-am_simulate-negative.R
git commit -m "cover negative assortment end to end in am_simulate"
```

---

### Task 4: int8 genotype layouts and sidecar metadata

**Files:**
- Create: `R/genotype_io.R`
- Test: `tests/testthat/test-genotype-io.R`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `write_genotypes(X, path, format = c("individual", "variant", "bed"))`, returns `path` invisibly.
  - `read_genotypes(path)`, returns an `n` by `m` integer matrix.
  - `.gt_data_path(path, format)`, `.gt_write_meta(path, n, m, format)`, `.gt_read_meta(path)` used by Task 7.
  - `path` is a prefix, not a filename. `individual` and `variant` write `<path>.int8`, `bed` writes `<path>.bed`. All write `<path>.meta`.

This task implements `individual` and `variant` only. Task 5 adds `bed`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-genotype-io.R`:

```r
test_that("individual and variant layouts round trip", {
  set.seed(31)
  X <- matrix(sample(0:2, 7 * 5, replace = TRUE), nrow = 7, ncol = 5)
  for (fmt in c("individual", "variant")) {
    p <- file.path(tempdir(), paste0("gt-", fmt))
    write_genotypes(X, p, format = fmt)
    expect_true(file.exists(paste0(p, ".int8")))
    expect_equal(read_genotypes(p), X)
  }
})

test_that("on-disk byte order matches the declared layout", {
  X <- matrix(c(0L, 1L, 2L,
                2L, 1L, 0L), nrow = 2, ncol = 3, byrow = TRUE)

  pv <- file.path(tempdir(), "gt-order-variant")
  write_genotypes(X, pv, format = "variant")
  expect_equal(readBin(paste0(pv, ".int8"), "integer", 6L, size = 1L),
               c(0L, 2L, 1L, 1L, 2L, 0L))

  pi <- file.path(tempdir(), "gt-order-individual")
  write_genotypes(X, pi, format = "individual")
  expect_equal(readBin(paste0(pi, ".int8"), "integer", 6L, size = 1L),
               c(0L, 1L, 2L, 2L, 1L, 0L))
})

test_that("file size is exactly one byte per genotype", {
  X <- matrix(1L, nrow = 11, ncol = 13)
  p <- file.path(tempdir(), "gt-size")
  write_genotypes(X, p, format = "variant")
  expect_equal(file.size(paste0(p, ".int8")), 11 * 13)
})

test_that("metadata round trips", {
  X <- matrix(0L, nrow = 4, ncol = 6)
  p <- file.path(tempdir(), "gt-meta")
  write_genotypes(X, p, format = "individual")
  meta <- rBahadur:::.gt_read_meta(p)
  expect_identical(meta$format, "individual")
  expect_identical(meta$n, 4L)
  expect_identical(meta$m, 6L)
  expect_identical(meta$dtype, "int8")
})

test_that("invalid genotypes and missing values are rejected", {
  p <- file.path(tempdir(), "gt-bad")
  expect_error(write_genotypes(matrix(3L, 2, 2), p, "variant"), "0, 1, or 2")
  expect_error(write_genotypes(matrix(NA_integer_, 2, 2), p, "variant"), "bed")
  expect_error(write_genotypes(1:4, p, "variant"), "matrix")
})

test_that("reading without metadata fails clearly", {
  expect_error(read_genotypes(file.path(tempdir(), "gt-does-not-exist")),
               "metadata file not found")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-genotype-io.R")'`

Expected: FAIL with `could not find function "write_genotypes"`.

- [ ] **Step 3: Create the file**

Create `R/genotype_io.R`:

```r
## ---- internal path and metadata helpers -------------------------------

.gt_data_path <- function(path, format) {
  if (format == "bed") paste0(path, ".bed") else paste0(path, ".int8")
}

.gt_write_meta <- function(path, n, m, format) {
  writeLines(c(
    "rBahadur_genotypes: 1",
    paste0("format: ", format),
    paste0("n: ", n),
    paste0("m: ", m),
    paste0("dtype: ", if (format == "bed") "bed2bit" else "int8")
  ), paste0(path, ".meta"))
  invisible(NULL)
}

.gt_read_meta <- function(path) {
  f <- paste0(path, ".meta")
  if (!file.exists(f)) stop("metadata file not found: ", f)
  kv <- strsplit(readLines(f), ":[[:space:]]*")
  keys <- vapply(kv, function(x) x[1], character(1))
  vals <- vapply(kv, function(x) x[2], character(1))
  names(vals) <- keys
  list(format = unname(vals["format"]),
       n = as.integer(vals["n"]),
       m = as.integer(vals["m"]),
       dtype = unname(vals["dtype"]))
}

.gt_check_matrix <- function(X, format) {
  if (!is.matrix(X)) stop("`X` must be a matrix")
  storage.mode(X) <- "integer"
  if (anyNA(X) && format != "bed") {
    stop("missing genotypes are only representable in the 'bed' format")
  }
  if (any(X < 0L | X > 2L, na.rm = TRUE)) stop("genotypes must be 0, 1, or 2")
  X
}

## ---- exported interface -----------------------------------------------

#' Write genotypes to a binary file
#'
#' @param X an integer matrix of genotypes with individuals in rows and
#'   variants in columns, taking values 0, 1, or 2
#' @param path file prefix. Layout `"individual"` and `"variant"` write
#'   `<path>.int8`; `"bed"` writes `<path>.bed` plus `<path>.bim` and
#'   `<path>.fam`. All layouts write `<path>.meta`.
#' @param format on-disk layout. `"individual"` (the default) stores each
#'   individual's variants contiguously, `"variant"` stores each variant's
#'   individuals contiguously, and `"bed"` writes a variant-major PLINK
#'   binary file at two bits per genotype.
#'
#' @return `path`, invisibly.
#' @export
#'
#' @examples
#' X <- matrix(sample(0:2, 20, replace = TRUE), nrow = 4)
#' p <- file.path(tempdir(), "example_genotypes")
#' write_genotypes(X, p)
#' identical(read_genotypes(p), X)
write_genotypes <- function(X, path, format = c("individual", "variant", "bed")) {
  format <- match.arg(format)
  X <- .gt_check_matrix(X, format)
  n <- nrow(X)
  m <- ncol(X)

  con <- file(.gt_data_path(path, format), "wb")
  on.exit(close(con))
  if (format == "variant") {
    writeBin(as.vector(X), con, size = 1L)
  } else if (format == "individual") {
    writeBin(as.vector(t(X)), con, size = 1L)
  } else {
    stop("the 'bed' format is not implemented yet")
  }
  .gt_write_meta(path, n, m, format)
  invisible(path)
}

#' Read genotypes from a binary file written by `write_genotypes()`
#'
#' @param path file prefix, the same value passed to [write_genotypes()]
#'
#' @return An integer matrix with individuals in rows and variants in columns,
#'   reconstructed to the same orientation regardless of the on-disk layout.
#' @export
#'
#' @examples
#' X <- matrix(sample(0:2, 20, replace = TRUE), nrow = 4)
#' p <- file.path(tempdir(), "example_genotypes")
#' write_genotypes(X, p, format = "variant")
#' read_genotypes(p)
read_genotypes <- function(path) {
  meta <- .gt_read_meta(path)
  n <- meta$n
  m <- meta$m
  con <- file(.gt_data_path(path, meta$format), "rb")
  on.exit(close(con))
  if (meta$format == "bed") {
    stop("the 'bed' format is not implemented yet")
  }
  v <- readBin(con, "integer", n = n * m, size = 1L, signed = TRUE)
  if (meta$format == "variant") {
    matrix(v, nrow = n, ncol = m)
  } else {
    t(matrix(v, nrow = m, ncol = n))
  }
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-genotype-io.R")'`

Expected: PASS, 6 tests, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/genotype_io.R tests/testthat/test-genotype-io.R
git commit -m "add int8 genotype layouts and sidecar metadata"
```

---

### Task 5: PLINK `.bed` format

**Files:**
- Modify: `R/genotype_io.R`
- Test: `tests/testthat/test-genotype-bed.R`

**Interfaces:**
- Consumes: `.gt_data_path()`, `.gt_write_meta()`, `.gt_read_meta()`, `.gt_check_matrix()` from Task 4.
- Produces: `.gt_pack_bed(g)` taking an integer dosage vector of length `n` and returning `ceiling(n/4)` raw bytes, `.gt_unpack_bed(bytes, n)` taking raw bytes and returning an integer vector of length `n`, and `.gt_write_plink_sidecars(path, n, m)`. Task 7 calls `.gt_pack_bed()` and `.gt_write_plink_sidecars()` directly.

Encoding contract, documented because A1/A2 handling is the usual source of error. The effect allele is written as A1 in the `.bim`, and PLINK codes count A1 copies, so dosage 2 maps to bits `00`, dosage 1 to `10`, dosage 0 to `11`, and missing to `01`. Four samples pack into one byte, lowest sample in the lowest bits.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-genotype-bed.R`:

```r
test_that("bed packing matches a hand-computed byte", {
  ## dosages 2,1,0,2 -> codes 0,2,3,0 -> 0*1 + 2*4 + 3*16 + 0*64 = 56
  expect_equal(as.integer(rBahadur:::.gt_pack_bed(c(2L, 1L, 0L, 2L))), 56L)
})

test_that("bed packing round trips including partial final bytes", {
  set.seed(41)
  for (n in c(1, 3, 4, 5, 17, 100)) {
    g <- sample(0:2, n, replace = TRUE)
    packed <- rBahadur:::.gt_pack_bed(as.integer(g))
    expect_equal(length(packed), ceiling(n / 4))
    expect_equal(rBahadur:::.gt_unpack_bed(packed, n), as.integer(g))
  }
})

test_that("missing genotypes round trip through bed", {
  g <- c(0L, NA_integer_, 2L, 1L, NA_integer_)
  expect_equal(rBahadur:::.gt_unpack_bed(rBahadur:::.gt_pack_bed(g), 5L), g)
})

test_that("bed files round trip and carry a valid header", {
  set.seed(42)
  X <- matrix(sample(0:2, 9 * 6, replace = TRUE), nrow = 9, ncol = 6)
  p <- file.path(tempdir(), "gt-bed")
  write_genotypes(X, p, format = "bed")

  hdr <- readBin(paste0(p, ".bed"), "raw", 3L)
  expect_equal(as.integer(hdr), c(0x6c, 0x1b, 0x01))
  expect_equal(file.size(paste0(p, ".bed")), 3 + 6 * ceiling(9 / 4))
  expect_equal(read_genotypes(p), X)
})

test_that("bed writes plink sidecars with correct dimensions", {
  X <- matrix(0L, nrow = 5, ncol = 4)
  p <- file.path(tempdir(), "gt-bed-sidecar")
  write_genotypes(X, p, format = "bed")
  expect_length(readLines(paste0(p, ".bim")), 4L)
  expect_length(readLines(paste0(p, ".fam")), 5L)
  expect_equal(length(strsplit(readLines(paste0(p, ".bim"))[1], "\t")[[1]]), 6L)
  expect_equal(length(strsplit(readLines(paste0(p, ".fam"))[1], "\t")[[1]]), 6L)
})

test_that("a bed file with a bad header is rejected", {
  p <- file.path(tempdir(), "gt-bed-bad")
  write_genotypes(matrix(0L, 4, 2), p, format = "bed")
  con <- file(paste0(p, ".bed"), "r+b")
  writeBin(as.raw(c(0x00, 0x00, 0x00)), con)
  close(con)
  expect_error(read_genotypes(p), "variant-major PLINK")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-genotype-bed.R")'`

Expected: FAIL with `could not find function ".gt_pack_bed"` and `the 'bed' format is not implemented yet`.

- [ ] **Step 3: Add the bed helpers**

Insert into `R/genotype_io.R`, after `.gt_check_matrix()`:

```r
## ---- PLINK bed helpers -------------------------------------------------
##
## Two bits per genotype, four samples per byte, lowest sample in the lowest
## bits. The effect allele is written as A1, and PLINK codes count A1 copies,
## so dosage 2 -> 00, dosage 1 -> 10, dosage 0 -> 11, and missing -> 01.

.gt_pack_bed <- function(g) {
  code <- integer(length(g))
  code[!is.na(g) & g == 2L] <- 0L
  code[!is.na(g) & g == 1L] <- 2L
  code[!is.na(g) & g == 0L] <- 3L
  code[is.na(g)] <- 1L
  pad <- (4L - (length(g) %% 4L)) %% 4L
  if (pad > 0L) code <- c(code, integer(pad))
  quad <- matrix(code, nrow = 4L)
  as.raw(quad[1, ] + quad[2, ] * 4L + quad[3, ] * 16L + quad[4, ] * 64L)
}

.gt_unpack_bed <- function(bytes, n) {
  b <- as.integer(bytes)
  codes <- as.vector(rbind(b %% 4L,
                           (b %/% 4L) %% 4L,
                           (b %/% 16L) %% 4L,
                           (b %/% 64L) %% 4L))[seq_len(n)]
  g <- integer(n)
  g[codes == 0L] <- 2L
  g[codes == 2L] <- 1L
  g[codes == 3L] <- 0L
  g[codes == 1L] <- NA_integer_
  g
}

.gt_write_plink_sidecars <- function(path, n, m) {
  writeLines(paste(1L, paste0("v", seq_len(m)), 0L, seq_len(m), "A", "G",
                   sep = "\t"), paste0(path, ".bim"))
  ids <- paste0("i", seq_len(n))
  writeLines(paste(ids, ids, 0L, 0L, 0L, -9L, sep = "\t"),
             paste0(path, ".fam"))
  invisible(NULL)
}
```

- [ ] **Step 4: Wire bed into the writer**

In `write_genotypes()`, replace the `stop("the 'bed' format is not implemented yet")` branch with:

```r
  } else {
    writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
    for (j in seq_len(m)) writeBin(.gt_pack_bed(X[, j]), con)
    .gt_write_plink_sidecars(path, n, m)
  }
```

- [ ] **Step 5: Wire bed into the reader**

In `read_genotypes()`, replace the `stop("the 'bed' format is not implemented yet")` branch with:

```r
  if (meta$format == "bed") {
    hdr <- readBin(con, "raw", n = 3L)
    if (!identical(as.integer(hdr), c(0x6cL, 0x1bL, 0x01L))) {
      stop("not a variant-major PLINK .bed file")
    }
    nb <- ceiling(n / 4)
    X <- matrix(NA_integer_, nrow = n, ncol = m)
    for (j in seq_len(m)) X[, j] <- .gt_unpack_bed(readBin(con, "raw", n = nb), n)
    return(X)
  }
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-genotype-bed.R")'`

Expected: PASS, 6 tests, 0 failures.

- [ ] **Step 7: Commit**

```bash
git add R/genotype_io.R tests/testthat/test-genotype-bed.R
git commit -m "add PLINK bed genotype output"
```

---

### Task 6: Locus-blocked streaming recursion

**Files:**
- Create: `R/am_stream.R`
- Test: `tests/testthat/test-am_stream.R`

**Interfaces:**
- Consumes: `.rb_sign()` from Task 2.
- Produces: `.rb_dplr_stream(n, mu, U, s = 1, block = 1024L, callback)`. The callback is invoked as `callback(B, col0)` where `B` is an `n` by `k` numeric 0/1 matrix and `col0` is the 1-based global index of `B`'s first column. Returns `NULL` invisibly.

The key property, verified before this plan was written, is that drawing `runif(n)` once per column consumes the random stream in exactly the order `matrix(runif(M*n), n, M)` does, because R fills matrices column-major. `.rb_dplr_stream()` is therefore bit-identical to `rb_dplr()` at any block size.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-am_stream.R`:

```r
arch <- function(m = 60, h2_0 = 0.5, seed = 11) {
  set.seed(seed)
  beta <- as.vector(scale(rnorm(m))) * sqrt(h2_0 / m)
  AF <- runif(m, 0.15, 0.85)
  list(beta = beta, AF = AF, mu = rep(AF, each = 2), m = m)
}

collect <- function(n, a, U, block) {
  acc <- matrix(0, nrow = n, ncol = 2 * a$m)
  rBahadur:::.rb_dplr_stream(
    n, a$mu, U, s = rBahadur:::.rb_sign(U), block = block,
    callback = function(B, col0) acc[, col0:(col0 + ncol(B) - 1L)] <<- B)
  acc
}

test_that("streaming is bit-identical to rb_dplr at every block size", {
  a <- arch()
  for (r in c(0.4, -0.4)) {
    U <- am_covariance_structure(a$beta, a$AF, r)
    set.seed(99); ref <- rb_dplr(40, a$mu, U)
    for (blk in c(2L, 7L, 16L, 2L * a$m, 5L * a$m)) {
      set.seed(99)
      expect_identical(collect(40, a, U, blk), ref)
    }
  }
})

test_that("callback receives contiguous blocks covering every column once", {
  a <- arch()
  U <- am_covariance_structure(a$beta, a$AF, 0.4)
  seen <- integer(0)
  set.seed(2)
  rBahadur:::.rb_dplr_stream(
    10, a$mu, U, s = 1, block = 7L,
    callback = function(B, col0) seen <<- c(seen, col0:(col0 + ncol(B) - 1L)))
  expect_identical(seen, seq_len(2L * a$m))
})

test_that("infeasible probabilities are reported with the offending column", {
  mu <- rep(0.5, 20)
  U <- rep(0.9, 20)
  expect_error(
    rBahadur:::.rb_dplr_stream(5, mu, U, s = 1, block = 4L,
                               callback = function(B, col0) NULL),
    "Infeasible probabilities")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-am_stream.R")'`

Expected: FAIL with `could not find function ".rb_dplr_stream"`.

- [ ] **Step 3: Create the file**

Create `R/am_stream.R`:

```r
## Locus-blocked form of the rb_dplr recursion.
##
## Emits column blocks through `callback` instead of allocating the full
## n-by-M matrix. Because R fills matrices column-major, drawing runif(n) once
## per column consumes the random stream in the same order that
## matrix(runif(M*n), n, M) does, so this is bit-identical to rb_dplr() at any
## block size.
##
## @param callback invoked as callback(B, col0); B is an n-by-k 0/1 matrix and
##   col0 is the 1-based global index of its first column.
.rb_dplr_stream <- function(n, mu, U, s = 1, block = 1024L, callback) {
  M <- length(mu)
  block <- max(1L, min(as.integer(block), M))
  buf <- matrix(0L, nrow = n, ncol = block)
  bi <- 0L
  col0 <- 1L
  x <- NULL
  cc <- 1

  for (m in seq_len(M)) {
    p <- if (m == 1L) rep(mu[1], n) else mu[m] + s * x * U[m]
    if (any(!is.finite(p) | p < 0 | p > 1)) {
      stop("Infeasible probabilities at column ", m,
           ". Are you sure the specified parameters correspond to a valid ",
           "Bahadur order-2 MVB distribution?")
    }
    km <- as.integer(runif(n) <= p)
    bi <- bi + 1L
    buf[, bi] <- km

    tmp_bool <- (km == 0L)
    pc <- tmp_bool * (1 - p) + (!tmp_bool) * p
    Bk0 <- tmp_bool * (1 - mu[m]) + (!tmp_bool) * mu[m]
    Bk1 <- tmp_bool * (-1) + (!tmp_bool) * 1
    if (m == 1L) {
      x <- Bk1 * U[1] / pc
      cc <- 1
    } else {
      x <- (x * Bk0 + cc * Bk1 * U[m]) / pc
      cc <- (Bk0 / pc) * cc
    }

    if (bi == block || m == M) {
      callback(buf[, seq_len(bi), drop = FALSE], col0)
      col0 <- m + 1L
      bi <- 0L
    }
  }
  invisible(NULL)
}
```

The block buffer uses `integer` storage, which is half the width of the
`matrix(NaN, ...)` double that `rb_dplr()` allocates. Bit-identity with
`rb_dplr()` is unaffected: the test accumulator in `collect()` is a double
matrix, and assigning an integer block into it coerces the values to double
without changing them, so `expect_identical()` against `rb_dplr()`'s double
output still holds. The dominant saving is the buffer being `n` by `block`
rather than `n` by `2*m`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-am_stream.R")'`

Expected: PASS, 3 tests, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/am_stream.R tests/testthat/test-am_stream.R
git commit -m "add locus-blocked streaming recursion"
```

---

### Task 7: Streaming output in `am_simulate()`

**Files:**
- Modify: `R/am_simulate.R`
- Test: `tests/testthat/test-am_simulate-stream.R`

**Interfaces:**
- Consumes: `.rb_dplr_stream()` (Task 6), `.rb_sign()` (Task 2), `.gt_data_path()`, `.gt_write_meta()`, `.gt_pack_bed()`, `.gt_write_plink_sidecars()` (Tasks 4 and 5).
- Produces: `am_simulate(h2_0, r, m, n, afs, min_MAF, haplotypes, path, format, batch_size)`. With `path` supplied it returns, invisibly, a list with `y`, `g`, `AF`, `beta_std`, `beta_raw`, `path`, `format`, `n`, `m`, and no `X`, and it also saves that list to `<path>.rds`.

Batching follows the layout. `individual` batches over people, because rows are independent draws and each batch's rows write contiguously. `variant` and `bed` batch over loci through `.rb_dplr_stream()`. Consequently `variant` and `bed` reproduce the in-memory genotypes at any `batch_size`, while `individual` does so only when `batch_size >= n`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-am_simulate-stream.R`:

```r
inmem <- function(seed, r = 0.4, m = 40, n = 24) {
  set.seed(seed)
  am_simulate(h2_0 = 0.5, r = r, m = m, n = n)
}
streamed <- function(seed, fmt, batch, r = 0.4, m = 40, n = 24) {
  p <- file.path(tempdir(), paste0("sim-", fmt, "-", batch, "-", seed))
  set.seed(seed)
  out <- am_simulate(h2_0 = 0.5, r = r, m = m, n = n,
                     path = p, format = fmt, batch_size = batch)
  list(out = out, X = read_genotypes(p), path = p)
}

test_that("locus-batched layouts reproduce the in-memory genotypes at any batch size", {
  ref <- inmem(51)
  for (fmt in c("variant", "bed")) {
    for (b in c(1L, 3L, 40L, 400L)) {
      s <- streamed(51, fmt, b)
      expect_equal(s$X, matrix(as.integer(ref$X), nrow = 24),
                   info = paste(fmt, b))
      expect_equal(as.vector(s$out$y), as.vector(ref$y), info = paste(fmt, b))
    }
  }
})

test_that("individual layout reproduces the in-memory result when batch_size >= n", {
  ref <- inmem(52)
  s <- streamed(52, "individual", 24L)
  expect_equal(s$X, matrix(as.integer(ref$X), nrow = 24))
  expect_equal(as.vector(s$out$y), as.vector(ref$y))
})

test_that("individual layout with a smaller batch is still a valid simulation", {
  s <- streamed(53, "individual", 5L, m = 40, n = 24)
  expect_equal(dim(s$X), c(24L, 40L))
  expect_true(all(s$X %in% c(0L, 1L, 2L)))
})

test_that("streaming works for negative r", {
  ref <- inmem(54, r = -0.5)
  s <- streamed(54, "variant", 7L, r = -0.5)
  expect_equal(s$X, matrix(as.integer(ref$X), nrow = 24))
})

test_that("streaming omits X, returns dimensions, and saves an rds", {
  s <- streamed(55, "variant", 8L)
  expect_null(s$out$X)
  expect_identical(s$out$n, 24L)
  expect_identical(s$out$m, 40L)
  expect_identical(s$out$format, "variant")
  expect_true(file.exists(paste0(s$path, ".rds")))
  expect_equal(readRDS(paste0(s$path, ".rds"))$AF, s$out$AF)
})

test_that("path = NULL leaves the in-memory result untouched", {
  a <- inmem(56)
  b <- inmem(56)
  expect_identical(a, b)
  expect_true(is.matrix(a$X))
})

test_that("haplotypes cannot be combined with streaming", {
  expect_error(
    am_simulate(0.5, 0.4, 20, 10, haplotypes = TRUE,
                path = file.path(tempdir(), "sim-hap")),
    "haplotypes")
})

test_that("an auto batch size is chosen when none is given", {
  s <- streamed(57, "variant", NULL)
  expect_equal(dim(s$X), c(24L, 40L))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-am_simulate-stream.R")'`

Expected: FAIL with `unused arguments (path = ..., format = ..., batch_size = ...)`.

- [ ] **Step 3: Replace the function body**

Replace `am_simulate()` in `R/am_simulate.R` (keep the roxygen block, which Task 8 updates):

```r
am_simulate <- function(h2_0, r, m, n, afs = NULL, min_MAF = .1,
                        haplotypes = FALSE, path = NULL,
                        format = c("individual", "variant", "bed"),
                        batch_size = NULL) {
  format <- match.arg(format)

  ## draw standardized diploid allele substitution effects
  beta <- scale(rnorm(m))*sqrt(h2_0 / m)
  ## draw allele frequencies if necessary
  if (is.null(afs)) {
    AF <- runif(m, min_MAF, 1 - min_MAF)
  } else if (length(afs) != m) {
    stop("`afs` must have length `m`")
  } else {
    AF <- afs
  }
  ## compute unstandardized effects
  beta_unscaled <- beta/sqrt(2*AF*(1-AF))
  ## generate corresponding haploid quantities
  AF_hap <- rep(AF, each=2)
  ## compute equilibrium outer product covariance component
  U <- am_covariance_structure(beta, AF, r)

  if (is.null(path)) {
    ## draw multivariate Bernoulli haplotypes
    H <- rb_dplr(n, AF_hap, U)
    ## convert haplotypes to diploid genotypes
    X <- (H[,seq(1,2*m,2)]+H[,seq(2,2*m,2)])
    ## compute genetic phenotypes
    g <- X %*% beta_unscaled
    ## compute full phenotype
    y <- g + rnorm(n, 0, sqrt(1 - h2_0))
    output <- list(
      y = y,
      g = g,
      X = X,
      AF = AF,
      beta_std = beta,
      beta_raw = beta_unscaled
      )
    if (haplotypes) {
      output$H <- H
    }
    return(output)
  }

  if (haplotypes) {
    stop("`haplotypes = TRUE` is not supported when streaming to `path`")
  }

  s <- .rb_sign(U)
  g <- numeric(n)
  con <- file(.gt_data_path(path, format), "wb")
  on.exit(close(con))

  if (format == "individual") {
    ## rows are independent draws, so batch over individuals; this writes each
    ## individual's variants contiguously
    if (is.null(batch_size)) {
      batch_size <- max(1L, min(n, as.integer(floor(128e6 / (16 * m)))))
    }
    for (start in seq(1L, n, by = batch_size)) {
      nb <- min(batch_size, n - start + 1L)
      Hb <- rb_dplr(nb, AF_hap, U)
      Xb <- Hb[, seq(1, 2*m, 2), drop = FALSE] + Hb[, seq(2, 2*m, 2), drop = FALSE]
      g[start:(start + nb - 1L)] <- as.vector(Xb %*% beta_unscaled)
      storage.mode(Xb) <- "integer"
      writeBin(as.vector(t(Xb)), con, size = 1L)
    }
  } else {
    ## variant-major layouts batch over loci and carry the recursion state
    if (is.null(batch_size)) {
      batch_size <- max(1L, min(m, as.integer(floor(128e6 / (8 * n)))))
    }
    if (format == "bed") writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
    .rb_dplr_stream(
      n, AF_hap, U, s = s, block = 2L * as.integer(batch_size),
      callback = function(B, col0) {
        nb <- ncol(B)
        Xb <- B[, seq(1, nb, 2), drop = FALSE] + B[, seq(2, nb, 2), drop = FALSE]
        loc0 <- (col0 + 1L) %/% 2L
        idx <- loc0:(loc0 + ncol(Xb) - 1L)
        g <<- g + as.vector(Xb %*% beta_unscaled[idx])
        storage.mode(Xb) <- "integer"
        if (format == "bed") {
          for (j in seq_len(ncol(Xb))) writeBin(.gt_pack_bed(Xb[, j]), con)
        } else {
          writeBin(as.vector(Xb), con, size = 1L)
        }
      })
  }

  g <- matrix(g, ncol = 1)
  y <- g + rnorm(n, 0, sqrt(1 - h2_0))
  .gt_write_meta(path, n, m, format)
  if (format == "bed") .gt_write_plink_sidecars(path, n, m)
  output <- list(
    y = y,
    g = g,
    AF = AF,
    beta_std = beta,
    beta_raw = beta_unscaled,
    path = path,
    format = format,
    n = as.integer(n),
    m = as.integer(m)
  )
  saveRDS(output, paste0(path, ".rds"))
  invisible(output)
}
```

Two ordering details matter and must not be changed. The `rnorm(m)`, `runif(m)`, sampler, `rnorm(n)` sequence is preserved in both branches, which is what makes the streamed and in-memory results comparable under one seed. The block passed to `.rb_dplr_stream()` is `2 * batch_size` haplotype columns, which stays even so every emitted block splits cleanly into diploid pairs.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-am_simulate-stream.R")'`

Expected: PASS, 8 tests, 0 failures.

- [ ] **Step 5: Run the whole suite**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`

Expected: 37 `test_that` blocks, 0 failures, 0 skips when run locally. Several
Monte Carlo tests carry `skip_on_cran()`, so a CRAN run reports skips instead.
If one of those tests fails marginally, raise `n` in that test rather than
loosening its tolerance, since a loosened tolerance stops detecting the bug it
was written for.

- [ ] **Step 6: Commit**

```bash
git add R/am_simulate.R tests/testthat/test-am_simulate-stream.R
git commit -m "stream genotypes to disk from am_simulate"
```

---

### Task 8: Documentation, namespace, and release metadata

**Files:**
- Modify: `R/am_covariance_structure.R`, `R/rb_dplr.R`, `R/am_simulate.R` (roxygen blocks only)
- Modify: `DESCRIPTION`, `NEWS.md`, `README.md`
- Regenerate: `NAMESPACE`, `man/`

- [ ] **Step 1: Update the roxygen for `am_covariance_structure()`**

Change the `@param r` line and the `@return` block to:

```r
#' @param r cross-mate phenotypic correlation, in the open interval (-1, 1).
#'   Negative values correspond to disassortative mating.
#'
#' @return Vector 'U' such that \eqn{D + s U U^T} corresponds to the expected
#' haploid LD-matrix given the specified genetic architecture (encoded by 'beta'
#' and 'AF') and cross-mate phenotypic correlation 'r', where the sign
#' \eqn{s} is `attr(U, "sign")`. It is assumed that the total phenotypic
#' variance at generation zero is one.
#'
#' @details For `r > 0` the low-rank term is added and `attr(U, "sign")` is 1.
#' For `r < 0` it is subtracted and the attribute is -1, because a positive
#' semidefinite rank-one term can only increase genetic variance whereas
#' disassortative mating reduces it. The returned vector in that case is the
#' analytic continuation of the positive branch, which is purely imaginary.
#' For `r = 0` the vector is zero, since panmixia induces no disequilibrium.
#'
#' @section Feasibility under negative assortment:
#' Disassortative mating leaves the Bahadur order-2 feasible region sooner
#' than assortative mating does. The returned vector can satisfy the
#' discriminant condition checked here and still drive [rb_dplr()] outside
#' \[0, 1\] during sampling, because feasibility there is a property of the
#' realized draws rather than of the parameters alone. Infeasibility becomes
#' more likely as `n` grows, since it only takes one individual to fall
#' outside the region. Positive `r` was feasible in every configuration
#' tested. Negative `r` was not: at `h2_0 = 0.5` and 1500 causal variants,
#' `r = -0.3` sampled reliably at 4000 individuals, while `r = -0.4` sampled
#' reliably at 2000 individuals but failed for some seeds at 4000. Raising
#' `min_MAF` widens the envelope. If [rb_dplr()] reports infeasible
#' probabilities, reduce the magnitude of `r`, raise `min_MAF`, or increase
#' the number of causal variants.
```

- [ ] **Step 2: Update the roxygen for `rb_dplr()`**

Add after the `@param U` line:

```r
#' @param sign either 1 or -1, selecting \eqn{C = D + U U^T} or
#'   \eqn{C = D - U U^T}. Defaults to `attr(U, "sign")` when present and to 1
#'   otherwise, so vectors from [am_covariance_structure()] carry the correct
#'   structure automatically.
```

- [ ] **Step 3: Update the roxygen for `am_simulate()`**

Add these parameter entries and extend the return block:

```r
#' @param r cross-mate phenotypic correlation, in the open interval (-1, 1).
#'   Negative values correspond to disassortative mating.
#' @param path (optional) file prefix. If supplied, genotypes are streamed to
#'   disk in batches rather than returned in memory, and `X` is omitted from
#'   the result. Also writes `<path>.meta` and `<path>.rds`.
#' @param format on-disk layout when `path` is supplied. `"individual"` (the
#'   default) stores each individual's variants contiguously in one byte per
#'   genotype, `"variant"` stores each variant's individuals contiguously, and
#'   `"bed"` writes a variant-major PLINK binary file at two bits per genotype
#'   alongside `.bim` and `.fam`.
#' @param batch_size (optional) number of individuals per batch for
#'   `"individual"`, or variants per block for `"variant"` and `"bed"`.
#'   Defaults to a value targeting a working buffer of roughly 128 MB.
#'
#' @return A list. Without `path` it contains `y`, `g`, `X`, `AF`, `beta_std`,
#' `beta_raw`, and `H` when `haplotypes` is TRUE. With `path` it is returned
#' invisibly, omits `X`, and adds `path`, `format`, `n`, and `m`.
#'
#' @details The `"variant"` and `"bed"` layouts stream over loci and reproduce
#' the in-memory genotypes exactly under a given seed at any `batch_size`. The
#' `"individual"` layout streams over people and matches only when
#' `batch_size >= n`, because a batch of rows is not contiguous in R's
#' column-major random draw. `haplotypes = TRUE` cannot be combined with
#' `path`.
```

Also add a streaming example to the existing `@examples` block:

```r
#' ## stream genotypes to disk instead of holding them in memory
#' p <- file.path(tempdir(), "am_sim")
#' meta <- am_simulate(h2_0, r, m, n, path = p, format = "variant")
#' dim(read_genotypes(p))
#'
#' ## disassortative mating
#' neg <- am_simulate(h2_0, -0.5, m, n)
#' var(neg$g) < var(sim_dat$g)
```

- [ ] **Step 4: Regenerate documentation**

Run: `Rscript -e 'roxygen2::roxygenise(".")'`

Expected: `NAMESPACE` gains `export(read_genotypes)` and `export(write_genotypes)`, and `man/` gains `read_genotypes.Rd` and `write_genotypes.Rd`.

- [ ] **Step 5: Update DESCRIPTION**

Set `Version: 1.1.0` and replace the `Description:` field with:

```
Description: Simulation of phenotype / genotype data under assortative and
    disassortative mating. Includes functions for generating Bahadur order-2
    multivariate Bernoulli variables with general and diagonal-plus-low-rank
    correlation structures, and for writing simulated genotypes to disk as
    binary int8 or PLINK bed files. Further details are provided in:
    Border and Malik (2022) <doi:10.1101/2022.10.13.512132>.
```

- [ ] **Step 6: Update NEWS.md**

Insert above the `## version 1.0.0` section:

```markdown
## version 1.1.0

---

- `am_simulate()` and `am_covariance_structure()` now support negative
  (disassortative) cross-mate correlations. The equilibrium covariance becomes
  diagonal minus low rank in that case, tracked by an `attr(U, "sign")` that
  `rb_dplr()` honors automatically.
- fixed `am_covariance_structure()` returning `NaN` at `r = 0`
- `rb_dplr()` gains a `sign` argument, and its infeasibility error now names
  the offending locus and suggests concrete remedies
- documented that negative `r` leaves the Bahadur feasible region sooner than
  positive `r`, reliably sampling to about `r = -0.4` at the default `min_MAF`
- `am_simulate()` gains `path`, `format`, and `batch_size` for streaming
  genotypes to disk in batches, removing the full-matrix allocation
- new `write_genotypes()` and `read_genotypes()` supporting individual-major
  int8, variant-major int8, and PLINK bed
- added a testthat suite
```

- [ ] **Step 7: Update README.md**

In the Features list, under the assortative mating tools, add:

```markdown
  * Genotype input / output
    * `write_genotypes`: write genotypes as int8 or PLINK bed
    * `read_genotypes`: read genotypes back into an R matrix
```

Then add this section after the existing Usage examples:

```markdown
Negative values of `r` correspond to disassortative mating, which reduces
genetic variance rather than inflating it:

```r
neg <- am_simulate(h2_0, r = -.5, m, n)
var(neg$g)
vg_eq(-.5, h2_0, h2_0)
```

For simulations too large to hold in memory, supply `path` to stream genotypes
to disk one batch at a time:

```r
p <- file.path(tempdir(), "am_sim")
meta <- am_simulate(h2_0, r, m = 2e4, n = 5e3, path = p, format = "variant")
X <- read_genotypes(p)
```
```

- [ ] **Step 8: Run the full check**

Run: `Rscript -e 'devtools::check(".", cran = TRUE)'`

Expected: 0 errors, 0 warnings, 0 notes. A note about the installed package size is acceptable. If examples are slow, reduce `m` and `n` in the examples rather than wrapping them in `\donttest{}`.

- [ ] **Step 9: Commit**

```bash
git add R/ man/ NAMESPACE DESCRIPTION NEWS.md README.md
git commit -m "document negative assortment and genotype io, bump to 1.1.0"
```

---

## Verification

After Task 8, confirm from a clean session:

```bash
Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'
Rscript -e 'devtools::check(".", cran = TRUE)'
```

Expected: 37 `test_that` blocks passing and a clean check. Do not claim
completion without pasting the actual output of both commands.

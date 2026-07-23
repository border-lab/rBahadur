.am_check_simulation_args <- function(h2_0, r, n) {
  if (!is.numeric(h2_0) || length(h2_0) != 1L || is.na(h2_0) ||
      !is.finite(h2_0) || h2_0 <= 0 || h2_0 >= 1) {
    stop("`h2_0` must be a finite number in the open interval (0, 1)")
  }
  if (!is.numeric(r) || length(r) != 1L || is.na(r) || !is.finite(r) ||
      r <= -1 || r >= 1) {
    stop("`r` must be a finite number in the open interval (-1, 1)")
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 1 || n != floor(n) || n > .Machine$integer.max) {
    stop("`n` must be a single positive whole number")
  }
  invisible(NULL)
}

.am_warn_small_m <- function(m, caller) {
  if (!identical(getOption("rBahadur.warn_small_m", TRUE), FALSE) &&
      is.numeric(m) && length(m) == 1L && is.finite(m) && m < 50L) {
    warning(
      sprintf(
        paste0(
          "`%s()` is using only %d causal variants. Small `m` is supported, ",
          "but `h2_0` and the equilibrium quantities are large-m targets, so ",
          "realized variances can differ materially below 50 variants. Set ",
          "`options(rBahadur.warn_small_m = FALSE)` to silence this warning."
        ),
        caller, as.integer(m)
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Simulate genotype/phenotype data under equilibrium univariate AM.
#'
#' @param h2_0 generation zero (panmictic) heritability, in `(0, 1)`
#' @param r cross-mate phenotypic correlation, in the open interval (-1, 1).
#'   Negative values correspond to disassortative mating.
#' @param m number of biallelic causal variants; a whole number of at least 2
#' @param n positive whole-number sample size
#' @param afs (optional). Allele frequencies to use. If not provided, `m` will be drawn
#'  uniformly from the interval \[`min_MAF`, 1-`min_MAF`\]
#' @param min_MAF (optional) minimum minor allele frequency for causal variants.
#' Ignored if if `afs` is not NULL. Defaults to 0.1
#' @param haplotypes logical. If TRUE, includes (phased) haploid genotypes in output.
#' Defaults to FALSE
#' @param path (optional) file prefix. If supplied, genotypes are streamed to
#'   disk in batches rather than returned in memory, and `X` is omitted from
#'   the result. Writes `<path>.int8` for `format` `"individual"` or
#'   `"variant"`, or `<path>.bed` plus `<path>.bim` and `<path>.fam` for
#'   `"bed"`, and in every case also writes `<path>.meta` and `<path>.rds`.
#' @param format on-disk layout when `path` is supplied. `"individual"` (the
#'   default) stores each individual's variants contiguously in one byte per
#'   genotype, `"variant"` stores each variant's individuals contiguously, and
#'   `"bed"` writes a variant-major PLINK binary file at two bits per genotype
#'   alongside `.bim` and `.fam`.
#' @param batch_size (optional) number of individuals per batch for
#'   `"individual"`, or variants per block for `"variant"` and `"bed"`.
#'   Defaults to a value targeting a working buffer of roughly 128 MB, but
#'   actual peak memory use is roughly twice that: the `"individual"` branch
#'   also holds the diploid `Xb` and its transposed copy, and the
#'   `"variant"`/`"bed"` branches also copy the buffer passed to the callback
#'   and build a double-precision `Xb` from it.
#'
#' @return A list. Without `path` it contains `y`, `g`, `X`, `AF`, `beta_std`,
#' `beta_raw`, and `H` when `haplotypes` is TRUE. With `path` it is returned
#' invisibly, omits `X`, and adds `path`, `format`, `n`, `m`, `h2_0`, `r`, and
#' `min_MAF`; these are also saved in `<path>.rds` as the call's provenance
#' record.
#'
#' @details The `"variant"` and `"bed"` layouts stream over loci and reproduce
#' the in-memory genotypes exactly under a given seed at any `batch_size`. The
#' `"individual"` layout streams over people and matches only when
#' `batch_size >= n`, because a batch of rows is not contiguous in R's
#' column-major random draw. `haplotypes = TRUE` cannot be combined with
#' `path`.
#'
#' The equilibrium formulas are large-locus results. Calls with fewer than 50
#' causal variants are allowed but warn because realized variances can differ
#' materially from the targets, especially when `m` is very small. Set
#' `options(rBahadur.warn_small_m = FALSE)` to silence this warning after
#' deciding that the finite-locus approximation is appropriate.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' h2_0 = .5; m = 200; n = 1000; r =.5
#'
#' ## simulate genotype/phenotype data
#' sim_dat <- am_simulate(h2_0, r, m, n)
#' str(sim_dat)
#'
#' ## empirical h2 vs expected equilibrium h2
#' (emp_h2 <- var(sim_dat$g)/var(sim_dat$y))
#' h2_eq(r, h2_0)
#'
#' ## stream genotypes to disk instead of holding them in memory
#' p <- file.path(tempdir(), "am_sim")
#' meta <- am_simulate(h2_0, r, m, n, path = p, format = "variant")
#' dim(read_genotypes(p))
#'
#' ## disassortative mating
#' neg <- am_simulate(h2_0, -0.3, m, n)
#' var(neg$g) < var(sim_dat$g)

am_simulate <- function(h2_0, r, m, n, afs = NULL, min_MAF = .1,
                        haplotypes = FALSE, path = NULL,
                        format = c("individual", "variant", "bed"),
                        batch_size = NULL) {
  format <- match.arg(format)
  .am_check_simulation_args(h2_0, r, n)
  n <- as.integer(n)
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m < 2 || m %% 1 != 0 || m > .Machine$integer.max) {
    stop("`m` must be a single whole number of at least 2")
  }
  m <- as.integer(m)
  if (!is.logical(haplotypes) || length(haplotypes) != 1L ||
      is.na(haplotypes)) {
    stop("`haplotypes` must be TRUE or FALSE")
  }
  if (!is.null(path)) .gt_check_path(path, write = TRUE)
  if (is.null(afs)) {
    if (!is.numeric(min_MAF) || length(min_MAF) != 1L || is.na(min_MAF) ||
        !is.finite(min_MAF) || min_MAF < 0 || min_MAF > 0.5) {
      stop("`min_MAF` must be a finite number in the closed interval [0, 0.5]")
    }
  } else if (length(afs) != m) {
    stop("`afs` must have length `m`")
  } else if (!is.numeric(afs) || anyNA(afs) || any(!is.finite(afs)) ||
             any(afs <= 0 | afs >= 1)) {
    stop("`afs` must contain finite allele frequencies strictly between 0 and 1")
  } else {
    AF <- afs
  }
  if (!is.null(batch_size) &&
      (!is.numeric(batch_size) || length(batch_size) != 1L ||
       is.na(batch_size) || !is.finite(batch_size) || batch_size < 1 ||
       batch_size %% 1 != 0 || batch_size > .Machine$integer.max)) {
    stop("`batch_size` must be a single positive whole number, or NULL")
  }
  if (!is.null(path) && haplotypes) {
    stop("`haplotypes = TRUE` is not supported when streaming to `path`")
  }
  .am_warn_small_m(m, "am_simulate")

  ## draw standardized diploid allele substitution effects
  beta <- scale(rnorm(m))*sqrt(h2_0 / m)
  ## draw allele frequencies if necessary
  if (is.null(afs)) AF <- runif(m, min_MAF, 1 - min_MAF)
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
    batch_size <- as.integer(min(batch_size, n))
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
    batch_size <- as.integer(min(batch_size, m))
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
          writeBin(.gt_pack_bed(Xb), con)
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
    m = as.integer(m),
    h2_0 = h2_0,
    r = r,
    min_MAF = min_MAF
  )
  saveRDS(output, paste0(path, ".rds"))
  invisible(output)
}

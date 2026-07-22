#' Simulate genotype/phenotype data under equilibrium univariate AM.
#'
#' @param h2_0 generation zero (panmictic) heritability
#' @param r cross-mate phenotypic correlation, in the open interval (-1, 1).
#'   Negative values correspond to disassortative mating.
#' @param m number of biallelic causal variants
#' @param n sample size
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
#'   actual peak memory use is roughly 3 times that: the `"individual"`
#'   branch also allocates [rb_dplr()]'s equally sized `rand_U` matrix plus
#'   the `Xb`/`t(Xb)` copies, and the `"variant"`/`"bed"` branches also copy
#'   the buffer passed to the callback and build a double-precision `Xb`
#'   from it.
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

  ## validate before opening any connection, so a bad value cannot leave a
  ## partial file behind. A non-integer batch_size would otherwise recycle
  ## silently and write duplicated rows.
  if (!is.null(batch_size)) {
    if (!is.numeric(batch_size) || length(batch_size) != 1L ||
        is.na(batch_size) || batch_size < 1 || batch_size %% 1 != 0) {
      stop("`batch_size` must be a single positive whole number, or NULL")
    }
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
    m = as.integer(m),
    h2_0 = h2_0,
    r = r,
    min_MAF = min_MAF
  )
  saveRDS(output, paste0(path, ".rds"))
  invisible(output)
}


#' Simulate genotype/phenotype data under equilibrium univariate AM.
#'
#' @param h2_0 generation zero (panmictic) heritability
#' @param r cross-mate phenotypic correlation
#' @param m number of biallelic causal variants
#' @param n sample size
#' @param afs (optional). Allele frequencies to use. If not provided, `m` will be drawn
#'  uniformly from the interval \[`min_MAF`, 1-`min_MAF`\]
#' @param min_MAF (optional) minimum minor allele frequency for causal variants. 
#' Ignored if if `afs` is not NULL. Defaults to 0.1
#' @param haplotypes logical. If TRUE, includes (phased) haploid genotypes in output. 
#' Defaults to FALSE
#'
#' @return A list including the following objects:
#' * `y`: phenotype vector
#' * `g`: heritable component of the phenotype vector
#' * `X`: matrix of diploid genotypes
#' * `AF`: vector of allele frequencies
#' * `beta_std`: standardized genetic effects
#' * `beta_raw`: unstandardized genetic effects
#' * `H`: matrix of haploid genotypes (returned only if `haplotypes`=TRUE)
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


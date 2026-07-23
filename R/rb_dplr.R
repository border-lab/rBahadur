## Resolve the low-rank sign for a loading vector, defaulting to +1 so that
## hand-built `U` vectors and pre-1.1.0 code keep working.
.rb_sign <- function(U) {
  s <- attr(U, "sign")
  if (is.null(s)) 1 else s
}

## Build a single, actionable infeasibility message naming the offending
## locus. "Infeasible probabilities" must stay the leading text: an existing
## test in test-am_stream.R matches on it and a later task depends on that
## wording.
.rb_infeasible_msg <- function(locus) {
  paste0(
    "Infeasible probabilities at locus ", locus, ". The specified ",
    "parameters do not correspond to a valid Bahadur order-2 MVB ",
    "distribution. Try reducing the magnitude of r, raising min_MAF, ",
    "or increasing the number of causal variants."
  )
}

.rb_check_bernoulli_inputs <- function(n, mu) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 1 || n != floor(n) || n > .Machine$integer.max) {
    stop("`n` must be a single positive whole number")
  }
  if (!is.numeric(mu) || !length(mu) || anyNA(mu) || any(!is.finite(mu)) ||
      any(mu <= 0 | mu >= 1)) {
    stop("`mu` must contain finite probabilities strictly between 0 and 1")
  }
  as.integer(n)
}

#' Binary random variates with Diagonal Plus Low Rank (dplr) correlations
#'
#' Generate second Bahadur order multivariate Bernoulli random variates with
#' Diagonal Plus Low Rank (dplr) correlation structures.
#'
#' @importFrom stats runif
#'
#' @param n positive whole-number count of observations
#' @param mu non-empty vector of means strictly between 0 and 1
#' @param U finite outer-product component vector, with the same length as `mu`
#' @param sign either 1 or -1, selecting \eqn{C = D + U U^T} or
#'   \eqn{C = D - U U^T}. Defaults to `attr(U, "sign")` when present and to 1
#'   otherwise, so vectors from [am_covariance_structure()] carry the correct
#'   structure automatically.
#'
#' @details This generates multivariate Bernoulli (MVB) random vectors with mean
#' vector 'mu' and correlation matrix \eqn{C = D + U U^T} where \eqn{D} is a diagonal
#'  matrix with values dictated by 'U'. 'mu' must take values in the open unit interval
#'  and 'U' must induce a valid second Bahadur order probability distribution. That is,
#'  there must exist an MVB probability distribution with first moments 'mu' and
#'  standardized central second moments \eqn{C} such that all higher order central
#'  moments are zero.
#'
#' @return An \eqn{n}-by-\eqn{m} matrix of binary random variates, where \eqn{m} is
#' the length of 'mu'.
#'
#' @section Warning, dropping the sign attribute:
#' `U` vectors produced by [am_covariance_structure()] carry
#' `attr(U, "sign")`. Subsetting or coercing `U`, including with
#' \code{as.vector()}, \code{c()}, or `U[i]`, drops that attribute. Because
#' `sign` here falls back to `+1` when the attribute is absent, doing so
#' silently converts a disassortative structure into an assortative one, with
#' no error at any point. If you manipulate `U` before calling `rb_dplr()`,
#' pass `sign = -1` explicitly rather than relying on the attribute to
#' survive.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' h2_0 = .5; m = 200; n = 1000; r =.5; min_MAF=.1
#'
#' ## draw standardized diploid allele substitution effects
#' beta <- scale(rnorm(m))*sqrt(h2_0 / m)
#'
#' ## draw allele frequencies
#' AF <- runif(m, min_MAF, 1 - min_MAF)
#'
#' ## compute unstandardized effects
#' beta_unscaled <- beta/sqrt(2*AF*(1-AF))
#'
#' ## generate corresponding haploid quantities
#' beta_hap <- rep(beta, each=2)
#' AF_hap <- rep(AF, each=2)
#'
#' ## compute equilibrium outer product covariance component
#' U <- am_covariance_structure(beta, AF, r)
#'
#' ## draw multivariate Bernoulli haplotypes
#' H <- rb_dplr(n, AF_hap, U)
#'
#' ## convert to diploid genotypes
#' G <- H[,seq(1,ncol(H),2)] + H[,seq(2,ncol(H),2)]
#'
#' ## empirical allele frequencies vs target frequencies
#' emp_afs <- colMeans(G)/2
#' plot(AF, emp_afs)
#'
#' ## construct phenotype
#' heritable_y <-  G%*%beta_unscaled
#' nonheritable_y <-  rnorm(n, 0, sqrt(1-h2_0))
#' y <- heritable_y + nonheritable_y
#'
#' ## empirical h2 vs expected equilibrium h2
#' (emp_h2 <- var(heritable_y)/var(y))
#' h2_eq(r, h2_0)

rb_dplr <- function(n, mu, U, sign = NULL) {

  n <- .rb_check_bernoulli_inputs(n, mu)
  M <- length(mu)
  if (!is.numeric(U) || length(U) != M || anyNA(U) || any(!is.finite(U))) {
    stop("`U` must be a finite numeric vector with the same length as `mu`")
  }
  if (is.null(sign)) sign <- .rb_sign(U)
  if (!is.numeric(sign) || length(sign) != 1L || is.na(sign) ||
      !is.finite(sign) || !sign %in% c(-1, 1)) {
    stop("`sign` must be either 1 or -1")
  }
  ## bind to a local name so the argument does not shadow base::sign()
  s <- sign

  k <- matrix(NaN, nrow=n, ncol=M)

  ## Uniforms are drawn one locus at a time rather than as an n by M matrix.
  ## R fills matrices column-major, so consuming the stream a column at a time
  ## is bit-identical to matrix(runif(M*n), n, M) while holding n values
  ## instead of n*M; .rb_dplr_stream() rests on the same equivalence.
  ##
  streamed <- M >= 3L
  rand_U <- if (streamed) NULL else matrix(runif(M*n), nrow=n, ncol=M)
  draw <- function(j) if (streamed) runif(n) else rand_U[, j]

  # initial step
  p <- rep(mu[1],n)
  k[ ,1] <- as.numeric(draw(1L) <= p)

  tmp_bool <- (k[ ,1]==0)
  p <- tmp_bool*(1-p) + (!tmp_bool)*p
  Bk0 <- tmp_bool*(1-mu[1]) + (!tmp_bool)*mu[1]
  Bk1 <- tmp_bool*(-1) + (!tmp_bool)*1

  x <- Bk1*U[1]/p
  c <- 1

  if (M == 1L) return(k)

  # recursive steps
  if (M > 2L) {
    for (m in 2:(M-1)) {
      p <- mu[m] + s * x * U[m]
      if (any(!is.finite(p) | p < 0 | p > 1)) {
        stop(.rb_infeasible_msg(m))
      }
      k[ ,m] <- (draw(m) <= p)

      tmp_bool <- (k[ ,m]==0)
      p <- tmp_bool*(1-p) + (!tmp_bool)*p
      Bk0 <- tmp_bool*(1-mu[m]) + (!tmp_bool)*mu[m]
      Bk1 <- tmp_bool*(-1) + (!tmp_bool)*1

      x <- (x*Bk0 + c*Bk1*U[m])/p
      c <- (Bk0/p)*c
    }
  }
  p <- mu[M] + s * x * U[M]
  if (any(!is.finite(p) | p < 0 | p > 1)) {
    stop(.rb_infeasible_msg(M))
  }
  k[ ,M] <-(draw(M) <= p)

  return(k)
}

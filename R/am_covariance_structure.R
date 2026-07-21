#' Compute Diagonal plus Low Rank equilibrium covariance structure
#'
#' @importFrom stats runif rnorm
#'
#' @param beta vector of standardized diploid allele-substitution effects
#' @param AF vector of allele frequencies
#' @param r cross-mate phenotypic correlation
#'
#' @return Vector 'U' such that $D + U U^T$ corresponds to the expected haploid 
#' LD-matrix given the specified genetic architecture (encoded by 'beta' and 'AF') 
#' and cross-mate phenotypic correlation 'r'. It is assumed that the total phenotypic
#' variance at generation zero is one.
#' @examples
#' set.seed(1)
#' h2_0 = .5; m = 200; n = 1000; r =.5; min_MAF=.1
#' betas <- rnorm(m,0,sqrt(h2_0/m))
#' afs <- runif(m, min_MAF, 1-min_MAF)
#' output <- am_covariance_structure(betas, afs, r)
#' @export
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


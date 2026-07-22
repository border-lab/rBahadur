#' Compute Diagonal plus Low Rank equilibrium covariance structure
#'
#' @importFrom stats runif rnorm
#'
#' @param beta vector of standardized diploid allele-substitution effects
#' @param AF vector of allele frequencies
#' @param r cross-mate phenotypic correlation, in the open interval (-1, 1).
#'   Negative values correspond to disassortative mating.
#'
#' @return Vector 'U' such that \eqn{D + s U U^T} corresponds to the expected
#' haploid LD-matrix given the specified genetic architecture (encoded by 'beta'
#' and 'AF') and cross-mate phenotypic correlation 'r', where the sign
#' \eqn{s} is `attr(U, "sign")`. It is assumed that the total phenotypic
#' variance at generation zero is one.
#'
#' @details For \code{r > 0} the low-rank term is added and `attr(U, "sign")` is 1.
#' For \code{r < 0} it is subtracted and the attribute is -1, because a positive
#' semidefinite rank-one term can only increase genetic variance whereas
#' disassortative mating reduces it. The returned vector in that case is the
#' analytic continuation of the positive branch, which is purely imaginary.
#' For \code{r = 0} the vector is zero, since panmixia induces no disequilibrium.
#'
#' @section Warning, dropping the sign attribute:
#' The returned vector carries `attr(U, "sign")`, which records whether the
#' rank-one term is added or subtracted. Subsetting or coercing the vector,
#' including with \code{as.vector()}, \code{c()}, or `U[i]`, drops that
#' attribute. Because [rb_dplr()] falls back to `sign = 1` when the attribute
#' is absent, this silently turns a disassortative structure into an
#' assortative one, with no error at any point. If you manipulate `U` before
#' passing it to [rb_dplr()], pass `sign = -1` explicitly for \code{r < 0}
#' rather than relying on the attribute to survive.
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
#' \code{r = -0.3} sampled reliably at 4000 individuals, while \code{r = -0.4} sampled
#' reliably at 2000 individuals but failed for some seeds at 4000. Raising
#' `min_MAF` widens the envelope. If [rb_dplr()] reports infeasible
#' probabilities, reduce the magnitude of `r`, raise `min_MAF`, or increase
#' the number of causal variants.
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


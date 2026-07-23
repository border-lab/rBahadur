#' Functions to compute equilibrium parameters under assortative mating
#'
#' Compute heritability ('h2_eq'), genetic variance ('vg_0'), and cross-mate genetic 
#' correlation ('rg_eq') at equilibrium under univariate primary-phenotypic assortative mating. 
#' These equations can be derived from Nagylaki's results (see below) under the assumption 
#' that number of causal variants is large (i.e., taking the limit as the number of causal 
#' variants approaches infinity).
#' @references Nagylaki, T. Assortative mating for a quantitative character. J. Math. Biology 
#' 16, 57–74 (1982). https://doi.org/10.1007/BF00275161
#'
#' @param r cross-mate phenotypic correlation, with values in `(-1, 1)`
#' @param h2_0 generation zero (panmictic) heritability, with values in
#'   `[0, 1)`
#' @param vg_0 non-negative generation zero (panmictic) additive genetic
#'   variance component
#' @return A single numerical quantity representing the equilibrium heritability (`h2_eq`),
#'  the equilibrium cross-mate genetic correlation (`rg_eq`), or the equilibrium genetic 
#'  variance (`vg_eq`).
#' @name am_equilibrium_parameters
NULL
#> NULL

.am_check_equilibrium_args <- function(r, h2_0) {
  if (!is.numeric(r) || !length(r) || anyNA(r) || any(!is.finite(r)) ||
      any(r <= -1 | r >= 1)) {
    stop("`r` must contain finite numbers in the open interval (-1, 1)")
  }
  if (!is.numeric(h2_0) || !length(h2_0) || anyNA(h2_0) ||
      any(!is.finite(h2_0)) || any(h2_0 < 0 | h2_0 >= 1)) {
    stop("`h2_0` must contain finite numbers in the interval [0, 1)")
  }
  invisible(NULL)
}
#' @examples
#' set.seed(1)
#' vg_0= .6; h2_0 = .5; r =.5
#' h2_eq(r, h2_0)
#' rg_eq(r, h2_0)
#' vg_eq(r, vg_0, h2_0)
#' @export
#' @rdname am_equilibrium_parameters
#' @export
h2_eq <- function(r, h2_0){
  .am_check_equilibrium_args(r, h2_0)
  ## Equilibrium heritability is equilibrium genetic variance over equilibrium
  ## total variance, and assortative mating leaves the environmental component
  ## at its generation zero value of 1 - h2_0. Heritability is scale free, so
  ## taking var(y) = 1 at generation zero, and hence vg_0 = h2_0, costs no
  ## generality.
  ##
  ## Written out directly this is 1/(2r) times a bracket that also vanishes at
  ## r = 0, which is 0/0 there: panmixia would report NaN rather than h2_0,
  ## and small r was ill conditioned by the same cancellation. Going through
  ## rg_eq() removes both, since nothing divides by r.
  vg <- vg_eq(r, h2_0, h2_0)
  vg / (vg + 1 - h2_0)
}

#' @rdname am_equilibrium_parameters
#' @export
rg_eq <- function(r, h2_0) {
  .am_check_equilibrium_args(r, h2_0)
  ## tmp <- 1/(1-h2_0)
  ## tmp*(tmp-sqrt(tmp^2-4*r*h2_0*tmp))/2
  ((1-h2_0)^-1-sqrt(1/(1-h2_0)^2-4*r*h2_0/(1-h2_0)))/2
}

#' @rdname am_equilibrium_parameters
#' @export
vg_eq <- function(r, vg_0, h2_0)  {
  .am_check_equilibrium_args(r, h2_0)
  if (!is.numeric(vg_0) || !length(vg_0) || anyNA(vg_0) ||
      any(!is.finite(vg_0)) || any(vg_0 < 0)) {
    stop("`vg_0` must contain finite non-negative numbers")
  }
  vg_0/(1-rg_eq(r, h2_0))
}

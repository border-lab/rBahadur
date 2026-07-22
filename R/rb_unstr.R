#' Binary random variates with unstructured correlations
#'
#' Generate Bahadur order-2 multivariate Bernoulli random variates with 
#' unstructured correlations.
#'
#' @importFrom stats runif
#'
#' @param n number of observations
#' @param mu vector of means
#' @param C correlation matrix
#'
#' @details This generates multivariate Bernoulli (MVB) random vectors with mean vector
#'  'mu' and correlation matrix 'C'. 'mu' must take values in the open unit interval and
#'   'C' must induce a valid second Bahadur order probability distribution. That is, 
#'   there must exist an MVB probability distribution with first moments 'mu' and 
#'   standardized central second moments 'C' such that all higher order central moments 
#'   are zero.
#'
#' @return An \eqn{n}-by-\eqn{m} matrix of binary random variates, where \eqn{m} is the 
#' length of 'mu'.
#' @export
#'
#' @examples
#' set.seed(1)
#' h2_0 = .5; m = 60; n = 200; r =.5; min_MAF=.1
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
#' ## construct Correlation matrix
#' S <- diag(1/sqrt(AF_hap*(1-AF_hap)))
#' DPLR <- U%o%U
#' diag(DPLR) <- 1
#' C <- cov2cor(S%*%DPLR%*%S)
#'
#' ## draw multivariate Bernoulli haplotypes
#' H <- rb_unstr(n, AF_hap, C)
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

rb_unstr <- function(n, mu, C) {

  M <- length(mu)
  k <- matrix(NaN, nrow=n, ncol=M)
  rand_U <- matrix(runif(M*n), nrow=n, ncol=M)

  # initial step
  p <- rep(mu[1],n)
  k[ ,1] <- as.numeric(rand_U[,1] <= p)


  Bk0_values_matrix <- rbind(1-mu,mu)
  Bk1_bvec <- c(-1,1)

  tmp_bool <- (k[ ,1]==0)
  p <- tmp_bool*(1-p) + (!tmp_bool)*p

  Bk0_mat <- Bk1_mat <- matrix(as.numeric(NA),
                               nrow=M,
                               ncol=n)
  Bk0 <- Bk0_mat[1,] <- Bk0_values_matrix[cbind(k[ ,1]+1,
                                                rep(1,n))]

  Bk1 <- Bk1_mat[1,] <- Bk1_bvec[k[ ,1]+1]

  v <- matrix(Bk1/p,ncol=1)
  c <- 1

  # recursive steps
  for (m in 2:(M-1)) {
    p <- drop(mu[m] + v%*%C[1:(m-1),m])
    if (any(p < 0 | p > 1)) {
      stop('Infeasible probabilities. Are you sure specified parameters 
           correspond to a valid Bahadur order-2 MVB distribution?')
    }
    k[ ,m] <- (rand_U[ ,m] <= p)

    tmp_bool <- (k[ ,m]==0)
    p <- tmp_bool*(1-p) + (!tmp_bool)*p
    Bk0_mat[m,] <- Bk0 <- Bk0_values_matrix[cbind(k[ ,m]+1,
                                                  rep(m,n))]
    Bk1_mat[m,] <- Bk1 <-  Bk1_bvec[k[ ,m]+1]
    v <- cbind(v*Bk0/p, + c*Bk1/p)
    c <- (Bk0/p)*c
  }
  p <- drop(mu[M] + v%*%C[1:(M-1),M])
  if (any(p < 0 | p > 1)) {
    stop(mu[M],' ',p,'Infeasible probabilities. Are you sure specified parameters 
         correspond to a valid Bahadur order-2 MVB distribution?')
  }
  k[ ,M] <-(rand_U[ ,M] <= p)

  return(k)
}

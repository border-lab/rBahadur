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

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

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

  ## Also test with U carrying the sign attribute from am_covariance_structure
  U_with_attr <- am_covariance_structure(beta, AF, 0.4)
  set.seed(8); no_sign_arg <- rb_dplr(30, mu, U_with_attr)
  set.seed(8); explicit_sign <- rb_dplr(30, mu, U_with_attr, sign = 1)
  expect_identical(no_sign_arg, explicit_sign)
})

test_that("an invalid sign is rejected", {
  expect_error(rb_dplr(5, rep(0.5, 6), rep(0.01, 6), sign = 0), "either 1 or -1")
})

test_that("one- and two-variable Bernoulli draws use the correct recursion", {
  set.seed(101)
  expected_one <- matrix(as.numeric(runif(12) <= 0.3), ncol = 1)
  set.seed(101)
  expect_identical(rb_dplr(12, 0.3, 0), expected_one)

  set.seed(102)
  uniforms <- matrix(runif(24), nrow = 12, ncol = 2)
  expected_two <- uniforms <= matrix(c(0.3, 0.7), 12, 2, byrow = TRUE)
  storage.mode(expected_two) <- "double"
  set.seed(102)
  expect_identical(rb_dplr(12, c(0.3, 0.7), c(0, 0)), expected_two)
})

test_that("rb_dplr validates dimensions and probabilities", {
  expect_error(rb_dplr(0, 0.5, 0), "`n`")
  expect_error(rb_dplr(2.5, 0.5, 0), "`n`")
  expect_error(rb_dplr(4, numeric(), numeric()), "`mu`")
  expect_error(rb_dplr(4, c(0, 0.5), c(0, 0)), "`mu`")
  expect_error(rb_dplr(4, c(0.2, 0.8), 0), "same length")
  expect_error(rb_dplr(4, c(0.2, 0.8), c(0, NA)), "finite numeric")
})

test_that("rb_unstr handles small dimensions and validates correlations", {
  set.seed(103)
  one <- rb_unstr(10, 0.4, matrix(1, 1, 1))
  expect_identical(dim(one), c(10L, 1L))
  expect_true(all(one %in% 0:1))

  set.seed(104)
  two <- rb_unstr(10, c(0.4, 0.6), diag(2))
  expect_identical(dim(two), c(10L, 2L))
  expect_true(all(two %in% 0:1))

  expect_error(rb_unstr(0, 0.5, matrix(1, 1, 1)), "`n`")
  expect_error(rb_unstr(4, c(0.5, 1), diag(2)), "`mu`")
  expect_error(rb_unstr(4, c(0.4, 0.6), matrix(1, 1, 1)), "matching `mu`")
  expect_error(rb_unstr(4, c(0.4, 0.6), matrix(c(1, 0.2, 0, 1), 2)),
               "symmetric")
})

test_that("infeasible probabilities trigger an actionable error message", {
  ## Use deterministic inputs that drive the recursion out of range quickly.
  ## The real invariant we want to guard is that the message is O(1) in the
  ## number of loci, not O(m): a reintroduced full-vector dump would make the
  ## message grow with the length of mu/U. An absolute length cap (e.g. under
  ## 500 chars) cannot detect that at a small m, since a dumped 20-element
  ## vector is already comfortably under such a cap. Instead, trigger the
  ## same error at two very different sizes and confirm message length does
  ## not grow with m.
  get_infeasible_message <- function(m) {
    mu <- rep(0.5, m)
    U <- rep(0.9, m)

    ## Guard the capture: if rb_dplr() unexpectedly succeeded, calling
    ## conditionMessage() on a non-error result would itself throw an
    ## obscure error instead of reporting a clear test failure.
    result <- tryCatch(
      rb_dplr(10, mu, U),
      error = function(e) e
    )
    if (!inherits(result, "error")) {
      fail(sprintf(
        "rb_dplr() unexpectedly succeeded for m = %d instead of raising the infeasibility error",
        m
      ))
      return(NA_character_)
    }
    conditionMessage(result)
  }

  msg_small <- get_infeasible_message(20)
  msg_large <- get_infeasible_message(400)

  ## Check both messages begin with the documented prefix
  expect_match(msg_small, "^Infeasible probabilities at locus")
  expect_match(msg_large, "^Infeasible probabilities at locus")

  ## Verify message length is essentially constant across a 20x difference in
  ## the number of loci: a dumped vector would make the m = 400 message
  ## dramatically longer than the m = 20 message. Allow a small allowance for
  ## extra locus-index digits.
  expect_lt(abs(nchar(msg_large) - nchar(msg_small)), 10)
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

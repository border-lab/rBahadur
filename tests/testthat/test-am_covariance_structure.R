make_arch <- function(m = 400, h2_0 = 0.5, seed = 1) {
  set.seed(seed)
  list(beta = as.vector(scale(rnorm(m))) * sqrt(h2_0 / m),
       AF = runif(m, 0.1, 0.9),
       h2_0 = h2_0)
}

test_that("positive r is unchanged and tagged with sign 1", {
  a <- make_arch()
  U <- am_covariance_structure(a$beta, a$AF, 0.5)
  expect_true(all(is.finite(U)))
  expect_identical(attr(U, "sign"), 1)
  expect_length(U, 2 * length(a$beta))
})

test_that("r = 0 returns zeros rather than NaN", {
  a <- make_arch()
  U <- am_covariance_structure(a$beta, a$AF, 0)
  expect_true(all(U == 0))
  expect_identical(attr(U, "sign"), 1)
})

test_that("negative r matches the imaginary part of the complex continuation", {
  a <- make_arch()
  for (r in c(-0.8, -0.5, -0.2)) {
    b <- rep(a$beta, each = 2)
    sdh <- rep(sqrt(a$AF * (1 - a$AF)), each = 2)
    h20 <- sum(a$beta^2)
    rg <- rg_eq(r, h20)
    vtot <- vg_eq(r, h20, h20) + (1 - h20)
    rc <- as.complex(r)
    Uc <- sqrt(vtot / 2) / (2 * b * sqrt(rc)) *
      (sqrt(4 * b^2 * rc / vtot + (1 - rg)^2) - (1 - rg)) * sdh

    U <- am_covariance_structure(a$beta, a$AF, r)
    expect_equal(max(abs(Re(Uc))), 0)
    expect_equal(as.vector(U), Im(Uc), tolerance = 1e-12)
    expect_identical(attr(U, "sign"), -1)
  }
})

test_that("negative r reproduces the equilibrium variance deficit", {
  a <- make_arch(m = 4000)
  r <- -0.5
  h20 <- sum(a$beta^2)
  a_hap <- rep(a$beta / sqrt(2 * a$AF * (1 - a$AF)), each = 2)
  V <- am_covariance_structure(a$beta, a$AF, r)
  expect_equal(-sum(a_hap * V)^2, vg_eq(r, h20, h20) - h20, tolerance = 1e-3)
})

test_that("r outside (-1, 1) is rejected", {
  a <- make_arch()
  expect_error(am_covariance_structure(a$beta, a$AF, 1), "open interval")
  expect_error(am_covariance_structure(a$beta, a$AF, -1), "open interval")
})

test_that("architecture inputs are validated and zero effects stay finite", {
  expect_error(am_covariance_structure(c(0.1, NA), c(0.2, 0.8), 0.2),
               "beta")
  expect_error(am_covariance_structure(c(0.1, 0.2), 0.5, 0.2), "AF")
  expect_error(am_covariance_structure(c(0.1, 0.2), c(0, 0.5), 0.2), "AF")
  expect_error(am_covariance_structure(c(1, 1), c(0.2, 0.8), 0.2), "h2_0")

  U <- am_covariance_structure(c(0, 0.3), c(0.2, 0.8), 0.4)
  expect_true(all(is.finite(U)))
  expect_identical(U[1:2], c(0, 0))
})

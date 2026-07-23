## The previous closed form for h2_eq, kept here as the reference. It is
## correct wherever it is finite, so it pins the rewrite everywhere except at
## the singularity that was the point of the rewrite.
previous_h2_eq <- function(r, h2_0) {
  1 / (2 * r) * (1 / (1 - h2_0) -
                   sqrt((1 - h2_0)^-2 - 4 * r * h2_0 / (1 - h2_0)))
}

test_that("h2_eq agrees with the closed form it replaced", {
  grid <- expand.grid(r = c(-0.9, -0.5, -0.3, -0.05, 0.05, 0.3, 0.5, 0.9),
                      h2_0 = c(0.05, 0.2, 0.5, 0.8, 0.95))
  expect_equal(mapply(h2_eq, grid$r, grid$h2_0),
               mapply(previous_h2_eq, grid$r, grid$h2_0))
})

test_that("h2_eq returns the panmictic heritability at r = 0", {
  ## the old form was 1/(2r) times a bracket that also vanishes at r = 0, so
  ## panmixia reported NaN instead of the one value that needs no derivation
  for (h2_0 in c(0.05, 0.2, 0.5, 0.8, 0.95)) {
    expect_equal(h2_eq(0, h2_0), h2_0)
    expect_false(is.na(h2_eq(0, h2_0)))
  }
})

test_that("h2_eq approaches r = 0 continuously from both sides", {
  for (h2_0 in c(0.2, 0.5, 0.8)) {
    for (eps in c(1e-4, 1e-8, 1e-12)) {
      expect_equal(h2_eq(eps, h2_0), h2_0, tolerance = 1e-3)
      expect_equal(h2_eq(-eps, h2_0), h2_0, tolerance = 1e-3)
    }
  }
})

test_that("h2_eq is vectorised over r, including through zero", {
  r <- c(-0.5, -1e-9, 0, 1e-9, 0.5)
  got <- h2_eq(r, 0.5)
  expect_length(got, length(r))
  expect_false(anyNA(got))
  expect_equal(got, vapply(r, h2_eq, numeric(1), 0.5))
  ## and it is monotone in r: more assortment, more heritability
  expect_false(is.unsorted(got))
})

test_that("equilibrium heritability is genetic over total variance", {
  ## the identity the definition now rests on, checked independently
  for (r in c(-0.5, -0.1, 0, 0.3, 0.7)) {
    for (h2_0 in c(0.2, 0.5, 0.8)) {
      vg <- vg_eq(r, h2_0, h2_0)
      expect_equal(h2_eq(r, h2_0), vg / (vg + 1 - h2_0))
    }
  }
})

test_that("rg_eq and vg_eq are already well behaved at r = 0", {
  ## no assortment means no cross-mate genetic correlation and no inflation
  for (h2_0 in c(0.2, 0.5, 0.8)) {
    expect_equal(rg_eq(0, h2_0), 0)
    expect_equal(vg_eq(0, h2_0, h2_0), h2_0)
  }
})

test_that("equilibrium helpers reject values outside their mathematical domain", {
  expect_error(h2_eq(1, 0.5), "open interval")
  expect_error(h2_eq(0.2, 1), "\\[0, 1\\)")
  expect_error(rg_eq(NA_real_, 0.5), "finite")
  expect_error(vg_eq(0.2, -0.1, 0.5), "non-negative")
})

test_that("panmictic simulation matches h2_eq at r = 0", {
  skip_on_cran()
  set.seed(11)
  sim <- am_simulate(0.5, 0, 400, 4000)
  emp <- var(as.vector(sim$g)) / var(as.vector(sim$y))
  expect_equal(emp, h2_eq(0, 0.5), tolerance = 0.05)
})

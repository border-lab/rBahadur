test_that("am_simulate runs for negative, zero, and positive r", {
  set.seed(21)
  for (r in c(-0.5, 0, 0.5)) {
    d <- am_simulate(h2_0 = 0.5, r = r, m = 200, n = 300)
    expect_true(all(is.finite(d$X)))
    expect_true(all(d$X %in% c(0, 1, 2)))
    expect_equal(dim(d$X), c(300L, 200L))
  }
})

test_that("empirical heritability tracks h2_eq across the sign range", {
  skip_on_cran()
  set.seed(22)
  h2_0 <- 0.5
  ## negative r is held within the reliable Bahadur feasible range; see the
  ## Global Constraints note on the negative-assortment envelope
  for (r in c(-0.3, -0.2, 0.3, 0.6)) {
    d <- am_simulate(h2_0 = h2_0, r = r, m = 1500, n = 4000)
    emp <- var(as.vector(d$g)) / var(as.vector(d$y))
    expect_equal(emp, h2_eq(r, h2_0), tolerance = 0.05)
  }
})

test_that("negative r reduces genetic variance and positive r inflates it", {
  skip_on_cran()
  set.seed(23)
  h2_0 <- 0.5
  vneg <- var(as.vector(am_simulate(h2_0, -0.3, 1500, 4000)$g))
  vpos <- var(as.vector(am_simulate(h2_0, 0.6, 1500, 4000)$g))
  expect_lt(vneg, h2_0)
  expect_gt(vpos, h2_0)
  expect_equal(vneg, vg_eq(-0.3, h2_0, h2_0), tolerance = 0.05)
  expect_equal(vpos, vg_eq(0.6, h2_0, h2_0), tolerance = 0.05)
})

test_that("allele frequencies are preserved under negative r", {
  skip_on_cran()
  set.seed(24)
  d <- am_simulate(0.5, -0.3, 1500, 4000)
  expect_gt(cor(d$AF, colMeans(d$X) / 2), 0.99)
})

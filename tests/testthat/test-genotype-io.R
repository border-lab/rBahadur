test_that("individual and variant layouts round trip", {
  set.seed(31)
  X <- matrix(sample(0:2, 7 * 5, replace = TRUE), nrow = 7, ncol = 5)
  for (fmt in c("individual", "variant")) {
    p <- file.path(tempdir(), paste0("gt-", fmt))
    write_genotypes(X, p, format = fmt)
    expect_true(file.exists(paste0(p, ".int8")))
    expect_equal(read_genotypes(p), X)
  }
})

test_that("on-disk byte order matches the declared layout", {
  X <- matrix(c(0L, 1L, 2L,
                2L, 1L, 0L), nrow = 2, ncol = 3, byrow = TRUE)

  pv <- file.path(tempdir(), "gt-order-variant")
  write_genotypes(X, pv, format = "variant")
  expect_equal(readBin(paste0(pv, ".int8"), "integer", 6L, size = 1L),
               c(0L, 2L, 1L, 1L, 2L, 0L))

  pi <- file.path(tempdir(), "gt-order-individual")
  write_genotypes(X, pi, format = "individual")
  expect_equal(readBin(paste0(pi, ".int8"), "integer", 6L, size = 1L),
               c(0L, 1L, 2L, 2L, 1L, 0L))
})

test_that("file size is exactly one byte per genotype", {
  X <- matrix(1L, nrow = 11, ncol = 13)
  p <- file.path(tempdir(), "gt-size")
  write_genotypes(X, p, format = "variant")
  expect_equal(file.size(paste0(p, ".int8")), 11 * 13)
})

test_that("metadata round trips", {
  X <- matrix(0L, nrow = 4, ncol = 6)
  p <- file.path(tempdir(), "gt-meta")
  write_genotypes(X, p, format = "individual")
  meta <- rBahadur:::.gt_read_meta(p)
  expect_identical(meta$format, "individual")
  expect_identical(meta$n, 4L)
  expect_identical(meta$m, 6L)
  expect_identical(meta$dtype, "int8")
})

test_that("invalid genotypes and missing values are rejected", {
  p <- file.path(tempdir(), "gt-bad")
  expect_error(write_genotypes(matrix(3L, 2, 2), p, "variant"), "0, 1, or 2")
  expect_error(write_genotypes(matrix(NA_integer_, 2, 2), p, "variant"), "bed")
  expect_error(write_genotypes(1:4, p, "variant"), "matrix")
})

test_that("reading without metadata fails clearly", {
  expect_error(read_genotypes(file.path(tempdir(), "gt-does-not-exist")),
               "metadata file not found")
})

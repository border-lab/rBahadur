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

test_that("corrupt metadata is rejected before its layout is interpreted", {
  p <- file.path(tempdir(), "gt-corrupt-meta")
  X <- matrix(0L, nrow = 4, ncol = 6)
  original <- c(
    "rBahadur_genotypes: 1", "format: individual", "n: 4", "m: 6",
    "dtype: int8"
  )
  cases <- list(
    c("rBahadur_genotypes: 2", original[-1]),
    sub("format: individual", "format: varaint", original, fixed = TRUE),
    sub("n: 4", "n: 0", original, fixed = TRUE),
    sub("m: 6", "m: 6.5", original, fixed = TRUE),
    sub("dtype: int8", "dtype: bed2bit", original, fixed = TRUE),
    original[-4],
    c(original, "m: 6"),
    c(original, "fromat: individual")
  )
  patterns <- c("version", "format", "positive whole", "positive whole",
                "incompatible", "missing", "duplicate", "unknown")

  for (i in seq_along(cases)) {
    write_genotypes(X, p, format = "individual")
    writeLines(cases[[i]], paste0(p, ".meta"))
    expect_error(read_genotypes(p), patterns[i], info = paste("case", i))
  }
})

test_that("invalid genotypes and missing values are rejected", {
  p <- file.path(tempdir(), "gt-bad")
  expect_error(write_genotypes(matrix(3L, 2, 2), p, "variant"), "0, 1, or 2")
  expect_error(write_genotypes(matrix(0.5, 2, 2), p, "variant"),
               "integer-valued")
  expect_error(write_genotypes(matrix(0.5, 2, 2), p, "bed"),
               "integer-valued")
  expect_error(write_genotypes(matrix(Inf, 2, 2), p, "variant"),
               "integer-valued")
  expect_error(write_genotypes(matrix("1", 2, 2), p, "variant"),
               "integer-valued")
  expect_error(write_genotypes(matrix(NA_integer_, 2, 2), p, "variant"), "bed")
  expect_error(write_genotypes(1:4, p, "variant"), "matrix")
  expect_error(write_genotypes(matrix(integer(), 0, 2), p, "variant"),
               "at least one individual")
  expect_error(write_genotypes(matrix(integer(), 2, 0), p, "variant"),
               "at least one individual")
  expect_error(write_genotypes(matrix(0L, 2, 2), character(), "variant"),
               "single non-empty")
})

test_that("PLINK sidecars reject invalid physical positions and marker text", {
  p <- file.path(tempdir(), "gt-bad-sidecar")
  expect_error(
    rBahadur:::.gt_write_plink_sidecars(p, 2, 2, pos = c(1, 2.5)),
    "positive whole numbers"
  )
  expect_error(
    rBahadur:::.gt_write_plink_sidecars(
      p, 2, 2, pos = c(1, .Machine$integer.max)
    ),
    "2\\^31 - 2"
  )
  expect_error(
    rBahadur:::.gt_write_plink_sidecars(p, 2, 2, chrom = "chr 1"),
    "no whitespace"
  )
})

test_that("reading without metadata fails clearly", {
  expect_error(read_genotypes(file.path(tempdir(), "gt-does-not-exist")),
               "metadata file not found")
})

test_that("a truncated int8 file is rejected instead of silently recycled", {
  ## 8 x 4 -> 32 bytes; truncating to 16 bytes (half) divides n*m exactly, so
  ## matrix(v, nrow = n, ncol = m) would otherwise recycle silently and
  ## produce a full 8x4 matrix whose columns 3-4 are copies of columns 1-2.
  set.seed(61)
  X <- matrix(sample(0:2, 8 * 4, replace = TRUE), nrow = 8, ncol = 4)
  p <- file.path(tempdir(), "gt-truncated-variant")
  write_genotypes(X, p, format = "variant")
  f <- paste0(p, ".int8")
  expect_equal(file.size(f), 32)

  ## sanity: the valid file still reads correctly before we corrupt it
  expect_equal(read_genotypes(p), X)

  raw <- readBin(f, "raw", n = 32)
  writeBin(raw[1:16], f)
  expect_equal(file.size(f), 16)
  expect_error(read_genotypes(p), "truncated|expected")

  ## a file that is too long should also be rejected
  writeBin(c(raw, as.raw(c(0L, 1L))), f)
  expect_error(read_genotypes(p), "truncated|expected")
})

test_that("a truncated individual-layout int8 file is rejected", {
  set.seed(62)
  X <- matrix(sample(0:2, 8 * 4, replace = TRUE), nrow = 8, ncol = 4)
  p <- file.path(tempdir(), "gt-truncated-individual")
  write_genotypes(X, p, format = "individual")
  f <- paste0(p, ".int8")

  expect_equal(read_genotypes(p), X)

  raw <- readBin(f, "raw", n = 32)
  writeBin(raw[1:16], f)
  expect_error(read_genotypes(p), "truncated|expected")

  writeBin(c(raw, as.raw(0L)), f)
  expect_error(read_genotypes(p), "truncated|expected")
})

test_that("invalid int8 dosage bytes are detected", {
  X <- matrix(0L, nrow = 3, ncol = 4)
  p <- file.path(tempdir(), "gt-invalid-byte")
  write_genotypes(X, p, format = "variant")
  f <- paste0(p, ".int8")
  bytes <- readBin(f, "raw", n = 12)
  bytes[5] <- as.raw(255)
  writeBin(bytes, f)
  expect_error(read_genotypes(p), "other than 0, 1, or 2")
})

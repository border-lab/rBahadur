test_that("bed packing matches a hand-computed byte", {
  ## dosages 2,1,0,2 -> codes 0,2,3,0 -> 0*1 + 2*4 + 3*16 + 0*64 = 56
  expect_equal(as.integer(rBahadur:::.gt_pack_bed(c(2L, 1L, 0L, 2L))), 56L)
})

test_that("bed packing round trips including partial final bytes", {
  set.seed(41)
  for (n in c(1, 3, 4, 5, 17, 100)) {
    g <- sample(0:2, n, replace = TRUE)
    packed <- rBahadur:::.gt_pack_bed(as.integer(g))
    expect_equal(length(packed), ceiling(n / 4))
    expect_equal(rBahadur:::.gt_unpack_bed(packed, n), as.integer(g))
  }
})

test_that("missing genotypes round trip through bed", {
  g <- c(0L, NA_integer_, 2L, 1L, NA_integer_)
  expect_equal(rBahadur:::.gt_unpack_bed(rBahadur:::.gt_pack_bed(g), 5L), g)
})

test_that("bed files round trip and carry a valid header", {
  set.seed(42)
  X <- matrix(sample(0:2, 9 * 6, replace = TRUE), nrow = 9, ncol = 6)
  p <- file.path(tempdir(), "gt-bed")
  write_genotypes(X, p, format = "bed")

  hdr <- readBin(paste0(p, ".bed"), "raw", 3L)
  expect_equal(as.integer(hdr), c(0x6c, 0x1b, 0x01))
  expect_equal(file.size(paste0(p, ".bed")), 3 + 6 * ceiling(9 / 4))
  expect_equal(read_genotypes(p), X)
})

test_that("bed writes plink sidecars with correct dimensions", {
  X <- matrix(0L, nrow = 5, ncol = 4)
  p <- file.path(tempdir(), "gt-bed-sidecar")
  write_genotypes(X, p, format = "bed")
  expect_length(readLines(paste0(p, ".bim")), 4L)
  expect_length(readLines(paste0(p, ".fam")), 5L)
  expect_equal(length(strsplit(readLines(paste0(p, ".bim"))[1], "\t")[[1]]), 6L)
  expect_equal(length(strsplit(readLines(paste0(p, ".fam"))[1], "\t")[[1]]), 6L)
})

test_that("a bed file with a bad header is rejected", {
  p <- file.path(tempdir(), "gt-bed-bad")
  write_genotypes(matrix(0L, 4, 2), p, format = "bed")
  con <- file(paste0(p, ".bed"), "r+b")
  writeBin(as.raw(c(0x00, 0x00, 0x00)), con)
  close(con)
  expect_error(read_genotypes(p), "variant-major PLINK")
})

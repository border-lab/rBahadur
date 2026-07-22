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

test_that("packing a block of variants matches packing them one at a time", {
  ## whole blocks are packed in one call, so the per-variant padding has to
  ## stay aligned to variant boundaries; sample sizes on either side of a
  ## multiple of four are what would expose a drift
  set.seed(43)
  for (n in c(4L, 5L, 7L, 8L, 41L)) {
    X <- matrix(sample(0:2, n * 6L, replace = TRUE), nrow = n)
    X[1L, 2L] <- NA_integer_
    one_at_a_time <- unlist(lapply(seq_len(ncol(X)),
                                   function(j) rBahadur:::.gt_pack_bed(X[, j])))
    expect_identical(rBahadur:::.gt_pack_bed(X), one_at_a_time)
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

test_that("a truncated bed file is rejected instead of silently reading as reference", {
  ## A short per-variant readBin() returns an empty raw vector, and
  ## .gt_unpack_bed() would otherwise fill that column with 0 (homozygous
  ## reference) for every individual instead of erroring, with anyNA() FALSE.
  set.seed(43)
  n <- 8; m <- 4
  X <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
  p <- file.path(tempdir(), "gt-bed-truncated")
  write_genotypes(X, p, format = "bed")
  f <- paste0(p, ".bed")
  expected_size <- 3 + m * ceiling(n / 4)
  expect_equal(file.size(f), expected_size)

  ## sanity: the valid file still reads correctly before we corrupt it
  expect_equal(read_genotypes(p), X)

  raw <- readBin(f, "raw", n = expected_size)

  ## drop the last variant's bytes entirely
  writeBin(raw[seq_len(expected_size - ceiling(n / 4))], f)
  expect_error(read_genotypes(p), "truncated|expected")

  ## a file that is too long should also be rejected
  writeBin(c(raw, as.raw(0L)), f)
  expect_error(read_genotypes(p), "truncated|expected")
})

inmem <- function(seed, r = 0.4, m = 40, n = 24) {
  set.seed(seed)
  am_simulate(h2_0 = 0.5, r = r, m = m, n = n)
}
streamed <- function(seed, fmt, batch, r = 0.4, m = 40, n = 24) {
  p <- file.path(tempdir(), paste0("sim-", fmt, "-", batch, "-", seed))
  set.seed(seed)
  out <- am_simulate(h2_0 = 0.5, r = r, m = m, n = n,
                     path = p, format = fmt, batch_size = batch)
  list(out = out, X = read_genotypes(p), path = p)
}

test_that("locus-batched layouts reproduce the in-memory genotypes at any batch size", {
  ref <- inmem(51)
  for (fmt in c("variant", "bed")) {
    for (b in c(1L, 3L, 40L, 400L)) {
      s <- streamed(51, fmt, b)
      expect_equal(s$X, matrix(as.integer(ref$X), nrow = 24),
                   info = paste(fmt, b))
      expect_equal(as.vector(s$out$y), as.vector(ref$y), info = paste(fmt, b))
    }
  }
})

test_that("individual layout reproduces the in-memory result when batch_size >= n", {
  ref <- inmem(52)
  s <- streamed(52, "individual", 24L)
  expect_equal(s$X, matrix(as.integer(ref$X), nrow = 24))
  expect_equal(as.vector(s$out$y), as.vector(ref$y))
})

test_that("individual layout with a smaller batch is still a valid simulation", {
  s <- streamed(53, "individual", 5L, m = 40, n = 24)
  expect_equal(dim(s$X), c(24L, 40L))
  expect_true(all(s$X %in% c(0L, 1L, 2L)))
})

test_that("streaming works for negative r", {
  ref <- inmem(54, r = -0.5)
  s <- streamed(54, "variant", 7L, r = -0.5)
  expect_equal(s$X, matrix(as.integer(ref$X), nrow = 24))
})

test_that("streaming omits X, returns dimensions, and saves an rds", {
  s <- streamed(55, "variant", 8L)
  expect_null(s$out$X)
  expect_identical(s$out$n, 24L)
  expect_identical(s$out$m, 40L)
  expect_identical(s$out$format, "variant")
  expect_true(file.exists(paste0(s$path, ".rds")))
  expect_equal(readRDS(paste0(s$path, ".rds"))$AF, s$out$AF)
})

test_that("the streamed rds records the call parameters as provenance", {
  p <- file.path(tempdir(), "sim-provenance")
  set.seed(58)
  out <- am_simulate(h2_0 = 0.5, r = 0.4, m = 40, n = 24,
                     min_MAF = 0.15, path = p, format = "variant",
                     batch_size = 8L)
  expect_identical(out$h2_0, 0.5)
  expect_identical(out$r, 0.4)
  expect_identical(out$min_MAF, 0.15)

  saved <- readRDS(paste0(p, ".rds"))
  expect_identical(saved$h2_0, 0.5)
  expect_identical(saved$r, 0.4)
  expect_identical(saved$min_MAF, 0.15)
})

test_that("path = NULL leaves the in-memory result untouched", {
  a <- inmem(56)
  b <- inmem(56)
  expect_identical(a, b)
  expect_true(is.matrix(a$X))
})

test_that("haplotypes cannot be combined with streaming", {
  expect_error(
    am_simulate(0.5, 0.4, 20, 10, haplotypes = TRUE,
                path = file.path(tempdir(), "sim-hap")),
    "haplotypes")
})

test_that("an auto batch size is chosen when none is given", {
  s <- streamed(57, "variant", NULL)
  expect_equal(dim(s$X), c(24L, 40L))
})

test_that("invalid batch_size values are rejected with a clear error", {
  bad <- list(
    "non-integer" = 2.7,
    "zero" = 0,
    "negative" = -3,
    "NA" = NA_real_,
    "vector" = c(1, 2),
    "non-numeric" = "5"
  )
  for (fmt in c("individual", "variant")) {
    for (nm in names(bad)) {
      p <- file.path(tempdir(), paste0("sim-badbatch-", fmt, "-", nm))
      expect_error(
        am_simulate(0.5, 0.4, 10, 8, path = p, format = fmt,
                    batch_size = bad[[nm]]),
        "positive whole number",
        info = paste(fmt, nm)
      )
    }
  }
})

test_that("the batch_size = 2.7 corruption case now errors instead of silently corrupting data", {
  p <- file.path(tempdir(), "sim-corruption-regression")
  expect_error(
    am_simulate(0.5, 0.4, 10, 8, path = p, format = "individual",
                batch_size = 2.7),
    "positive whole number"
  )
})

test_that("a rejected batch_size leaves no file behind", {
  p <- file.path(tempdir(), paste0("sim-nofile-", as.integer(Sys.time())))
  expect_false(any(file.exists(
    c(paste0(p, ".int8"), paste0(p, ".bed"), paste0(p, ".meta"),
      paste0(p, ".bim"), paste0(p, ".fam"), paste0(p, ".rds"))
  )))
  expect_error(
    am_simulate(0.5, 0.4, 10, 8, path = p, format = "bed",
                batch_size = 0),
    "positive whole number"
  )
  leftover <- list.files(dirname(p), pattern = paste0("^", basename(p)),
                          full.names = TRUE)
  expect_equal(leftover, character(0))
})

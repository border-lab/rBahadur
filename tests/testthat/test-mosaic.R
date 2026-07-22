## A small synthetic panel with a known map, so map-dependent behaviour can be
## asserted exactly rather than inferred from real data.
toy_panel <- function(N = 40, p = 60, hotspot = NULL, seed = 1) {
  set.seed(seed)
  H <- matrix(rbinom(N * p, 1, 0.5), nrow = N, ncol = p)
  ## guarantee both alleles are present at every marker
  H[1, ] <- 0L
  H[2, ] <- 1L
  pos <- seq(1000L, by = 1000L, length.out = p)
  cM <- if (is.null(hotspot)) {
    seq(0, 1, length.out = p)
  } else {
    ## flat everywhere except one interval carrying nearly all recombination
    step <- rep(1e-6, p - 1)
    step[hotspot] <- 1
    c(0, cumsum(step))
  }
  list(haplotypes = matrix(as.raw(H), nrow = N), pos = pos, cM = cM)
}

test_that("the bundled 1000 Genomes panel loads and looks right", {
  panel <- kg_reference()
  expect_true(is.raw(panel$haplotypes))
  expect_identical(dim(panel$haplotypes), c(520L, 2500L))
  expect_length(panel$pos, 2500L)
  expect_length(panel$cM, 2500L)
  expect_false(is.unsorted(panel$pos, strictly = TRUE))
  expect_false(is.unsorted(panel$cM))
  H <- matrix(as.integer(panel$haplotypes), nrow = 520)
  expect_true(all(H %in% c(0L, 1L)))
  ## it is a real panel, so it must carry real local LD
  expect_gt(cor(H[, 1], H[, 2])^2, cor(H[, 1], H[, 2500])^2)
})

test_that("malformed panels are rejected", {
  good <- toy_panel()
  expect_error(rBahadur:::.mosaic_check_panel(list(pos = 1:3)), "must be a list with")
  bad <- good; bad$pos <- rev(good$pos)
  expect_error(rBahadur:::.mosaic_check_panel(bad), "strictly increasing")
  bad <- good; bad$pos <- good$pos[-1]
  expect_error(rBahadur:::.mosaic_check_panel(bad), "one entry per column")
  bad <- good
  bad$haplotypes <- matrix(as.raw(rep(2L, length(good$pos) * 4)), nrow = 4)
  expect_error(rBahadur:::.mosaic_check_panel(bad), "only 0 and 1")
  bad <- good; bad$cM <- rev(good$cM)
  expect_error(rBahadur:::.mosaic_check_panel(bad), "non-decreasing")
})

test_that("a panel with no map falls back to physical distance", {
  panel <- toy_panel()
  panel$cM <- NULL
  checked <- rBahadur:::.mosaic_check_panel(panel)
  expect_length(checked$cM, length(panel$pos))
  expect_false(is.unsorted(checked$cM))
})

test_that("causal indices are validated", {
  panel <- toy_panel()
  expect_error(am_mosaic(0.5, 0.3, 5, panel), "either `causal_idx` or `m`")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, m = 1), "between 2 and")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, m = 1e6), "between 2 and")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, causal_idx = c(5L, 3L)),
               "strictly increasing")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, causal_idx = c(1L, 10000L)),
               "strictly increasing|within")
})

test_that("block boundaries follow the genetic map, not physical distance", {
  ## all recombination sits in interval 30, so essentially every boundary
  ## between causal markers 1 and 60 should land there
  panel <- rBahadur:::.mosaic_check_panel(toy_panel(hotspot = 30L))
  set.seed(1)
  ends <- rBahadur:::.mosaic_boundaries(c(1L, 60L), panel$cM, panel$p, 500)
  expect_gt(mean(ends[1, ] == 30L), 0.99)

  ## with a flat map the same call spreads boundaries out instead
  flat <- rBahadur:::.mosaic_check_panel(toy_panel())
  set.seed(1)
  ends_flat <- rBahadur:::.mosaic_boundaries(c(1L, 60L), flat$cM, flat$p, 500)
  expect_gt(length(unique(ends_flat[1, ])), 20L)
})

test_that("adjacent causal markers with no room between them are handled", {
  panel <- rBahadur:::.mosaic_check_panel(toy_panel())
  ends <- rBahadur:::.mosaic_boundaries(c(10L, 11L, 60L), panel$cM, panel$p, 20)
  expect_true(all(ends[1, ] == 10L))
  expect_true(all(ends[3, ] == panel$p))
})

test_that("causal loci carry exactly the alleles that were drawn", {
  ## the defining property: the mosaic must not perturb the causal variants,
  ## or the assortative mating structure is lost
  panel <- kg_reference()
  set.seed(3)
  sim <- am_mosaic(0.5, 0.5, n = 120, panel = panel, m = 25)
  g_from_X <- as.vector(sim$X[, sim$causal_idx, drop = FALSE] %*% sim$beta_raw)
  expect_equal(g_from_X, as.vector(sim$g))
})

test_that("simulated data carries local LD that am_simulate cannot produce", {
  panel <- kg_reference()
  set.seed(4)
  ## few causal variants means long blocks, so LD should reach well beyond them
  sim <- am_mosaic(0.5, 0.3, n = 400, panel = panel, m = 8)

  ## aggregate over many pairs: a single arbitrary pair says nothing, since
  ## two neighbouring markers may simply not be correlated in the panel
  r2 <- function(a, b) {
    v <- suppressWarnings(mapply(function(i, j) cor(sim$X[, i], sim$X[, j]),
                                 a, b))
    mean(v^2, na.rm = TRUE)
  }
  i <- seq(50, 2000, by = 25)
  near <- r2(i, i + 1)        # immediate neighbours
  far <- r2(i, i + 400)       # far apart, mostly across block boundaries
  expect_gt(near, 0.1)
  expect_lt(far, near / 4)

  ## and the panel itself should show the same ordering
  H <- matrix(as.integer(panel$haplotypes), nrow = nrow(panel$haplotypes))
  pn <- mean(suppressWarnings(mapply(function(a, b) cor(H[, a], H[, b]),
                                     i, i + 1))^2, na.rm = TRUE)
  expect_gt(pn, 0.1)
})

test_that("every output value is a valid diploid dosage", {
  panel <- kg_reference()
  set.seed(5)
  sim <- am_mosaic(0.5, -0.3, n = 60, panel = panel, m = 15)
  expect_true(all(sim$X %in% c(0L, 1L, 2L)))
  expect_identical(dim(sim$X), c(60L, 2500L))
  expect_length(sim$AF, 15L)
})

test_that("streaming reproduces the in-memory matrix for every format", {
  panel <- kg_reference()
  for (fmt in c("individual", "variant", "bed")) {
    for (bs in list(NULL, 9L)) {
      set.seed(31)
      ref <- am_mosaic(0.5, 0.3, n = 40, panel = panel, m = 12)$X
      p <- file.path(tempdir(), paste0("mos_", fmt, "_",
                                       if (is.null(bs)) "auto" else bs))
      set.seed(31)
      out <- am_mosaic(0.5, 0.3, n = 40, panel = panel, m = 12, path = p,
                       format = fmt, batch_size = bs)
      expect_null(out$X)
      expect_identical(read_genotypes(p), matrix(as.integer(ref), nrow = 40))
    }
  }
})

test_that("a bad batch_size is rejected before anything is written", {
  panel <- toy_panel()
  p <- file.path(tempdir(), "mos_badbatch")
  unlink(list.files(dirname(p), pattern = "mos_badbatch", full.names = TRUE))
  expect_error(am_mosaic(0.5, 0.3, 5, panel, m = 3, path = p,
                         batch_size = 2.5), "positive whole number")
  expect_length(list.files(dirname(p), pattern = "mos_badbatch"), 0L)
})

test_that("vcf_to_panel reads a phased VCF and rejects an unphased one", {
  hdr <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                 "INFO", "FORMAT", "s1", "s2", "s3", "s4"), collapse = "\t")
  rec <- function(pos, gts) paste(c("22", pos, ".", "A", "G", ".", ".", ".",
                                    "GT", gts), collapse = "\t")
  phased <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", hdr,
               rec(100, c("0|1", "1|1", "0|0", "1|0")),
               rec(200, c("1|0", "0|1", "1|1", "0|0"))), phased)
  panel <- vcf_to_panel(phased, min_maf = 0)
  expect_identical(dim(panel$haplotypes), c(8L, 2L))
  expect_identical(panel$pos, c(100L, 200L))
  expect_null(panel$cM)

  unphased <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", hdr,
               rec(100, c("0/1", "1/1", "0/0", "1/0"))), unphased)
  expect_error(vcf_to_panel(unphased, min_maf = 0), "not phased")

  expect_error(vcf_to_panel(tempfile()), "VCF not found")
})

test_that("vcf_to_panel applies the maf filter and attaches a map", {
  hdr <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                 "INFO", "FORMAT", "s1", "s2", "s3", "s4"), collapse = "\t")
  rec <- function(pos, gts) paste(c("22", pos, ".", "A", "G", ".", ".", ".",
                                    "GT", gts), collapse = "\t")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", hdr,
               rec(100, c("0|1", "1|1", "0|0", "1|0")),   # common
               rec(200, c("0|0", "0|0", "0|0", "1|0"))),  # rare
             vcf)
  expect_identical(ncol(vcf_to_panel(vcf, min_maf = 0.2)$haplotypes), 1L)
  expect_error(vcf_to_panel(vcf, min_maf = 0.9), "no markers survived")

  map <- tempfile(fileext = ".map")
  writeLines(c("22\t.\t0.0\t50", "22\t.\t1.0\t250"), map)
  panel <- vcf_to_panel(vcf, map = map, min_maf = 0)
  expect_length(panel$cM, 2L)
  expect_false(is.unsorted(panel$cM))
})

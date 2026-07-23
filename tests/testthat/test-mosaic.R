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
  bad <- good
  bad$haplotypes <- matrix(as.numeric(bad$haplotypes), nrow = nrow(bad$haplotypes))
  bad$haplotypes[1, 1] <- 0.5
  expect_error(rBahadur:::.mosaic_check_panel(bad), "only 0 and 1")
  bad <- good; bad$cM <- rev(good$cM)
  expect_error(rBahadur:::.mosaic_check_panel(bad), "non-decreasing")
  bad <- good; bad$pos[2] <- bad$pos[2] + 0.5
  expect_error(rBahadur:::.mosaic_check_panel(bad), "whole-number")
  bad <- good; bad$haplotypes <- bad$haplotypes[, 0, drop = FALSE]; bad$pos <- numeric()
  expect_error(rBahadur:::.mosaic_check_panel(bad), "two markers")
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
  expect_error(am_mosaic(0.5, 0.3, 5, panel, m = 3.5), "whole number")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, causal_idx = c(5L, 3L)),
               "strictly increasing")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, causal_idx = c(1L, 10000L)),
               "strictly increasing|within")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, causal_idx = c(1.9, 20)),
               "whole-number")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, causal_idx = c(1, NA)),
               "whole-number")
  expect_error(am_mosaic(0.5, 0.3, 5, panel, causal_idx = c("1", "20")),
               "whole-number")
  expect_error(am_mosaic(0, 0.3, 5, panel, m = 3), "h2_0")
  expect_error(am_mosaic(0.5, 1, 5, panel, m = 3), "`r`")
  expect_error(am_mosaic(0.5, 0.3, 0, panel, m = 3), "`n`")
})

test_that("am_mosaic warns but continues with few causal variants", {
  old <- getOption("rBahadur.warn_small_m")
  on.exit(options(rBahadur.warn_small_m = old), add = TRUE)
  options(rBahadur.warn_small_m = TRUE)

  expect_warning(
    sim <- am_mosaic(0.5, 0, n = 4, panel = toy_panel(), m = 3),
    "am_mosaic.*only 3 causal variants"
  )
  expect_identical(dim(sim$X), c(4L, 60L))
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

test_that("mosaic PLINK output preserves panel marker metadata", {
  panel <- toy_panel(p = 60)
  panel$chrom <- "22"
  panel$id <- paste0("rs", seq_len(60))
  panel$ref <- rep(c("A", "C"), 30)
  panel$alt <- rep(c("G", "T"), 30)
  p <- file.path(tempdir(), "mos_metadata")

  set.seed(32)
  am_mosaic(0.5, 0, n = 8, panel = panel, m = 10, path = p,
            format = "bed")
  bim <- utils::read.table(paste0(p, ".bim"), header = FALSE,
                           stringsAsFactors = FALSE)

  expect_identical(as.character(bim[[1]]), rep("22", 60))
  expect_identical(bim[[2]], panel$id)
  expect_equal(bim[[3]], panel$cM, tolerance = 1e-12)
  expect_equal(bim[[4]], panel$pos)
  expect_identical(bim[[5]], panel$alt)
  expect_identical(bim[[6]], panel$ref)
})

test_that("mosaic PLINK metadata is preflighted before simulation or writing", {
  panel <- toy_panel()
  panel$pos <- as.double(.Machine$integer.max) + seq_len(ncol(panel$haplotypes))
  p <- file.path(tempdir(), "mos_invalid_plink_position")
  unlink(paste0(p, c(".bed", ".bim", ".fam", ".meta", ".rds")))
  expect_error(
    am_mosaic(0.5, 0, 5, panel, m = 3, path = p, format = "bed"),
    "2\\^31 - 2"
  )
  expect_false(any(file.exists(paste0(
    p, c(".bed", ".bim", ".fam", ".meta", ".rds")
  ))))

  panel <- toy_panel()
  panel$chrom <- "chr 22"
  expect_error(
    am_mosaic(0.5, 0, 5, panel, m = 3, path = p, format = "bed"),
    "no whitespace"
  )
})

test_that("a haplotype is the concatenation of its donors' blocks", {
  ## .mosaic_hap switches to double indices once markers times haplotypes
  ## overflows integer. Such a panel is far too large to build here, so both
  ## branches are exercised on a small one and required to agree.
  panel <- rBahadur:::.mosaic_check_panel(toy_panel(N = 12, p = 40))
  Ht <- t(panel$H)
  donor_i <- c(3L, 7L, 1L, 12L)
  ends_i <- c(9L, 20L, 33L, 40L)
  narrow <- rBahadur:::.mosaic_hap(donor_i, ends_i, Ht, 40L, FALSE)
  wide <- rBahadur:::.mosaic_hap(donor_i, ends_i, Ht, 40L, TRUE)
  expect_identical(narrow, wide)

  H <- matrix(as.integer(panel$H), nrow = 12)
  expect_identical(narrow,
                   c(H[3, 1:9], H[7, 10:20], H[1, 21:33], H[12, 34:40]))
})

test_that("streaming matches memory when causal loci sit next to each other", {
  ## Adjacent causal markers leave no room for a boundary, so a block spans a
  ## single marker and the sweep has to advance pointers at neighbouring
  ## markers. batch_size 1 additionally puts a chunk boundary between every
  ## pair, so the carried-over pointer state is exercised at every step.
  panel <- toy_panel(N = 20, p = 50)
  idx <- c(1L, 2L, 3L, 25L, 49L, 50L)
  set.seed(21)
  ref <- am_mosaic(0.5, 0.3, n = 12, panel = panel, causal_idx = idx)$X
  for (fmt in c("individual", "variant", "bed")) {
    f <- file.path(tempdir(), paste0("mos_adj_", fmt))
    set.seed(21)
    am_mosaic(0.5, 0.3, n = 12, panel = panel, causal_idx = idx, path = f,
              format = fmt, batch_size = 1L)
    expect_identical(read_genotypes(f), ref)
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

test_that("vcf_to_panel reads phased VCFs and warns on unphased calls", {
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
  expect_identical(panel$chrom, c("22", "22"))
  expect_identical(panel$id, c(".", "."))
  expect_identical(panel$ref, c("A", "A"))
  expect_identical(panel$alt, c("G", "G"))
  expect_null(panel$cM)

  unphased <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", hdr,
               rec(100, c("0/1", "1/1", "0/0", "1/0"))), unphased)
  expect_warning(
    unphased_panel <- vcf_to_panel(unphased, min_maf = 0),
    "unphased.*treated as phase"
  )
  expect_identical(dim(unphased_panel$haplotypes), c(8L, 1L))

  mixed <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", hdr,
               rec(100, c("0|1", "1/1", "0|0", "1|0"))), mixed)
  expect_warning(vcf_to_panel(mixed, min_maf = 0), "unphased")

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
  expect_error(vcf_to_panel(vcf, min_maf = 0.9), "\\[0, 0.5\\]")

  map <- tempfile(fileext = ".map")
  writeLines(c("22\t.\t0.0\t50", "22\t.\t1.0\t250"), map)
  panel <- vcf_to_panel(vcf, map = map, min_maf = 0)
  expect_length(panel$cM, 2L)
  expect_false(is.unsorted(panel$cM))
})

test_that("VCF records are checked independently and symbolic alleles are excluded", {
  header <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
              "INFO", "FORMAT", "s1")
  hdr <- paste(header, collapse = "\t")
  vcf <- tempfile(fileext = ".vcf")

  ## The two bad widths sum to exactly two full records. Flattening all fields
  ## before checking would therefore accept this malformed file.
  short <- paste(c("22", "100", ".", "A", "G", ".", ".", ".", "GT"),
                 collapse = "\t")
  long <- paste(c("22", "200", ".", "A", "G", ".", ".", ".", "GT",
                  "0|1", "extra"), collapse = "\t")
  writeLines(c("##fileformat=VCFv4.2", hdr, short, long), vcf)
  expect_error(vcf_to_panel(vcf, min_maf = 0), "record on line 2 has 9 fields")

  rec <- function(pos, ref, alt) {
    paste(c("22", pos, ".", ref, alt, ".", ".", ".", "GT", "0|1"),
          collapse = "\t")
  }
  writeLines(c("##fileformat=VCFv4.2", hdr,
               rec(100, "A", "*"), rec(200, "C", "T")), vcf)
  panel <- vcf_to_panel(vcf, min_maf = 0)
  expect_identical(panel$pos, 200L)
  expect_identical(panel$ref, "C")
  expect_identical(panel$alt, "T")
})

test_that("VCFs and genetic maps cannot mix chromosomes", {
  hdr <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                 "INFO", "FORMAT", "s1"), collapse = "\t")
  rec <- function(chrom, pos) {
    paste(c(chrom, pos, ".", "A", "G", ".", ".", ".", "GT", "0|1"),
          collapse = "\t")
  }
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", hdr, rec("22", 100), rec("1", 200)),
             vcf)
  expect_error(vcf_to_panel(vcf, min_maf = 0), "exactly one chromosome")

  writeLines(c("##fileformat=VCFv4.2", hdr, rec("22", 100), rec("22", 200)),
             vcf)
  map <- tempfile(fileext = ".map")
  writeLines(c("1\t.\t50\t100", "1\t.\t100\t200",
               "chr22\t.\t1\t100", "chr22\t.\t2\t200"), map)
  panel <- vcf_to_panel(vcf, map = map, min_maf = 0)
  expect_identical(panel$cM, c(1, 2))

  writeLines(c("1\t.\t50\t100", "1\t.\t100\t200"), map)
  expect_error(vcf_to_panel(vcf, map = map, min_maf = 0),
               "at least two rows for chromosome 22")
})

test_that("vcf_to_panel rejects malformed genotype fields", {
  hdr <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                 "INFO", "FORMAT", "s1"), collapse = "\t")
  rec <- function(format, gt) paste(c("22", "100", "rs1", "A", "G", ".",
                                      ".", ".", format, gt), collapse = "\t")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", hdr, rec("DP:GT", "8:0|1")), vcf)
  expect_error(vcf_to_panel(vcf, min_maf = 0), "FORMAT.*begin with GT")

  writeLines(c("##fileformat=VCFv4.2", hdr, rec("GT", "0|2")), vcf)
  expect_error(vcf_to_panel(vcf, min_maf = 0), "alleles 0, 1, or")

  expect_error(vcf_to_panel(character()), "single non-empty path")
})

test_that("missing external panel commands are named clearly", {
  expect_error(
    rBahadur:::.panel_require_commands("rbahadur-command-that-does-not-exist"),
    "requires external command.*rbahadur-command-that-does-not-exist"
  )
})

test_that("download coordinates and cache paths are validated before networking", {
  expect_error(download_1kg_panel("../22", 100, 200), "autosome")
  expect_error(download_1kg_panel("X", 100, 200), "autosome")
  expect_error(download_1kg_panel("22", 100.5, 200), "whole-number")
  expect_error(download_1kg_panel("22", 200, 100), "start <= end")
  expect_error(download_1kg_panel("22", 100, 200, dest = character()),
               "single non-empty")
})

test_that("streamed regions are complete and published atomically", {
  skip_on_os("windows")
  skip_if(any(!nzchar(Sys.which(c("curl", "zcat", "awk")))))

  make_gzip <- function(pos) {
    path <- tempfile(fileext = ".vcf.gz")
    con <- gzfile(path, "wt")
    on.exit(close(con))
    header <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                      "FILTER", "INFO", "FORMAT", "s1"), collapse = "\t")
    records <- vapply(
      pos,
      function(bp) paste(c("22", bp, ".", "A", "G", ".", ".", ".",
                           "GT", "0|1"), collapse = "\t"),
      character(1)
    )
    writeLines(c("##fileformat=VCFv4.2", header, records), con)
    path
  }

  complete <- make_gzip(c(100, 200, 300))
  region <- tempfile(fileext = ".vcf")
  unlink(region)
  rBahadur:::.panel_stream_region(
    paste0("file://", normalizePath(complete, winslash = "/")),
    100, 200, region
  )
  out <- readLines(region)
  expect_length(out, 4L)
  expect_match(out[1], "^##fileformat=")
  expect_match(out[3], "^22\\t100\\t")
  expect_match(out[4], "^22\\t200\\t")

  incomplete <- make_gzip(c(100, 200))
  failed_region <- tempfile(fileext = ".vcf")
  unlink(failed_region)
  expect_error(
    rBahadur:::.panel_stream_region(
      paste0("file://", normalizePath(incomplete, winslash = "/")),
      100, 200, failed_region
    ),
    "did not reach"
  )
  expect_false(file.exists(failed_region))
})

tmp <- function(x) file.path(tempdir(), x)

## Run the CLI and return only its exit status, swallowing both the stderr
## diagnostics and the stdout report so the test log stays pristine.
quiet_main <- function(args) {
  status <- NULL
  invisible(utils::capture.output(
    status <- suppressMessages(rbahadur_main(args))))
  status
}

## Same, but leaves stdout intact for the tests that assert on what is printed.
loud_main <- function(args) suppressMessages(rbahadur_main(args))

test_that("the argument parser handles values, flags, and negative numbers", {
  p <- rBahadur:::.cli_parse(c("simulate", "--r", "-0.3", "--csv",
                               "--batch-size", "40", "--out", "x"))
  expect_identical(p$positional, "simulate")
  expect_identical(p$opts$r, "-0.3")          # a leading dash is a value
  expect_true(p$opts$csv)                     # flags take no value
  expect_identical(p$opts$batch_size, "40")   # --batch-size maps to batch_size
  expect_identical(p$opts$out, "x")
})

test_that("an option missing its value is rejected rather than silently eaten", {
  expect_error(rBahadur:::.cli_parse(c("--h2", "--r", "0.3")),
               "requires a value", class = "rbahadur_usage_error")
  expect_error(rBahadur:::.cli_parse(c("--r", "0.2", "--r", "0.3")),
               "more than once", class = "rbahadur_usage_error")
})

test_that("exit statuses distinguish success, usage errors, and failures", {
  p <- tmp("cli_status")
  expect_identical(quiet_main("--help"), 0L)
  expect_identical(quiet_main("--version"), 0L)
  expect_identical(quiet_main(character(0)), 1L)
  expect_identical(quiet_main("frobnicate"), 1L)
  expect_identical(quiet_main(c("simulate", "--h2", ".5", "--r", ".3",
                                "--m", "10", "--n", "10")), 1L)
  expect_identical(quiet_main(c("simulate", "--h2", "abc", "--r", ".3",
                                "--m", "10", "--n", "10", "--out", p)), 1L)
  expect_identical(quiet_main(c("simulate", "--h2", ".5", "--r", ".3",
                                "--m", "10", "--n", "10", "--out", p,
                                "--format", "vcf")), 1L)
  expect_identical(quiet_main("info"), 1L)
  expect_identical(quiet_main(c("info", tmp("does_not_exist"))), 2L)
})

test_that("unknown options and extra positional arguments are usage errors", {
  p <- tmp("cli_unknown_option")
  old <- list.files(dirname(p), pattern = "^cli_unknown_option", full.names = TRUE)
  if (length(old)) unlink(old)

  expect_identical(
    quiet_main(c("simulate", "--h2", ".5", "--r", ".3", "--m", "50",
                 "--n", "10", "--out", p, "--formt", "bed")),
    1L
  )
  expect_length(list.files(dirname(p), pattern = "^cli_unknown_option"), 0L)
  expect_identical(quiet_main(c("simulate", "extra")), 1L)
  expect_identical(quiet_main(c("info", "prefix", "extra")), 1L)
  expect_identical(quiet_main(c("info", "prefix", "--quiet")), 1L)
})

test_that("simulate writes the expected files for every format", {
  expected <- list(individual = ".int8", variant = ".int8", bed = ".bed")
  for (fmt in names(expected)) {
    p <- tmp(paste0("cli_fmt_", fmt))
    expect_identical(
      quiet_main(c("simulate", "--h2", "0.5", "--r", "-0.3", "--m", "40",
                   "--n", "37", "--out", p, "--format", fmt, "--quiet")), 0L)
    expect_true(file.exists(paste0(p, expected[[fmt]])))
    expect_true(file.exists(paste0(p, ".meta")))
    expect_true(file.exists(paste0(p, ".rds")))
    if (fmt == "bed") {
      expect_true(file.exists(paste0(p, ".bim")))
      expect_true(file.exists(paste0(p, ".fam")))
    }
    expect_equal(dim(read_genotypes(p)), c(37L, 40L))
  }
})

test_that("--seed makes runs reproducible and differing seeds differ", {
  a <- tmp("cli_seed_a"); b <- tmp("cli_seed_b"); c3 <- tmp("cli_seed_c")
  args <- function(p, s) c("simulate", "--h2", "0.5", "--r", "0.3", "--m", "40",
                           "--n", "30", "--out", p, "--seed", s, "--quiet")
  quiet_main(args(a, "99")); quiet_main(args(b, "99")); quiet_main(args(c3, "100"))
  expect_identical(read_genotypes(a), read_genotypes(b))
  expect_false(identical(read_genotypes(a), read_genotypes(c3)))

  zero <- tmp("cli_seed_zero")
  expect_identical(quiet_main(c("simulate", "--h2", "0.5", "--r", "0",
                                "--m", "40", "--n", "10", "--out", zero,
                                "--seed", "0", "--quiet")), 0L)
})

test_that("--csv writes portable sidecars with the right shape", {
  p <- tmp("cli_csv")
  quiet_main(c("simulate", "--h2", "0.5", "--r", "0.3", "--m", "25",
               "--n", "20", "--out", p, "--csv", "--quiet"))
  pheno <- read.csv(paste0(p, "_pheno.csv"))
  vars <- read.csv(paste0(p, "_variants.csv"))
  expect_identical(nrow(pheno), 20L)
  expect_identical(names(pheno), c("y", "g"))
  expect_identical(nrow(vars), 25L)
  expect_identical(names(vars), c("AF", "beta_std", "beta_raw"))
})

test_that("--quiet suppresses the summary and its absence prints one", {
  p <- tmp("cli_quiet")
  args <- c("simulate", "--h2", "0.5", "--r", "0.3", "--m", "20", "--n", "20",
            "--out", p)
  expect_silent(loud_main(c(args, "--quiet")))
  expect_output(loud_main(args), "individuals x")
})

test_that("info reports ok for an intact run and flags a truncated one", {
  p <- tmp("cli_info")
  quiet_main(c("simulate", "--h2", "0.5", "--r", "0.3", "--m", "30",
               "--n", "24", "--out", p, "--quiet"))
  expect_output(expect_identical(loud_main(c("info", p)), 0L), "status    ok")

  ## truncate to a length that still divides n*m, the silent-recycling case
  f <- paste0(p, ".int8")
  all_bytes <- readBin(f, "raw", file.size(f))
  writeBin(all_bytes[seq_len(length(all_bytes) / 2)], f)
  expect_output(expect_identical(loud_main(c("info", p)), 2L), "CORRUPT")
})

test_that("info flags malformed metadata instead of guessing a layout", {
  p <- tmp("cli_bad_meta")
  quiet_main(c("simulate", "--h2", "0.5", "--r", "0.3", "--m", "30",
               "--n", "24", "--out", p, "--quiet"))
  meta <- readLines(paste0(p, ".meta"))
  meta <- sub("format: individual", "format: varaint", meta, fixed = TRUE)
  writeLines(meta, paste0(p, ".meta"))
  expect_output(expect_identical(loud_main(c("info", p)), 2L), "CORRUPT")
})

test_that("info detects invalid int8 dosage bytes", {
  p <- tmp("cli_bad_byte")
  quiet_main(c("simulate", "--h2", "0.5", "--r", "0.3", "--m", "30",
               "--n", "24", "--out", p, "--quiet"))
  f <- paste0(p, ".int8")
  bytes <- readBin(f, "raw", n = file.size(f))
  bytes[10] <- as.raw(255)
  writeBin(bytes, f)
  expect_output(expect_identical(loud_main(c("info", p)), 2L), "CORRUPT")
})

test_that("the shipped script exists and is a thin shim", {
  path <- rbahadur_cli_path()
  expect_true(nzchar(path))
  expect_true(file.exists(path))
  script <- readLines(path)
  expect_match(script[1], "^#!")
  expect_true(any(grepl("rbahadur_main", script)))
})

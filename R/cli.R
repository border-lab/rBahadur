## Command line interface for rBahadur.
##
## The logic lives here rather than in the exec/ script so that it can be unit
## tested: `rbahadur_main()` returns an exit status instead of calling quit(),
## and the shipped script is a one-line shim around it.
##
## Argument parsing is hand rolled to keep the package free of dependencies.

.cli_flags <- c("csv", "quiet", "help", "version")

## Signal a bad invocation, as distinct from a run that started and then
## failed. The two exit with different statuses so calling scripts can tell a
## typo from, say, an infeasible set of simulation parameters.
.cli_abort <- function(...) {
  stop(structure(
    class = c("rbahadur_usage_error", "error", "condition"),
    list(message = paste0(...), call = NULL)))
}

## Split `args` into named options and positional arguments. Options are
## `--name value`, except for the flags above which take no value. A value is
## allowed to start with a single dash so that negative numbers work, as in
## `--r -0.3`.
.cli_parse <- function(args) {
  opts <- list()
  positional <- character(0)
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (grepl("^--", a)) {
      key <- gsub("-", "_", sub("^--", "", a))
      if (!is.null(opts[[key]])) {
        .cli_abort("option `", a, "` was supplied more than once")
      }
      if (key %in% .cli_flags) {
        opts[[key]] <- TRUE
        i <- i + 1L
      } else {
        if (i == length(args) || grepl("^--", args[i + 1L])) {
          .cli_abort("option `", a, "` requires a value")
        }
        opts[[key]] <- args[i + 1L]
        i <- i + 2L
      }
    } else {
      positional <- c(positional, a)
      i <- i + 1L
    }
  }
  list(opts = opts, positional = positional)
}

## Command-specific allowlists keep misspelled options from being silently
## ignored (for example, `--formt bed` previously wrote the default layout).
.cli_validate_opts <- function(opts, allowed, command = NULL) {
  unknown <- setdiff(names(opts), allowed)
  if (!length(unknown)) return(invisible(NULL))
  shown <- paste0("--", gsub("_", "-", unknown))
  where <- if (is.null(command)) "" else paste0(" for `", command, "`")
  .cli_abort(
    "unknown option", if (length(shown) > 1L) "s" else "",
    where, ": ", paste(shown, collapse = ", ")
  )
}

## Fetch and validate a numeric option.
.cli_num <- function(opts, name, required = TRUE, default = NULL,
                     integer = FALSE, positive = TRUE) {
  if (is.null(opts[[name]])) {
    if (required) {
      .cli_abort("missing required option `--", gsub("_", "-", name), "`")
    }
    return(default)
  }
  value <- suppressWarnings(as.numeric(opts[[name]]))
  if (is.na(value) || !is.finite(value)) {
    .cli_abort("option `--", gsub("_", "-", name), "` must be numeric, got '",
               opts[[name]], "'")
  }
  if (integer) {
    lower <- if (positive) 1 else 0
    if (value < lower || value %% 1 != 0 || value > .Machine$integer.max) {
      .cli_abort("option `--", gsub("_", "-", name),
                "` must be a ", if (positive) "positive" else "non-negative",
                " whole number, got '", opts[[name]], "'")
    }
    value <- as.integer(value)
  }
  value
}

.cli_usage <- function() {
  c("rbahadur: simulate genotype/phenotype data under assortative mating",
    "",
    "Usage:",
    "  rbahadur simulate --h2 <num> --r <num> --m <int> --n <int> --out <prefix>",
    "                    [--format individual|variant|bed] [--batch-size <int>]",
    "                    [--min-maf <num>] [--seed <int>] [--csv] [--quiet]",
    "  rbahadur info <prefix>",
    "  rbahadur --help | --version",
    "",
    "simulate options:",
    "  --h2          generation zero (panmictic) heritability, in (0, 1)",
    "  --r           cross-mate phenotypic correlation, in (-1, 1).",
    "                Negative values give disassortative mating.",
    "  --m           number of biallelic causal variants",
    "  --n           number of individuals",
    "  --out         output file prefix",
    "  --format      on-disk layout, default 'individual'. 'individual' and",
    "                'variant' write <prefix>.int8 at one byte per genotype;",
    "                'bed' writes variant-major PLINK at two bits per genotype",
    "                plus <prefix>.bim and <prefix>.fam",
    "  --batch-size  individuals per batch for 'individual', variants per block",
    "                for 'variant' and 'bed'. Defaults to a memory-based guess.",
    "  --min-maf     minimum minor allele frequency, default 0.1",
    "  --seed        integer seed, for reproducible output",
    "  --csv         also write <prefix>_pheno.csv and <prefix>_variants.csv,",
    "                which are readable outside R",
    "  --quiet       suppress the summary written to stdout",
    "",
    "Every run writes <prefix>.meta describing the layout, and <prefix>.rds",
    "holding allele frequencies, effect sizes, and phenotypes.",
    "",
    "Note that disassortative mating leaves the Bahadur order-2 feasible region",
    "sooner than assortative mating does, and more readily as n grows. If the",
    "run reports infeasible probabilities, reduce the magnitude of --r, raise",
    "--min-maf, or increase --m.")
}

.cli_simulate <- function(opts) {
  .cli_validate_opts(
    opts,
    c("h2", "r", "m", "n", "out", "format", "batch_size", "min_maf",
      "seed", "csv", "quiet", "help", "version"),
    "simulate"
  )
  h2_0 <- .cli_num(opts, "h2")
  r <- .cli_num(opts, "r")
  m <- .cli_num(opts, "m", integer = TRUE)
  n <- .cli_num(opts, "n", integer = TRUE)
  min_maf <- .cli_num(opts, "min_maf", required = FALSE, default = 0.1)
  batch_size <- .cli_num(opts, "batch_size", required = FALSE, integer = TRUE)
  seed <- .cli_num(opts, "seed", required = FALSE, integer = TRUE,
                   positive = FALSE)

  if (is.null(opts$out)) .cli_abort("missing required option `--out`")
  out <- opts$out
  outdir <- dirname(out)
  if (!dir.exists(outdir)) {
    .cli_abort("output directory does not exist: ", outdir)
  }

  format <- if (is.null(opts$format)) "individual" else opts$format
  if (!format %in% c("individual", "variant", "bed")) {
    .cli_abort("`--format` must be one of individual, variant, or bed, got '",
               format, "'")
  }
  if (h2_0 <= 0 || h2_0 >= 1) {
    .cli_abort("`--h2` must lie in the open interval (0, 1)")
  }
  if (r <= -1 || r >= 1) {
    .cli_abort("`--r` must lie in the open interval (-1, 1)")
  }
  if (m < 2L) .cli_abort("`--m` must be at least 2")
  if (min_maf < 0 || min_maf > 0.5) {
    .cli_abort("`--min-maf` must lie in the closed interval [0, 0.5]")
  }

  if (!is.null(seed)) set.seed(seed)
  res <- am_simulate(h2_0 = h2_0, r = r, m = m, n = n, min_MAF = min_maf,
                     path = out, format = format, batch_size = batch_size)

  if (isTRUE(opts$csv)) {
    utils::write.csv(data.frame(y = as.vector(res$y), g = as.vector(res$g)),
                     paste0(out, "_pheno.csv"), row.names = FALSE)
    utils::write.csv(data.frame(AF = res$AF,
                                beta_std = as.vector(res$beta_std),
                                beta_raw = as.vector(res$beta_raw)),
                     paste0(out, "_variants.csv"), row.names = FALSE)
  }

  if (!isTRUE(opts$quiet)) {
    data_file <- .gt_data_path(out, format)
    cat(sprintf("wrote %s (%s)\n", data_file,
                format(structure(file.size(data_file), class = "object_size"),
                       units = "auto")))
    cat(sprintf("  %d individuals x %d variants, format %s\n", n, m, format))
    cat(sprintf("  equilibrium h2 %.4f, genetic variance %.4f\n",
                h2_eq(r, h2_0), vg_eq(r, h2_0, h2_0)))
    if (isTRUE(opts$csv)) {
      cat(sprintf("  also wrote %s_pheno.csv and %s_variants.csv\n", out, out))
    }
  }
  0L
}

.cli_info <- function(prefix) {
  meta_error <- NULL
  meta <- tryCatch(
    .gt_read_meta(prefix),
    error = function(e) {
      meta_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(meta)) {
    cat(sprintf("prefix    %s\n", prefix))
    cat(sprintf("status    CORRUPT: %s\n", meta_error))
    return(2L)
  }
  data_file <- .gt_data_path(prefix, meta$format)
  expected <- if (meta$format == "bed") {
    3 + as.double(meta$m) * ceiling(meta$n / 4)
  } else {
    as.double(meta$n) * meta$m
  }
  actual <- if (file.exists(data_file)) file.size(data_file) else NA_real_

  cat(sprintf("prefix    %s\n", prefix))
  cat(sprintf("format    %s\n", meta$format))
  cat(sprintf("dtype     %s\n", meta$dtype))
  cat(sprintf("n         %d individuals\n", meta$n))
  cat(sprintf("m         %d variants\n", meta$m))
  cat(sprintf("data      %s\n", data_file))
  cat(sprintf("size      %s bytes actual, %s expected\n",
              format(actual, scientific = FALSE),
              format(expected, scientific = FALSE)))
  if (is.na(actual)) {
    cat("status    MISSING: data file not found\n")
    return(2L)
  }
  if (actual != expected) {
    cat("status    CORRUPT: size does not match the metadata\n")
    return(2L)
  }
  integrity_error <- tryCatch({
    if (meta$format == "bed") {
      .gt_check_bed_header(data_file)
    } else {
      .gt_check_int8_values(data_file)
    }
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(integrity_error)) {
    cat(sprintf("status    CORRUPT: %s\n", integrity_error))
    return(2L)
  }
  cat("status    ok\n")
  0L
}

#' Path to the `rbahadur` command line script
#'
#' Returns the location of the executable shipped with the package, so it can
#' be placed on the search path.
#'
#' @return A single string giving the path to the script, or `""` if the
#'   package was installed without it.
#' @export
#'
#' @examples
#' rbahadur_cli_path()
#'
#' ## to put it on your PATH:
#' ## ln -s $(Rscript -e 'cat(rBahadur::rbahadur_cli_path())') ~/bin/rbahadur
rbahadur_cli_path <- function() {
  system.file("exec", "rbahadur", package = "rBahadur")
}

#' Run the `rbahadur` command line interface
#'
#' Dispatches a vector of command line arguments. This is what the shipped
#' `rbahadur` script calls; it is exported so the interface can be driven from
#' R and tested. It returns an exit status rather than terminating the session.
#'
#' @param args character vector of command line arguments, as returned by
#'   `commandArgs(trailingOnly = TRUE)`
#'
#' @return Invisibly, an integer exit status: `0` on success, `1` for a usage
#'   error, and `2` for a runtime failure. Usage text goes to stdout when
#'   requested and to stderr when the invocation was wrong.
#' @export
#'
#' @examples
#' p <- file.path(tempdir(), "cli_example")
#' rbahadur_main(c("simulate", "--h2", "0.5", "--r", "0.3",
#'                 "--m", "50", "--n", "40", "--out", p, "--quiet"))
#' rbahadur_main(c("info", p))
rbahadur_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  status <- tryCatch({
    parsed <- .cli_parse(args)
    opts <- parsed$opts
    positional <- parsed$positional

    if (isTRUE(opts$version)) {
      .cli_validate_opts(opts, c("version"))
      cat(as.character(utils::packageVersion("rBahadur")), "\n", sep = "")
      return(invisible(0L))
    }
    if (isTRUE(opts$help)) {
      .cli_validate_opts(opts, c("help"))
      cat(.cli_usage(), sep = "\n")
      return(invisible(0L))
    }
    if (length(positional) == 0L) {
      .cli_validate_opts(opts, character(0))
      cat(.cli_usage(), sep = "\n")
      return(invisible(1L))
    }

    command <- positional[1]
    switch(command,
      simulate = {
        if (length(positional) != 1L) {
          .cli_abort("`simulate` does not take positional arguments")
        }
        .cli_simulate(opts)
      },
      info = {
        .cli_validate_opts(opts, character(0), "info")
        if (length(positional) < 2L) {
          .cli_abort("`info` requires a prefix, as in `rbahadur info my_run`")
        }
        if (length(positional) > 2L) {
          .cli_abort("`info` takes exactly one prefix")
        }
        .cli_info(positional[2])
      },
      {
        message("unknown command '", command, "'")
        message(paste(.cli_usage(), collapse = "\n"))
        1L
      })
  },
  rbahadur_usage_error = function(e) {
    message("rbahadur: ", conditionMessage(e))
    message("try `rbahadur --help`")
    1L
  },
  error = function(e) {
    message("rbahadur: ", conditionMessage(e))
    2L
  })
  invisible(status)
}

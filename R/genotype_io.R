## ---- internal path and metadata helpers -------------------------------

.gt_data_path <- function(path, format) {
  if (format == "bed") paste0(path, ".bed") else paste0(path, ".int8")
}

.gt_write_meta <- function(path, n, m, format) {
  writeLines(c(
    "rBahadur_genotypes: 1",
    paste0("format: ", format),
    paste0("n: ", n),
    paste0("m: ", m),
    paste0("dtype: ", if (format == "bed") "bed2bit" else "int8")
  ), paste0(path, ".meta"))
  invisible(NULL)
}

.gt_read_meta <- function(path) {
  f <- paste0(path, ".meta")
  if (!file.exists(f)) stop("metadata file not found: ", f)
  kv <- strsplit(readLines(f), ":[[:space:]]*")
  keys <- vapply(kv, function(x) x[1], character(1))
  vals <- vapply(kv, function(x) x[2], character(1))
  names(vals) <- keys
  list(format = unname(vals["format"]),
       n = as.integer(vals["n"]),
       m = as.integer(vals["m"]),
       dtype = unname(vals["dtype"]))
}

.gt_check_matrix <- function(X, format) {
  if (!is.matrix(X)) stop("`X` must be a matrix")
  storage.mode(X) <- "integer"
  if (anyNA(X) && format != "bed") {
    stop("missing genotypes are only representable in the 'bed' format")
  }
  if (any(X < 0L | X > 2L, na.rm = TRUE)) stop("genotypes must be 0, 1, or 2")
  X
}

## ---- PLINK bed helpers -------------------------------------------------
##
## Two bits per genotype, four samples per byte, lowest sample in the lowest
## bits. The effect allele is written as A1, and PLINK codes count A1 copies,
## so dosage 2 -> 00, dosage 1 -> 10, dosage 0 -> 11, and missing -> 01.

.gt_pack_bed <- function(g) {
  code <- integer(length(g))
  code[!is.na(g) & g == 2L] <- 0L
  code[!is.na(g) & g == 1L] <- 2L
  code[!is.na(g) & g == 0L] <- 3L
  code[is.na(g)] <- 1L
  pad <- (4L - (length(g) %% 4L)) %% 4L
  if (pad > 0L) code <- c(code, integer(pad))
  quad <- matrix(code, nrow = 4L)
  as.raw(quad[1, ] + quad[2, ] * 4L + quad[3, ] * 16L + quad[4, ] * 64L)
}

.gt_unpack_bed <- function(bytes, n) {
  b <- as.integer(bytes)
  codes <- as.vector(rbind(b %% 4L,
                           (b %/% 4L) %% 4L,
                           (b %/% 16L) %% 4L,
                           (b %/% 64L) %% 4L))[seq_len(n)]
  g <- integer(n)
  g[codes == 0L] <- 2L
  g[codes == 2L] <- 1L
  g[codes == 3L] <- 0L
  g[codes == 1L] <- NA_integer_
  g
}

.gt_write_plink_sidecars <- function(path, n, m) {
  writeLines(paste(1L, paste0("v", seq_len(m)), 0L, seq_len(m), "A", "G",
                   sep = "\t"), paste0(path, ".bim"))
  ids <- paste0("i", seq_len(n))
  writeLines(paste(ids, ids, 0L, 0L, 0L, -9L, sep = "\t"),
             paste0(path, ".fam"))
  invisible(NULL)
}

## ---- exported interface -----------------------------------------------

#' Write genotypes to a binary file
#'
#' @param X an integer matrix of genotypes with individuals in rows and
#'   variants in columns, taking values 0, 1, or 2
#' @param path file prefix. Layout `"individual"` and `"variant"` write
#'   `<path>.int8`; `"bed"` writes `<path>.bed` plus `<path>.bim` and
#'   `<path>.fam`. All layouts write `<path>.meta`.
#' @param format on-disk layout. `"individual"` (the default) stores each
#'   individual's variants contiguously, `"variant"` stores each variant's
#'   individuals contiguously, and `"bed"` writes a variant-major PLINK
#'   binary file at two bits per genotype.
#'
#' @return `path`, invisibly.
#' @export
#'
#' @examples
#' X <- matrix(sample(0:2, 20, replace = TRUE), nrow = 4)
#' p <- file.path(tempdir(), "example_genotypes")
#' write_genotypes(X, p)
#' identical(read_genotypes(p), X)
write_genotypes <- function(X, path, format = c("individual", "variant", "bed")) {
  format <- match.arg(format)
  X <- .gt_check_matrix(X, format)
  n <- nrow(X)
  m <- ncol(X)

  con <- file(.gt_data_path(path, format), "wb")
  on.exit(close(con))
  if (format == "variant") {
    writeBin(as.vector(X), con, size = 1L)
  } else if (format == "individual") {
    writeBin(as.vector(t(X)), con, size = 1L)
  } else {
    writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
    for (j in seq_len(m)) writeBin(.gt_pack_bed(X[, j]), con)
    .gt_write_plink_sidecars(path, n, m)
  }
  .gt_write_meta(path, n, m, format)
  invisible(path)
}

#' Read genotypes from a binary file written by `write_genotypes()`
#'
#' @param path file prefix, the same value passed to [write_genotypes()]
#'
#' @return An integer matrix with individuals in rows and variants in columns,
#'   reconstructed to the same orientation regardless of the on-disk layout.
#' @export
#'
#' @examples
#' X <- matrix(sample(0:2, 20, replace = TRUE), nrow = 4)
#' p <- file.path(tempdir(), "example_genotypes")
#' write_genotypes(X, p, format = "variant")
#' read_genotypes(p)
read_genotypes <- function(path) {
  meta <- .gt_read_meta(path)
  n <- meta$n
  m <- meta$m
  con <- file(.gt_data_path(path, meta$format), "rb")
  on.exit(close(con))
  if (meta$format == "bed") {
    hdr <- readBin(con, "raw", n = 3L)
    if (!identical(as.integer(hdr), c(0x6cL, 0x1bL, 0x01L))) {
      stop("not a variant-major PLINK .bed file")
    }
    nb <- ceiling(n / 4)
    X <- matrix(NA_integer_, nrow = n, ncol = m)
    for (j in seq_len(m)) X[, j] <- .gt_unpack_bed(readBin(con, "raw", n = nb), n)
    return(X)
  }
  v <- readBin(con, "integer", n = n * m, size = 1L, signed = TRUE)
  if (meta$format == "variant") {
    matrix(v, nrow = n, ncol = m)
  } else {
    t(matrix(v, nrow = m, ncol = n))
  }
}

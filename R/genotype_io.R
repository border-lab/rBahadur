## ---- internal path and metadata helpers -------------------------------

.gt_check_path <- function(path, write = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    stop("`path` must be a single non-empty character string")
  }
  if (write && !dir.exists(dirname(path))) {
    stop("output directory does not exist: ", dirname(path))
  }
  invisible(path)
}

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
  lines <- readLines(f, warn = FALSE)
  colon <- regexpr(":", lines, fixed = TRUE)
  if (!length(lines) || any(!nzchar(lines)) || any(colon < 2L)) {
    stop("malformed genotype metadata in ", f, call. = FALSE)
  }
  keys <- trimws(substring(lines, 1L, colon - 1L))
  vals <- trimws(substring(lines, colon + 1L))
  required <- c("rBahadur_genotypes", "format", "n", "m", "dtype")
  if (any(!nzchar(keys)) || any(!nzchar(vals)) || anyDuplicated(keys)) {
    stop("malformed or duplicate fields in genotype metadata: ", f,
         call. = FALSE)
  }
  unknown <- setdiff(keys, required)
  missing <- setdiff(required, keys)
  if (length(unknown) || length(missing)) {
    details <- c(
      if (length(missing)) paste0("missing ", paste(missing, collapse = ", ")),
      if (length(unknown)) paste0("unknown ", paste(unknown, collapse = ", "))
    )
    stop("invalid genotype metadata fields (", paste(details, collapse = "; "),
         "): ", f, call. = FALSE)
  }
  names(vals) <- keys
  if (!identical(unname(vals["rBahadur_genotypes"]), "1")) {
    stop("unsupported genotype metadata version in ", f, call. = FALSE)
  }
  format <- unname(vals["format"])
  if (!format %in% c("individual", "variant", "bed")) {
    stop("unrecognized genotype format in ", f, ": ", format,
         call. = FALSE)
  }
  dimension <- function(name) {
    value <- unname(vals[name])
    if (!grepl("^[1-9][0-9]*$", value)) {
      stop("metadata dimension `", name, "` must be a positive whole number",
           call. = FALSE)
    }
    number <- suppressWarnings(as.numeric(value))
    if (!is.finite(number) || number > .Machine$integer.max) {
      stop("metadata dimension `", name, "` is too large", call. = FALSE)
    }
    as.integer(number)
  }
  n <- dimension("n")
  m <- dimension("m")
  dtype <- unname(vals["dtype"])
  expected_dtype <- if (format == "bed") "bed2bit" else "int8"
  if (!identical(dtype, expected_dtype)) {
    stop("metadata dtype `", dtype, "` is incompatible with format `", format,
         "`", call. = FALSE)
  }
  list(format = format, n = n, m = m, dtype = dtype)
}

.gt_check_matrix <- function(X, format) {
  if (!is.matrix(X)) stop("`X` must be a matrix")
  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop("`X` must have at least one individual and one variant")
  }
  if ((!is.numeric(X) && !is.logical(X)) || any(is.nan(X)) ||
      any(is.infinite(X)) ||
      any(X != trunc(X), na.rm = TRUE)) {
    stop("genotypes must be integer-valued 0, 1, or 2")
  }
  if (anyNA(X) && format != "bed") {
    stop("missing genotypes are only representable in the 'bed' format")
  }
  if (any(X < 0 | X > 2, na.rm = TRUE)) stop("genotypes must be 0, 1, or 2")
  storage.mode(X) <- "integer"
  X
}

.gt_check_int8_values <- function(path, block = 16e6) {
  con <- file(path, "rb")
  on.exit(close(con))
  repeat {
    bytes <- readBin(con, "raw", n = block)
    if (!length(bytes)) break
    if (any(as.integer(bytes) > 2L)) {
      stop("genotype file contains values other than 0, 1, or 2: ", path)
    }
  }
  invisible(NULL)
}

.gt_check_bed_header <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  hdr <- readBin(con, "raw", n = 3L)
  if (!identical(as.integer(hdr), c(0x6cL, 0x1bL, 0x01L))) {
    stop("not a variant-major PLINK .bed file")
  }
  invisible(NULL)
}

## ---- PLINK bed helpers -------------------------------------------------
##
## Two bits per genotype, four samples per byte, lowest sample in the lowest
## bits. The effect allele is written as A1, and PLINK codes count A1 copies,
## so dosage 2 -> 00, dosage 1 -> 10, dosage 0 -> 11, and missing -> 01.

## `g` is one variant's genotypes, or a matrix of individuals by variants, in
## which case every column is packed and padded independently and the bytes are
## returned variant after variant, exactly as writing the columns one at a time
## would. Packing a whole block at once matters: at genome scale this is called
## once per variant otherwise, and the per-call vector overhead dominates.
.gt_pack_bed <- function(g) {
  nr <- if (is.matrix(g)) nrow(g) else length(g)
  ## a lookup recode rather than four masked assignments, which cost several
  ## full passes each
  code <- c(3L, 2L, 0L)[g + 1L]
  if (anyNA(code)) code[is.na(code)] <- 1L
  pad <- (4L - (nr %% 4L)) %% 4L
  if (pad > 0L) {
    padded <- matrix(0L, nrow = nr + pad, ncol = length(code) %/% nr)
    padded[seq_len(nr), ] <- code
    code <- padded
  }
  dim(code) <- c(4L, length(code) %/% 4L)
  as.raw(code[1L, ] + code[2L, ] * 4L + code[3L, ] * 16L + code[4L, ] * 64L)
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

.gt_prepare_plink_sidecars <- function(m, chrom = "1", id = NULL,
                                       cM = 0, pos = seq_len(m),
                                       a1 = "A", a2 = "G") {
  marker_field <- function(value, name) {
    if (length(value) == 1L) value <- rep(value, m)
    if (length(value) != m || anyNA(value)) {
      stop("`", name, "` must have length 1 or `m`")
    }
    value
  }

  chrom <- as.character(marker_field(chrom, "chrom"))
  cM <- as.numeric(marker_field(cM, "cM"))
  pos <- as.numeric(marker_field(pos, "pos"))
  a1 <- as.character(marker_field(a1, "a1"))
  a2 <- as.character(marker_field(a2, "a2"))
  if (any(!is.finite(cM)) || any(!is.finite(pos))) {
    stop("PLINK genetic and physical positions must be finite")
  }
  if (any(pos < 1) || any(pos != floor(pos)) ||
      any(pos > .Machine$integer.max - 1)) {
    stop(paste0("PLINK physical positions must be positive whole numbers no ",
                "larger than 2^31 - 2"))
  }

  text_fields <- list(chrom = chrom, a1 = a1, a2 = a2)
  bad_text <- vapply(
    text_fields,
    function(x) any(!nzchar(x) | grepl("[[:space:]]", x)),
    logical(1)
  )
  if (any(bad_text)) {
    stop("PLINK marker fields must be non-empty and contain no whitespace: `",
         names(bad_text)[which(bad_text)[1L]], "`")
  }

  if (is.null(id)) id <- paste0("v", seq_len(m))
  id <- as.character(marker_field(id, "id"))
  missing_id <- !nzchar(id) | id == "."
  id[missing_id] <- paste0("v", which(missing_id))
  id <- make.unique(id, sep = "_")
  if (any(grepl("[[:space:]]", id))) {
    stop("PLINK marker IDs must contain no whitespace")
  }

  list(chrom = chrom, id = id, cM = cM, pos = pos, a1 = a1, a2 = a2)
}

.gt_write_plink_sidecars <- function(path, n, m, chrom = "1", id = NULL,
                                     cM = 0, pos = seq_len(m),
                                     a1 = "A", a2 = "G") {
  fields <- .gt_prepare_plink_sidecars(m, chrom, id, cM, pos, a1, a2)

  writeLines(
    paste(
      fields$chrom,
      fields$id,
      format(fields$cM, scientific = FALSE, trim = TRUE, digits = 15),
      format(fields$pos, scientific = FALSE, trim = TRUE, digits = 15),
      fields$a1,
      fields$a2,
      sep = "\t"
    ),
    paste0(path, ".bim")
  )
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
#' @param path non-empty file prefix in an existing directory. Layout
#'   `"individual"` and `"variant"` write
#'   `<path>.int8`; `"bed"` writes `<path>.bed` plus `<path>.bim` and
#'   `<path>.fam`. All layouts write `<path>.meta`.
#' @param format on-disk layout. `"individual"` (the default) stores each
#'   individual's variants contiguously, `"variant"` stores each variant's
#'   individuals contiguously, and `"bed"` writes a variant-major PLINK
#'   binary file at two bits per genotype.
#'
#' @details Because this function receives no marker annotations, PLINK `.bim`
#' output uses placeholder chromosome, position, and allele values. PLINK files
#' written by [am_mosaic()] instead preserve the reference panel's available
#' chromosome, physical position, genetic map, ID, and allele metadata.
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
  .gt_check_path(path, write = TRUE)
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
    ## pack in blocks, sized so the intermediate stays around 64 MB
    step <- max(1L, min(m, as.integer(floor(16e6 / max(n, 1L)))))
    for (start in seq(1L, m, by = step)) {
      cols <- start:min(start + step - 1L, m)
      writeBin(.gt_pack_bed(X[, cols, drop = FALSE]), con)
    }
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
  .gt_check_path(path)
  meta <- .gt_read_meta(path)
  n <- meta$n
  m <- meta$m
  data_path <- .gt_data_path(path, meta$format)

  ## Validate the file size up front. A truncated or padded file must be
  ## rejected here: matrix(v, nrow, ncol) recycles silently whenever a short
  ## read divides n*m, and a bed file with too few or too many bytes would
  ## otherwise be read as if it were complete.
  expected_size <- if (meta$format == "bed") {
    3 + as.double(m) * ceiling(n / 4)
  } else {
    as.double(n) * m
  }
  actual_size <- file.size(data_path)
  if (is.na(actual_size)) {
    stop("genotype file not found: ", data_path)
  }
  if (actual_size != expected_size) {
    stop(sprintf(
      "genotype file '%s' has size %.0f bytes, but %.0f bytes were expected for n = %d, m = %d, format = '%s'. The file may be truncated or corrupted.",
      data_path, actual_size, expected_size, n, m, meta$format
    ))
  }

  con <- file(data_path, "rb")
  on.exit(close(con))
  if (meta$format == "bed") {
    hdr <- readBin(con, "raw", n = 3L)
    if (!identical(as.integer(hdr), c(0x6cL, 0x1bL, 0x01L))) {
      stop("not a variant-major PLINK .bed file")
    }
    nb <- ceiling(n / 4)
    X <- matrix(NA_integer_, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      bytes <- readBin(con, "raw", n = nb)
      if (length(bytes) != nb) {
        stop(sprintf(
          "genotype file '%s' is truncated: variant %d expected %d bytes but only %d were read.",
          data_path, j, nb, length(bytes)
        ))
      }
      X[, j] <- .gt_unpack_bed(bytes, n)
    }
    return(X)
  }
  v <- readBin(con, "integer", n = n * m, size = 1L, signed = TRUE)
  if (length(v) != n * m) {
    stop(sprintf(
      "genotype file '%s' is truncated: expected %d values but only %d were read.",
      data_path, n * m, length(v)
    ))
  }
  if (any(!v %in% 0:2)) {
    stop("genotype file contains values other than 0, 1, or 2: ", data_path)
  }
  if (meta$format == "variant") {
    matrix(v, nrow = n, ncol = m)
  } else {
    t(matrix(v, nrow = m, ncol = n))
  }
}

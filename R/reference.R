## Reference haplotype panels for am_mosaic().
##
## A panel is a plain list: `haplotypes` (haplotypes by markers, 0/1, stored as
## raw to keep it compact), `pos` (base pair positions), and `cM` (genetic map
## positions). Keeping it a list rather than a class means users can assemble
## one from any source without going through this package.

#' Load the bundled 1000 Genomes reference panel
#'
#' A small real reference panel, included so that examples, tests, and the
#' vignette run offline. It covers a 1 Mb window of chromosome 22 and is
#' intended for demonstration rather than analysis; use [vcf_to_panel()] or
#' [download_1kg_panel()] to build a panel at realistic scale.
#'
#' @return A list with `haplotypes` (520 by 2500 raw matrix of 0/1), `pos`,
#'   `cM`, `chrom`, `build`, and `source`.
#'
#' @source Phase 3 of the 1000 Genomes Project, GRCh38, chromosome 22
#'   positions 20,000,000 to 21,000,000, restricted to biallelic single
#'   nucleotide variants with minor allele frequency at least 0.05 in the first
#'   260 samples. Genetic map positions are interpolated from the GRCh38 PLINK
#'   maps distributed with BEAGLE.
#' @export
#'
#' @examples
#' panel <- kg_reference()
#' dim(panel$haplotypes)
#' range(panel$pos)
#'
#' ## realistic local LD is the point: correlation decays with distance
#' H <- matrix(as.integer(panel$haplotypes), nrow = nrow(panel$haplotypes))
#' cor(H[, 1], H[, 2])^2
#' cor(H[, 1], H[, 2000])^2
kg_reference <- function() {
  readRDS(system.file("extdata", "kg_chr22_panel.rds", package = "rBahadur"))
}

## Interpolate genetic map positions onto marker positions.
.panel_interpolate_cM <- function(pos, map_pos, map_cM) {
  if (!is.numeric(pos) || anyNA(pos) || any(!is.finite(pos))) {
    stop("marker positions must be finite numbers")
  }
  if (!is.numeric(map_pos) || !is.numeric(map_cM) ||
      length(map_pos) != length(map_cM) || length(map_pos) < 2L ||
      anyNA(map_pos) || anyNA(map_cM) ||
      any(!is.finite(map_pos)) || any(!is.finite(map_cM))) {
    stop("the genetic map must contain at least two finite positions")
  }
  if (any(map_pos < 1) || any(map_pos != floor(map_pos))) {
    stop("genetic-map base-pair positions must be positive whole numbers")
  }
  if (any(map_cM < 0)) {
    stop("genetic-map centimorgan positions must be non-negative")
  }
  ord <- order(map_pos)
  map_pos <- map_pos[ord]
  map_cM <- map_cM[ord]
  if (anyDuplicated(map_pos)) {
    stop("the genetic map contains duplicate base-pair positions")
  }
  if (is.unsorted(map_cM)) {
    stop("genetic-map centimorgan positions must be non-decreasing")
  }
  stats::approx(map_pos, map_cM, xout = pos, rule = 2)$y
}

## Treat the conventional "chr22" and "22" spellings as the same chromosome
## while otherwise preserving labels exactly for panel metadata and PLINK
## output.
.panel_chrom_key <- function(x) sub("^chr", "", tolower(as.character(x)))

#' Build a reference panel from a VCF
#'
#' Reads a biallelic VCF and returns a panel in the form [am_mosaic()] expects.
#' Phased genotypes such as `0|1` are strongly recommended because the mosaic
#' copies haplotypes rather than genotypes. Unphased genotypes such as `0/1`
#' are accepted with a warning and their written allele order is treated as if
#' it were phase.
#'
#' @param vcf path to a VCF, optionally gzipped
#' @param map optional path to a PLINK-format genetic map with columns
#'   chromosome, marker, centimorgans, and base pair position. Without one, the
#'   panel carries no `cM` and [am_mosaic()] places breakpoints uniformly in
#'   physical distance.
#' @param min_maf drop markers whose minor allele frequency in the panel falls
#'   below this. Must lie in `[0, 0.5]`; defaults to 0.05.
#' @param max_markers if given, thin evenly to at most this many markers
#' @param max_samples if given, keep only the first this many samples
#'
#' @return A panel list suitable for [am_mosaic()], including chromosome, ID,
#'   reference allele, and alternate allele metadata for each retained marker.
#' @export
#'
#' @examples
#' ## a tiny phased VCF written on the fly
#' vcf <- tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
#'           "FORMAT", "s1", "s2", "s3", "s4"), collapse = "\t"),
#'   paste(c("22", "100", ".", "A", "G", ".", ".", ".", "GT",
#'           "0|1", "1|1", "0|0", "1|0"), collapse = "\t"),
#'   paste(c("22", "200", ".", "C", "T", ".", ".", ".", "GT",
#'           "1|0", "0|1", "1|1", "0|0"), collapse = "\t")), vcf)
#' panel <- vcf_to_panel(vcf, min_maf = 0)
#' dim(panel$haplotypes)
vcf_to_panel <- function(vcf, map = NULL, min_maf = 0.05,
                         max_markers = NULL, max_samples = NULL) {
  if (!is.character(vcf) || length(vcf) != 1L || is.na(vcf) || !nzchar(vcf)) {
    stop("`vcf` must be a single non-empty path")
  }
  if (!file.exists(vcf)) stop("VCF not found: ", vcf)
  if (!is.null(map) && (!is.character(map) || length(map) != 1L ||
      is.na(map) || !nzchar(map))) {
    stop("`map` must be NULL or a single non-empty path")
  }
  if (!is.numeric(min_maf) || length(min_maf) != 1L || is.na(min_maf) ||
      !is.finite(min_maf) || min_maf < 0 || min_maf > 0.5) {
    stop("`min_maf` must be a single number in the interval [0, 0.5]")
  }
  positive_integer <- function(x, name) {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 1 || x %% 1 != 0 || x > .Machine$integer.max) {
      stop("`", name, "` must be a positive whole number")
    }
    as.integer(x)
  }
  max_markers <- positive_integer(max_markers, "max_markers")
  max_samples <- positive_integer(max_samples, "max_samples")

  con <- if (grepl("\\.gz$", vcf)) gzfile(vcf, "rt") else file(vcf, "rt")
  on.exit(close(con))
  lines <- readLines(con)
  lines <- lines[!startsWith(lines, "##")]
  if (!length(lines) || !startsWith(lines[1], "#CHROM")) {
    stop("no #CHROM header found in ", vcf)
  }
  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  body <- lines[-1]
  if (!length(body)) stop("no variant records found in ", vcf)

  nf <- length(header)
  if (nf < 10L) stop("VCF contains no sample genotype columns")
  fixed_header <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT")
  if (!identical(header[seq_along(fixed_header)], fixed_header)) {
    stop("malformed VCF: the first nine #CHROM header fields are invalid")
  }
  samples <- header[-seq_along(fixed_header)]
  if (any(!nzchar(samples)) || anyDuplicated(samples)) {
    stop("malformed VCF: sample names must be non-empty and unique")
  }
  records <- strsplit(body, "\t", fixed = TRUE)
  widths <- lengths(records)
  if (any(widths != nf)) {
    bad_line <- which(widths != nf)[1L] + 1L
    stop("malformed VCF: record on line ", bad_line, " has ",
         widths[bad_line - 1L], " fields; expected ", nf)
  }
  mat <- matrix(unlist(records, use.names = FALSE), nrow = nf)
  chrom <- mat[1, ]
  if (any(!nzchar(chrom) | chrom == "." | grepl("[[:space:]]", chrom))) {
    stop("VCF contains an invalid chromosome field")
  }
  pos_num <- suppressWarnings(as.numeric(mat[2, ]))
  if (anyNA(pos_num) || any(!is.finite(pos_num)) || any(pos_num < 1) ||
      any(pos_num != floor(pos_num)) || any(pos_num > .Machine$integer.max)) {
    stop("VCF contains an invalid marker position")
  }
  pos <- as.integer(pos_num)
  id <- mat[3, ]
  ref <- mat[4, ]
  alt <- mat[5, ]
  ## A one-character symbolic allele such as `*` is not an SNV. Restricting
  ## both alleles to the DNA alphabet also excludes missing and ambiguous
  ## alleles without letting them reach genotype parsing.
  snv <- grepl("^[ACGTacgt]$", ref) & grepl("^[ACGTacgt]$", alt)

  sample_rows <- 10:nf
  if (!is.null(max_samples)) {
    sample_rows <- sample_rows[seq_len(min(length(sample_rows), max_samples))]
  }
  if (any(sub(":.*$", "", mat[9, snv]) != "GT")) {
    stop("the FORMAT field must begin with GT for every retained SNV")
  }
  gt <- sub(":.*$", "", mat[sample_rows, , drop = FALSE])
  valid_gt <- grepl("^[01.][|/][01.]$", gt[, snv, drop = FALSE])
  if (any(!valid_gt)) {
    stop("SNV genotypes must be diploid calls using alleles 0, 1, or ., ",
         "such as 0|1 or 0/1")
  }
  if (any(grepl("/", gt[, snv, drop = FALSE], fixed = TRUE))) {
    warning(
      paste0(
        "unphased VCF genotypes were found. Their supplied allele order is ",
        "being treated as phase; this can create artificial haplotypes and ",
        "local LD. Use a phased VCF for scientifically reliable mosaics."
      ),
      call. = FALSE
    )
  }
  ## Records that are not biallelic SNVs are filtered below. Give them a valid
  ## placeholder here so their unrelated FORMAT fields cannot trigger parsing
  ## warnings while the complete rectangular matrix is converted.
  gt[, !snv] <- ".|."
  h1 <- matrix(suppressWarnings(as.integer(substr(gt, 1, 1))),
               nrow = nrow(gt))
  h2 <- matrix(suppressWarnings(as.integer(substr(gt, 3, 3))),
               nrow = nrow(gt))
  H <- rbind(h1, h2)

  keep <- snv & !apply(is.na(H), 2, any)
  af <- colMeans(H[, keep, drop = FALSE])
  maf <- pmin(af, 1 - af)
  idx <- which(keep)[maf >= min_maf]
  if (!length(idx)) {
    stop("no markers survived the min_maf = ", min_maf, " filter")
  }
  if (!is.null(max_markers) && length(idx) > max_markers) {
    idx <- idx[round(seq(1, length(idx), length.out = max_markers))]
    idx <- unique(idx)
  }
  H <- H[, idx, drop = FALSE]
  chrom <- chrom[idx]
  pos <- pos[idx]
  id <- id[idx]
  ref <- ref[idx]
  alt <- alt[idx]

  chrom_key <- unique(.panel_chrom_key(chrom))
  if (length(chrom_key) != 1L || !nzchar(chrom_key)) {
    stop("a reference panel must contain retained markers from exactly one chromosome")
  }
  if (is.unsorted(pos, strictly = TRUE)) {
    stop("retained VCF marker positions must be strictly increasing")
  }

  panel <- list(
    haplotypes = matrix(as.raw(H), nrow = nrow(H)),
    pos = pos,
    chrom = chrom,
    id = id,
    ref = ref,
    alt = alt
  )
  if (!is.null(map)) {
    if (!file.exists(map)) stop("genetic map not found: ", map)
    mp <- tryCatch(
      utils::read.table(
        map, header = FALSE, col.names = c("chr", "id", "cM", "bp"),
        colClasses = c("character", "character", "numeric", "numeric"),
        comment.char = "", quote = "", fill = FALSE,
        na.strings = character()
      ),
      error = function(e) {
        stop("could not read the four-column genetic map: ", conditionMessage(e),
             call. = FALSE)
      }
    )
    if (ncol(mp) != 4L) stop("the genetic map must have exactly four columns")
    mp_key <- .panel_chrom_key(mp$chr)
    mp <- mp[!is.na(mp_key) & mp_key == chrom_key, , drop = FALSE]
    if (nrow(mp) < 2L) {
      stop("the genetic map must contain at least two rows for chromosome ",
           unique(chrom)[1L])
    }
    panel$cM <- .panel_interpolate_cM(pos, mp$bp, mp$cM)
  }
  panel
}

## Fail before starting a long shell pipeline, and name every missing program
## in one diagnostic rather than reporting a generic download failure later.
.panel_require_commands <- function(commands) {
  found <- Sys.which(commands)
  missing <- commands[!nzchar(found)]
  if (length(missing)) {
    stop(
      "`download_1kg_panel()` requires external command",
      if (length(missing) > 1L) "s" else "",
      ": ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(found)
}

## Stream just enough of a chromosome VCF to prove that the requested end
## coordinate was reached. The destination is published only after that proof
## and a minimal structural check, so an interrupted command can never become
## a seemingly valid cache entry. A non-zero curl/zcat status cannot be used
## directly here: stopping after `end` deliberately closes their pipe early.
.panel_stream_region <- function(vcf_url, start, end, region) {
  .panel_require_commands(c("curl", "zcat", "awk"))
  tmp <- tempfile(paste0(".", basename(region), "-"), dirname(region))
  done <- tempfile(paste0(".", basename(region), "-complete-"), dirname(region))
  on.exit(unlink(c(tmp, done)), add = TRUE)

  awk <- paste0(
    "'/^##/ {print; next} /^#CHROM/ {print; next} ",
    "$2 > end {print \"complete\" > done; close(done); exit} ",
    "$2 >= start {print}'"
  )
  cmd <- sprintf(
    "curl -fsSL %s | zcat | awk -v start=%.0f -v end=%.0f -v done=%s -F'\\t' %s > %s",
    shQuote(vcf_url), start, end, shQuote(done), awk, shQuote(tmp)
  )
  status <- system(cmd)
  complete <- file.exists(done) &&
    identical(readLines(done, n = 1L, warn = FALSE), "complete")
  preview <- if (file.exists(tmp)) {
    readLines(tmp, n = 10000L, warn = FALSE)
  } else {
    character()
  }
  if (status != 0L || !complete) {
    stop(
      "download did not reach a VCF record beyond the requested end; ",
      "the network stream may be incomplete",
      call. = FALSE
    )
  }
  header_at <- which(startsWith(preview, "#CHROM"))[1L]
  if (is.na(header_at) || length(preview) <= header_at) {
    stop("the downloaded region contains no VCF records", call. = FALSE)
  }
  if (file.exists(region)) {
    stop("refusing to overwrite an existing region cache: ", region,
         call. = FALSE)
  }
  if (!file.rename(tmp, region)) {
    stop("could not publish the downloaded region cache: ", region,
         call. = FALSE)
  }
  invisible(region)
}

.panel_zip_listing <- function(path) {
  tryCatch(
    utils::unzip(path, list = TRUE),
    warning = function(w) {
      stop("invalid genetic-map archive: ", conditionMessage(w), call. = FALSE)
    },
    error = function(e) {
      stop("invalid genetic-map archive: ", conditionMessage(e), call. = FALSE)
    }
  )
}

#' Download a 1000 Genomes region and build a reference panel
#'
#' Convenience wrapper that fetches a window of phased 1000 Genomes data along
#' with a genetic map and hands both to [vcf_to_panel()]. It exists so the
#' vignette's workflow can be reproduced at realistic scale; it requires
#' network access and downloads a large file, so it is not run in examples or
#' tests.
#'
#' @param chrom autosome, as a string from `"1"` through `"22"`
#' @param start,end base pair bounds of the region to keep
#' @param dest directory in which to cache downloads. Defaults to a temporary
#'   directory.
#' @param ... further arguments passed to [vcf_to_panel()], such as `min_maf`
#'   and `max_samples`
#'
#' @return A panel list suitable for [am_mosaic()].
#'
#' @details Streams the VCF and stops reading once `end` is passed, so a small
#'   window costs far less than the whole chromosome. The data are the GRCh38
#'   phased biallelic release, and the genetic map is the GRCh38 PLINK map
#'   distributed with BEAGLE. A new cache file is published only after the
#'   stream reaches a variant beyond `end`, preventing an interrupted download
#'   from being reused as a complete region.
#' @export
#'
#' @examples
#' \donttest{
#' ## requires network access
#' if (interactive()) {
#'   panel <- download_1kg_panel("22", 20e6, 20.5e6, max_samples = 200)
#'   dim(panel$haplotypes)
#' }
#' }
download_1kg_panel <- function(chrom = "22", start = 20e6, end = 21e6,
                               dest = tempdir(), ...) {
  if (!is.character(chrom) || length(chrom) != 1L || is.na(chrom) ||
      !grepl("^([1-9]|1[0-9]|2[0-2])$", chrom)) {
    stop("`chrom` must name one autosome from '1' through '22'")
  }
  if (!is.numeric(start) || length(start) != 1L || is.na(start) ||
      !is.finite(start) || !is.numeric(end) || length(end) != 1L ||
      is.na(end) || !is.finite(end) || start < 1 || end < start ||
      start != floor(start) || end != floor(end) ||
      end > .Machine$integer.max) {
    stop(paste0("`start` and `end` must be finite positive whole-number ",
                "bounds with `start <= end`"))
  }
  if (!is.character(dest) || length(dest) != 1L || is.na(dest) ||
      !nzchar(dest)) {
    stop("`dest` must be a single non-empty directory path")
  }
  if (!dir.exists(dest) && !dir.create(dest, recursive = TRUE)) {
    stop("could not create destination directory: ", dest)
  }
  base <- paste0("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/",
                 "1000_genomes_project/release/20181203_biallelic_SNV/")
  vcf_url <- paste0(base, "ALL.chr", chrom,
                    ".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz")
  map_url <- "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"

  start_label <- format(start, scientific = FALSE, trim = TRUE)
  end_label <- format(end, scientific = FALSE, trim = TRUE)
  region <- file.path(
    dest, paste0("chr", chrom, "_", start_label, "_", end_label, ".vcf")
  )
  if (!file.exists(region)) {
    message("streaming ", vcf_url)
    .panel_stream_region(vcf_url, start, end, region)
  }

  map_file <- file.path(dest, paste0("plink.chr", chrom, ".GRCh38.map"))
  if (!file.exists(map_file)) {
    zipf <- file.path(dest, "plink.GRCh38.map.zip")
    if (!file.exists(zipf)) {
      zip_tmp <- tempfile(".plink-map-", dest, fileext = ".zip")
      on.exit(unlink(zip_tmp), add = TRUE)
      status <- utils::download.file(map_url, zip_tmp, mode = "wb", quiet = TRUE)
      if (length(status) != 1L || is.na(status) || status != 0L ||
          !file.exists(zip_tmp) ||
          file.size(zip_tmp) == 0) {
        stop("genetic-map download failed", call. = FALSE)
      }
      .panel_zip_listing(zip_tmp)
      if (!file.rename(zip_tmp, zipf)) {
        stop("could not publish the genetic-map archive cache", call. = FALSE)
      }
    }
    listing <- .panel_zip_listing(zipf)
    inner <- grep(paste0("plink\\.chr", chrom, "\\.GRCh38\\.map$"),
                  listing$Name, value = TRUE)
    inner <- inner[!grepl("chrchr", inner)][1]
    if (is.na(inner)) stop("no genetic map for chromosome ", chrom)
    map_tmp_dir <- tempfile(".plink-map-extract-", dest)
    if (!dir.create(map_tmp_dir)) {
      stop("could not create a temporary map extraction directory")
    }
    on.exit(unlink(map_tmp_dir, recursive = TRUE), add = TRUE)
    utils::unzip(zipf, files = inner, exdir = map_tmp_dir, junkpaths = TRUE)
    extracted <- file.path(map_tmp_dir, basename(inner))
    if (!file.exists(extracted) || file.size(extracted) == 0) {
      stop("the genetic-map archive did not yield a usable map")
    }
    if (!file.rename(extracted, map_file)) {
      stop("could not publish the extracted genetic-map cache")
    }
  }

  vcf_to_panel(region, map = map_file, ...)
}

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
  ord <- order(map_pos)
  stats::approx(map_pos[ord], map_cM[ord], xout = pos, rule = 2)$y
}

#' Build a reference panel from a phased VCF
#'
#' Reads a phased, biallelic VCF and returns a panel in the form
#' [am_mosaic()] expects. The VCF must be phased, meaning genotypes look like
#' `0|1` rather than `0/1`, because the mosaic copies haplotypes rather than
#' genotypes.
#'
#' @param vcf path to a VCF, optionally gzipped
#' @param map optional path to a PLINK-format genetic map with columns
#'   chromosome, marker, centimorgans, and base pair position. Without one, the
#'   panel carries no `cM` and [am_mosaic()] places breakpoints uniformly in
#'   physical distance.
#' @param min_maf drop markers whose minor allele frequency in the panel falls
#'   below this. Defaults to 0.05.
#' @param max_markers if given, thin evenly to at most this many markers
#' @param max_samples if given, keep only the first this many samples
#'
#' @return A panel list suitable for [am_mosaic()].
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
  if (!file.exists(vcf)) stop("VCF not found: ", vcf)
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
  fields <- scan(text = body, what = "", sep = "\t", quiet = TRUE)
  if (length(fields) %% nf != 0L) {
    stop("malformed VCF: records do not all have ", nf, " fields")
  }
  mat <- matrix(fields, nrow = nf)
  pos <- as.integer(mat[2, ])
  ref <- mat[4, ]
  alt <- mat[5, ]

  sample_rows <- 10:nf
  if (!is.null(max_samples)) {
    sample_rows <- sample_rows[seq_len(min(length(sample_rows), max_samples))]
  }
  gt <- mat[sample_rows, , drop = FALSE]
  if (!any(grepl("|", gt[1, ], fixed = TRUE))) {
    stop("genotypes are not phased; am_mosaic() copies haplotypes and needs ",
         "phased data such as 0|1")
  }
  h1 <- matrix(as.integer(substr(gt, 1, 1)), nrow = nrow(gt))
  h2 <- matrix(as.integer(substr(gt, 3, 3)), nrow = nrow(gt))
  H <- rbind(h1, h2)

  keep <- nchar(ref) == 1L & nchar(alt) == 1L & !apply(is.na(H), 2, any)
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
  pos <- pos[idx]

  panel <- list(haplotypes = matrix(as.raw(H), nrow = nrow(H)), pos = pos)
  if (!is.null(map)) {
    if (!file.exists(map)) stop("genetic map not found: ", map)
    mp <- utils::read.table(map, header = FALSE,
                            col.names = c("chr", "id", "cM", "bp"))
    panel$cM <- .panel_interpolate_cM(pos, mp$bp, mp$cM)
  }
  panel
}

#' Download a 1000 Genomes region and build a reference panel
#'
#' Convenience wrapper that fetches a window of phased 1000 Genomes data along
#' with a genetic map and hands both to [vcf_to_panel()]. It exists so the
#' vignette's workflow can be reproduced at realistic scale; it requires
#' network access and downloads a large file, so it is not run in examples or
#' tests.
#'
#' @param chrom chromosome, as a string, for example `"22"`
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
#'   distributed with BEAGLE.
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
  if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)
  base <- paste0("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/",
                 "1000_genomes_project/release/20181203_biallelic_SNV/")
  vcf_url <- paste0(base, "ALL.chr", chrom,
                    ".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz")
  map_url <- "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"

  region <- file.path(dest, paste0("chr", chrom, "_", start, "_", end, ".vcf"))
  if (!file.exists(region)) {
    message("streaming ", vcf_url)
    cmd <- sprintf(
      "curl -sL %s | zcat | awk -F'\\t' '/^##/ {next} /^#CHROM/ {print; next} $2 > %.0f {exit} $2 >= %.0f {print}' > %s",
      shQuote(vcf_url), end, start, shQuote(region))
    if (system(cmd) != 0L || !file.exists(region)) {
      stop("download failed; check network access and that curl, zcat, and ",
           "awk are available")
    }
  }

  map_file <- file.path(dest, paste0("plink.chr", chrom, ".GRCh38.map"))
  if (!file.exists(map_file)) {
    zipf <- file.path(dest, "plink.GRCh38.map.zip")
    if (!file.exists(zipf)) {
      utils::download.file(map_url, zipf, mode = "wb", quiet = TRUE)
    }
    inner <- grep(paste0("plink\\.chr", chrom, "\\.GRCh38\\.map$"),
                  utils::unzip(zipf, list = TRUE)$Name, value = TRUE)
    inner <- inner[!grepl("chrchr", inner)][1]
    if (is.na(inner)) stop("no genetic map for chromosome ", chrom)
    utils::unzip(zipf, files = inner, exdir = dest, junkpaths = TRUE)
    map_file <- file.path(dest, basename(inner))
  }

  vcf_to_panel(region, map = map_file, ...)
}

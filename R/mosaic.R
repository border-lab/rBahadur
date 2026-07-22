## Simulating equilibrium AM with realistic local LD.
##
## Assortative mating induces dense, genome-wide covariance among causal
## variants, while recombination induces banded covariance among physically
## close markers. A single Bahadur order-2 distribution cannot represent both:
## strong local LD pushes the target correlation matrix outside the feasible
## region described in `?am_covariance_structure`.
##
## The way around it, following Algorithm S4 of the rBahadur supplementary
## note, is to split the problem. Causal variants are drawn with rb_dplr(),
## which supplies the global AM structure and is well behaved because causal
## loci are far apart. The intervening markers are then filled in by copying
## contiguous blocks from a real reference panel, which supplies local LD for
## free because the blocks are real human haplotypes. Each block contains
## exactly one causal locus and is copied from a panel haplotype carrying the
## allele already drawn there, so the AM structure survives intact.
##
## Block boundaries are drawn from the genetic map rather than uniformly, so
## breakpoints concentrate where recombination actually happens.
##
## The per-individual state is small: one donor haplotype and one block end per
## causal locus, so O(n * m) rather than O(n * p). That is what lets the
## variant-major layouts stream over markers while the individual-major layout
## streams over people.

## Validate a reference panel and return it in canonical form.
.mosaic_check_panel <- function(panel) {
  if (!is.list(panel) || !all(c("haplotypes", "pos") %in% names(panel))) {
    stop("`panel` must be a list with at least `haplotypes` and `pos`")
  }
  H <- panel$haplotypes
  if (is.raw(H)) H <- matrix(as.integer(H), nrow = nrow(H))
  if (!is.matrix(H)) stop("`panel$haplotypes` must be a matrix")
  storage.mode(H) <- "integer"
  if (anyNA(H) || any(H != 0L & H != 1L)) {
    stop("`panel$haplotypes` must contain only 0 and 1, with no missing values")
  }
  p <- ncol(H)
  if (length(panel$pos) != p) {
    stop("`panel$pos` must have one entry per column of `panel$haplotypes`")
  }
  if (is.unsorted(panel$pos, strictly = TRUE)) {
    stop("`panel$pos` must be strictly increasing")
  }
  cM <- panel$cM
  if (is.null(cM)) {
    ## with no genetic map, fall back to physical distance, which makes the
    ## breakpoint distribution uniform in base pairs
    cM <- (panel$pos - panel$pos[1]) / 1e6
  }
  if (length(cM) != p) stop("`panel$cM` must have one entry per marker")
  if (is.unsorted(cM)) stop("`panel$cM` must be non-decreasing")
  list(H = H, pos = panel$pos, cM = cM, N = nrow(H), p = p)
}

## Resolve and validate the causal marker indices.
.mosaic_causal_idx <- function(causal_idx, m, p) {
  if (is.null(causal_idx)) {
    if (is.null(m)) stop("supply either `causal_idx` or `m`")
    if (m < 2 || m > p) {
      stop("`m` must be between 2 and the number of markers (", p, ")")
    }
    causal_idx <- unique(round(seq(1, p, length.out = m)))
  }
  causal_idx <- as.integer(causal_idx)
  if (any(causal_idx < 1L | causal_idx > p) ||
      is.unsorted(causal_idx, strictly = TRUE)) {
    stop("`causal_idx` must be strictly increasing and within 1:", p)
  }
  if (length(causal_idx) < 2L) {
    stop("at least two causal variants are required")
  }
  causal_idx
}

## Precompute, for each causal locus, which panel haplotypes carry each allele.
.mosaic_match_lists <- function(panel, causal_idx) {
  lapply(causal_idx, function(j) {
    col <- panel$H[, j]
    list(which(col == 0L), which(col == 1L))
  })
}

## Draw `n` block boundaries in each gap between consecutive causal loci.
##
## Within a gap the boundary falls between markers j and j+1 with probability
## proportional to the recombination distance cM[j+1] - cM[j], which is what
## makes breakpoints cluster in hotspots. Returns an m by n matrix of block end
## indices, where row k holds the last marker belonging to block k.
.mosaic_boundaries <- function(causal_idx, cM, p, n) {
  m <- length(causal_idx)
  ends <- matrix(0L, nrow = m, ncol = n)
  ends[m, ] <- p                      # the final block runs to the last marker
  for (k in seq_len(m - 1L)) {
    lo <- causal_idx[k]               # boundary falls at or after this marker
    hi <- causal_idx[k + 1L] - 1L     # and strictly before the next causal one
    if (hi <= lo) {
      ends[k, ] <- lo                 # adjacent causal loci leave no room
      next
    }
    w <- diff(cM[lo:(hi + 1L)])
    if (anyNA(w) || all(w <= 0)) w <- rep(1, length(w))
    ends[k, ] <- sample(lo:hi, n, replace = TRUE, prob = w)
  }
  ends
}

## Choose a donor panel haplotype per individual and block, matching the drawn
## causal allele. `alleles` is n by m of 0/1.
.mosaic_donors <- function(alleles, match_lists, causal_idx) {
  n <- nrow(alleles)
  m <- ncol(alleles)
  donor <- matrix(0L, nrow = n, ncol = m)
  for (k in seq_len(m)) {
    for (allele in 0:1) {
      who <- which(alleles[, k] == allele)
      if (!length(who)) next
      cand <- match_lists[[k]][[allele + 1L]]
      if (!length(cand)) {
        stop("no reference haplotype carries allele ", allele,
             " at causal marker ", causal_idx[k],
             "; causal markers must be polymorphic in the panel")
      }
      donor[who, k] <- cand[sample.int(length(cand), length(who),
                                       replace = TRUE)]
    }
  }
  donor
}

## Expand one individual's block donors into a full length-p haplotype.
.mosaic_one <- function(donor_i, ends_i, H, p) {
  starts <- c(1L, ends_i[-length(ends_i)] + 1L)
  lens <- ends_i - starts + 1L
  H[cbind(rep.int(donor_i, lens), seq_len(p))]
}

## Gather markers `cols` for every individual, advancing a per-individual block
## pointer. `cols` must be increasing, and `state` carries the pointer between
## calls so successive chunks stay O(n) per marker.
.mosaic_gather_cols <- function(cols, donor, ends, H, state) {
  n <- nrow(donor)
  m <- ncol(donor)
  out <- matrix(0L, nrow = n, ncol = length(cols))
  cur <- state$cur
  rows <- seq_len(n)
  for (ci in seq_along(cols)) {
    j <- cols[ci]
    repeat {
      behind <- which(ends[cbind(cur, rows)] < j)
      if (!length(behind)) break
      cur[behind] <- pmin(cur[behind] + 1L, m)
    }
    out[, ci] <- H[cbind(donor[cbind(rows, cur)], j)]
  }
  list(values = out, state = list(cur = cur))
}

#' Simulate equilibrium assortative mating with realistic local LD
#'
#' Combines the global linkage disequilibrium induced by assortative mating
#' with the local linkage disequilibrium induced by limited recombination.
#' Causal variants are drawn with [rb_dplr()], giving the dense genome-wide
#' structure assortative mating produces, and the remaining markers are filled
#' in by copying contiguous haplotype blocks from a reference panel, giving
#' realistic short-range structure. Block boundaries are sampled from the
#' panel's genetic map, so breakpoints concentrate where recombination is high.
#'
#' @param h2_0 generation zero (panmictic) heritability
#' @param r cross-mate phenotypic correlation, in the open interval (-1, 1).
#'   Negative values give disassortative mating.
#' @param n number of individuals to simulate
#' @param panel reference panel: a list with `haplotypes` (a haplotypes by
#'   markers matrix of 0 and 1, optionally stored as `raw`), `pos` (strictly
#'   increasing base pair positions), and optionally `cM` (genetic map position
#'   of each marker). Without `cM`, breakpoints are drawn uniformly in physical
#'   distance. See [kg_reference()] for the bundled example.
#' @param causal_idx integer indices of the markers to treat as causal. If
#'   `NULL`, `m` evenly spaced markers are used.
#' @param m number of causal variants, used only when `causal_idx` is `NULL`
#' @param path,format,batch_size streaming options, exactly as in
#'   [am_simulate()]. With `path = NULL` the genotype matrix is returned in
#'   memory; otherwise it is streamed to disk and omitted from the result.
#'
#' @return A list with `y`, `g`, `AF` (allele frequencies at the causal loci),
#'   `beta_std`, `beta_raw`, `causal_idx`, and `pos`. With `path = NULL` it also
#'   carries `X`, an `n` by `p` integer matrix of diploid genotypes at every
#'   panel marker. With `path` supplied, `X` is omitted and `path`, `format`,
#'   `n`, `m`, `h2_0`, and `r` are added, where `m` counts all panel markers
#'   written rather than only the causal ones.
#'
#' @details Genotypes at the causal loci are exactly what [rb_dplr()] drew, so
#'   the equilibrium relationships in [h2_eq()] and [vg_eq()] hold there, while
#'   surrounding markers inherit the panel's correlation structure. Combining
#'   the two rather than approximating them jointly sidesteps the feasibility
#'   limit in [am_covariance_structure()]: representing strong local LD and
#'   genome-wide assortative mating in a single Bahadur order-2 distribution is
#'   generally not possible.
#'
#'   Each block is copied from one panel haplotype, so within a block the
#'   simulated data reproduces panel LD exactly, while correlation across a
#'   breakpoint is broken apart from what the causal variants carry. More
#'   causal variants therefore means more, shorter blocks.
#'
#'   Because the panel is finite, the simulated data cannot contain haplotypes
#'   the panel does not, and a small panel will show inflated identity by
#'   descent between simulated individuals.
#'
#' @seealso [am_simulate()] for the unlinked-loci case, and [kg_reference()]
#'   for the bundled 1000 Genomes panel.
#' @export
#'
#' @examples
#' panel <- kg_reference()
#' sim <- am_mosaic(h2_0 = 0.5, r = 0.4, n = 50, panel = panel, m = 20)
#' dim(sim$X)
#'
#' ## neighbouring markers are correlated because they are copied together,
#' ## which is the local LD that am_simulate() cannot produce
#' j <- sim$causal_idx[10]
#' cor(sim$X[, j], sim$X[, j + 1])
am_mosaic <- function(h2_0, r, n, panel, causal_idx = NULL, m = NULL,
                      path = NULL,
                      format = c("individual", "variant", "bed"),
                      batch_size = NULL) {
  format <- match.arg(format)
  panel <- .mosaic_check_panel(panel)
  p <- panel$p
  causal_idx <- .mosaic_causal_idx(causal_idx, m, p)
  m <- length(causal_idx)

  AF <- colMeans(panel$H[, causal_idx, drop = FALSE])
  if (any(AF <= 0 | AF >= 1)) {
    stop("every causal marker must be polymorphic in the panel; ",
         sum(AF <= 0 | AF >= 1), " of ", m, " are not")
  }
  if (!is.null(batch_size)) {
    if (!is.numeric(batch_size) || length(batch_size) != 1L ||
        is.na(batch_size) || batch_size < 1 || batch_size %% 1 != 0) {
      stop("`batch_size` must be a single positive whole number, or NULL")
    }
  }

  beta <- scale(rnorm(m)) * sqrt(h2_0 / m)
  beta_unscaled <- beta / sqrt(2 * AF * (1 - AF))
  U <- am_covariance_structure(beta, AF, r)
  match_lists <- .mosaic_match_lists(panel, causal_idx)

  ## Draw the causal variants for all n individuals and reduce them to the
  ## compact per-individual state: a donor haplotype and a block end per causal
  ## locus, for each of the two parental copies.
  H <- rb_dplr(n, rep(AF, each = 2), U)
  copies <- lapply(1:2, function(cp) {
    alleles <- H[, seq(cp, 2 * m, 2), drop = FALSE]
    list(donor = .mosaic_donors(alleles, match_lists, causal_idx),
         ends = .mosaic_boundaries(causal_idx, panel$cM, p, n),
         alleles = alleles)
  })

  ## genetic values come from the causal loci only, so they need no expansion
  causal_dosage <- copies[[1]]$alleles + copies[[2]]$alleles
  g <- matrix(as.vector(causal_dosage %*% beta_unscaled), ncol = 1)
  y <- g + rnorm(n, 0, sqrt(1 - h2_0))

  base <- list(y = y, g = g, AF = AF, beta_std = beta,
               beta_raw = beta_unscaled, causal_idx = causal_idx,
               pos = panel$pos)

  if (is.null(path)) {
    X <- matrix(0L, nrow = n, ncol = p)
    for (i in seq_len(n)) {
      X[i, ] <- .mosaic_one(copies[[1]]$donor[i, ], copies[[1]]$ends[, i],
                            panel$H, p) +
                .mosaic_one(copies[[2]]$donor[i, ], copies[[2]]$ends[, i],
                            panel$H, p)
    }
    base$X <- X
    return(base)
  }

  con <- file(.gt_data_path(path, format), "wb")
  on.exit(close(con))

  if (format == "individual") {
    ## batch over people; each individual's row is written contiguously
    if (is.null(batch_size)) {
      batch_size <- max(1L, min(n, as.integer(floor(64e6 / (4 * p)))))
    }
    batch_size <- as.integer(min(batch_size, n))
    for (start in seq(1L, n, by = batch_size)) {
      nb <- min(batch_size, n - start + 1L)
      Xb <- matrix(0L, nrow = nb, ncol = p)
      for (b in seq_len(nb)) {
        i <- start + b - 1L
        Xb[b, ] <- .mosaic_one(copies[[1]]$donor[i, ], copies[[1]]$ends[, i],
                               panel$H, p) +
                   .mosaic_one(copies[[2]]$donor[i, ], copies[[2]]$ends[, i],
                               panel$H, p)
      }
      writeBin(as.vector(t(Xb)), con, size = 1L)
    }
  } else {
    ## variant-major layouts batch over markers, walking a per-individual block
    ## pointer so no full matrix is ever held
    if (is.null(batch_size)) {
      batch_size <- max(1L, min(p, as.integer(floor(64e6 / (4 * n)))))
    }
    batch_size <- as.integer(min(batch_size, p))
    if (format == "bed") writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
    state <- lapply(1:2, function(cp) list(cur = rep(1L, n)))
    for (start in seq(1L, p, by = batch_size)) {
      cols <- start:min(start + batch_size - 1L, p)
      Xb <- matrix(0L, nrow = n, ncol = length(cols))
      for (cp in 1:2) {
        got <- .mosaic_gather_cols(cols, copies[[cp]]$donor, copies[[cp]]$ends,
                                   panel$H, state[[cp]])
        state[[cp]] <- got$state
        Xb <- Xb + got$values
      }
      if (format == "variant") {
        writeBin(as.vector(Xb), con, size = 1L)
      } else {
        for (jj in seq_along(cols)) writeBin(.gt_pack_bed(Xb[, jj]), con)
      }
    }
  }

  .gt_write_meta(path, n, p, format)
  if (format == "bed") .gt_write_plink_sidecars(path, n, p)
  out <- c(base, list(path = path, format = format, n = as.integer(n),
                      m = as.integer(p), h2_0 = h2_0, r = r))
  saveRDS(out, paste0(path, ".rds"))
  invisible(out)
}

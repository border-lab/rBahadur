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
##
## Expanding that state into genotypes is where the time goes, and the two
## layouts want opposite things from memory, so there are two kernels rather
## than one.
##
## Layouts that emit whole markers at a time (in memory, variant-major, and
## bed) sweep markers left to right through .mosaic_sweep(), carrying one
## donor pointer per individual. The pointers advance from a schedule compiled
## once by .mosaic_schedule(), which lists the individuals whose block boundary
## falls at each marker, so a marker costs one gather per parental copy and
## nothing proportional to the number of blocks behind it.
##
## The individual-major layout has to emit whole people, so it expands one
## haplotype at a time through .mosaic_hap() against a transposed panel, where
## a block is a contiguous run of memory rather than a stride.
##
## Throughout, the panel stays `raw` at one byte per allele and is converted to
## integer a block of markers at a time, where it is used. A genome-scale panel
## expanded to integer up front is several gigabytes that nothing needs.

## Validate a reference panel and return it in canonical form.
##
## The panel is held as `raw`, one byte per allele, and never expanded to
## integer. A panel at genome scale runs to hundreds of millions of alleles, so
## the four-fold expansion costs more memory than the simulation itself and
## more time than the expansion it was meant to serve. Markers are converted a
## block at a time instead, where they are actually used.
.mosaic_check_panel <- function(panel) {
  if (!is.list(panel) || !all(c("haplotypes", "pos") %in% names(panel))) {
    stop("`panel` must be a list with at least `haplotypes` and `pos`")
  }
  H <- panel$haplotypes
  if (!is.matrix(H)) stop("`panel$haplotypes` must be a matrix")
  if (nrow(H) < 1L || ncol(H) < 2L) {
    stop("`panel$haplotypes` must have at least one haplotype and two markers")
  }
  bad <- "`panel$haplotypes` must contain only 0 and 1, with no missing values"
  if (is.raw(H)) {
    ## raw cannot be NA, so a single range check settles it
    if (any(H > as.raw(1L))) stop(bad)
  } else {
    if ((!is.numeric(H) && !is.logical(H)) || anyNA(H) ||
        any(!is.finite(H)) || any(H != 0 & H != 1)) {
      stop(bad)
    }
    storage.mode(H) <- "integer"
    H <- matrix(as.raw(H), nrow = nrow(H))
  }
  p <- ncol(H)
  if (length(panel$pos) != p) {
    stop("`panel$pos` must have one entry per column of `panel$haplotypes`")
  }
  if (!is.numeric(panel$pos) || anyNA(panel$pos) ||
      any(!is.finite(panel$pos)) || any(panel$pos < 1) ||
      any(panel$pos != floor(panel$pos))) {
    stop("`panel$pos` must contain positive whole-number base-pair positions")
  }
  if (is.unsorted(panel$pos, strictly = TRUE)) {
    stop("`panel$pos` must be strictly increasing")
  }
  cM <- panel$cM
  has_map <- !is.null(cM)
  if (is.null(cM)) {
    ## with no genetic map, fall back to physical distance, which makes the
    ## breakpoint distribution uniform in base pairs
    cM <- (panel$pos - panel$pos[1]) / 1e6
  }
  if (length(cM) != p) stop("`panel$cM` must have one entry per marker")
  if (!is.numeric(cM) || anyNA(cM) || any(!is.finite(cM))) {
    stop("`panel$cM` must contain finite numeric positions")
  }
  if (is.unsorted(cM)) stop("`panel$cM` must be non-decreasing")

  marker_field <- function(name, default) {
    value <- panel[[name]]
    if (is.null(value)) value <- default
    if (length(value) == 1L) value <- rep(value, p)
    if (length(value) != p || anyNA(value)) {
      stop("`panel$", name, "` must have length 1 or one entry per marker")
    }
    as.character(value)
  }

  list(
    H = H,
    pos = panel$pos,
    cM = cM,
    plink_cM = if (has_map) cM else rep(0, p),
    chrom = marker_field("chrom", "1"),
    id = marker_field("id", paste0("v", seq_len(p))),
    ref = marker_field("ref", "G"),
    alt = marker_field("alt", "A"),
    N = nrow(H),
    p = p
  )
}

## Resolve and validate the causal marker indices.
.mosaic_causal_idx <- function(causal_idx, m, p) {
  if (is.null(causal_idx)) {
    if (is.null(m)) stop("supply either `causal_idx` or `m`")
    if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
        m != floor(m)) {
      stop("`m` must be a single whole number")
    }
    if (m < 2 || m > p) {
      stop("`m` must be between 2 and the number of markers (", p, ")")
    }
    m <- as.integer(m)
    causal_idx <- unique(round(seq(1, p, length.out = m)))
  }
  if (!is.numeric(causal_idx) || !length(causal_idx) || anyNA(causal_idx) ||
      any(!is.finite(causal_idx)) || any(causal_idx != floor(causal_idx))) {
    stop("`causal_idx` must contain finite whole-number marker indices")
  }
  if (any(causal_idx < 1 | causal_idx > p)) {
    stop("`causal_idx` must be strictly increasing and within 1:", p)
  }
  causal_idx <- as.integer(causal_idx)
  if (is.unsorted(causal_idx, strictly = TRUE)) {
    stop("`causal_idx` must be strictly increasing and within 1:", p)
  }
  if (length(causal_idx) < 2L) {
    stop("at least two causal variants are required")
  }
  causal_idx
}

## Allele frequency at each causal locus, read one panel column at a time so
## the causal columns are never gathered into a block of their own.
.mosaic_causal_af <- function(H, causal_idx, N) {
  vapply(causal_idx, function(j) sum(as.integer(H[, j])) / N, numeric(1))
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
##
## The carrier lists are rebuilt per locus rather than precomputed for all of
## them at once: holding them all costs one integer per haplotype per causal
## variant, which at genome scale is larger than everything else here put
## together, and rebuilding one panel column is cheap.
.mosaic_donors <- function(alleles, H, causal_idx) {
  n <- nrow(alleles)
  m <- ncol(alleles)
  donor <- matrix(0L, nrow = n, ncol = m)
  for (k in seq_len(m)) {
    col <- as.integer(H[, causal_idx[k]])
    carriers <- list(which(col == 0L), which(col == 1L))
    for (allele in 0:1) {
      who <- which(alleles[, k] == allele)
      if (!length(who)) next
      cand <- carriers[[allele + 1L]]
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

## Compile per-individual block ends into a marker-ordered schedule of donor
## changes: the individuals whose block boundary falls at marker j, and the
## donor each of them switches to there.
##
## This is what lets the sweep below advance the block pointers with one
## vectorised assignment per marker, touching only the individuals that
## actually change, rather than rescanning all n of them at every marker.
## Block ends are strictly increasing within an individual, because block k
## ends before causal locus k+1 and block k+1 ends at or after it, so no
## individual has two changes at the same marker and the order within a marker
## does not matter.
.mosaic_schedule <- function(donor, ends, p) {
  m <- nrow(ends)
  n <- ncol(ends)
  if (m < 2L) {
    return(list(ind = integer(0), don = integer(0), off = integer(p + 1L)))
  }
  ## column-major over (block, individual) in both, so the three line up
  keep <- seq_len(m - 1L)
  at <- as.integer(ends[keep, , drop = FALSE]) + 1L
  ind <- rep.int(seq_len(n), rep.int(m - 1L, n))
  don <- as.integer(t(donor[, keep + 1L, drop = FALSE]))
  o <- sort.list(at, method = "radix")
  list(ind = ind[o], don = don[o],
       off = c(0L, cumsum(tabulate(at, nbins = p))))
}

## Expand markers `cols` into an n by length(cols) dosage block, summing both
## parental copies. `cols` must be a contiguous increasing run continuing where
## the previous call stopped: `state` carries the block pointers across calls,
## so a marker costs one gather per copy however many blocks precede it.
.mosaic_sweep <- function(cols, sched, state, H, N, n) {
  nc <- length(cols)
  out <- matrix(0L, nrow = n, ncol = nc)
  d1 <- state$d1
  d2 <- state$d2
  o1 <- sched[[1L]]$off; i1 <- sched[[1L]]$ind; c1 <- sched[[1L]]$don
  o2 <- sched[[2L]]$off; i2 <- sched[[2L]]$ind; c2 <- sched[[2L]]$don
  ## Convert the panel a block of markers at a time. One bulk raw to integer
  ## pass costs a fraction of one pass per marker, and bounding the block keeps
  ## the converted copy at a few megabytes however many markers were asked for.
  step <- max(1L, as.integer(2e6 %/% max(N, 1L)))
  for (s in seq(1L, nc, by = step)) {
    ci <- s:min(s + step - 1L, nc)
    Hc <- matrix(as.integer(H[, cols[ci], drop = FALSE]), nrow = N)
    for (k in seq_along(ci)) {
      j <- cols[ci[k]]
      a <- o1[j]; b <- o1[j + 1L]
      if (b > a) { e <- (a + 1L):b; d1[i1[e]] <- c1[e] }
      a <- o2[j]; b <- o2[j + 1L]
      if (b > a) { e <- (a + 1L):b; d2[i2[e]] <- c2[e] }
      hj <- Hc[, k]
      out[, ci[k]] <- hj[d1] + hj[d2]
    }
  }
  list(values = out, state = list(d1 = d1, d2 = d2))
}

## Expand one individual's block donors into a full length-p haplotype.
##
## `Ht` is the panel transposed, markers by haplotypes, so each block is a
## contiguous run of memory and the whole haplotype is one sequential read.
## Against the untransposed panel the same gather strides by the number of
## haplotypes and misses cache on every marker.
.mosaic_hap <- function(donor_i, ends_i, Ht, p, wide) {
  lens <- diff(c(0L, ends_i))
  base <- if (wide) (as.double(donor_i) - 1) * p else (donor_i - 1L) * p
  as.integer(Ht[rep.int(base, lens) + seq_len(p)])
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
#' @param h2_0 generation zero (panmictic) heritability, in `(0, 1)`
#' @param r cross-mate phenotypic correlation, in the open interval (-1, 1).
#'   Negative values give disassortative mating.
#' @param n positive whole-number count of individuals to simulate
#' @param panel reference panel: a list with `haplotypes` (a haplotypes by
#'   markers matrix of 0 and 1, optionally stored as `raw`), `pos` (strictly
#'   increasing base pair positions), and optionally `cM` (genetic map position
#'   of each marker), `chrom`, `id`, `ref`, and `alt`. Marker metadata can have
#'   length one or one entry per marker and is carried into PLINK `.bim` output.
#'   Haplotype value 0 denotes the reference allele and value 1 the alternate
#'   allele. Without `cM`, breakpoints are drawn uniformly in physical distance.
#'   See [kg_reference()] for the bundled example.
#' @param causal_idx integer indices of the markers to treat as causal. If
#'   `NULL`, `m` evenly spaced markers are used.
#' @param m whole-number count of causal variants, used only when `causal_idx`
#'   is `NULL`
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
#'   The equilibrium formulas are large-locus results. Calls with fewer than 50
#'   causal variants are allowed but warn because realized variances can differ
#'   materially from the targets. Set
#'   `options(rBahadur.warn_small_m = FALSE)` to silence this warning after
#'   deciding that the finite-locus approximation is appropriate.
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
#' sim <- am_mosaic(h2_0 = 0.5, r = 0.4, n = 50, panel = panel, m = 50)
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
  .am_check_simulation_args(h2_0, r, n)
  n <- as.integer(n)
  if (!is.null(path)) .gt_check_path(path, write = TRUE)
  panel <- .mosaic_check_panel(panel)
  if (!is.null(path) && format == "bed") {
    .gt_prepare_plink_sidecars(
      panel$p, panel$chrom, panel$id, panel$plink_cM, panel$pos,
      panel$alt, panel$ref
    )
  }
  p <- panel$p
  causal_idx <- .mosaic_causal_idx(causal_idx, m, p)
  m <- length(causal_idx)
  .am_warn_small_m(m, "am_mosaic")

  AF <- .mosaic_causal_af(panel$H, causal_idx, panel$N)
  if (any(AF <= 0 | AF >= 1)) {
    stop("every causal marker must be polymorphic in the panel; ",
         sum(AF <= 0 | AF >= 1), " of ", m, " are not")
  }
  if (!is.null(batch_size)) {
    if (!is.numeric(batch_size) || length(batch_size) != 1L ||
        is.na(batch_size) || !is.finite(batch_size) || batch_size < 1 ||
        batch_size %% 1 != 0 || batch_size > .Machine$integer.max) {
      stop("`batch_size` must be a single positive whole number, or NULL")
    }
  }

  beta <- scale(rnorm(m)) * sqrt(h2_0 / m)
  beta_unscaled <- beta / sqrt(2 * AF * (1 - AF))
  U <- am_covariance_structure(beta, AF, r)

  ## Draw the causal variants for all n individuals and reduce them to the
  ## compact per-individual state: a donor haplotype and a block end per causal
  ## locus, for each of the two parental copies.
  H <- rb_dplr(n, rep(AF, each = 2), U)
  copies <- lapply(1:2, function(cp) {
    alleles <- H[, seq(cp, 2 * m, 2), drop = FALSE]
    list(donor = .mosaic_donors(alleles, panel$H, causal_idx),
         ends = .mosaic_boundaries(causal_idx, panel$cM, p, n),
         alleles = alleles)
  })
  rm(H)

  ## genetic values come from the causal loci only, so they need no expansion
  causal_dosage <- copies[[1]]$alleles + copies[[2]]$alleles
  g <- matrix(as.vector(causal_dosage %*% beta_unscaled), ncol = 1)
  y <- g + rnorm(n, 0, sqrt(1 - h2_0))
  rm(causal_dosage)
  copies[[1]]$alleles <- NULL
  copies[[2]]$alleles <- NULL

  base <- list(y = y, g = g, AF = AF, beta_std = beta,
               beta_raw = beta_unscaled, causal_idx = causal_idx,
               pos = panel$pos)

  ## Every layout except individual-major emits whole markers at a time, so it
  ## goes through the marker-ordered sweep. Compiling the schedule lets the
  ## block ends and donors go, leaving only the two pointer vectors alive.
  if (is.null(path) || format != "individual") {
    sched <- lapply(copies, function(cp) .mosaic_schedule(cp$donor, cp$ends, p))
    state <- list(d1 = copies[[1]]$donor[, 1L], d2 = copies[[2]]$donor[, 1L])
    copies <- NULL
  }

  if (is.null(path)) {
    ## one sweep over every marker, so the block it fills is X itself rather
    ## than a chunk that then has to be copied into it
    base$X <- .mosaic_sweep(seq_len(p), sched, state, panel$H, panel$N, n)$values
    return(base)
  }

  con <- file(.gt_data_path(path, format), "wb")
  on.exit(close(con))

  if (format == "individual") {
    ## batch over people; each individual's row is written contiguously, so the
    ## buffer is filled in place and never transposed
    if (is.null(batch_size)) {
      batch_size <- max(1L, min(n, as.integer(floor(64e6 / (4 * p)))))
    }
    batch_size <- as.integer(min(batch_size, n))
    Ht <- t(panel$H)
    wide <- as.double(p) * panel$N > .Machine$integer.max
    d1 <- copies[[1]]$donor; e1 <- copies[[1]]$ends
    d2 <- copies[[2]]$donor; e2 <- copies[[2]]$ends
    copies <- NULL
    for (start in seq(1L, n, by = batch_size)) {
      nb <- min(batch_size, n - start + 1L)
      buf <- integer(nb * p)
      for (b in seq_len(nb)) {
        i <- start + b - 1L
        buf[((b - 1L) * p + 1L):(b * p)] <-
          .mosaic_hap(d1[i, ], e1[, i], Ht, p, wide) +
          .mosaic_hap(d2[i, ], e2[, i], Ht, p, wide)
      }
      writeBin(buf, con, size = 1L)
    }
  } else {
    ## variant-major layouts batch over markers, walking a per-individual block
    ## pointer so no full matrix is ever held
    if (is.null(batch_size)) {
      batch_size <- max(1L, min(p, as.integer(floor(64e6 / (4 * n)))))
    }
    batch_size <- as.integer(min(batch_size, p))
    if (format == "bed") writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
    for (start in seq(1L, p, by = batch_size)) {
      cols <- start:min(start + batch_size - 1L, p)
      got <- .mosaic_sweep(cols, sched, state, panel$H, panel$N, n)
      state <- got$state
      if (format == "variant") {
        writeBin(as.vector(got$values), con, size = 1L)
      } else {
        writeBin(.gt_pack_bed(got$values), con)
      }
    }
  }

  .gt_write_meta(path, n, p, format)
  if (format == "bed") {
    .gt_write_plink_sidecars(
      path, n, p,
      chrom = panel$chrom,
      id = panel$id,
      cM = panel$plink_cM,
      pos = panel$pos,
      a1 = panel$alt,
      a2 = panel$ref
    )
  }
  out <- c(base, list(path = path, format = format, n = as.integer(n),
                      m = as.integer(p), h2_0 = h2_0, r = r))
  saveRDS(out, paste0(path, ".rds"))
  invisible(out)
}

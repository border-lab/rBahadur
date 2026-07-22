## Locus-blocked form of the rb_dplr recursion.
##
## Emits column blocks through `callback` instead of allocating the full
## n-by-M matrix. Because R fills matrices column-major, drawing runif(n) once
## per column consumes the random stream in the same order that
## matrix(runif(M*n), n, M) does, so this is bit-identical to rb_dplr() at any
## block size.
##
## @param callback invoked as callback(B, col0); B is an n-by-k 0/1 matrix and
##   col0 is the 1-based global index of its first column.
.rb_dplr_stream <- function(n, mu, U, s = 1, block = 1024L, callback) {
  M <- length(mu)
  block <- max(1L, min(as.integer(block), M))
  buf <- matrix(0L, nrow = n, ncol = block)
  bi <- 0L
  col0 <- 1L
  x <- NULL
  cc <- 1

  for (m in seq_len(M)) {
    p <- if (m == 1L) rep(mu[1], n) else mu[m] + s * x * U[m]
    if (any(!is.finite(p) | p < 0 | p > 1)) {
      stop(.rb_infeasible_msg(m))
    }
    km <- as.integer(runif(n) <= p)
    bi <- bi + 1L
    buf[, bi] <- km

    tmp_bool <- (km == 0L)
    pc <- tmp_bool * (1 - p) + (!tmp_bool) * p
    Bk0 <- tmp_bool * (1 - mu[m]) + (!tmp_bool) * mu[m]
    Bk1 <- tmp_bool * (-1) + (!tmp_bool) * 1
    if (m == 1L) {
      x <- Bk1 * U[1] / pc
      cc <- 1
    } else {
      x <- (x * Bk0 + cc * Bk1 * U[m]) / pc
      cc <- (Bk0 / pc) * cc
    }

    if (bi == block || m == M) {
      callback(buf[, seq_len(bi), drop = FALSE], col0)
      col0 <- m + 1L
      bi <- 0L
    }
  }
  invisible(NULL)
}

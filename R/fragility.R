# Neural fragility (Li et al. 2021) per recording, with EZFragility.
#
# For each window, EZFragility fits a linear model x(t+1) = A x(t) by ridge regression,
# raising the ridge penalty until A is stable (ridgeSearch), then finds for each channel i
# the smallest perturbation of row i of A that puts an eigenvalue on the unit circle
# (fragilityRow); fragility = (max - norm) / max over the channels of that window, so the
# most fragile channel of each window is 1 and the least is 0. ENGINE=stable (default)
# computes that norm without the cancellation of EZFragility's and v0.1's formula; fast is
# v0.1.x, ezfragility is EZFragility's calcAdjFrag.

suppressPackageStartupMessages({
  library(EZFragility)
  library(gsignal)
})

# The trial eigenvalues of EZFragility's search: lambda = sqrt(1 - w^2) + i w, w in (0, 1].
trial_lambdas <- function(nSearch) {
  om <- seq(0, 1, length.out = nSearch + 1L)[-1L]
  complex(real = sqrt(1 - om^2), imaginary = om)
}

# (max - norm) / max over the channels of a window; NA where no trial eigenvalue gave a norm.
normalise_norms <- function(best) {
  best[!is.finite(best)] <- NA
  m <- max(best, na.rm = TRUE)
  if (!is.finite(m)) return(rep(NA_real_, length(best)))
  (m - best) / m
}

# ---- ENGINE=fast (v0.1.x) ---------------------------------------------------------------
# fragilityRow, vectorised over the channels. EZFragility loops in R over every channel and
# every trial eigenvalue lambda, computing for row r of (A - lambda I)^-1 the minimum-norm
# Delta = B' (B B')^-1 b with B = [Im r; Re r], b = (0, -1). Its norm is
# sqrt(b' (B B')^-1 b) = sqrt(|Im r|^2 / (|Im r|^2 |Re r|^2 - (Im r . Re r)^2)), so all rows come
# from one inverse per trial eigenvalue. Same result as fragilityRow to rounding, but that Gram
# determinant cancels badly when A is small (see fragility_norms_stable): kept to reproduce
# v0.1.x, not the default.
fragility_row_fast <- function(A, nSearch = 100L) {
  n <- nrow(A)
  I <- diag(n)
  best <- rep(Inf, n)
  for (lam in trial_lambdas(nSearch)) {
    M <- solve(A - lam * I)
    re <- Re(M); im <- Im(M)
    ii <- rowSums(im * im); rr <- rowSums(re * re); ir <- rowSums(im * re)
    # When A - lambda I is nearly singular (a model that only just passed the stability
    # check) the Gram determinant can cancel to <= 0 in floating point: skip that lambda
    # for that channel rather than produce NaN (EZFragility's loop gives an arbitrary
    # value there, or fails in solve()).
    d <- ii * rr - ir * ir
    v <- rep(NA_real_, n)
    ok <- is.finite(d) & d > 0
    v[ok] <- sqrt(ii[ok] / d[ok])
    best <- pmin(best, v, na.rm = TRUE)
  }
  normalise_norms(best)
}

# ---- ENGINE=stable (default from v0.2.0) -------------------------------------------------
# The same minimum-norm perturbation, computed without cancellation.
#
# Row i of (A - lambda I)^-1 is r = (N_i - e_i) / lambda with N = (A - lambda I)^-1 A: when A is
# small (all eigenvalues ~0.01-0.1, typical with more channels than samples per window), r is
# e_i/lambda plus a small part, Re r and Im r are nearly parallel, and the Gram determinant
# |Im r|^2 |Re r|^2 - (Im r . Re r)^2 of the fast engine loses up to 7 digits. Instead:
#  - N is computed directly (one complex solve with A as right-hand side), so its small entries
#    carry full relative accuracy, rather than as the difference of the inverse and -I/lambda.
#  - Each row is rotated by the unit complex number that makes its diagonal entry real,
#    r' = r conj(r_ii)/|r_ii|. A rotation of (Re r, Im r) leaves their Gram determinant unchanged,
#    and since Im r'_ii = 0 the determinant splits exactly into two non-negative terms over the
#    off-diagonal part w = r'_{j != i}:
#        det = |r_ii|^2 |Im w|^2 + (|Re w|^2 |Im w|^2 - (Re w . Im w)^2)
#    The first dominates and is a product of sums of squares; the second (which still cancels)
#    is smaller by |w|^2/|r_ii|^2, so its rounding no longer matters.
#  - |Delta|^2 = |Im r|^2 / det, with |Im r|^2 a sum of squares.
# Checked against 256-bit ball arithmetic (README "Validation", tests/testthat/test-fragility.R):
# relative error <= ~1e-13 where the fast engine's reaches 8e-3.
#
# Neumann series: when ||A||_2 = q <= NEUMANN_Q, N = A (A - lambda I)^-1 = -sum_k lambda^-k A^k,
# truncated where the tail (<= q^(K+1) / (1 - q)) is below NEUMANN_TOL ||A||: K real matrix
# products per window replace the 100 complex solves (as accurate, ~1.5-2x faster for those
# windows). Larger models use the solve.
NEUMANN_Q <- 0.7   # K <= 120 terms; beyond that the solve is as fast
NEUMANN_TOL <- 2^-60
NEUMANN_CHUNK <- 25L   # trial eigenvalues per block (bounds memory to 2 n^2 x 25 doubles)

# Norms for every channel from N = (A - lambda I)^-1 A = Nr + i Ni at one trial eigenvalue.
# Real arithmetic, one pass per quantity (this runs 100 times per window).
norms_from_N <- function(Nr, Ni, lam) {
  n <- nrow(Nr)
  dg <- seq.int(1L, n * n, by = n + 1L)
  # diagonal of (A - lambda I)^-1 = (N - I) / lambda, then N's off-diagonal part only
  il <- 1 / lam
  z <- complex(real = Nr[dg] - 1, imaginary = Ni[dg]) * il
  Nr[dg] <- 0; Ni[dg] <- 0
  # |Im r|^2: off-diagonal part of Im((N - I) / lambda), plus the diagonal entry
  pr <- Re(il); pi <- Im(il)
  ri <- Nr * pi + Ni * pr
  ii <- .rowSums(ri * ri, n, n) + Im(z)^2
  # off-diagonal part of the rotated row, r * conj(z)/|z| = N * g with g = conj(z) / (|z| lambda)
  mz <- Mod(z)
  g <- Conj(z) / mz * il
  gr <- Re(g); gi <- Im(g)
  a <- Nr * gr - Ni * gi               # gr, gi recycle down the columns: row i scaled by g_i
  b <- Nr * gi + Ni * gr
  aa <- .rowSums(a * a, n, n); bb <- .rowSums(b * b, n, n); ab <- .rowSums(a * b, n, n)
  d <- mz^2 * bb + pmax(aa * bb - ab * ab, 0)
  v <- sqrt(ii / d)
  # d = 0: Im r = 0 after rotation, no real perturbation of row i reaches lambda (EZFragility's
  # B B' is singular there); a non-finite N: A - lambda I singular. Skip that lambda.
  v[!is.finite(v) | !(d > 0)] <- NA
  v
}

neumann_terms <- function(A) {
  q <- norm(A, "2")
  if (!is.finite(q) || q > NEUMANN_Q || q == 0) return(NA_integer_)
  max(1L, as.integer(ceiling(log(NEUMANN_TOL * (1 - q)) / log(q))))
}

# channels x nSearch matrix of the minimum perturbation norms (NA where undefined).
# series = NA: the Neumann series when it applies; TRUE/FALSE force it on (if it applies) / off.
fragility_norms_stable <- function(A, nSearch = 100L, series = NA) {
  n <- nrow(A)
  lam <- trial_lambdas(nSearch)
  out <- matrix(NA_real_, n, nSearch)
  K <- if (isFALSE(series)) NA_integer_ else neumann_terms(A)
  if (is.na(K)) {
    I <- diag(n)
    for (k in seq_along(lam)) {
      N <- solve(A - lam[k] * I, A)
      out[, k] <- norms_from_N(Re(N), Im(N), lam[k])
    }
    attr(out, "terms") <- 0L
    return(out)
  }
  P <- matrix(0, n * n, K)                  # columns: A, A^2, ..., A^K
  Ak <- A
  P[, 1] <- Ak
  for (k in seq_len(K)[-1]) { Ak <- Ak %*% A; P[, k] <- Ak }
  for (blk in split(seq_along(lam), ceiling(seq_along(lam) / NEUMANN_CHUNK))) {
    C <- -outer(seq_len(K), lam[blk], function(k, l) l^(-k))
    Nr <- P %*% Re(C); Ni <- P %*% Im(C)
    for (j in seq_along(blk))
      out[, blk[j]] <- norms_from_N(matrix(Nr[, j], n), matrix(Ni[, j], n), lam[blk[j]])
  }
  attr(out, "terms") <- K
  out
}

fragility_row_stable <- function(A, nSearch = 100L, series = NA) {
  v <- fragility_norms_stable(A, nSearch, series)
  best <- suppressWarnings(apply(v, 1, min, na.rm = TRUE))
  structure(normalise_norms(best), terms = attr(v, "terms"))
}

ez_ridge_search <- utils::getFromNamespace("ridgeSearch", "EZFragility")
ez_standardize <- utils::getFromNamespace("standardizeIEEG", "EZFragility")

# x: channels x samples (rows named). Returns list(frag = channels x windows, start = window
# start in samples (0-based), lambda, stable, terms = Neumann terms used per window, 0 = solve).
fragility_matrix <- function(x, win, step, lambda = NULL, nSearch = 100L,
                             engine = "stable", threads = 1L) {
  nW <- floor((ncol(x) - win) / step) + 1L
  if (nW < 1) stop(sprintf("segment of %d samples is shorter than one window (%d)", ncol(x), win))
  if (engine == "ezfragility") {
    if (threads > 1) {
      doParallel::registerDoParallel(cores = threads)
      on.exit(doParallel::stopImplicitCluster(), add = TRUE)
    }
    ep <- Epoch::Epoch(x, startTime = 0, samplingRate = 1)
    fr <- calcAdjFrag(ep, window = win, step = step, lambda = lambda, nSearch = nSearch,
                      parallel = threads > 1)
    lam <- as.numeric(fr$lambdas)
    st <- rep(NA, nW)   # calcAdjFrag does not return the matrices (unless in debug mode)
    f <- fr$frag
    dimnames(f) <- list(rownames(x), NULL)
    return(list(frag = f, start = (seq_len(nW) - 1L) * step, lambda = lam, stable = st))
  }
  x <- ez_standardize(x)
  one <- function(iw) {
    si <- (iw - 1L) * step + seq_len(win - 1L)
    A <- ez_ridge_search(x[, si, drop = FALSE], x[, si + 1L, drop = FALSE], lambda)
    f <- if (engine == "fast") fragility_row_fast(A, nSearch) else fragility_row_stable(A, nSearch)
    list(f = as.vector(f), lambda = attr(A, "lambda"), terms = attr(f, "terms") %||% NA_integer_,
         stable = isTRUE(attr(A, "stable")) || (!is.null(lambda) &&
                   max(Mod(eigen(A, only.values = TRUE)$values)) < 1))
  }
  res <- if (threads > 1) parallel::mclapply(seq_len(nW), one, mc.cores = threads, mc.preschedule = TRUE)
         else lapply(seq_len(nW), one)
  bad <- vapply(res, function(r) !is.list(r) || is.null(r$f), TRUE)
  if (any(bad)) stop(sprintf("%d window(s) failed: %s", sum(bad),
                             paste(unique(vapply(res[bad], function(r) paste(as.character(r), collapse = " "), "")), collapse = "; ")))
  f <- do.call(cbind, lapply(res, `[[`, "f"))
  dimnames(f) <- list(rownames(x), NULL)
  list(frag = f, start = (seq_len(nW) - 1L) * step,
       lambda = vapply(res, `[[`, 0, "lambda"), stable = vapply(res, `[[`, TRUE, "stable"),
       terms = vapply(res, function(r) as.integer(r$terms), 0L))
}

# Smallest p/q (q <= max_q) with p/q ~= r, for rational resampling.
rational <- function(r, max_q = 1000L) {
  best <- c(round(r), 1L); err <- abs(r - best[1])
  for (q in seq_len(max_q)) {
    p <- round(r * q)
    e <- abs(r - p / q)
    if (e < err - 1e-12) { best <- c(p, q); err <- e }
    if (err < 1e-9) break
  }
  as.integer(best)
}

# Band-pass (4th-order Butterworth, zero phase, in second-order sections) then resample, one channel at a time (so
# memory holds the input and the output, not several filtered copies of the whole matrix).
preprocess <- function(x, fs, bandpass, resample_hz, threads = 1L) {
  notes <- character()
  target <- fs
  pq <- c(1L, 1L)
  if (resample_hz > 0 && fs > resample_hz * (1 + 1e-6)) {
    pq <- rational(resample_hz / fs)
    target <- fs * pq[1] / pq[2]
  }
  bf <- NULL
  if (!is.null(bandpass)) {
    lo <- bandpass[1]; hi <- bandpass[2]
    nyq <- min(fs, target) / 2
    if (hi >= 0.95 * nyq) {
      notes <- c(notes, sprintf("band-pass upper edge %g Hz lowered to %g Hz (Nyquist %g Hz)", hi, 0.9 * nyq, nyq))
      hi <- 0.9 * nyq
    }
    # Second-order sections, not the transfer function (b, a): with a 0.5 Hz edge at 1 kHz the
    # 8th-order (b, a) polynomial loses ~6 digits, so the filtered signal differed by ~2e-6
    # (relative) between amd64 and arm64 and fragility_pre by up to ~1e-4 (v0.1.x, and v0.2.0
    # before this change). In sections the same filter is accurate to ~1e-12 on both.
    bf <- butter(4, c(lo, hi) / (fs / 2), "pass", output = "Sos")
  }
  one <- function(v) {
    if (!is.null(bf)) v <- filtfilt(bf, v)
    if (!identical(pq, c(1L, 1L))) v <- resample(v, pq[1], pq[2])
    v
  }
  first <- one(x[1, ])
  out <- matrix(0, nrow(x), length(first), dimnames = list(rownames(x), NULL))
  out[1, ] <- first
  rest <- seq_len(nrow(x))[-1]
  if (threads > 1 && length(rest) > 1) {
    # channels in forked workers (the filter in sections is ~4x slower than in (b, a) form)
    res <- parallel::mclapply(rest, function(i) one(x[i, ]), mc.cores = threads)
    bad <- vapply(res, function(r) !is.numeric(r), TRUE)
    if (any(bad)) stop(sprintf("preprocessing failed for %d channel(s)", sum(bad)))
    for (k in seq_along(rest)) out[rest[k], ] <- res[[k]]
  } else for (i in rest) out[i, ] <- one(x[i, ])
  list(x = out, fs = target, notes = notes)
}

# One recording -> list(rows = data.frame for fragility.csv, info = provenance entry).
fragility_recording <- function(rec, onset, onset_source, p) {
  t_start <- Sys.time()
  h <- edf_header(rec$edf)
  ch <- recording_channels(rec$edf, h)
  notes <- character()
  if (length(ch$missing)) notes <- c(notes, sprintf("%d channel(s) of channels.tsv not in the EDF: %s",
                                                    length(ch$missing), paste(utils::head(ch$missing, 5), collapse = ", ")))
  if (length(ch$names) < 2) stop("fewer than 2 good channels")
  if (h$discontinuous) notes <- c(notes, "EDF+D (discontinuous): read as continuous")
  w0 <- onset + p$WINDOW_S[1]; w1 <- onset + p$WINDOW_S[2]
  if (w0 < 0) { notes <- c(notes, sprintf("window starts %.1f s before the recording; cut", -w0)); w0 <- 0 }
  if (w1 > h$duration) { notes <- c(notes, sprintf("window ends %.1f s after the recording; cut", w1 - h$duration)); w1 <- h$duration }
  if (w1 - w0 < p$WIN_MS / 1000) stop(sprintf("analysis window %.2f-%.2f s is shorter than one window", w0, w1))
  pad <- 5
  r0 <- max(0, w0 - pad); r1 <- min(h$duration, w1 + pad)
  x <- edf_read(h, ch$edf_labels, r0, r1)
  rownames(x) <- ch$names
  fs <- attr(x, "sfreq")
  x_start <- attr(x, "start_s")
  x <- x * to_uv(h$unit[match(ch$edf_labels[1], h$label)])
  pp <- preprocess(x, fs, p$BANDPASS_HZ, p$RESAMPLE_HZ, p$THREADS)
  rm(x)
  notes <- c(notes, pp$notes)
  fs2 <- pp$fs
  i0 <- max(1L, round((w0 - x_start) * fs2) + 1L)
  i1 <- min(ncol(pp$x), round((w1 - x_start) * fs2))
  seg <- pp$x[, i0:i1, drop = FALSE]
  seg_start <- x_start + (i0 - 1) / fs2
  rm(pp)
  # Workers are forked from this process: collect garbage first, or each worker's own
  # garbage collection touches (and so copies) every page of the parent's dead objects.
  invisible(gc(verbose = FALSE))
  win <- max(2L, round(p$WIN_MS / 1000 * fs2)); step <- max(1L, round(p$STEP_MS / 1000 * fs2))
  fm <- fragility_matrix(seg, win, step, p$LAMBDA, p$N_SEARCH, p$ENGINE, p$THREADS)
  ws <- seg_start + fm$start / fs2 - onset          # window start, s from onset
  we <- ws + win / fs2
  pre <- we <= 1e-9
  post <- ws >= -1e-9
  mean_rows <- function(cols) if (any(cols)) rowMeans(fm$frag[, cols, drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(fm$frag))
  f_pre <- mean_rows(pre); f_post <- mean_rows(post)
  if (!any(pre)) notes <- c(notes, "no window before the onset: fragility_pre empty")
  if (!any(post)) notes <- c(notes, "no window after the onset: fragility_post empty")
  f_pre[is.nan(f_pre)] <- NA; f_post[is.nan(f_post)] <- NA
  n_nan <- sum(colSums(is.na(fm$frag)) > 0)
  if (n_nan) notes <- c(notes, sprintf("%d window(s) with numerically undefined fragility for some channels (left out of the means)", n_nan))
  rk <- rank(-f_post, ties.method = "min", na.last = "keep")
  rows <- data.frame(recording_id = rec$recording_id, channel = rownames(fm$frag),
                     fragility_pre = round(f_pre, 6), fragility_post = round(f_post, 6),
                     rank = as.integer(rk), stringsAsFactors = FALSE)
  nun <- sum(!fm$stable, na.rm = TRUE)
  if (nun) notes <- c(notes, sprintf("%d of %d windows not stable at the largest ridge penalty", nun, length(fm$stable)))
  info <- list(recording_id = rec$recording_id, status = "ok",
               path = rec$path, channels = nrow(fm$frag), channels_source = ch$source,
               sfreq = fs, analysis_sfreq = fs2, onset_s = onset, onset_source = onset_source,
               window_s = c(round(w0, 3), round(w1, 3)),
               windows = length(ws), windows_pre = sum(pre), windows_post = sum(post),
               window_samples = win, step_samples = step,
               lambda = list(min = min(fm$lambda), median = stats::median(fm$lambda), max = max(fm$lambda),
                             raised = sum(fm$lambda > 1e-4 + 1e-12)),
               unstable_windows = nun, undefined_windows = n_nan,
               series_windows = sum(fm$terms > 0, na.rm = TRUE),
               seconds = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
               notes = paste(notes, collapse = "; "))
  list(rows = rows, info = info)
}

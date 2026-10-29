# Neural fragility (Li et al. 2021) per recording, with EZFragility.
#
# For each window, EZFragility fits a linear model x(t+1) = A x(t) by ridge regression,
# raising the ridge penalty until A is stable (ridgeSearch), then finds for each channel i
# the smallest perturbation of row i of A that puts an eigenvalue on the unit circle
# (fragilityRow); fragility = (max - norm) / max over the channels of that window, so the
# most fragile channel of each window is 1 and the least is 0.

suppressPackageStartupMessages({
  library(EZFragility)
  library(gsignal)
})

# fragilityRow, vectorised over the channels. EZFragility loops in R over every channel and
# every trial eigenvalue lambda = sqrt(1 - w^2) + i w, w in (0, 1], computing for row r of
# (A - lambda I)^-1 the minimum-norm Delta = B' (B B')^-1 b with B = [Im r; Re r], b = (0, -1).
# Its norm is sqrt(b' (B B')^-1 b) = sqrt(|Im r|^2 / (|Im r|^2 |Re r|^2 - (Im r . Re r)^2)),
# so all rows come from one inverse per trial eigenvalue. Same result (to rounding, see
# tests/testthat/test-fragility.R), about 3x faster, and the R loop overhead is gone.
fragility_row_fast <- function(A, nSearch = 100L) {
  n <- nrow(A)
  I <- diag(n)
  om <- seq(0, 1, length.out = nSearch + 1L)[-1L]
  best <- rep(Inf, n)
  for (w in om) {
    M <- solve(A - complex(real = sqrt(1 - w^2), imaginary = w) * I)
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
  best[!is.finite(best)] <- NA
  m <- max(best, na.rm = TRUE)
  if (!is.finite(m)) return(rep(NA_real_, n))
  (m - best) / m
}

ez_ridge_search <- utils::getFromNamespace("ridgeSearch", "EZFragility")
ez_standardize <- utils::getFromNamespace("standardizeIEEG", "EZFragility")

# x: channels x samples (rows named). Returns list(frag = channels x windows, start = window
# start in samples (0-based), lambda, stable).
fragility_matrix <- function(x, win, step, lambda = NULL, nSearch = 100L,
                             engine = "fast", threads = 1L) {
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
    list(f = fragility_row_fast(A, nSearch), lambda = attr(A, "lambda"),
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
       lambda = vapply(res, `[[`, 0, "lambda"), stable = vapply(res, `[[`, TRUE, "stable"))
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

# Band-pass (4th-order Butterworth, zero phase) then resample, one channel at a time (so
# memory holds the input and the output, not several filtered copies of the whole matrix).
preprocess <- function(x, fs, bandpass, resample_hz) {
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
    bf <- butter(4, c(lo, hi) / (fs / 2), "pass")
  }
  one <- function(v) {
    if (!is.null(bf)) v <- filtfilt(bf, v)
    if (!identical(pq, c(1L, 1L))) v <- resample(v, pq[1], pq[2])
    v
  }
  first <- one(x[1, ])
  out <- matrix(0, nrow(x), length(first), dimnames = list(rownames(x), NULL))
  out[1, ] <- first
  for (i in seq_len(nrow(x))[-1]) out[i, ] <- one(x[i, ])
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
  pp <- preprocess(x, fs, p$BANDPASS_HZ, p$RESAMPLE_HZ)
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
               seconds = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
               notes = paste(notes, collapse = "; "))
  list(rows = rows, info = info)
}

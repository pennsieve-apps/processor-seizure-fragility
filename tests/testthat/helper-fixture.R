# Test helpers: source the app, write tiny EDF+ files and a BIDS-iEEG derivative fixture
# that follows conventions.md (recordings.csv, _channels.tsv with onset labels, _events.tsv).

root <- Sys.getenv("APP_ROOT", normalizePath(file.path(testthat::test_path(), "..", "..")))
source(file.path(root, "R", "main.R"))

# Write signals (channels x samples, physical units) as EDF+C with 1-s records.
write_edf <- function(path, x, fs, unit = "uV", annotations = TRUE) {
  stopifnot(fs == round(fs))
  ns_data <- nrow(x)
  nrec <- ceiling(ncol(x) / fs)
  if (ncol(x) < nrec * fs) x <- cbind(x, matrix(0, nrow(x), nrec * fs - ncol(x)))
  labels <- c(rownames(x), if (annotations) "EDF Annotations")
  ns <- length(labels)
  spr <- c(rep(fs, ns_data), if (annotations) 30)
  pmax <- c(apply(abs(x), 1, function(v) ceiling(max(v, 1) * 1.1)), if (annotations) 1)
  pmin <- -pmax
  dmin <- rep(-32768, ns); dmax <- rep(32767, ns)
  pad <- function(s, n) formatC(substr(as.character(s), 1, n), width = -n)
  con <- file(path, "wb")
  on.exit(close(con))
  w <- function(s, n) writeChar(pad(s, n), con, eos = NULL)
  w("0", 8); w("X X X X", 80); w("Startdate 01-JAN-2020 X X X", 80)
  w("01.01.20", 8); w("00.00.00", 8); w(256 * (ns + 1), 8)
  w("EDF+C", 44); w(nrec, 8); w(1, 8); w(ns, 4)
  for (l in labels) w(l, 16)
  for (i in seq_len(ns)) w("", 80)
  for (i in seq_len(ns)) w(if (i <= ns_data) unit else "", 8)
  for (v in pmin) w(v, 8)
  for (v in pmax) w(v, 8)
  for (v in dmin) w(v, 8)
  for (v in dmax) w(v, 8)
  for (i in seq_len(ns)) w("", 80)
  for (v in spr) w(v, 8)
  for (i in seq_len(ns)) w("", 32)
  dig <- round((x - pmin[seq_len(ns_data)]) / (pmax[seq_len(ns_data)] - pmin[seq_len(ns_data)]) * 65535 - 32768)
  for (r in seq_len(nrec)) {
    cols <- (r - 1) * fs + seq_len(fs)
    writeBin(as.integer(t(dig[, cols, drop = FALSE])), con, size = 2, endian = "little")
    if (annotations) {
      tal <- c(charToRaw(sprintf("+%d", r - 1)), as.raw(c(20, 20, 0)))
      writeBin(c(tal, as.raw(rep(0, 60 - length(tal)))), con)
    }
  }
  invisible(path)
}

# A small ictal recording: n_ch bipolar channels of noise with a shared slow rhythm; from
# `onset` on, channels `soz` carry a growing 8 Hz oscillation that spreads to their neighbours.
synth_recording <- function(n_ch = 8, fs = 256, dur = 40, onset = 20, soz = 1:2, seed = 1) {
  set.seed(seed)
  t <- seq(0, dur - 1 / fs, by = 1 / fs)
  base <- sin(2 * pi * 1.5 * t)
  x <- t(vapply(seq_len(n_ch), function(i) 20 * stats::filter(stats::rnorm(length(t)), rep(1 / 3, 3), circular = TRUE) + 15 * base * (i %% 3), numeric(length(t))))
  g <- pmax(0, pmin(1, (t - onset) / 5))
  for (i in soz) x[i, ] <- x[i, ] + 200 * g * sin(2 * pi * 8 * t + i)
  for (i in setdiff(seq_len(n_ch), soz)) x[i, ] <- x[i, ] + 40 * g * sin(2 * pi * 8 * t - 0.3 * i) * (abs(i - max(soz)) <= 1)
  rownames(x) <- sprintf("A%02d-A%02d", seq_len(n_ch), seq_len(n_ch) + 1)
  x
}

# BIDS-iEEG derivative with `n_rec` recordings under <dir>/ieeg (or <dir> when bare).
make_fixture <- function(dir, n_rec = 2, fs = 256, dur = 40, onset = 20, n_ch = 8,
                         soz = 1:2, recordings_csv = TRUE, bad = NULL, subdir = "ieeg") {
  out <- if (nzchar(subdir)) file.path(dir, subdir) else dir
  rows <- list()
  for (k in seq_len(n_rec)) {
    rid <- sprintf("sub-01_ses-01_task-ictal%d_run-01", k)
    d <- file.path(out, "sub-01", "ses-01", "ieeg")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    x <- synth_recording(n_ch, fs, dur, onset, soz, seed = k)
    edf <- file.path(d, paste0(rid, "_desc-preproc_ieeg.edf"))
    write_edf(edf, x, fs)
    ch <- data.frame(name = rownames(x), type = "SEEG", units = "µV", status = "good",
                     status_description = "n/a", contacts = sub("-", "|", rownames(x)),
                     soz = "", resected = "", onset = ifelse(seq_len(n_ch) %in% soz, "true", "false"),
                     spread = "false", stringsAsFactors = FALSE)
    if (!is.null(bad)) ch$status[bad] <- "bad"
    utils::write.table(ch, file.path(d, paste0(rid, "_desc-preproc_channels.tsv")), sep = "\t",
                       quote = FALSE, row.names = FALSE, na = "")
    ev <- data.frame(onset = c(onset, 2), duration = c(10, 10), trial_type = c("seizure", "baseline"),
                     source = c("test", "test"))
    utils::write.table(ev, file.path(d, paste0(rid, "_desc-preproc_events.tsv")), sep = "\t",
                       quote = FALSE, row.names = FALSE)
    rows[[k]] <- data.frame(recording_id = rid, subject = "01", session = "01", task = paste0("ictal", k),
                            run = "01", path = file.path("sub-01", "ses-01", "ieeg", basename(edf)),
                            sfreq = fs, n_channels = n_ch, n_good = n_ch - length(bad), duration_s = dur,
                            montage = "bipolar", kind = "ictal", onset_s = onset,
                            baseline_start_s = 2, baseline_end_s = 12, source_format = "bids-edf")
  }
  if (recordings_csv) utils::write.csv(do.call(rbind, rows), file.path(out, "recordings.csv"), row.names = FALSE)
  out
}

# Parameters for fast tests (short windows, few trial eigenvalues).
fast_env <- function(...) {
  withr::local_envvar(c(WINDOW_S = "-4,4", RESAMPLE_HZ = "0", BANDPASS_HZ = "0.5,60",
                        N_SEARCH = "10", THREADS = "1", ...), .local_envir = parent.frame())
}

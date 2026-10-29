# A small EDF/EDF+ reader in base R: the header, and the samples of chosen signals over a
# time range (only the data records that cover it are read from disk).
#
# EDF (Kemp et al. 1992) and EDF+ (Kemp & Olivan 2003): a 256-byte header, 256 bytes per
# signal, then data records of little-endian 16-bit integers, each holding
# samples_per_record[i] samples of every signal in turn. Physical value =
# (digital - dmin) * (pmax - pmin) / (dmax - dmin) + pmin. The "EDF Annotations" signal of
# EDF+ is skipped. BDF (24-bit) is not supported.

edf_header <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  field <- function(n) trimws(rawToChar(readBin(con, "raw", n)))
  num <- function(n) suppressWarnings(as.numeric(field(n)))
  version <- field(8)
  if (!identical(version, "0")) stop(sprintf("%s: not an EDF file (version field '%s')", path, version))
  patient <- field(80); recording <- field(80)
  startdate <- field(8); starttime <- field(8)
  header_bytes <- num(8)
  reserved <- field(44)
  n_records <- num(8)
  record_duration <- num(8)
  ns <- num(4)
  if (is.na(ns) || ns < 1) stop(sprintf("%s: bad number of signals", path))
  sig <- function(n) vapply(seq_len(ns), function(i) field(n), "")
  label <- sig(16); transducer <- sig(80); unit <- sig(8)
  pmin <- as.numeric(sig(8)); pmax <- as.numeric(sig(8))
  dmin <- as.numeric(sig(8)); dmax <- as.numeric(sig(8))
  prefilter <- sig(80); spr <- as.numeric(sig(8))
  size <- file.size(path)
  rec_bytes <- 2 * sum(spr)
  # n_records may be -1 while a file is being written: derive it from the size.
  if (is.na(n_records) || n_records < 0) n_records <- floor((size - header_bytes) / rec_bytes)
  list(path = path, header_bytes = header_bytes, reserved = reserved, n_records = n_records,
       record_duration = record_duration, label = label, unit = unit, pmin = pmin, pmax = pmax,
       dmin = dmin, dmax = dmax, prefilter = prefilter, spr = spr,
       sfreq = spr / record_duration, duration = n_records * record_duration,
       rec_bytes = rec_bytes, annotation = label == "EDF Annotations",
       discontinuous = startsWith(reserved, "EDF+D"))
}

# Signals `labels` (all data signals when NULL) from start_s to end_s (seconds from the
# start of the file), as a channels x samples numeric matrix in physical units. All chosen
# signals must share one sampling rate. Attributes: sfreq, start_s (time of column 1).
edf_read <- function(h, labels = NULL, start_s = 0, end_s = Inf) {
  data_idx <- which(!h$annotation)
  idx <- if (is.null(labels)) data_idx else match(labels, h$label)
  if (anyNA(idx)) stop(sprintf("%s: signal(s) not found: %s", h$path,
                               paste(labels[is.na(idx)], collapse = ", ")))
  fs <- unique(h$sfreq[idx])
  if (length(fs) != 1) stop(sprintf("%s: the chosen signals have different sampling rates (%s)",
                                    h$path, paste(unique(h$sfreq[idx]), collapse = ", ")))
  rd <- h$record_duration
  start_s <- max(0, start_s); end_s <- min(h$duration, end_s)
  if (end_s <= start_s) stop(sprintf("%s: empty time range %.3f-%.3f s", h$path, start_s, end_s))
  r0 <- floor(start_s / rd); r1 <- min(h$n_records - 1, ceiling(end_s / rd) - 1)
  nrec <- r1 - r0 + 1
  con <- file(h$path, "rb")
  on.exit(close(con))
  seek(con, h$header_bytes + r0 * h$rec_bytes)
  words <- h$rec_bytes / 2
  raw <- readBin(con, "integer", n = nrec * words, size = 2, signed = TRUE, endian = "little")
  nrec <- length(raw) %/% words
  if (nrec < 1) stop(sprintf("%s: no data records in range", h$path))
  raw <- matrix(raw[seq_len(nrec * words)], nrow = words)
  offs <- c(0, cumsum(h$spr))[seq_along(h$spr)]
  spr <- h$spr[idx[1]]
  out <- matrix(0, nrow = length(idx), ncol = nrec * spr)
  for (k in seq_along(idx)) {
    i <- idx[k]
    d <- as.vector(raw[offs[i] + seq_len(spr), , drop = FALSE])
    gain <- (h$pmax[i] - h$pmin[i]) / (h$dmax[i] - h$dmin[i])
    out[k, ] <- (d - h$dmin[i]) * gain + h$pmin[i]
  }
  rm(raw)
  t0 <- r0 * rd
  first <- max(1, round((start_s - t0) * fs) + 1)
  last <- min(ncol(out), round((end_s - t0) * fs))
  out <- out[, first:last, drop = FALSE]
  rownames(out) <- h$label[idx]
  attr(out, "sfreq") <- fs
  attr(out, "start_s") <- t0 + (first - 1) / fs
  out
}

# Scale factor from an EDF physical dimension to microvolts (1 when unknown).
to_uv <- function(unit) {
  u <- tolower(gsub("µ|μ", "u", trimws(unit)))
  switch(u, v = 1e6, mv = 1e3, uv = 1, nv = 1e-3, 1)
}

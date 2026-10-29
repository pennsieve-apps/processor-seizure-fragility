# processor-seizure-fragility entry point: Rscript /app/R/main.R
# Reads $INPUT_DIR/{ieeg,onsets} (or everything side by side), writes $OUTPUT_DIR/fragility.csv
# and $OUTPUT_DIR/provenance.json.

app_dir <- function() {
  if (nzchar(Sys.getenv("APP_ROOT"))) return(file.path(Sys.getenv("APP_ROOT"), "R"))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "R"
}
source(file.path(app_dir(), "common.R"))
source(file.path(app_dir(), "edf.R"))
source(file.path(app_dir(), "fragility.R"))

COLUMNS <- c("recording_id", "channel", "fragility_pre", "fragility_post", "rank")

run <- function(input_dir = Sys.getenv("INPUT_DIR", "/data/input"),
                output_dir = Sys.getenv("OUTPUT_DIR", "/data/output")) {
  options(warn = 1)
  started <- Sys.time()
  p <- tryCatch(fragility_params(), error = function(e) { log_line("ERROR: %s", conditionMessage(e)); NULL })
  if (is.null(p)) return(1L)
  log_line("processor-seizure-fragility %s (EZFragility %s): WINDOW_S=%g,%g RESAMPLE_HZ=%g BANDPASS_HZ=%s WIN_MS=%g STEP_MS=%g LAMBDA=%s N_SEARCH=%d ENGINE=%s ONSET_SOURCE=%s THREADS=%d",
           app_version(), as.character(utils::packageVersion("EZFragility")), p$WINDOW_S[1], p$WINDOW_S[2],
           p$RESAMPLE_HZ, if (is.null(p$BANDPASS_HZ)) "none" else paste(p$BANDPASS_HZ, collapse = ","),
           p$WIN_MS, p$STEP_MS, if (is.null(p$LAMBDA)) "auto" else format(p$LAMBDA), p$N_SEARCH,
           p$ENGINE, p$ONSET_SOURCE, p$THREADS)

  ieeg_dir <- port_dir(input_dir, "ieeg")
  recs <- tryCatch(find_recordings(ieeg_dir), error = function(e) { log_line("ERROR: %s", conditionMessage(e)); NULL })
  if (is.null(recs) || !nrow(recs)) {
    log_line("ERROR: no recordings (recordings.csv or *.edf) under %s", ieeg_dir)
    return(1L)
  }
  onsets <- if (p$ONSET_SOURCE != "expert") read_onsets(port_dir(input_dir, "onsets")) else NULL
  log_line("%d recording(s) in %s; detected onsets: %s", nrow(recs), rel_path(ieeg_dir, input_dir) %||% ".",
           if (is.null(onsets)) "none" else sprintf("%d from %s", length(onsets),
             paste(rel_path(attr(onsets, "files"), input_dir), collapse = ", ")))
  if (p$ONSET_SOURCE == "onsets" && is.null(onsets)) {
    log_line("ERROR: ONSET_SOURCE=onsets but no seizure_onsets.csv in the input")
    return(1L)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv <- file.path(output_dir, "fragility.csv")
  part <- paste0(out_csv, ".part")
  utils::write.csv(as.data.frame(stats::setNames(replicate(length(COLUMNS), character(), simplify = FALSE), COLUMNS)),
                   part, row.names = FALSE)

  infos <- list()
  for (i in seq_len(nrow(recs))) {
    rec <- as.list(recs[i, , drop = FALSE])
    rid <- rec$recording_id
    kind <- tolower(as.character(rec$kind %||% NA))
    det <- if (!is.null(onsets) && rid %in% names(onsets)) onsets[[rid]] else NA_real_
    ex <- if (p$ONSET_SOURCE == "onsets") list(onset = NA_real_, source = NA) else expert_onset(rec)
    if (!is.na(det)) { onset <- det; src <- "onsets" }
    else { onset <- ex$onset; src <- ex$source }
    if (is.na(onset)) {
      why <- if (identical(kind, "interictal")) "interictal recording" else "no seizure onset known"
      log_line("%s: skipped (%s)", rid, why)
      infos[[length(infos) + 1]] <- list(recording_id = rid, status = "skipped", notes = why)
      next
    }
    r <- tryCatch(fragility_recording(rec, onset, src, p), error = function(e) e)
    if (inherits(r, "error")) {
      log_line("%s: FAILED: %s", rid, conditionMessage(r))
      infos[[length(infos) + 1]] <- list(recording_id = rid, status = "failed", path = rec$path,
                                         notes = conditionMessage(r))
      next
    }
    utils::write.table(r$rows, part, sep = ",", row.names = FALSE, col.names = FALSE,
                       append = TRUE, na = "", qmethod = "double")
    inf <- r$info
    infos[[length(infos) + 1]] <- inf
    top <- utils::head(r$rows$channel[order(r$rows$rank)], 3)
    log_line("%s: ok %d channels, onset %.2f s (%s), %d windows (%d pre, %d post) at %g Hz, lambda raised in %d, top %s, %.1fs%s",
             rid, inf$channels, inf$onset_s, inf$onset_source, inf$windows, inf$windows_pre, inf$windows_post,
             inf$analysis_sfreq, inf$lambda$raised, paste(top, collapse = " "), inf$seconds,
             if (nzchar(inf$notes)) paste0(" [", inf$notes, "]") else "")
    gc(verbose = FALSE)
  }
  file.rename(part, out_csv)
  status <- vapply(infos, function(x) x$status, "")
  params <- p
  params$BANDPASS_HZ <- if (is.null(p$BANDPASS_HZ)) "none" else p$BANDPASS_HZ
  params$LAMBDA <- if (is.null(p$LAMBDA)) "auto" else p$LAMBDA
  prov <- list(
    app = "processor-seizure-fragility", version = app_version(),
    tools = c(list(R = as.character(getRversion())), pkg_versions(c("EZFragility", "Epoch", "gsignal", "jsonlite")),
              list(BLAS = basename(extSoftVersion()[["BLAS"]] %||% ""))),
    parameters = params,
    inputs = list(recordings = rel_path(unique(stats::na.omit(recs$index)), input_dir),
                  edf = rel_path(recs$edf, input_dir),
                  onsets = if (is.null(onsets)) NULL else rel_path(attr(onsets, "files"), input_dir)),
    started = format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    recordings = infos
  )
  write_json_file(prov, file.path(output_dir, "provenance.json"))
  log_line("done: %d ok, %d skipped, %d failed; %.0fs", sum(status == "ok"), sum(status == "skipped"),
           sum(status == "failed"), prov$seconds)
  if (!any(status == "ok")) {
    log_line("ERROR: no recording was analysed")
    return(1L)
  }
  0L
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

if (sys.nframe() == 0L) quit(status = run(), save = "no")

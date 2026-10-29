# Shared helpers: logging, parameters, input discovery, BIDS sidecars, provenance.

suppressPackageStartupMessages(library(jsonlite))

log_line <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")
  flush(stdout())
}

app_version <- function() {
  for (f in c("/app/VERSION", file.path(dirname(app_dir()), "VERSION"), "VERSION"))
    if (file.exists(f)) return(trimws(readLines(f, n = 1, warn = FALSE)))
  "0.0.0"
}

# ---- parameters ---------------------------------------------------------------------

env_string <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (identical(trimws(v), "")) default else trimws(v)
}

env_number <- function(name, default, min = -Inf, max = Inf, integer = FALSE) {
  raw <- env_string(name, as.character(default))
  v <- suppressWarnings(as.numeric(raw))
  if (is.na(v) || v < min || v > max || (integer && v != round(v)))
    stop(sprintf("%s must be a %s in [%s, %s], got '%s'", name,
                 if (integer) "whole number" else "number", min, max, raw))
  if (integer) as.integer(v) else v
}

# "a,b" -> c(a, b); "none"/"off"/"" -> NULL.
env_pair <- function(name, default, allow_none = FALSE) {
  raw <- env_string(name, default)
  if (allow_none && tolower(raw) %in% c("none", "off", "false", "0")) return(NULL)
  v <- suppressWarnings(as.numeric(strsplit(gsub("[][ ]", "", raw), "[,;]")[[1]]))
  if (length(v) != 2 || anyNA(v) || v[1] >= v[2])
    stop(sprintf("%s must be two increasing numbers 'a,b'%s, got '%s'", name,
                 if (allow_none) " or none" else "", raw))
  v
}

env_choice <- function(name, default, choices) {
  v <- tolower(env_string(name, default))
  if (!v %in% choices) stop(sprintf("%s must be one of %s, got '%s'", name,
                                    paste(choices, collapse = ", "), v))
  v
}

fragility_params <- function() {
  lam <- env_string("LAMBDA", "auto")
  if (!identical(tolower(lam), "auto")) {
    l <- suppressWarnings(as.numeric(lam))
    if (is.na(l) || l <= 0) stop(sprintf("LAMBDA must be auto or a positive number, got '%s'", lam))
    lam <- l
  } else lam <- NULL
  list(
    WINDOW_S = env_pair("WINDOW_S", "-60,30"),
    RESAMPLE_HZ = env_number("RESAMPLE_HZ", 512, min = 0),
    BANDPASS_HZ = env_pair("BANDPASS_HZ", "0.5,150", allow_none = TRUE),
    WIN_MS = env_number("WIN_MS", 250, min = 10),
    STEP_MS = env_number("STEP_MS", 125, min = 1),
    LAMBDA = lam,
    N_SEARCH = env_number("N_SEARCH", 100, min = 1, max = 10000, integer = TRUE),
    ENGINE = env_choice("ENGINE", "fast", c("fast", "ezfragility")),
    ONSET_SOURCE = env_choice("ONSET_SOURCE", "auto", c("auto", "onsets", "expert")),
    THREADS = env_number("THREADS", 0, min = 0, max = 128, integer = TRUE)
  )
}

# ---- inputs -------------------------------------------------------------------------

# Files under dir matching pattern (recursive, symlinks followed), skipping hidden files.
list_files <- function(dir, pattern) {
  if (!length(dir) || !dir.exists(dir)) return(character())
  f <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  sort(f[!startsWith(basename(f), ".")])
}

# A port's folder: INPUT_DIR/<port>/ when it exists (linked input, placement subfolder),
# else INPUT_DIR itself (no links: everything side by side).
port_dir <- function(input_dir, port) {
  sub <- file.path(input_dir, port)
  if (dir.exists(sub)) sub else input_dir
}

read_table_any <- function(path) {
  sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
  utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                    na.strings = c("", "n/a", "NA"), encoding = "UTF-8", fileEncoding = "UTF-8-BOM")
}

as_num <- function(x) suppressWarnings(as.numeric(x))

rel_path <- function(paths, root) {
  if (!length(paths)) return(character())
  root <- paste0(normalizePath(root), "/")
  p <- normalizePath(paths, mustWork = FALSE)
  ifelse(startsWith(p, root), substring(p, nchar(root) + 1), p)
}

# BIDS stem of an EDF: name without _ieeg.edf, and without a trailing _desc-<x> entity.
bids_stem <- function(path) {
  s <- sub("\\.edf$", "", basename(path), ignore.case = TRUE)
  s <- sub("_ieeg$", "", s)
  sub("_desc-[A-Za-z0-9]+$", "", s)
}

# The sidecar <prefix>_<suffix> next to an EDF named <prefix>_ieeg.edf, if it exists.
sidecar <- function(edf, suffix) {
  prefix <- sub("_ieeg\\.edf$", "", edf, ignore.case = TRUE)
  f <- paste0(prefix, "_", suffix)
  if (file.exists(f)) f else NA_character_
}

# Recordings to analyse: rows of recordings.csv (ieeg-prep's index) when there is one, else
# every *.edf under the port (a plain BIDS-iEEG dataset or loose EDF files).
find_recordings <- function(ieeg_dir) {
  idx <- list_files(ieeg_dir, "^recordings\\.csv$")
  if (length(idx)) {
    rows <- lapply(idx, function(f) {
      r <- read_table_any(f)
      need <- c("recording_id", "path")
      if (!all(need %in% names(r)))
        stop(sprintf("%s lacks column(s) %s", f, paste(setdiff(need, names(r)), collapse = ", ")))
      r$edf <- normalizePath(file.path(dirname(f), r$path), mustWork = FALSE)
      r$index <- f
      r
    })
    cols <- unique(unlist(lapply(rows, names)))
    rows <- do.call(rbind, lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
    return(rows[!duplicated(rows$recording_id), , drop = FALSE])
  }
  edfs <- list_files(ieeg_dir, "\\.edf$")
  if (!length(edfs)) return(NULL)
  r <- data.frame(recording_id = bids_stem(edfs), path = rel_path(edfs, ieeg_dir),
                  edf = normalizePath(edfs), onset_s = NA_character_, kind = NA_character_,
                  index = NA_character_, stringsAsFactors = FALSE)
  r[!duplicated(r$recording_id), , drop = FALSE]
}

# Channels to analyse: the good iEEG channels of _channels.tsv (status != bad; types SEEG,
# ECOG, DBS, EEG when a type column exists), else every data signal of the EDF.
recording_channels <- function(edf, header) {
  tsv <- sidecar(edf, "channels.tsv")
  data_labels <- header$label[!header$annotation]
  if (is.na(tsv)) return(list(names = data_labels, edf_labels = data_labels, missing = character(), source = "edf"))
  ch <- read_table_any(tsv)
  keep <- rep(TRUE, nrow(ch))
  if ("status" %in% names(ch)) keep <- keep & !(tolower(trimws(ch$status)) %in% "bad")
  if ("type" %in% names(ch)) keep <- keep & toupper(trimws(ch$type)) %in% c("SEEG", "ECOG", "DBS", "EEG", "IEEG")
  names <- ch$name[keep]
  # Match to EDF labels exactly, else ignoring case and surrounding blanks.
  m <- match(names, data_labels)
  miss <- is.na(m)
  if (any(miss)) m[miss] <- match(toupper(trimws(names[miss])), toupper(trimws(data_labels)))
  list(names = names[!is.na(m)], edf_labels = data_labels[m[!is.na(m)]],
       missing = names[is.na(m)], source = basename(tsv))
}

# Seizure onsets from the onsets port (seizure-dynasd's seizure_onsets.csv): a named vector
# recording_id -> onset_s (the earliest when a recording has several rows).
read_onsets <- function(onsets_dir) {
  f <- list_files(onsets_dir, "^seizure_onsets\\.csv$")
  if (!length(f)) return(NULL)
  t <- do.call(rbind, lapply(f, function(x) {
    r <- read_table_any(x)
    if (!all(c("recording_id", "onset_s") %in% names(r))) return(NULL)
    r[, c("recording_id", "onset_s")]
  }))
  if (is.null(t) || !nrow(t)) return(NULL)
  t$onset_s <- as_num(t$onset_s)
  t <- t[!is.na(t$onset_s), , drop = FALSE]
  v <- tapply(t$onset_s, t$recording_id, min)
  structure(as.numeric(v), names = names(v), files = f)
}

# The expert/clinical onset: recordings.csv onset_s, else the first `seizure` event of the
# recording's _events.tsv (also accepting trial_type "ictal", as in some raw BIDS datasets).
expert_onset <- function(rec) {
  o <- as_num(rec$onset_s)
  if (length(o) && !is.na(o)) return(list(onset = o, source = "recordings.csv"))
  ev <- sidecar(rec$edf, "events.tsv")
  if (!is.na(ev)) {
    e <- read_table_any(ev)
    if (all(c("onset", "trial_type") %in% names(e))) {
      tt <- tolower(trimws(e$trial_type))
      hit <- which(tt %in% c("seizure", "ictal", "sz", "seizure_onset", "sz_onset"))
      if (length(hit)) return(list(onset = min(as_num(e$onset[hit]), na.rm = TRUE),
                                   source = basename(ev)))
    }
  }
  list(onset = NA_real_, source = NA_character_)
}

write_json_file <- function(x, path) {
  writeLines(toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA), path)
}

pkg_versions <- function(pkgs) {
  stats::setNames(lapply(pkgs, function(p) as.character(utils::packageVersion(p))), pkgs)
}

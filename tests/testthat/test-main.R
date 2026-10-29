read_out <- function(out) utils::read.csv(file.path(out, "fragility.csv"), stringsAsFactors = FALSE)

test_that("linked inputs: fragility.csv per conventions and provenance.json", {
  fast_env()
  d <- withr::local_tempdir()
  make_fixture(d, n_rec = 2, bad = 8)
  out <- file.path(d, "out")
  expect_equal(run(d, out), 0L)
  f <- read_out(out)
  expect_equal(names(f), c("recording_id", "channel", "fragility_pre", "fragility_post", "rank"))
  expect_equal(length(unique(f$recording_id)), 2)
  expect_equal(nrow(f), 2 * 7)                       # the bad channel is left out
  expect_false("A08-A09" %in% f$channel)
  expect_true(all(f$fragility_post >= 0 & f$fragility_post <= 1))
  r1 <- f[f$recording_id == f$recording_id[1], ]
  expect_equal(sort(r1$rank), rank(-r1$fragility_post, ties.method = "min")[order(rank(-r1$fragility_post, ties.method = "min"))])
  expect_equal(min(r1$rank), 1)
  prov <- jsonlite::fromJSON(file.path(out, "provenance.json"), simplifyVector = FALSE)
  expect_equal(prov$app, "processor-seizure-fragility")
  expect_equal(prov$parameters$N_SEARCH, 10)
  rec <- prov$recordings[[1]]
  expect_equal(rec$status, "ok")
  expect_equal(rec$onset_source, "recordings.csv")
  # the window that straddles the onset belongs to neither
  expect_true((rec$windows - rec$windows_pre - rec$windows_post) %in% 0:1)
  expect_equal(rec$channels, 7)
})

test_that("onsets from the onsets port win over recordings.csv; ONSET_SOURCE=expert ignores them", {
  fast_env()
  d <- withr::local_tempdir()
  make_fixture(d, n_rec = 1)
  dir.create(file.path(d, "onsets"))
  utils::write.csv(data.frame(recording_id = "sub-01_ses-01_task-ictal1_run-01", onset_s = 22.5, method = "dynasd"),
                   file.path(d, "onsets", "seizure_onsets.csv"), row.names = FALSE)
  out <- file.path(d, "out")
  expect_equal(run(d, out), 0L)
  prov <- jsonlite::fromJSON(file.path(out, "provenance.json"), simplifyVector = FALSE)
  expect_equal(prov$recordings[[1]]$onset_source, "onsets")
  expect_equal(prov$recordings[[1]]$onset_s, 22.5)
  withr::local_envvar(ONSET_SOURCE = "expert")
  out2 <- file.path(d, "out2")
  expect_equal(run(d, out2), 0L)
  prov2 <- jsonlite::fromJSON(file.path(out2, "provenance.json"), simplifyVector = FALSE)
  expect_equal(prov2$recordings[[1]]$onset_s, 20)
})

test_that("no links, plain BIDS without recordings.csv: onset from _events.tsv", {
  fast_env()
  d <- withr::local_tempdir()
  make_fixture(d, n_rec = 1, recordings_csv = FALSE, subdir = "")
  out <- file.path(d, "out")
  expect_equal(run(d, out), 0L)
  f <- read_out(out)
  expect_equal(unique(f$recording_id), "sub-01_ses-01_task-ictal1_run-01")
  prov <- jsonlite::fromJSON(file.path(out, "provenance.json"), simplifyVector = FALSE)
  expect_match(prov$recordings[[1]]$onset_source, "events.tsv")
})

test_that("a recording without an onset is skipped, a broken one fails, the step goes on", {
  fast_env()
  d <- withr::local_tempdir()
  ieeg <- make_fixture(d, n_rec = 3)
  rc <- utils::read.csv(file.path(ieeg, "recordings.csv"), stringsAsFactors = FALSE)
  rc$onset_s[2] <- NA
  utils::write.csv(rc, file.path(ieeg, "recordings.csv"), row.names = FALSE, na = "")
  ev <- list.files(ieeg, "ictal2_run-01_desc-preproc_events.tsv", recursive = TRUE, full.names = TRUE)
  unlink(ev)
  writeBin(as.raw(1:10), file.path(ieeg, rc$path[3]))
  out <- file.path(d, "out")
  expect_equal(run(d, out), 0L)
  prov <- jsonlite::fromJSON(file.path(out, "provenance.json"), simplifyVector = FALSE)
  st <- vapply(prov$recordings, function(r) r$status, "")
  expect_equal(st, c("ok", "skipped", "failed"))
  expect_equal(unique(read_out(out)$recording_id), rc$recording_id[1])
})

test_that("the step fails without input, and with every recording failing", {
  fast_env()
  d <- withr::local_tempdir()
  expect_equal(run(d, file.path(d, "out")), 1L)
  make_fixture(d, n_rec = 1, onset = 1000)   # onset after the end: no window fits
  expect_equal(run(d, file.path(d, "out2")), 1L)
  withr::local_envvar(WINDOW_S = "5,1")
  expect_equal(run(d, file.path(d, "out3")), 1L)   # bad parameter
})

test_that("windows that run past the recording are cut, and noted", {
  fast_env(WINDOW_S = "-30,4")
  d <- withr::local_tempdir()
  make_fixture(d, n_rec = 1)
  out <- file.path(d, "out")
  expect_equal(run(d, out), 0L)
  prov <- jsonlite::fromJSON(file.path(out, "provenance.json"), simplifyVector = FALSE)
  expect_match(prov$recordings[[1]]$notes, "before the recording")
  expect_equal(prov$recordings[[1]]$window_s[[1]], 0)
})

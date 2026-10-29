test_that("the EDF reader returns what was written, over any time range", {
  d <- withr::local_tempdir()
  x <- synth_recording(n_ch = 4, fs = 128, dur = 10)
  f <- write_edf(file.path(d, "a_ieeg.edf"), x, 128)
  h <- edf_header(f)
  expect_equal(h$label, c(rownames(x), "EDF Annotations"))
  expect_equal(h$duration, 10)
  y <- edf_read(h)
  expect_equal(dim(y), dim(x))
  step <- 2 * max(abs(x)) * 1.1 / 65535
  expect_lt(max(abs(y - x)), step)
  z <- edf_read(h, rownames(x)[3:4], 2.5, 4.25)
  expect_equal(ncol(z), round(1.75 * 128))
  expect_equal(attr(z, "start_s"), 2.5)
  expect_lt(max(abs(z - x[3:4, 2.5 * 128 + seq_len(ncol(z))])), step)
  expect_error(edf_read(h, "nope"), "not found")
})

test_that("the vectorised fragility row equals EZFragility's fragilityRow", {
  set.seed(3)
  ez_row <- utils::getFromNamespace("fragilityRow", "EZFragility")
  for (n in c(5, 12, 30)) {
    xt <- matrix(stats::rnorm(n * 80), n)
    A <- ez_ridge_search(xt[, -80], xt[, -1])
    expect_equal(fragility_row_fast(A, 40L), ez_row(A, 40L), tolerance = 1e-10)
  }
})

test_that("the fast engine gives EZFragility's calcAdjFrag result", {
  x <- synth_recording(n_ch = 6, fs = 128, dur = 6)
  a <- fragility_matrix(x, win = 32, step = 16, nSearch = 20L, engine = "fast")
  b <- fragility_matrix(x, win = 32, step = 16, nSearch = 20L, engine = "ezfragility")
  expect_equal(a$frag, b$frag, tolerance = 1e-10)
  expect_equal(a$lambda, b$lambda)
  expect_equal(a$start, b$start)
  expect_true(all(a$frag >= 0 & a$frag <= 1))
  # parallel workers give the same numbers
  p <- fragility_matrix(x, win = 32, step = 16, nSearch = 20L, engine = "fast", threads = 2L)
  expect_equal(p$frag, a$frag)
})

test_that("rational resampling ratios and preprocessing", {
  expect_equal(rational(512 / 1024), c(1L, 2L))
  expect_equal(rational(512 / 1000), c(64L, 125L))
  x <- synth_recording(n_ch = 3, fs = 256, dur = 4)
  p <- preprocess(x, 256, c(0.5, 150), 128)
  expect_equal(p$fs, 128)
  expect_equal(ncol(p$x), 512)
  expect_match(p$notes, "lowered")
  expect_equal(rownames(p$x), rownames(x))
  expect_equal(preprocess(x, 256, NULL, 512)$fs, 256)   # never upsampled
})

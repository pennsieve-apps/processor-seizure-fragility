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

test_that("the stable engine equals EZFragility's fragilityRow where that is well conditioned", {
  set.seed(3)
  ez_row <- utils::getFromNamespace("fragilityRow", "EZFragility")
  for (n in c(5, 12, 30)) {
    xt <- matrix(stats::rnorm(n * 80), n)
    A <- ez_ridge_search(xt[, -80], xt[, -1])
    expect_equal(as.vector(fragility_row_stable(A, 40L)), ez_row(A, 40L), tolerance = 1e-10)
    # and so does its Neumann series, on the same model scaled into its range
    As <- A * (0.5 / norm(A, "2"))
    s <- fragility_row_stable(As, 40L)
    expect_gt(attr(s, "terms"), 0L)
    expect_equal(as.vector(s), ez_row(As, 40L), tolerance = 1e-10)
  }
})

# A case written by tests/reference/make_cases.R with its 256-bit reference norms
# (tests/reference/reference.py): list(A, nSearch, norms = channels x nSearch).
read_case <- function(name) {
  f <- testthat::test_path("reference", paste0(name, "_A.bin"))
  con <- file(f, "rb"); on.exit(close(con))
  h <- readBin(con, "integer", 2L, size = 4, endian = "little")
  A <- matrix(readBin(con, "double", h[1]^2, size = 8, endian = "little"), h[1])
  lam <- readBin(con, "double", 2L * h[2], size = 8, endian = "little")
  ref <- readBin(testthat::test_path("reference", paste0(name, "_norms.bin")), "double",
                 h[1] * h[2], size = 8, endian = "little")
  list(A = A, nSearch = h[2], lam = complex(real = lam[seq_len(h[2])], imaginary = lam[-seq_len(h[2])]),
       norms = matrix(ref, h[1]))
}
rel_err <- function(v, ref) max(abs(v - ref) / abs(ref))
frag_of <- function(v) normalise_norms(suppressWarnings(apply(v, 1, min, na.rm = TRUE)))

test_that("the stable engine matches a 256-bit reference to 1e-12 (small and large models)", {
  for (case in c("wide", "wide06", "narrow")) {
    r <- read_case(case)
    expect_identical(trial_lambdas(r$nSearch), r$lam)   # the reference used these exact values
    expect_false(anyNA(r$norms))
    solve_only <- fragility_norms_stable(r$A, r$nSearch, series = FALSE)
    hybrid <- fragility_norms_stable(r$A, r$nSearch)
    expect_equal(attr(hybrid, "terms") > 0, case != "narrow", info = case)
    expect_lt(rel_err(solve_only, r$norms), 1e-12)
    expect_lt(rel_err(hybrid, r$norms), 1e-12)
    expect_lt(max(abs(frag_of(hybrid) - frag_of(r$norms))), 1e-12)
    expect_equal(as.vector(fragility_row_stable(r$A, r$nSearch)), frag_of(r$norms), tolerance = 1e-12)
  }
})

test_that("the fast engine's Gram determinant cancels on a small model (why v0.2.0 changed engine)", {
  r <- read_case("wide")
  I <- diag(nrow(r$A))
  fast <- vapply(r$lam, function(l) {
    M <- solve(r$A - l * I)
    ii <- rowSums(Im(M)^2); rr <- rowSums(Re(M)^2); ir <- rowSums(Im(M) * Re(M))
    sqrt(ii / (ii * rr - ir * ir))
  }, numeric(nrow(I)))
  expect_gt(rel_err(fast, r$norms), 1e-8)                    # 5e-7 here, up to 8e-3 on real windows
  expect_lt(rel_err(fragility_norms_stable(r$A, r$nSearch, series = FALSE), r$norms), 1e-12)
  # where the model is not small the two agree
  n <- read_case("narrow")
  expect_equal(fragility_row_fast(n$A, n$nSearch), frag_of(n$norms), tolerance = 1e-10)
})

test_that("a nearly singular A - lambda I gives no NaN", {
  set.seed(5)
  n <- 6
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  A <- Q %*% diag(c(1 - 1e-12, 0.5, 0.4, 0.3, 0.2, 0.1)) %*% t(Q)   # eigenvalue on the circle
  for (f in list(fragility_row_fast(A, 50L), fragility_row_stable(A, 50L))) {
    expect_true(all(is.finite(f)))
    expect_true(all(f >= 0 & f <= 1))
  }
})

test_that("the stable and fast engines give EZFragility's calcAdjFrag result", {
  x <- synth_recording(n_ch = 6, fs = 128, dur = 6)
  b <- fragility_matrix(x, win = 32, step = 16, nSearch = 20L, engine = "ezfragility")
  for (e in c("stable", "fast")) {
    a <- fragility_matrix(x, win = 32, step = 16, nSearch = 20L, engine = e)
    expect_equal(a$frag, b$frag, tolerance = 1e-10)
    expect_equal(a$lambda, b$lambda)
    expect_equal(a$start, b$start)
    expect_true(all(a$frag >= 0 & a$frag <= 1))
  }
  # parallel workers give the same numbers
  a <- fragility_matrix(x, win = 32, step = 16, nSearch = 20L)
  p <- fragility_matrix(x, win = 32, step = 16, nSearch = 20L, threads = 2L)
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
  expect_identical(preprocess(x, 256, c(0.5, 150), 128, threads = 2L), p)   # channels in workers
})

test_that("ENGINE=fast reproduces the fragility v0.1.0 computed (stored)", {
  # v0.1.0-fragility.rds: fragility_matrix() of v0.1.0 (amd64 image) on this synthetic recording,
  # with the penalty searched (raised in 18 of 95 windows) and fixed. Only the fast engine (v0.1.x)
  # is held to it; the stable engine (v0.2.0) is held to the 256-bit reference above.
  ref <- readRDS(testthat::test_path("v0.1.0-fragility.rds"))
  x <- synth_recording(n_ch = 12, fs = 128, dur = 12, onset = 6)
  a <- fragility_matrix(x, win = 32, step = 16, nSearch = 30L, engine = "fast")
  f <- fragility_matrix(x, win = 32, step = 16, lambda = 0.01, nSearch = 30L, engine = "fast")
  expect_equal(a$frag, ref$auto$frag, tolerance = 1e-10)
  expect_identical(a$lambda, ref$auto$lambda)
  expect_equal(f$frag, ref$fixed$frag, tolerance = 1e-10)
  expect_identical(a$stable, ref$auto$stable)
})

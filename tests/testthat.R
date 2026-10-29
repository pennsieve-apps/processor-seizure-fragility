# Run the unit tests: Rscript tests/testthat.R (inside the image: /app/tests/testthat.R).
library(testthat)
root <- if (dir.exists("/app/R")) "/app" else normalizePath(".")
Sys.setenv(APP_ROOT = root)
res <- test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)

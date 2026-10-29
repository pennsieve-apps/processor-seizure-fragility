# processor-seizure-fragility: neural fragility (Li et al. 2021) with the EZFragility R package.
#
# rocker/r-ver:4.6.0 pins CRAN to the Posit Package Manager snapshot of 2026-06-23 and installs
# Linux binaries from it, so every package version below is fixed by the base image digest.
# EZFragility 2.1.1 is in that snapshot (CRAN archived the package on 2026-10-09 because a
# vignette could no longer download its data; the package code is unchanged). The versions are
# checked so a silent change fails the build. The EDF reader is our own (R/edf.R, base R).
FROM rocker/r-ver:4.6.0@sha256:6f05a1a8b8c52328f181593923909b01cbfd14c9caea93bf75ddc65e806d8eac

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    TZ=UTC

# libuv: needed by the fs package's binary (pulled in by EZFragility's dependencies).
# OpenBLAS (Ubuntu 24.04's 0.3.26) instead of the reference BLAS/LAPACK: the complex solves
# of the fragility search and the eigenvalues of the stability check run ~1.3-4x faster.
RUN apt-get update \
 && apt-get install -y --no-install-recommends libuv1t64 libopenblas0-pthread \
 && rm -rf /var/lib/apt/lists/* \
 && Rscript -e ' \
  options(Ncpus = 4L); \
  install.packages(c("EZFragility", "Epoch", "gsignal", "doParallel", "foreach", "jsonlite", "testthat")); \
  want <- c(EZFragility = "2.1.1", Epoch = "1.0.7", gsignal = "0.3.7", doParallel = "1.0.17", \
            foreach = "1.5.2", jsonlite = "2.0.0"); \
  got <- vapply(names(want), function(p) as.character(packageVersion(p)), ""); \
  print(got); stopifnot(identical(unname(got), unname(want)))' > /tmp/install.log 2>&1 \
  || { grep -n -E "ERROR|error:|cannot|not available|Error" /tmp/install.log | tail -40; tail -40 /tmp/install.log; exit 1; } \
 && tail -4 /tmp/install.log \
 && rm -rf /usr/local/lib/R/site-library/*/doc /usr/local/lib/R/site-library/*/html \
           /usr/local/lib/R/site-library/*/help /tmp/Rtmp* /tmp/downloaded_packages /tmp/install.log

WORKDIR /app
COPY VERSION /app/VERSION
COPY R/ /app/R/
COPY tests/testthat/ /app/tests/testthat/
COPY tests/testthat.R /app/tests/testthat.R

# One BLAS thread per worker: parallelism is over time windows (THREADS worker processes).
ENV HOME=/tmp \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    INPUT_DIR=/data/input \
    OUTPUT_DIR=/data/output

USER 1000:1000
ENTRYPOINT ["Rscript", "--vanilla", "/app/R/main.R"]

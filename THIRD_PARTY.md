# Third-party software

**This image is distributed under the GNU General Public License, version 3 or later.** The
app's R code (`R/`) runs in the same R process as EZFragility, Epoch and gsignal, which are
GPL-3.0(-or-later), and it calls EZFragility's internal functions (`ridgeSearch`,
`standardizeIEEG`) directly. Because of that tight, in-process linking, this repository's own
code is also licensed **GPL-3.0-or-later** (see `LICENSE`), not Apache-2.0 like most Pennsieve
apps. The fixture script `tests/make_fixture.py` is part of the same repository and under the
same licence.

All R packages are installed unmodified as Linux binaries from the Posit Package Manager CRAN
snapshot of 2026-06-23 that `rocker/r-ver:4.6.0` is pinned to; the Dockerfile checks the
versions of the ones the app calls. Each package's licence is inside the image under
`/usr/local/lib/R/site-library/<package>/`.

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| EZFragility | 2.1.1 | Ridge regression of the linear network model with the stability search (`ridgeSearch`), scaling (`standardizeIEEG`), and `calcAdjFrag` with `ENGINE=ezfragility`; its `fragilityRow` is the reference the fast engine is tested against | GPL (>= 3) | https://CRAN.R-project.org/package=EZFragility (archived by CRAN on 2026-10-09, see README), https://github.com/Jiefei-Wang/EZFragility |
| Epoch | 1.0.7 | iEEG container class EZFragility depends on | GPL (>= 3) | CRAN |
| gsignal | 0.3-7 | Butterworth band-pass (`butter`, `filtfilt`) and polyphase resampling (`resample`) | GPL-3 | CRAN |
| doParallel, foreach, iterators | 1.0.17, 1.5.2, 1.0.14 | Parallel windows for `ENGINE=ezfragility` | GPL-2; Apache-2.0; Apache-2.0 | CRAN |
| jsonlite | 2.0.0 | provenance.json | MIT | CRAN |
| ggplot2, ggtext, gridtext, viridis, glue, rlang, progress, ramify, reshape2, pracma, … | snapshot 2026-06-23 | Imported by EZFragility / Epoch / gsignal (plots and helpers, not used by the app) | MIT; GPL-2 (ggtext); GPL (>= 2) (ramify); GPL (>= 3) (pracma) | CRAN |
| testthat (+ deps) | 3.3.2 | Unit tests (`/app/tests`) | MIT | CRAN |
| R | 4.6.0 | Interpreter | GPL-2 / GPL-3 | https://www.r-project.org |
| OpenBLAS | 0.3.26 (Ubuntu 24.04 `libopenblas0-pthread`) | BLAS/LAPACK for the matrix solves | BSD-3-Clause | https://www.openblas.net |
| libuv | Ubuntu 24.04 `libuv1t64` | Needed by the `fs` package binary | MIT | https://libuv.org |
| rocker/r-ver:4.6.0 | sha256:6f05a1a8… | Base image (Ubuntu 24.04 + R) | GPL-2+ (rocker scripts); Ubuntu packages under their own free licences | https://hub.docker.com/r/rocker/r-ver |

The EDF reader (`R/edf.R`) is the app's own (base R), written from the EDF/EDF+ specification.

The complete list of R packages in the image with versions and licences:
`docker run --rm --entrypoint Rscript processor-seizure-fragility -e 'installed.packages(fields="License")[, c("Version","License")]'`.

## Data

The demo data used to test, measure and validate the app (not in the image or the repository):
Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. *Multi-Expert Seizure Annotation*
(version 3). Pennsieve Discover (2026). doi:10.26275/dely-fzwx. CC BY 4.0.

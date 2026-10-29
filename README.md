# processor-seizure-fragility

**Neural fragility** (Li et al., *Nature Neuroscience* 2021) of every iEEG channel around a
seizure onset, computed with the **EZFragility** R package (Wang, Karas et al.). A step of the
Seizure Onset Localisation workflow:
`ieeg-prep → seizure-dynasd → {seizure-ei, seizure-fragility} → seizure-consensus → ieeg-report`.

**For the epileptologist.** Fragility asks, for each contact, how easily the network around
it could be tipped into a runaway (unstable) state by a small change in that contact's
connections. Li et al. found that, in the seconds around seizure onset, contacts in the
clinically identified seizure onset zone are the most fragile, and that this was more
pronounced in patients with good surgical outcome. The app gives each channel a fragility
score before and after the onset and ranks the channels; `seizure-consensus` fuses that
ranking with DynaSD's onset times and the Epileptogenicity Index. It is a localisation aid
for review alongside the EEG, not a decision.

**For the researcher.** In windows of 250 ms (step 125 ms) from 60 s before to 30 s after the
onset, the app fits the discrete linear model x(t+1) = A x(t) by ridge regression, raising the
ridge penalty from 1e-4 until A is stable (all eigenvalues inside the unit circle). For each
channel *i* it finds the minimum 2-norm perturbation Δ of row *i* that places an eigenvalue of
A + Δ on the unit circle, searching 100 points λ = √(1−ω²) + iω, ω ∈ (0, 1]. Fragility of a
channel in a window is (max − ‖Δᵢ‖) / max over the window's channels: 1 for the most fragile
channel, 0 for the least. The table reports the mean over the windows before and after the
onset. This is EZFragility's reading of the method; see *How faithful is this to Li et al.?*

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `ieeg` | input (subfolder), required | `recordings.csv`, `*_desc-preproc_ieeg.edf`, `_channels.tsv`, `_events.tsv` | ieeg-prep's BIDS-iEEG derivative. Also a plain BIDS-iEEG dataset (`*_ieeg.edf` + sidecars) or EDF files. Found in `INPUT_DIR/ieeg/`, else anywhere under `INPUT_DIR` (symlinks followed). |
| `onsets` | input (subfolder), optional | `seizure_onsets.csv` (`recording_id, onset_s, method`) | seizure-dynasd's detected onsets (it writes them under `dynasd/`); used in place of the expert onset for the recordings it lists. Found anywhere under `INPUT_DIR/onsets/`, else anywhere under `INPUT_DIR`. |
| `fragility` | output | `fragility/fragility.csv`, `fragility/provenance.json` | One output port: the whole `OUTPUT_DIR`, everything under `fragility/` so it can be one of several links into a fan-in port (seizure-consensus `methods`) without name collisions. |

Channels used: those of `_channels.tsv` with `status` not `bad` and `type` SEEG, ECOG, DBS or
EEG (all EDF signals except EDF+ annotations when there is no `_channels.tsv`).

**Onset.** `ONSET_SOURCE=auto` (default) takes, per recording, the onset from `seizure_onsets.csv`
when it lists the recording, else `recordings.csv` `onset_s` (the expert/clinical onset), else the
first `seizure` (or `ictal`) row of `_events.tsv`. Recordings without any onset (interictal ones)
are skipped and listed in `provenance.json`. The source used is recorded per recording.

## Parameters

| Parameter | Default | What it does |
|---|---|---|
| `WINDOW_S` | `-60,30` | Analysed segment, seconds relative to the onset. Cut to the recording when it runs past either end (noted in `provenance.json`). |
| `ONSET_SOURCE` | `auto` | `auto`, `onsets` (detected onsets only; recordings without one are skipped), `expert` (ignore the onsets port). |
| `RESAMPLE_HZ` | `512` | Faster signals are resampled (polyphase, gsignal `resample`) to this rate before the model is fitted; `0` keeps the recorded rate; never upsampled. |
| `BANDPASS_HZ` | `0.5,150` | Zero-phase 4th-order Butterworth band-pass before resampling, as in EZFragility's reproduction of the paper; `none` skips it. The upper edge is lowered to 0.9 × Nyquist when needed. |
| `WIN_MS` | `250` | Model window. |
| `STEP_MS` | `125` | Step between windows. |
| `LAMBDA` | `auto` | Ridge penalty. `auto`: EZFragility's search (1e-4, raised by bisection up to 10 until A is stable). A number fixes it for every window (EZFragility's manuscript compares 1e-4, 1e-3, 1e-2). |
| `N_SEARCH` | `100` | Trial eigenvalues on the unit circle per channel (EZFragility's `nSearch`). |
| `ENGINE` | `fast` | `fast`: EZFragility's ridge search, with the fragility row computed for all channels at once (below). `ezfragility`: EZFragility's `calcAdjFrag` unchanged (~3x slower). |
| `THREADS` | `0` | Worker processes over windows (fork). `0` = one per vCPU the container is given (cgroup CPU quota, else CPU affinity, else the CPU count), capped so that ~550 MB per worker plus ~500 MB for the parent fit the memory limit; the number chosen and why are logged and kept in `provenance.json` (`threads`). One recording is processed at a time. |

## Output: `fragility.csv`

One row per recording and channel (conventions.md):

| Column | Meaning |
|---|---|
| `recording_id` | BIDS stem of the recording (`sub-HUP224_ses-postimplant_task-ictal339135_run-00`) |
| `channel` | `_channels.tsv` name, e.g. `LB03-LB04` |
| `fragility_pre` | Mean fragility (0–1) over the windows that end at or before the onset (−60 s to 0 by default) |
| `fragility_post` | Mean fragility over the windows that start at or after the onset (0 to +30 s) |
| `rank` | Rank of `fragility_post` within the recording, 1 = most fragile (ties share the lower rank) |

The one window that straddles the onset belongs to neither mean. `provenance.json` holds per
recording: the onset used and its source, the analysed interval, the analysis sampling rate,
number of windows before/after, window and step in samples, the ridge penalties (min, median,
max and in how many windows the penalty had to be raised), windows still unstable at the
largest penalty, seconds, notes; plus app version, R and package versions, parameters and
inputs.

`seizure-consensus` ranks channels by `fragility_post` (or `fragility_post − fragility_pre`).

## The fast engine

EZFragility's `fragilityRow` loops in R over every channel and every trial eigenvalue λ,
taking row *i* of (A − λI)⁻¹, r, building B = [Im r; Re r] and the minimum-norm solution
Δ = Bᵀ(BBᵀ)⁻¹b with b = (0, −1). Its norm has a closed form, ‖Δ‖² = bᵀ(BBᵀ)⁻¹b =
‖Im r‖² / (‖Im r‖²‖Re r‖² − (Im r · Re r)²), so all channels follow from one inverse per λ.
`ENGINE=fast` uses EZFragility's own `ridgeSearch` and `standardizeIEEG` and this vectorised
row; the unit tests check it against EZFragility's `fragilityRow` and against `calcAdjFrag`
end to end (equal to 1e-10). On a 208-channel window it takes 0.8 s instead of 2.5 s (amd64
emulation), and it does not keep the per-window matrices, so memory does not grow with the
number of windows.

## How faithful is this to Li et al.?

EZFragility's README says plainly that the paper leaves out details needed to reproduce it and
documents the choices it makes instead:

- **Ridge penalty.** The paper recommends λ = 1e-4; on the authors' pt01 data that does not
  always give a stable A, so EZFragility searches upward until it does. On the demo clips the
  penalty had to be raised in 4.5 % (777 of 17,256; at most 0.12), and every window ended up stable of the windows.
- **The perturbation's eigenvalue.** The paper places the eigenvalue at some |λ| = 1 without
  saying which; EZFragility takes the λ that minimises the perturbation norm over a grid
  (`N_SEARCH`).
- **Normalisation.** Fragility is normalised per window ((max − norm)/max), so it ranks
  channels within a window; values are not comparable across recordings in absolute terms.
- **Preprocessing.** Li et al. analysed data around 1 kHz; EZFragility's reproduction band-passes
  0.5–150 Hz and works on clips of ±10–30 s. The app uses the same band and resamples to 512 Hz
  by default. On four demo clips, 512 Hz and 1024 Hz gave the same median AUC (0.970 vs
  0.976); 512 Hz halves the filtering work. Runtime is set by the number of channels and
  windows, not the sampling rate.
- **Summary.** The paper's figures show fragility heat maps and compare SOZ with non-SOZ
  distributions over time; collapsing to a mean before and after onset is this app's choice,
  so that the result fits a channel table.
- The CRAN package was **archived on 2026-10-09** because its reproduction vignette could not
  download its data during CRAN checks; the code is unchanged. The image installs 2.1.1 from
  the Posit Package Manager snapshot of 2026-06-23.

## Results on the demo data

All 24 SEEG clips of the Multi-Expert Seizure Annotation dataset (6 HUP patients, 4 seizures
each, 1 of them stimulation-induced), after the real ieeg-prep (bipolar, 72–204 good channels,
1024 Hz), default parameters, expert onset (`recordings.csv` `onset_s`). Validation: per clip,
the AUC of `fragility_post` for telling the expert consensus onset channels (`onset` in
`_channels.tsv`, 1–13 per clip) from all other channels.

| Score | Median AUC | Range | Mean |
|---|---|---|---|
| `fragility_post` (0 to +30 s) | **0.986** | 0.551–1.000 | 0.921 |
| `fragility_pre` (−60 to 0 s) | 0.935 | 0.549–0.994 | 0.885 |
| `fragility_post − fragility_pre` | 0.643 | 0.014–1.000 | 0.610 |

- Per patient: HUP235, HUP238 and HUP250 0.93–1.00 in every seizure; HUP224 0.95–0.99;
  HUP257 0.90–0.99; **HUP230 0.55–0.72** in all four seizures (the most fragile channels, LO/LM
  and LB, are not the experts' onset channels).
- Spontaneous seizures median 0.992 (n = 18), stimulation-induced 0.936 (n = 6); with low-voltage
  fast activity 0.946 (n = 9), without 0.995 (n = 15). (These two flags are not in ieeg-prep's
  `recordings.csv` yet; they were taken from the dataset's `annotations.tsv` for this check only.)
- **With seizure-dynasd's detected onsets** (the workflow default when `onsets` is linked;
  median absolute onset error 5.0 s, up to 59.5 s): `fragility_post` median AUC 0.981 (range
  0.384–1.000), `fragility_pre` 0.924, post − pre 0.795. The worst clip falls to 0.38 when the
  detected onset is far from the expert one. These outputs are in
  `workflow-epilepsy-ieeg/data/work/fragility/`.
- Through seizure-consensus with fragility alone: precision 0.35, recall 0.87, F1 0.46 (medians)
  for the top 10 channels.
- **Read with care.** The most fragile channels are largely the same in every seizure of a
  patient, and already before onset (pre-onset AUC 0.935), so fragility here behaves much like a
  stable property of the channel (its signal and its place in the network) that coincides with
  the expert onset channels, rather than a change at onset: the post − pre difference is barely
  better than chance. 24 seizures from 6 patients are 6 independent tests, not 24.
- 512 Hz vs 1024 Hz analysis rate on four clips (stand-in preprocessing): median AUC 0.970 vs 0.976,
  same top channels.

## Runtime

`computeTypes: [standard]`, 4 vCPU / 8 GiB (4 GiB is not a valid Fargate size for 4 vCPU), 4 h
timeout. Cost per window grows with the cube
of the number of channels (100 complex matrix inverses per window); the number of windows is
fixed by `WINDOW_S`/`STEP_MS` (719 by default), not by the sampling rate. One recording at a time,
windows in `THREADS` forked workers.

Measured on the 24 demo clips after ieeg-prep, amd64 image under emulation on an M-series Mac
(Docker Desktop), `--cpus 4 --memory 4g`, `THREADS=4`; anonymous memory of the container's cgroup
sampled every 0.5 s:

| | channels | wall time per clip | peak memory |
|---|---|---|---|
| HUP238 (smallest) | 72–113 | 22–55 s | |
| HUP250 | 148–152 | 92–107 s | |
| HUP224 (largest) | 186–204 | 177–200 s | 2.29 GB anon (memory.peak 2.34 GB) |
| all 24 clips | 72–204 | 39 min 25 s | 2.29 GB |

Memory is ~0.5 GB per worker for a 200-channel recording (each forked worker's own heap), plus
the parent; ~1.3 GB with 4 workers for a 15 s window. One clip (HUP257, 105 channels) took 682 s
in that run, against 43–49 s for its siblings, while other containers shared the Docker VM; we
read that as contention, not the app. Before the vectorised fragility row, `ENGINE=ezfragility`
measured 2.5 s per 208-channel window instead of 0.8 s. The image is 1.44 GB unpacked, 394 MB
compressed (`rocker/r-ver:4.6.0` 1.36 GB plus ~80 MB of R packages).

Re-measured for v0.1.1 on 3 clips (HUP224 186 ch, HUP230 95 ch, HUP250 151 ch; same emulation,
`--cpus 4 --memory 8g`, `THREADS=4`): 279 s wall, 1,030 CPU-s, cgroup `memory.peak` 2.33 GB
(148, 37 and 93 s per clip); natively (arm64, same Mac, 4 workers) 115 s. Where the time goes
(Rprof, natively, 186 channels): 90% in the fragility row, of which two thirds are the 100
complex inverses `solve(A - lambda I)` (already at BLAS speed) and a fifth the row sums; the
ridge search is 6%, band-pass and resampling 3%. Bad channels never enter the model (only
`good` channels are read). Parallelism is over windows in forked processes with one OpenBLAS
thread each, which keeps every vCPU busy; recordings in parallel would only add the ~3% serial
part. v0.1.1 is therefore as fast as v0.1.0: it changes no computation.

**Why it was not made faster.** For the recordings with more channels than samples per window
(about 128 at 512 Hz and 250 ms), the ridge model `A` is tiny (all eigenvalues below ~0.1) and
`A - lambda I` is close to `-lambda I`. The fragility row then rests on a nearly rank-1 Gram
determinant `|Im r|^2 |Re r|^2 - (Im r . Re r)^2`, which cancels by a factor up to ~10^7, so the
per-window fragility carries rounding noise of up to ~10^-2 (8e-3 against the same quantity
computed without the cancellation, as `|Re r|^2 |Im r_perp|^2`). Any reformulation (one
eigen-decomposition per window with closed forms per trial eigenvalue, a Neumann series, BLAS
threads, a different BLAS kernel) changes that noise: the same v0.1.0 image differs by up to
6e-3 per window between two OpenBLAS kernels (forced Haswell vs the one picked under emulation), and its `fragility.csv`
by up to 7e-5 between amd64 and arm64 (ranks identical). A ~2x faster series form for these
windows moved the CSV means by up to 1.4e-3 (ranks identical, but less accurate than v0.1.0
against the cancellation-free reference), so it was not adopted. Making the computation
accurate first, then faster, changes results and is left for a minor version.

## Running it locally

```sh
make test                 # unit tests in the image (synthetic EDF fixtures written in R)
make fixture              # synthetic derivative + seizure_onsets.csv -> data/fixture/ (needs NumPy)
make run                  # build, run over data/fixture/ -> data/output/
```

On real data see `sample-input/README.md`.

## Limitations

- Fragility ranks channels **within** a recording; it does not say whether a recording has a
  focal onset at all, and the most fragile channel of a diffuse or poorly sampled seizure is
  still ranked first.
- It needs a seizure onset time; the window is placed around it, so an onset off by several
  seconds moves the before/after split. With seizure-dynasd linked, the detected onset is used.
- High fragility before onset is common in the onset zone (also on the demo clips), so
  `fragility_pre` is not a "baseline" to subtract blindly.
- Linear model in short windows: artefacts, flat or bridged channels and strong common
  signals distort A. Re-reference first (ieeg-prep's bipolar montage) and mark bad channels.
- Stimulation-induced seizures carry stimulation artefacts in the pre-onset window. On such
  windows the fitted model can sit right at the stability boundary, where the perturbation norm
  is numerically undefined for some trial eigenvalues; those are skipped (EZFragility's loop
  would give an arbitrary value or fail) and windows with any undefined channel are counted in
  `provenance.json` (`undefined_windows`).
- The time-resolved fragility matrix (channels × windows) is not written, only its means.
- Runtime grows with channels³ per window: a 200-channel montage takes ~4–5 min per seizure
  on 4 vCPUs; very large montages (> 300 channels) are slow.

## Licences

**GPL-3.0-or-later.** The app's R code runs in one R process with EZFragility, Epoch and
gsignal (GPL-3.0(-or-later)) and calls EZFragility's internal functions, so the repository
(`LICENSE`) and the image are GPL-3.0-or-later, unlike most Pennsieve apps (Apache-2.0).
See `THIRD_PARTY.md`.

## Citations

- Li A, Huynh C, Fitzgerald Z, Cajigas I, Brusko D, Jagid J, Claudio AO, Kanner AM, Hopp J,
  Chen S, Haagensen J, Johnson E, Anderson W, Crone N, Inati S, Zaghloul KA, Bulacio J,
  Gonzalez-Martinez J, Sarma SV. Neural fragility as an EEG marker of the seizure onset zone.
  *Nat Neurosci* 24(10):1465–1474 (2021). doi:10.1038/s41593-021-00901-w. PMID 34354282.
- Wang J, Karas P, Lesage A-C, Malagaris I, Zhou O, Camarillo Rodriguez L, O'Leary S, Zhang Y.
  *EZFragility: Compute Neural Fragility for Ictal iEEG Time Series*. R package version 2.1.1,
  CRAN (2025). doi:10.32614/CRAN.package.EZFragility. https://CRAN.R-project.org/package=EZFragility
  (GPL-3.0-or-later).
- Wang J, Lesage A-C, Zhou O, Malagaris I, O'Leary S, Camarillo Rodriguez L, Silveira DC,
  Wang Z, Zhang Y, Karas PJ. Building an Ecosystem of Seizure Localization Methods: Neural
  Fragility as the First Step. *eNeuro* 13(3):ENEURO.0340-25.2026 (2026).
  doi:10.1523/ENEURO.0340-25.2026. (EZFragility's companion paper.)
- Demo data: Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. *Multi-Expert Seizure
  Annotation* (version 3). Pennsieve Discover (2026). doi:10.26275/dely-fzwx. CC BY 4.0.

DOIs, titles, author lists and years were checked against Crossref (articles, the CRAN DOI)
and DataCite (the dataset); the PMID against PubMed E-utilities (the eNeuro paper had no PMID
at the time of writing).

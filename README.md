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
| `BANDPASS_HZ` | `0.5,150` | Zero-phase 4th-order Butterworth band-pass before resampling, as in EZFragility's reproduction of the paper, applied in second-order sections (from v0.2.0; see *Method*); `none` skips it. The upper edge is lowered to 0.9 × Nyquist when needed. |
| `WIN_MS` | `250` | Model window. |
| `STEP_MS` | `125` | Step between windows. |
| `LAMBDA` | `auto` | Ridge penalty. `auto`: EZFragility's search (1e-4, raised by bisection up to 10 until A is stable). A number fixes it for every window (EZFragility's manuscript compares 1e-4, 1e-3, 1e-2). |
| `N_SEARCH` | `100` | Trial eigenvalues on the unit circle per channel (EZFragility's `nSearch`). |
| `ENGINE` | `stable` | `stable`: EZFragility's ridge search and minimum-norm perturbation, computed without rounding cancellation (*Method*). `fast`: the v0.1.x computation, kept for comparison. `ezfragility`: EZFragility's `calcAdjFrag` unchanged (~3x slower; fails on some windows of recordings with many channels). |
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

## Method: the minimum-norm perturbation, computed accurately (v0.2.0)

For each window the model A is EZFragility's (`ridgeSearch` on `standardizeIEEG`d data). For
channel *i* and trial eigenvalue λ, let r be row *i* of (A − λI)⁻¹. EZFragility's `fragilityRow`
solves for the minimum-norm real Δ with [Im r; Re r] Δ = (0, −1); its norm has the closed form

‖Δ‖² = ‖Im r‖² / det G,  G = Gram matrix of (Re r, Im r),  det G = ‖Im r‖²‖Re r‖² − (Im r · Re r)².

That is the quantity of Li et al. 2021 and EZFragility, and every engine computes it. They differ
in how:

- **`ENGINE=stable` (default from v0.2.0).** When A is small (all eigenvalues ≲ 0.1, which is
  what the ridge fit gives on recordings with more channels than samples per window, ≥ ~128
  channels at 512 Hz and 250 ms), r ≈ e*ᵢ*/λ plus a small part, Re r and Im r are nearly
  parallel, and det G as written cancels by up to ~10⁷. The stable engine instead
  1. computes N = (A − λI)⁻¹A with one complex solve (A as right-hand side), so the small
     off-diagonal entries of r = (N*ᵢ* − e*ᵢ*)/λ carry full relative accuracy rather than being
     the difference of the inverse and −I/λ;
  2. rotates r by the unit complex number that makes its diagonal entry real,
     r′ = r · conj(r*ᵢᵢ*)/|r*ᵢᵢ*|. A rotation of (Re r, Im r) leaves det G unchanged, and with
     Im r′*ᵢᵢ* = 0 it splits exactly into two non-negative terms over the off-diagonal part w of r′:
     det G = |r*ᵢᵢ*|² ‖Im w‖² + (‖Re w‖²‖Im w‖² − (Re w · Im w)²).
     The first term dominates and is a product of sums of squares; the second (which still
     cancels internally) is smaller by ‖w‖²/|r*ᵢᵢ*|², so its rounding no longer matters;
  3. ‖Im r‖² is a sum of squares.
  When ‖A‖₂ ≤ 0.7, N = −Σₖ λ⁻ᵏAᵏ (Neumann series) truncated where the tail bound
  q^(K+1)/(1 − q) falls below 2⁻⁶⁰‖A‖ (K ≤ 120 terms): K real matrix products per window
  replace the 100 complex solves, with the same accuracy (checked against the same reference).
  `provenance.json` counts these windows (`series_windows`).
- **`ENGINE=fast`** (v0.1.x): the closed form above evaluated directly from (A − λI)⁻¹ for all
  channels at once. Kept to reproduce v0.1.x results.
- **`ENGINE=ezfragility`**: EZFragility's `calcAdjFrag` unchanged. It solves the 2 × 2 system
  B Bᵀ per channel and λ and has the same cancellation; on 6 of 20 sampled demo windows with
  ≥ 151 channels its `solve()` stops with "system is computationally singular", and on another
  its normalised fragility is off by 0.49.

**Band-pass in second-order sections (v0.2.0).** The 0.5–150 Hz Butterworth band-pass at
1 kHz, applied as one 8th-order (b, a) polynomial as in v0.1.x, loses ~6 digits: the filtered
signal differed by 2e-6 (relative) between amd64 and arm64, enough to move `fragility_pre` by up
to 9e-5. The same filter in second-order sections (gsignal `butter(..., output = "Sos")`) agrees
to 2e-12 between the two, and differs from the (b, a) result by ~7e-6. It is ~4x slower, so
channels are now filtered in the `THREADS` workers.

## Validation

**Against a high-precision reference.** For 20 windows of 5 demo clips (72, 95, 151, 186 and
204 channels; 10 %, 50 %, 70 % and 90 % into each clip; model norms ‖A‖₂ from 0.02 to 3.8) the
norm of every channel at all 100 trial eigenvalues was computed in 256-bit ball arithmetic
(python-flint, `tests/reference/reference.py`; the balls bound the reference's own error below
1e-60) from the exact double-precision A, and compared with each engine:

| Engine | max relative error of ‖Δ‖, ≥ 151 channels | same, ≤ 95 channels | max error of the window's fragility (0–1) |
|---|---|---|---|
| `stable` (solve) | 4.8e-14 | 7.2e-14 | 1.5e-14 |
| `stable` (Neumann series, the 6 windows with ‖A‖₂ ≤ 0.7) | 1.6e-15 | – | 3.5e-16 |
| `fast` (v0.1.x) | **7.9e-3** | 8.6e-10 | **1.9e-3** |
| `ezfragility` (`fragilityRow`) | – | – | fails in 6 of 12 windows; 0.49 in another |

`tests/testthat/test-fragility.R` holds the engines to stored references of the same kind
(`tests/testthat/reference/`: a 140-channel model with ‖A‖₂ = 0.004, the same scaled to 0.6, and
a 16-channel model with an eigenvalue at 0.994): the stable engine, by solve and by series,
within 1e-12; the fast engine misses by 5e-7 on the small model. `tests/reference/make_cases.R`
and `reference.py` regenerate them.

**Across BLAS kernels and architectures.** Every 2nd window of two clips (186 and 151
channels) on arm64 (native) and amd64 (emulated) with OpenBLAS's automatic kernel,
`OPENBLAS_CORETYPE=Haswell` and `Prescott`: the stable engine's per-window fragility agrees to
6e-12 in every pair; the fast engine's differs by up to 1.1e-2. With the band-pass in sections
as well, v0.2.0's `fragility.csv` for all 24 demo clips is identical on amd64 and arm64 (to its
6 decimals; v0.1.1: up to 8.8e-5, 20 channels in 6 clips ranked differently).

**v0.1.1 vs v0.2.0 on the 24 demo clips** (amd64 images, defaults, expert onset; AUC of the
expert onset channels as in *Results*):

| | all 24 clips | 8 clips with ≥ 128 channels | 16 clips with < 128 channels |
|---|---|---|---|
| max \|Δ `fragility_pre`\| | 9.8e-4 | 9.8e-4 | 8.4e-4 |
| max \|Δ `fragility_post`\| | 4.4e-4 | 2.7e-4 | 4.4e-4 |
| Spearman of `fragility_post`, min (median) | 0.99992 (0.999998) | 0.99996 (0.999995) | 0.99992 (1) |
| top-10 channels the same | 24 of 24 | 8 of 8 | 16 of 16 |
| channels whose rank changed (clips) | 69 (14), by ≤ 3 places | 43 (8) | 26 (6) |
| AUC `fragility_post`, median (range) | 0.986 (0.551–1.000) → 0.986 (0.551–1.000) | 0.986 (0.946–1.000) → 0.986 (0.945–1.000) | 0.981 (0.551–1.000) → 0.981 (0.551–1.000) |
| AUC `fragility_pre`, median | 0.935 → 0.934 | 0.931 → 0.930 | 0.936 → 0.936 |
| largest change of a clip's AUC (post / pre / post − pre) | 0.004 / 0.001 / 0.007 | 0.001 / 0.001 / 0.007 | 0.004 / 0.000 / 0.002 |

The engine alone (v0.2.0 with the old (b, a) band-pass) changed `fragility_pre` by up to 9.5e-4
and `fragility_post` by up to 1.7e-4, mostly in the clips with ≥ 128 channels or stimulation
artefacts, and two ranks in one clip; the band-pass in sections accounts for the rest. Both
are below what matters for the ranking: no clip's top 10 changed and no AUC moved by more than
0.007.

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
1024 Hz), default parameters (v0.2.0), expert onset (`recordings.csv` `onset_s`). Validation: per clip,
the AUC of `fragility_post` for telling the expert consensus onset channels (`onset` in
`_channels.tsv`, 1–13 per clip) from all other channels.

| Score | Median AUC | Range | Mean |
|---|---|---|---|
| `fragility_post` (0 to +30 s) | **0.986** | 0.551–1.000 | 0.922 |
| `fragility_pre` (−60 to 0 s) | 0.934 | 0.549–0.994 | 0.885 |
| `fragility_post − fragility_pre` | 0.646 | 0.014–1.000 | 0.610 |

- Per patient: HUP235, HUP238 and HUP250 0.93–1.00 in every seizure; HUP224 0.95–0.99;
  HUP257 0.89–0.99; **HUP230 0.55–0.72** in all four seizures (the most fragile channels, LO/LM
  and LB, are not the experts' onset channels).
- Spontaneous seizures median 0.992 (n = 18), stimulation-induced 0.935 (n = 6); with low-voltage
  fast activity 0.945 (n = 9), without 0.995 (n = 15) (`stim_induced`, `lvfa` of `recordings.csv`).
- **With seizure-dynasd's detected onsets** (the workflow default when `onsets` is linked;
  median absolute onset error 5.0 s, up to 59.5 s): `fragility_post` median AUC 0.981 (range
  0.384–1.000), `fragility_pre` 0.924, post − pre 0.794. The worst clip falls to 0.38 when the
  detected onset is far from the expert one. Against v0.1.1 with the same onsets, one clip
  (HUP250, 151 channels) moves from 0.966 to 0.926: its post-onset means are within ~1e-3 of each
  other around the onset channels, which is the size of v0.1.x's rounding bias; every other
  clip's AUC moves by ≤ 0.004.
- Through seizure-consensus with fragility alone (v0.1.0): precision 0.35, recall 0.87, F1 0.46
  (medians) for the top 10 channels; v0.2.0 changes no clip's top 10.
- **Read with care.** The most fragile channels are largely the same in every seizure of a
  patient, and already before onset (pre-onset AUC 0.934), so fragility here behaves much like a
  stable property of the channel (its signal and its place in the network) that coincides with
  the expert onset channels, rather than a change at onset: the post − pre difference is barely
  better than chance. 24 seizures from 6 patients are 6 independent tests, not 24.
- 512 Hz vs 1024 Hz analysis rate on four clips (stand-in preprocessing): median AUC 0.970 vs 0.976,
  same top channels.

## Runtime

`computeTypes: [standard]`, 4 vCPU / 8 GiB (4 GiB is not a valid Fargate size for 4 vCPU), 4 h
timeout. Cost per window grows with the cube of the number of channels (100 complex solves, or
for small models K ≤ 120 real matrix products, per window); the number of windows is fixed by
`WINDOW_S`/`STEP_MS` (719 by default), not by the sampling rate. One recording at a time;
channels (band-pass, resampling) and then windows in `THREADS` forked workers with one OpenBLAS
thread each.

**v0.2.0 against v0.1.1**, 3 clips (HUP224 186 ch, HUP230 95 ch, HUP250 151 ch), one run at a
time on an M-series Mac, `--cpus 4 --memory 8g`, `THREADS=4`:

| | wall (per clip) | CPU | cgroup `memory.peak` |
|---|---|---|---|
| amd64 image under emulation, v0.1.1 | 266 s (141, 35, 91 s) | 987 s | 2.41 GB |
| amd64 image under emulation, v0.2.0 | 246 s (126, 41, 80 s) | 966 s | 2.36 GB |
| arm64 natively, v0.1.1 | 116 s (62, 16, 38 s) | 425 s | 2.11 GB |
| arm64 natively, v0.2.0 | 90 s (44, 18, 27 s) | 343 s | 2.15 GB |

The 204-channel clip alone (amd64, v0.2.0): 217 s, `memory.peak` 2.64 GB.

Time per window, natively, one thread, every 3rd window of a clip (240 windows):

| clip (channels) | `fast` (v0.1.x) | `stable`, solve only | `stable` (default: series where ‖A‖₂ ≤ 0.7) | windows on the series |
|---|---|---|---|---|
| HUP238 (72) | 0.031 s | 0.038 s | 0.038 s | 0 % |
| HUP230 (95) | 0.061 s | 0.072 s | 0.072 s | 0 % |
| HUP250 (151) | 0.174 s | 0.200 s | **0.107 s** | 78 % (K median 13) |
| HUP224 (186) | 0.294 s | 0.335 s | **0.201 s** | 68 % (K median 15) |
| HUP224 stim (204) | 0.384 s | 0.436 s | **0.313 s** | 49 % (K median 27) |

The accurate formulation by itself costs ~15–25 % per window (more element-wise work per trial
eigenvalue); the Neumann series wins that back and more on the large recordings, where nearly
all the time goes, and is not used on the small ones (their models have eigenvalues near the
unit circle). The band-pass in second-order sections is ~4x slower per channel than in (b, a)
form; filtering the channels in the workers more than pays for it. Overall: ~20 % faster
natively on the 3 clips, ~8 % under emulation, small recordings (< 100 channels) ~10–20 %
slower, large ones 10–35 % faster.

**What else was tried for speed.** With the cancellation gone, the rejected approaches were
re-checked against the 256-bit reference:

- *Neumann series* (adopted, above): as accurate as the solve (≤ 1.6e-15 relative) when
  truncated by the bound q^(K+1)/(1 − q); it only converges usefully for ‖A‖₂ well below 1.
  Beyond q = 0.7 (K > 120) the solve is as fast.
- *One eigen-decomposition per window with closed forms per trial eigenvalue*: not adopted.
  The per-(channel, λ) sums become quadratic forms in the eigen-basis, which brings back the
  cancellation this version removes (‖Im w‖² as a difference of rotated Gram terms loses
  ~1/ω² ≈ 10⁴ at ω = 0.01), and A is non-normal, so its eigenvectors can be ill-conditioned.
  It also saves little: the per-λ work stays O(n²) per channel row, as now.
- *One SVD per ridge search* (EZFragility's `ridge` recomputes `svd(xt)` at every bisection
  step): the ridge search is ~6 % of the time and only ~5 % of windows raise the penalty, so it
  was left as EZFragility's code (which keeps the penalties identical to EZFragility's).
- Where the time goes now (natively, 186 channels, solve path): about two thirds in the 100
  complex solves, one third in the element-wise row sums per trial eigenvalue (~1 ms each).

The image is ~1.44 GB unpacked, ~396 MB compressed (`rocker/r-ver:4.6.0` 1.36 GB plus ~80 MB of R
packages). Memory is ~0.5 GB per worker for a 200-channel recording (each forked worker's own
heap) plus the parent; the series adds ≤ ~60 MB per worker (the powers of A).

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
- Runtime grows with channels³ per window: a 200-channel montage takes ~2–4 min per seizure
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

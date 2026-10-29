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
| `onsets` | input (subfolder), optional | `seizure_onsets.csv` (`recording_id, onset_s, method`) | seizure-dynasd's detected onsets; used in place of the expert onset for the recordings it lists. Found in `INPUT_DIR/onsets/`, else anywhere under `INPUT_DIR`. |
| `fragility` | output | `fragility.csv`, `provenance.json` | One output port: the whole `OUTPUT_DIR`. |

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
| `THREADS` | `4` | Worker processes over windows (fork); match the vCPUs. One recording is processed at a time. |

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
  penalty had to be raised in {{RAISED}} of the windows.
- **The perturbation's eigenvalue.** The paper places the eigenvalue at some |λ| = 1 without
  saying which; EZFragility takes the λ that minimises the perturbation norm over a grid
  (`N_SEARCH`).
- **Normalisation.** Fragility is normalised per window ((max − norm)/max), so it ranks
  channels within a window; values are not comparable across recordings in absolute terms.
- **Preprocessing.** Li et al. analysed data around 1 kHz; EZFragility's reproduction band-passes
  0.5–150 Hz and works on clips of ±10–30 s. The app uses the same band and resamples to 512 Hz
  by default. On four demo clips, 512 Hz and 1024 Hz gave the same median AUC ({{AUC512}} vs
  {{AUC1024}}); 512 Hz halves the filtering work. Runtime is set by the number of channels and
  windows, not the sampling rate.
- **Summary.** The paper's figures show fragility heat maps and compare SOZ with non-SOZ
  distributions over time; collapsing to a mean before and after onset is this app's choice,
  so that the result fits a channel table.
- The CRAN package was **archived on 2026-10-09** because its reproduction vignette could not
  download its data during CRAN checks; the code is unchanged. The image installs 2.1.1 from
  the Posit Package Manager snapshot of 2026-06-23.

## Results on the demo data

{{RESULTS}}

## Runtime

{{RUNTIME}}

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

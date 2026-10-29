# Sample input

The app reads the `ieeg` output of **ieeg-prep** (a BIDS-iEEG derivative, see
`workflow-epilepsy-ieeg/docs/conventions.md`) and, optionally, the `onsets` output of
**seizure-dynasd**. Real data are not kept in this repository.

## Demo data (the 24 seizure clips used for the README numbers)

1. Fetch the Multi-Expert Seizure Annotation dataset (Pennsieve Discover 530 v3,
   doi:10.26275/dely-fzwx, CC BY 4.0): in `workflow-epilepsy-ieeg`, run
   `scripts/fetch-sample-data.sh` (6 HUP patients, 24 SEEG clips at 1024 Hz, ~2.4 GB).
2. Run ieeg-prep on `data/seizures/` (bipolar montage, notch, expert labels). Its `ieeg`
   output lands in `workflow-epilepsy-ieeg/data/work/prep/ieeg/`.
3. Optionally run seizure-dynasd; its `onsets` output lands in `data/work/dynasd/`.
4. Run this app with the two folders mounted as the input ports:

```sh
W=../workflow-epilepsy-ieeg/data/work
docker run --rm --platform linux/amd64 --memory 4g --cpus 4 --env-file dev.env \
  -v "$(realpath $W/prep/ieeg)":/data/input/ieeg:ro \
  -v "$(realpath $W/dynasd)":/data/input/onsets:ro \
  -v "$PWD/data/output":/data/output processor-seizure-fragility:dev
```

Leave out the `onsets` mount to use the expert onset (`recordings.csv` `onset_s`).

## Your own data

Any BIDS-iEEG dataset works without ieeg-prep: the app finds `*_ieeg.edf` files, reads
`_channels.tsv` (good SEEG/ECOG channels) and takes the onset from the first `seizure` (or
`ictal`) row of `_events.tsv`. Loose EDF files work too when their onsets come in through the
`onsets` port (`seizure_onsets.csv`: `recording_id` = the file name without `.edf`,
`onset_s`). Re-reference first (bipolar or common average): ieeg-prep does this.

## Synthetic fixture

`make fixture` writes a small synthetic derivative (2 recordings x 12 channels, 512 Hz, 2 min)
and a `seizure_onsets.csv` to `data/fixture/`; `make run` runs the image over it.

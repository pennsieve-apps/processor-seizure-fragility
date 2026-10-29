"""Write a small synthetic BIDS-iEEG derivative, as ieeg-prep's `ieeg` port would (see
workflow-epilepsy-ieeg/docs/conventions.md), for `make run`. Standard library + NumPy.

    python3 tests/make_fixture.py data/fixture [--recordings 2] [--channels 12] [--seconds 120]

Writes <out>/ieeg/{dataset_description.json, recordings.csv, sub-01/ses-01/ieeg/...} with
EDF+ signals (16 bit), _channels.tsv (onset labels on the first two channels), _events.tsv
(seizure, baseline) and <out>/onsets/seizure_onsets.csv (a detected onset 0.5 s after the
expert one, as seizure-dynasd would write).

Each recording: bipolar channels of filtered noise; from the onset on, channels 1-2 carry a
growing 8 Hz rhythm (low-voltage fast activity is not modelled) that spreads to channel 3.
"""
import argparse
import csv
import json
import os

import numpy as np


def write_edf(path, x, fs, labels):
    n, total = x.shape
    nrec = int(np.ceil(total / fs))
    x = np.pad(x, ((0, 0), (0, nrec * fs - total)))
    pmax = np.ceil(np.maximum(np.abs(x).max(axis=1), 1.0) * 1.1)
    ns = n + 1
    labels = list(labels) + ["EDF Annotations"]
    spr = [fs] * n + [30]

    def f(v, w):
        return str(v)[:w].ljust(w).encode("ascii")

    h = b"".join([f(0, 8), f("X X X X", 80), f("Startdate 01-JAN-2020 X X X", 80), f("01.01.20", 8),
                  f("00.00.00", 8), f(256 * (ns + 1), 8), f("EDF+C", 44), f(nrec, 8), f(1, 8), f(ns, 4)])
    h += b"".join(f(l, 16) for l in labels) + f("", 80) * ns
    h += b"".join(f("uV", 8) for _ in range(n)) + f("", 8)
    h += b"".join(f(-int(v), 8) for v in pmax) + f(-1, 8)
    h += b"".join(f(int(v), 8) for v in pmax) + f(1, 8)
    h += f(-32768, 8) * ns + f(32767, 8) * ns + f("", 80) * ns
    h += b"".join(f(v, 8) for v in spr) + f("", 32) * ns
    dig = np.round((x + pmax[:, None]) / (2 * pmax[:, None]) * 65535 - 32768).astype("<i2")
    with open(path, "wb") as fo:
        fo.write(h)
        for r in range(nrec):
            fo.write(dig[:, r * fs:(r + 1) * fs].tobytes())
            tal = f"+{r}".encode() + bytes([20, 20, 0])
            fo.write(tal.ljust(60, b"\0"))


def recording(rng, n_ch, fs, dur, onset):
    t = np.arange(int(dur * fs)) / fs
    noise = rng.standard_normal((n_ch, t.size))
    k = np.ones(5) / 5
    x = np.array([np.convolve(v, k, mode="same") for v in noise]) * 30
    x += 10 * np.sin(2 * np.pi * 1.5 * t)[None, :] * (np.arange(n_ch)[:, None] % 3)
    g = np.clip((t - onset) / 5, 0, 1)
    for i, a in ((0, 200), (1, 200), (2, 60)):
        x[i] += a * g * np.sin(2 * np.pi * 8 * t + i)
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--recordings", type=int, default=2)
    ap.add_argument("--channels", type=int, default=12)
    ap.add_argument("--seconds", type=float, default=120)
    ap.add_argument("--fs", type=int, default=512)
    a = ap.parse_args()
    rng = np.random.default_rng(42)
    root = os.path.join(a.out, "ieeg")
    os.makedirs(root, exist_ok=True)
    names = [f"LA{i + 1:02d}-LA{i + 2:02d}" for i in range(a.channels)]
    rows, onsets = [], []
    for k in range(1, a.recordings + 1):
        rid = f"sub-01_ses-01_task-ictal{k}_run-01"
        d = os.path.join(root, "sub-01", "ses-01", "ieeg")
        os.makedirs(d, exist_ok=True)
        onset = a.seconds * 0.6 + k
        x = recording(rng, a.channels, a.fs, a.seconds, onset)
        edf = os.path.join(d, f"{rid}_desc-preproc_ieeg.edf")
        write_edf(edf, x, a.fs, names)
        with open(os.path.join(d, f"{rid}_desc-preproc_channels.tsv"), "w", encoding="utf-8") as fo:
            fo.write("name\ttype\tunits\tlow_cutoff\thigh_cutoff\tsampling_frequency\tstatus\tstatus_description\tcontacts\tsoz\tresected\tonset\tspread\n")
            for i, nm in enumerate(names):
                fo.write(f"{nm}\tSEEG\tµV\tn/a\tn/a\t{a.fs}\tgood\tn/a\t{nm.replace('-', '|')}\t\t\t"
                         f"{'true' if i < 2 else 'false'}\t{'true' if i == 2 else 'false'}\n")
        with open(os.path.join(d, f"{rid}_desc-preproc_events.tsv"), "w") as fo:
            fo.write("onset\tduration\ttrial_type\tsource\n")
            fo.write(f"{onset:.3f}\t20.000\tseizure\tfixture\n{onset - 60:.3f}\t30.000\tbaseline\tfixture\n")
        with open(os.path.join(d, f"{rid}_desc-preproc_ieeg.json"), "w") as fo:
            json.dump({"SamplingFrequency": a.fs, "iEEGReference": "bipolar", "PowerLineFrequency": 60}, fo, indent=2)
        rows.append({"recording_id": rid, "subject": "01", "session": "01", "task": f"ictal{k}", "run": "01",
                     "path": os.path.relpath(edf, root), "sfreq": a.fs, "n_channels": a.channels,
                     "n_good": a.channels, "duration_s": a.seconds, "montage": "bipolar", "kind": "ictal",
                     "onset_s": round(onset, 3), "baseline_start_s": round(onset - 60, 3),
                     "baseline_end_s": round(onset - 30, 3), "source_format": "bids-edf"})
        onsets.append({"recording_id": rid, "onset_s": round(onset + 0.5, 3), "method": "fixture"})
    with open(os.path.join(root, "recordings.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(root, "dataset_description.json"), "w") as fo:
        json.dump({"Name": "synthetic fixture", "BIDSVersion": "1.9.0", "DatasetType": "derivative",
                   "GeneratedBy": [{"Name": "tests/make_fixture.py"}]}, fo, indent=2)
    os.makedirs(os.path.join(a.out, "onsets"), exist_ok=True)
    with open(os.path.join(a.out, "onsets", "seizure_onsets.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["recording_id", "onset_s", "method"])
        w.writeheader()
        w.writerows(onsets)
    print(f"wrote {len(rows)} recording(s) x {a.channels} channels to {a.out}/ieeg and onsets to {a.out}/onsets")


if __name__ == "__main__":
    main()

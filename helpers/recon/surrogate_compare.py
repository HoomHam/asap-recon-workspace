"""Cross-check the respiratory surrogates: do signal / pneumotach / diaphragm
agree on the breathing period? (Stage-1 verification.)

Plots each available surrogate's volume-vs-time (zoomed to a few breaths), the
derived phase, and the dominant breathing period from the FFT of each curve.

Usage: ../.venv/bin/python surrogate_compare.py <dump_dir> [--diaphragm]
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import binning


def dominant_period(t, v):
    """Breathing period (s) from the FFT peak of the (detrended) curve."""
    good = np.isfinite(v)
    v = v[good] - np.mean(v[good])
    tt = t[good]
    dt = np.median(np.diff(tt))
    f = np.fft.rfftfreq(len(v), dt)
    P = np.abs(np.fft.rfft(v))
    P[f < 0.05] = 0          # ignore < 0.05 Hz drift
    fpk = f[np.argmax(P)]
    return 1.0 / fpk if fpk > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--diaphragm", action="store_true",
                    help="also compute + include the diaphragm surrogate (slow)")
    args = ap.parse_args()

    meta = json.load(open(os.path.join(args.dump_dir, "meta.json")))
    ilvtime = np.load(os.path.join(args.dump_dir, "ilvtime.npy"))
    exclude = np.load(os.path.join(args.dump_dir, "exclude_ilv.npy"))

    curves = {}
    for name, fn in [("signal", "ilvvol_signal.npy"), ("pneumo", "ilvvol_pneumo.npy")]:
        p = os.path.join(args.dump_dir, fn)
        if os.path.exists(p):
            curves[name] = np.load(p)
    if args.diaphragm:
        import surrogates
        curves["diaphragm"] = surrogates.diaphragm_surrogate(args.dump_dir)

    B = meta["nbins"]
    fig, axes = plt.subplots(3, 1, figsize=(11, 9))
    zoom = ilvtime < min(40, ilvtime[-1])           # first ~40 s
    for name, v in curves.items():
        per = dominant_period(ilvtime, v)
        axes[0].plot(ilvtime[zoom], v[zoom], label=f"{name} (T~{per:.1f}s)")
        phi, _ = binning.membership_from_surrogate(v, n_bins=B, exclude=exclude)
        axes[1].plot(ilvtime[zoom], phi[zoom], ".", ms=2, label=name)
        occ = binning.soft_membership(phi, n_bins=B).sum(1)
        axes[2].plot(np.arange(B), occ, "o-", label=name)
        print(f"{name:10s} breathing period ~ {per:.2f} s")

    axes[0].set(title="surrogate volume vs time (first 40 s)", xlabel="t (s)",
                ylabel="normalized volume"); axes[0].legend(fontsize=8)
    axes[1].set(title="derived respiratory phase", xlabel="t (s)", ylabel="phase")
    axes[1].legend(fontsize=8)
    axes[2].set(title="bin occupancy (effective interleaves)", xlabel="bin",
                ylabel="sum membership"); axes[2].legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(args.dump_dir, "surrogate_compare.png")
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

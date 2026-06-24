"""Validate the diaphragm surrogate + dump its position curves.

Runs the CS-nav diaphragm extractor (surrogates.diaphragm_curve), saves the raw
edge positions and the per-interleave volume, and plots:
  1. diaphragm edge position vs time (raw + smoothed), full scan
  2. same, zoomed to the first ~40 s (individual breaths visible)
  3. normalized diaphragm volume overlaid with signal (+ pneumo) for cross-check

Usage:
    ../.venv/bin/python diaphragm_check.py <dump_dir> [--method cg|wavelet]
        [--win-ilv 26] [--nav-n 64]
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import surrogates


def period(t, v):
    good = np.isfinite(v); v = v[good] - np.mean(v[good]); tt = t[good]
    dt = np.median(np.diff(tt)); f = np.fft.rfftfreq(len(v), dt)
    P = np.abs(np.fft.rfft(v)); P[f < 0.05] = 0
    fpk = f[np.argmax(P)]
    return 1.0 / fpk if fpk > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--method", choices=["cg", "wavelet"], default="cg")
    ap.add_argument("--metric", choices=["centroid", "edge"], default="edge")
    ap.add_argument("--prefer", choices=["hi", "lo", "auto"], default="auto",
                    help="which lung boundary to use as binning surrogate: "
                         "auto=corr-based (default, lo wins on 025JC, period 3.52s clean); "
                         "hi=dome (nav_movie display line -- clips FOV at inspiration); "
                         "lo=apex (force)")
    ap.add_argument("--win-ilv", type=int, default=20)
    ap.add_argument("--smooth-win", type=int, default=5)
    ap.add_argument("--nav-n", type=int, default=64)
    args = ap.parse_args()

    d = surrogates.diaphragm_curve(args.dump_dir, win_ilv=args.win_ilv,
                                   nav_N=args.nav_n, method=args.method,
                                   metric=args.metric, prefer=args.prefer,
                                   smooth_win=args.smooth_win)
    np.save(os.path.join(args.dump_dir, "diaphragm_times.npy"), d["times"])
    np.save(os.path.join(args.dump_dir, "diaphragm_pos_raw.npy"), d["pos_raw"])
    np.save(os.path.join(args.dump_dir, "diaphragm_pos_smooth.npy"), d["pos_smooth"])
    np.save(os.path.join(args.dump_dir, "ilvvol_diaphragm.npy"), d["ilvvol"])

    t, ilvt = d["times"], d["ilvtime"]
    per = period(t, d["pos_smooth"])
    print(f"diaphragm breathing period ~ {per:.2f} s   "
          f"edge range {d['pos_raw'].min():.1f}-{d['pos_raw'].max():.1f} nav-vox")

    fig, ax = plt.subplots(3, 1, figsize=(11, 9))
    ax[0].plot(t, d["pos_raw"], ".", ms=3, alpha=.5, label="raw edge")
    ax[0].plot(t, d["pos_smooth"], "-", lw=1.2, label="smoothed")
    ax[0].set(title=f"diaphragm edge position vs time (nav={args.method}, "
              f"T~{per:.1f}s)", xlabel="t (s)", ylabel="position (nav-vox)")
    ax[0].legend(fontsize=8); ax[0].invert_yaxis()

    z = t < min(40, t[-1])
    ax[1].plot(t[z], d["pos_raw"][z], ".-", ms=4, label="raw")
    ax[1].plot(t[z], d["pos_smooth"][z], "-", lw=1.5, label="smoothed")
    ax[1].set(title="zoom: first 40 s (individual breaths)", xlabel="t (s)",
              ylabel="position"); ax[1].legend(fontsize=8); ax[1].invert_yaxis()

    ax[2].plot(ilvt, d["ilvvol"], "-", lw=1, label=f"diaphragm (T~{per:.1f}s)")
    for nm, fn in [("signal", "ilvvol_signal.npy"), ("pneumo", "ilvvol_pneumo.npy")]:
        p = os.path.join(args.dump_dir, fn)
        if os.path.exists(p):
            v = np.load(p)
            ax[2].plot(ilvt, v, "-", lw=.8, alpha=.7,
                       label=f"{nm} (T~{period(ilvt, v):.1f}s)")
    ax[2].set(title="normalized surrogate cross-check (full scan)", xlabel="t (s)",
              ylabel="norm volume", xlim=(0, min(40, ilvt[-1])))
    ax[2].legend(fontsize=8)

    fig.tight_layout()
    out = os.path.join(args.dump_dir, "diaphragm_check.png")
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

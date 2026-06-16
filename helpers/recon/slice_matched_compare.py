"""Slice-MATCHED ours/BART/Lustig at best-detail lambda — same physical planes.

Correction (2026-06-15, Hooman): the cross-hatch is a REAL ACR feature (the
resolution / grid insert), not a BART artifact. Earlier montages compared
[:, :, 50] of each volume; after orientation those ARE ~the same plane (slicewise
corr ~0.95), so the right question is which pipeline RESOLVES the real features
at MATCHED slices — over-smoothing (low lfCV) that erases the grid is WORSE, not
better.

This NEW script reconstructs each pipeline at a detail-preserving lambda, orients
BART and Lustig onto our frame (full 48-orientation, whole-volume correlation),
and lays out the SAME z-slices as columns, pipelines as rows, so matched planes
sit above each other. Pick lambda per pipeline with --our-t / --bart-l / --lustig
(reads the existing Lustig TV sweep .mat).

    ../.venv/bin/python slice_matched_compare.py ../../data/v3_fov250/recon_io \
        --bart /Users/hoomham/bin/bart-src/bart \
        --our-t 0.003 --bart-l 0.0003 --lustig-tv 0.001 \
        --slices 38 44 50 56 62

Outputs: recon_io/slice_matched_montage.png
"""

import argparse
import itertools
import os
import subprocess

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sigpy as sp

import asap_recon as ar
import cs_recon as cs
from bartio import readcfl
from bart_compare import our_setup

N = cs.N


def best_orient_full(vol, ref):
    """48-orientation search maximizing WHOLE-VOLUME correlation -> aligns z too."""
    a = np.abs(ref); a = (a - a.mean()) / (a.std() + 1e-12)
    best, best_c = np.abs(vol), -np.inf
    for perm in itertools.permutations(range(3)):
        vp = np.transpose(np.abs(vol), perm)
        for fl in itertools.product([False, True], repeat=3):
            v = vp
            for ax, do in enumerate(fl):
                if do:
                    v = np.flip(v, axis=ax)
            b = (v - v.mean()) / (v.std() + 1e-12)
            cc = float(np.mean(a * b))
            if cc > best_c:
                best_c, best = cc, v
    return best, best_c


def norm01(v):
    p = np.percentile(np.abs(v), 99.5)
    return np.clip(np.abs(v) / (p + 1e-12), 0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("--bart", default="/Users/hoomham/bin/bart-src/bart")
    ap.add_argument("--our-t", type=float, default=0.003)
    ap.add_argument("--bart-l", type=float, default=0.0003)
    ap.add_argument("--lustig-tv", type=float, default=0.001)
    ap.add_argument("--slices", type=int, nargs="+", default=[38, 44, 50, 56, 62])
    ap.add_argument("--iters", type=int, default=100)
    a = ap.parse_args()
    D = a.recon_io

    o = our_setup(D)
    ref = np.abs(o["img_cg"])

    # ours wavelet at detail-preserving t
    W = sp.linop.Wavelet(o["A"].ishape, wave_name="db4")
    tref = float(np.percentile(np.abs(W(o["img_cg"])), 99))
    ours = np.abs(cs.wavelet_recon(o["A"], o["y_n"], a.our_t * tref, a.iters))
    print(f"ours wavelet t={a.our_t}")

    # BART wavelet at bart_l (ADMM, as established)
    bout = os.path.join(D, "bart", f"bart_W_l{a.bart_l:g}")
    if not os.path.exists(bout + ".cfl"):
        subprocess.run([a.bart, "pics", "-i", str(a.iters), "-m",
                        "-t", os.path.join(D, "bart", "traj"),
                        "-R", f"W:7:0:{a.bart_l}",
                        os.path.join(D, "bart", "ksp"),
                        os.path.join(D, "bart", "sens"), bout],
                       check=True, capture_output=True, text=True)
    bart, ccb = best_orient_full(np.squeeze(readcfl(bout)), ref)
    print(f"BART wavelet l={a.bart_l}  orient cc {ccb:.3f}")

    # Lustig TV from the sweep .mat
    mat = os.path.join(D, "lustig", "lustig_tv_sweep.mat")
    g = sio.loadmat(mat); gf = np.asarray(g["gas_final"], float)
    tvW = np.asarray(g["tvW"], float).ravel()
    li = int(np.argmin(np.abs(tvW - a.lustig_tv)))
    lustig, ccl = best_orient_full(gf[li], ref)
    print(f"Lustig TV w={tvW[li]:g}  orient cc {ccl:.3f}")

    rows = [(f"ours wav t{a.our_t:g}", ours),
            (f"BART wav l{a.bart_l:g}", bart),
            (f"Lustig TV w{tvW[li]:g}", lustig)]
    nz = len(a.slices)
    fig, axes = plt.subplots(3, nz, figsize=(3 * nz, 9.5))
    for ri, (label, vol) in enumerate(rows):
        vn = norm01(vol)
        for ci, z in enumerate(a.slices):
            ax = axes[ri, ci]
            ax.imshow(vn[:, :, z], cmap="gray", vmin=0, vmax=1)
            ax.axis("off")
            if ri == 0:
                ax.set_title(f"z={z}", fontsize=10)
            if ci == 0:
                ax.text(-0.12, 0.5, label, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=11)
    fig.suptitle("Matched physical slices (same z, oriented to ours). "
                 "Resolution-insert grid is REAL ACR structure — judge who RESOLVES it.",
                 fontsize=12)
    fig.tight_layout()
    png = os.path.join(D, "slice_matched_montage.png")
    fig.savefig(png, dpi=110)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()

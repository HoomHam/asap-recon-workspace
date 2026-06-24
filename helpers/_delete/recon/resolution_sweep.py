"""Per-pipeline lambda sweep at the RESOLUTION-INSERT slice — pick best by eye.

Correction context (2026-06-15): the cross-hatch is the REAL ACR resolution
insert, not an artifact. So 'best lambda' = the one that RESOLVES the insert
without drowning it in noise — NOT the lowest-lfCV (that just blurs it away).

Rows = pipeline (ours wavelet / BART wavelet / Lustig TV), cols = increasing
regularization. Plus ours CG-20 (unsmoothed) as the max-resolution anchor. One
matched slice z (oriented to ours) so every panel is the SAME physical plane.

    ../.venv/bin/python resolution_sweep.py ../../data/v3_fov250/recon_io \
        --bart /Users/hoomham/bin/bart-src/bart --z 50

Outputs: recon_io/resolution_sweep_z<z>.png
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
OUR_T = [0.001, 0.003, 0.01, 0.03]
BART_L = [0.0001, 0.0003, 0.001, 0.003, 0.01]


def best_orient_full(vol, ref):
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
    return best


def n01(sl):
    p = np.percentile(sl, 99.5)
    return np.clip(sl / (p + 1e-12), 0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("--bart", default="/Users/hoomham/bin/bart-src/bart")
    ap.add_argument("--z", type=int, default=50)
    ap.add_argument("--iters", type=int, default=100)
    a = ap.parse_args()
    D = a.recon_io
    z = a.z

    o = our_setup(D)
    ref = np.abs(o["img_cg"])
    W = sp.linop.Wavelet(o["A"].ishape, wave_name="db4")
    tref = float(np.percentile(np.abs(W(o["img_cg"])), 99))

    # build rows: list of (row_label, [(col_label, slice2d), ...])
    rows = []

    # ours: CG anchor + wavelet sweep
    cols = [("CG-20 (no reg)", ref[:, :, z])]
    for t in OUR_T:
        img = np.abs(cs.wavelet_recon(o["A"], o["y_n"], t * tref, a.iters))
        cols.append((f"t={t:g}", img[:, :, z]))
    rows.append(("ours wavelet", cols))
    print("ours done")

    # BART wavelet sweep (ADMM); reuse cfls if present
    cols = [("", np.zeros((N, N)))]  # pad col0 to align with ours CG anchor
    for lam in BART_L:
        out = os.path.join(D, "bart", f"bart_W_l{lam:g}")
        if not os.path.exists(out + ".cfl"):
            subprocess.run([a.bart, "pics", "-i", str(a.iters), "-m",
                            "-t", os.path.join(D, "bart", "traj"),
                            "-R", f"W:7:0:{lam}",
                            os.path.join(D, "bart", "ksp"),
                            os.path.join(D, "bart", "sens"), out],
                           check=True, capture_output=True, text=True)
        v = best_orient_full(np.squeeze(readcfl(out)), ref)
        cols.append((f"l={lam:g}", v[:, :, z]))
    rows.append(("BART wavelet", cols))
    print("BART done")

    # Lustig TV sweep from .mat
    g = sio.loadmat(os.path.join(D, "lustig", "lustig_tv_sweep.mat"))
    gf = np.asarray(g["gas_final"], float); tvW = np.asarray(g["tvW"], float).ravel()
    cols = [("", np.zeros((N, N)))]
    for i in range(len(tvW)):
        v = best_orient_full(gf[i], ref)
        cols.append((f"w={tvW[i]:g}", v[:, :, z]))
    rows.append(("Lustig TV", cols))
    print("Lustig done")

    ncol = max(len(c) for _, c in rows)
    fig, axes = plt.subplots(len(rows), ncol, figsize=(2.6 * ncol, 2.6 * len(rows)))
    for ri, (rlabel, cols) in enumerate(rows):
        for ci in range(ncol):
            ax = axes[ri, ci]; ax.axis("off")
            if ci < len(cols):
                clabel, sl = cols[ci]
                if sl.max() > 0:
                    ax.imshow(n01(sl), cmap="gray", vmin=0, vmax=1)
                if ri == 0 and ci > 0:
                    ax.set_title(clabel, fontsize=9)
                elif ci > 0:
                    ax.set_title(clabel, fontsize=9)
            if ci == 0:
                ax.text(-0.15, 0.5, rlabel, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=11)
    fig.suptitle(f"Resolution-insert slice z={z} (matched). Best lambda = RESOLVES "
                 f"the grid without noise. Col1 = ours CG-20 anchor.", fontsize=12)
    fig.tight_layout()
    png = os.path.join(D, f"resolution_sweep_z{z}.png")
    fig.savefig(png, dpi=115)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()

"""Render our recons in Faraz's montage format, next to his result.

Replicates his MATLAB figure exactly:
    imagesc([rot90(squeeze(img_gp(:,16,:))), ..., rot90(squeeze(img_gp(:,63,:)))])
    6 rows x 8 cols, y-slices 16..63 (1-based), each panel rot90 of (:,j,:).

Our volumes (100^3 over the same FOV) are auto-oriented against his 80^3
volume (search over axis permutations + flips, correlation on a resampled
copy), then sliced at fraction-matched positions.

Usage:
    .venv/bin/python faraz_montage.py <recon_io folder> [faraz_recon.mat]
Writes montage_faraz.png, montage_cg.png, montage_steve.png (+ combined).
"""

import itertools
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.ndimage import zoom

import asap_recon as ar
from steve_kernel_numpy import steve_recon

ROWS, COLS = 6, 8
SLICES_1B = list(range(16, 64))               # his 22-6 .. 69-6, 1-based


def montage(vol, slice_idx_0b):
    """His layout: panels rot90(vol[:, j, :]), ROWS x COLS."""
    panels = [np.rot90(vol[:, j, :]) for j in slice_idx_0b]
    rows = [np.hstack(panels[r * COLS:(r + 1) * COLS]) for r in range(ROWS)]
    return np.vstack(rows)


def orient_to_match(vol, ref):
    """Find axis permutation + flips of |vol| best matching |ref| (volume
    correlation on a resampled copy). Returns transformed full-res volume."""
    a = np.abs(vol)
    r = np.abs(ref)
    a_small = zoom(a, np.array(r.shape) / np.array(a.shape), order=1)
    best = (-2, None)
    for perm in itertools.permutations(range(3)):
        for flips in itertools.product((False, True), repeat=3):
            t = np.transpose(a_small, perm)
            for ax, f in enumerate(flips):
                if f:
                    t = np.flip(t, axis=ax)
            c = np.corrcoef(t.ravel(), r.ravel())[0, 1]
            if c > best[0]:
                best = (c, (perm, flips))
    c, (perm, flips) = best
    print(f"  orientation: perm {perm} flips {flips}  (vol corr {c:.4f})")
    out = np.transpose(np.abs(vol), perm)
    for ax, f in enumerate(flips):
        if f:
            out = np.flip(out, axis=ax)
    return out


def fraction_matched_slices(n_ours, n_ref=80):
    """Map his 1-based slices 16..63 (center 40.5) to our grid (center n/2+.5)."""
    return [int(round((j - 0.5 - n_ref / 2) * n_ours / n_ref + n_ours / 2 - 0.5))
            for j in SLICES_1B]


def save_montage(img, path, title):
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.imshow(img, cmap="gray")
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"  wrote {path}")


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    matpath = sys.argv[2] if len(sys.argv) > 2 else os.path.join(folder, "faraz", "faraz_recon.mat")

    vol_f = np.abs(loadmat(matpath, squeeze_me=True)["img_gp"])
    print(f"faraz img_gp {vol_f.shape}")

    d = ar.load_steve_npy(folder)
    meta = json.load(open(os.path.join(folder, "meta.json")))

    print("CG recon ...")
    vol_cg = ar.recon(d["traj"], d["acq"], method="cg", cg_iters=20)
    print("steve-equiv recon ...")
    vol_st = steve_recon(d["traj"], d["acq"], npts=meta["npts"], MS=meta["MS"],
                         IS=meta["IS"], smoothing=meta["gplb"], axes="xyz")

    print("orienting CG to faraz:")
    vol_cg = orient_to_match(vol_cg, vol_f)
    print("orienting steve-equiv to faraz:")
    vol_st = orient_to_match(vol_st, vol_f)

    idx_f = [j - 1 for j in SLICES_1B]
    idx_o = fraction_matched_slices(vol_cg.shape[1], vol_f.shape[1])

    m_f = montage(vol_f, idx_f)
    m_cg = montage(vol_cg, idx_o)
    m_st = montage(vol_st, idx_o)

    save_montage(m_f, os.path.join(folder, "montage_faraz.png"), "Faraz (MATLAB KB+DCF)")
    save_montage(m_cg, os.path.join(folder, "montage_cg.png"), "Ours (FINUFFT CG)")
    save_montage(m_st, os.path.join(folder, "montage_steve.png"), "Steve-equiv (numpy kernel)")

    fig, axes = plt.subplots(3, 1, figsize=(16, 34))
    for ax, (img, t) in zip(axes, [(m_f, "Faraz (MATLAB KB+DCF)"),
                                   (m_cg, "Ours (FINUFFT CG)"),
                                   (m_st, "Steve-equiv (numpy kernel)")]):
        ax.imshow(img, cmap="gray")
        ax.set_title(t, fontsize=14)
        ax.axis("off")
    fig.tight_layout()
    out = os.path.join(folder, "montage_all.png")
    fig.savefig(out, dpi=100)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

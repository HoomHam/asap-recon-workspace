"""Quick visual+metric check: Lustig MATLAB CS (tv01g01.mat) vs our CS.

NEW standalone script — imports existing helpers read-only, modifies nothing.
Reproduces our wavelet-CS volume at the requested t_rel(s) using cs_recon
primitives, loads the final fnlCg iteration from the MATLAB .mat, puts both
on a shared p99.5 / [0,1] axis, prints cg_tune metrics, writes a montage.

Run:
    ../.venv/bin/python lustig_compare.py \
        ../../data/v3_fov250/recon_io \
        ../../codes/2025-09-24_ACR/tv01g01.mat \
        --t 0.003 0.01 --max-iter 100
"""

import argparse
import os

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sigpy as sp

import asap_recon as ar
from cg_tune import metrics
import cs_recon as cs


def our_cs_volume(folder, t_rel, max_iter):
    """Full 3D wavelet-CS volume at one t_rel — mirrors cs_recon.main() setup
    exactly (same DCF precond, normalization, threshold reference) but returns
    the whole |volume| instead of one slice."""
    d = ar.load_steve_npy(folder)
    traj_rad = ar.grid_to_radians(
        ar.tile_traj(d["traj"].astype(float), len(d["acq"])))
    y = np.ascontiguousarray(d["acq"], dtype=np.complex128)

    A_raw = cs.FinufftForward(traj_rad, cs.N)
    M = len(y)
    dens = np.abs(A_raw(A_raw.H(np.ones(M, dtype=complex))))
    w = 1.0 / np.clip(dens, dens.max() * 1e-4, None)
    w /= w.mean()
    Aw_raw = sp.linop.Multiply((M,), np.sqrt(w)) * A_raw
    L = sp.app.MaxEig(Aw_raw.H * Aw_raw, dtype=np.complex128,
                      show_pbar=False).run()
    c = 1.0 / np.sqrt(L)
    A = sp.linop.Multiply((M,), c * np.sqrt(w)) * A_raw
    y_n = c * np.sqrt(w) * y

    img_cg = ar.recon(d["traj"], y, method="cg", cg_iters=20)
    W = sp.linop.Wavelet(A.ishape, wave_name="db4")
    t_ref = float(np.percentile(np.abs(W(img_cg)), 99))
    lam = t_rel * t_ref
    img = cs.wavelet_recon(A, y_n, lam, max_iter)
    return np.abs(img)


def load_lustig_final(mat_path):
    """Final fnlCg iteration from gas(15,100,100,100). Already |.|,
    max-normalized in the MATLAB script. The script also rot90'd each axis-0
    slice in the (dim1,dim2) plane for display — undo it (k=-1) so the volume
    sits in the raw recon frame, comparable to ours."""
    try:
        g = sio.loadmat(mat_path)["gas"]
    except NotImplementedError:  # v7.3 / HDF5
        import h5py
        g = np.array(h5py.File(mat_path, "r")["gas"]).transpose(3, 2, 1, 0)
    v = np.asarray(g[-1], dtype=float)
    return np.rot90(v, k=-1, axes=(1, 2))


def norm01(vol):
    p = np.percentile(vol, 99.5)
    return np.clip(vol / (p + 1e-12), 0, 1)


def montage_slices(vol, idxs, axis=2):
    out = []
    for i in idxs:
        sl = np.take(vol, i, axis=axis)
        out.append(sl)
    return np.concatenate(out, axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("lustig_mat")
    ap.add_argument("--t", type=float, nargs="+", default=[0.003, 0.01])
    ap.add_argument("--max-iter", type=int, default=100)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    lus = load_lustig_final(args.lustig_mat)
    print(f"Lustig tv01g01 final iter: {lus.shape}  metrics: {metrics(lus)}")

    ours = {}
    for t in args.t:
        v = our_cs_volume(args.recon_io, t, args.max_iter)
        ours[t] = v
        print(f"our CS wavelet t={t:g}: {v.shape}  metrics: {metrics(v)}")

    # central slices along each axis for the eye
    N = lus.shape[0]
    idxs = [N // 2 - 12, N // 2, N // 2 + 12]
    rows = [("Lustig tv01 (matlab)", norm01(lus))]
    for t in args.t:
        rows.append((f"ours wavelet t={t:g}", norm01(ours[t])))

    fig, axes = plt.subplots(len(rows), 3, figsize=(12, 4 * len(rows)))
    if len(rows) == 1:
        axes = axes[None, :]
    for r, (label, vol) in enumerate(rows):
        for c, ax_idx in enumerate([0, 1, 2]):
            sl = norm01(vol)
            img = montage_slices(vol, [idxs[c]], axis=ax_idx)
            axes[r, c].imshow(img, cmap="gray", vmin=0, vmax=1)
            axes[r, c].axis("off")
            axes[r, c].set_title(f"{label}  axis{ax_idx} s{idxs[c]}",
                                 fontsize=9)
    fig.tight_layout()
    out = args.out or os.path.join(
        os.path.dirname(args.lustig_mat), "compare_lustig_vs_ours.png")
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

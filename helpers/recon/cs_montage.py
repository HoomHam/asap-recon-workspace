"""Multi-slice montage: CS wavelet recon alongside Faraz (zoom-corrected) and
Steve-equiv — same slice layout as montage_zoomfix.png, separate output file.

Reuses maxeig and t_ref from cs_sweep_metrics.json (written by cs_recon.py)
so only the DCF weights and the one CS recon are recomputed.

Usage: .venv/bin/python cs_montage.py <recon_io folder> [--t-rel 0.003]
Output: <folder>/montage_cs_t<t_rel>.png
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat
import sigpy as sp

import asap_recon as ar
from steve_kernel_numpy import steve_recon
from faraz_montage import montage, orient_to_match, fraction_matched_slices, SLICES_1B
from faraz_zoom_check import center_zoom, FOV_MM
from cs_recon import FinufftForward, wavelet_recon, N


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", nargs="?", default=".")
    ap.add_argument("--t-rel", type=float, default=0.003)
    ap.add_argument("--max-iter", type=int, default=100)
    args = ap.parse_args()

    d = ar.load_steve_npy(args.folder)
    meta = json.load(open(os.path.join(args.folder, "meta.json")))
    cs_meta = json.load(open(os.path.join(args.folder, "cs_sweep_metrics.json")))
    y = np.ascontiguousarray(d["acq"], dtype=np.complex128)
    M = len(y)

    # Faraz, zoom-corrected to true geometry (alpha from his resizing formula)
    vol_f = np.abs(loadmat(os.path.join(args.folder, "faraz", "faraz_recon.mat"),
                           squeeze_me=True)["img_gp"])
    MS, IS = meta["MS"], meta["IS"]
    kdk = d["traj"] - MS / 2
    kmax_meas_invmm = np.linalg.norm(kdk, axis=1).max() / (MS / IS) / FOV_MM
    alpha = ((80 / 2 - 1) / FOV_MM) / kmax_meas_invmm
    vol_f_fix = center_zoom(vol_f, alpha, out_shape=vol_f.shape)
    print(f"alpha = {alpha:.4f} (x{1/alpha:.3f} magnification corrected)")

    print("steve-equiv recon (numpy kernel) ...")
    vol_st = np.abs(steve_recon(d["traj"], y, npts=meta["npts"], MS=MS, IS=IS,
                                smoothing=meta["gplb"], axes="xyz"))

    # CS wavelet recon, same construction as cs_recon.py main
    traj_rad = ar.grid_to_radians(ar.tile_traj(d["traj"].astype(float), M))
    A_raw = FinufftForward(traj_rad, N)
    dens = np.abs(A_raw(A_raw.H(np.ones(M, dtype=complex))))
    w = 1.0 / np.clip(dens, dens.max() * 1e-4, None)
    w /= w.mean()
    c = 1.0 / np.sqrt(cs_meta["maxeig"])
    A = sp.linop.Multiply((M,), c * np.sqrt(w)) * A_raw
    lam = args.t_rel * cs_meta["t_ref"]["wavelet"]
    print(f"CS wavelet recon (t_rel={args.t_rel:g}, lam={lam:.3e}, "
          f"{args.max_iter} iters) ...")
    vol_cs = np.abs(wavelet_recon(A, c * np.sqrt(w) * y, lam, args.max_iter))

    print("orienting CS volume to Faraz layout:")
    cs_d = orient_to_match(vol_cs, vol_f_fix)
    st_d = orient_to_match(vol_st, vol_f_fix)

    # volumes are in unrelated absolute units (mat / kernel / CS): normalize
    # each by its p99.5 and display on a shared [0, 1] color axis
    def nrm(v):
        return np.clip(v / np.percentile(v, 99.5), 0, 1)
    vol_f_fix, st_d, cs_d = nrm(vol_f_fix), nrm(st_d), nrm(cs_d)

    idx_f = [j - 1 for j in SLICES_1B]
    idx_o = fraction_matched_slices(N, 80)
    panels = [
        (montage(vol_f_fix, idx_f), "Faraz CORRECTED (true geometry)"),
        (montage(st_d, idx_o), "Steve-equiv (numpy kernel)"),
        (montage(cs_d, idx_o), f"CS L1-wavelet t={args.t_rel:g}·p99 (FINUFFT, DCF-weighted FISTA)"),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(16, 34))
    for ax, (img, t) in zip(axes, panels):
        ax.imshow(img, cmap="gray", vmin=0, vmax=1)
        ax.set_title(t, fontsize=14)
        ax.axis("off")
    fig.tight_layout()
    out = os.path.join(args.folder, f"montage_cs_t{args.t_rel:g}.png")
    fig.savefig(out, dpi=100)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

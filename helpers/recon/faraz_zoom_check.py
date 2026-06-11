"""Faraz zoom-bug verification + geometry-corrected comparison.

Documents the finding of 2026-06-11: his analysis script's `resizing` block
(spiral_gpt_ACR_20250915.m:445-452, flag at :16) rescales KSpaceCoor by
alpha = ((matsize/2-1)/fov) / kmax_meas = 0.8298, magnifying the image x1.205.

This script:
  1. recomputes alpha from the actual trajectory (verifies 0.8298)
  2. reruns our CG with alpha-emulated trajectory -> scale match vs his volume
     must peak at s = 1.00 (the proof the bug explains everything)
  3. builds a geometry-CORRECTED Faraz volume (zoom by alpha about center)
  4. renders a 4-row montage: faraz original / faraz corrected / our CG /
     steve-equiv  — all at true geometry except row 1
  5. measures object extents + uniformity/SNR metrics, saves to JSON

Usage: .venv/bin/python faraz_zoom_check.py <recon_io folder>
Outputs (in folder): montage_zoomfix.png, zoom_check_metrics.json
"""

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.ndimage import zoom, binary_erosion, binary_dilation

import asap_recon as ar
from steve_kernel_numpy import steve_recon
from faraz_montage import montage, orient_to_match, fraction_matched_slices, SLICES_1B

ACR_TRUE_MM = (190.0, 190.0, 148.0)
FOV_MM = 250.0


def center_zoom(vol, factor, out_shape=None):
    """Zoom about the volume center, crop/pad back to out_shape."""
    if out_shape is None:
        out_shape = vol.shape
    z = zoom(vol, factor, order=1)
    out = np.zeros(out_shape, dtype=vol.dtype)
    src, dst = [], []
    for n_z, n_o in zip(z.shape, out_shape):
        if n_z >= n_o:
            a = (n_z - n_o) // 2; src.append(slice(a, a + n_o)); dst.append(slice(0, n_o))
        else:
            a = (n_o - n_z) // 2; src.append(slice(0, n_z)); dst.append(slice(a, a + n_z))
    out[tuple(dst)] = z[tuple(src)]
    return out


def extent_mm(vol, name, fov=FOV_MM):
    n = vol.shape[0]
    thr = 0.25 * np.percentile(vol, 99)
    mask = vol > thr
    ext = []
    for ax in range(3):
        proj = mask.any(axis=tuple(i for i in range(3) if i != ax))
        idx = np.where(proj)[0]
        ext.append(float((idx[-1] - idx[0] + 1) / n * fov))
    print(f"  {name:24s} extent x/y/z = {ext[0]:.0f}/{ext[1]:.0f}/{ext[2]:.0f} mm"
          f"   (ACR true {ACR_TRUE_MM[0]:.0f}/{ACR_TRUE_MM[1]:.0f}/{ACR_TRUE_MM[2]:.0f})")
    return ext


def quality(vol, name):
    thr = 0.4 * np.percentile(vol, 99.5)
    obj = binary_erosion(vol > thr, iterations=3)
    bg = ~binary_dilation(vol > thr, iterations=5)
    cv = float(vol[obj].std() / vol[obj].mean())
    snr = float(vol[obj].mean() / vol[bg].std())
    print(f"  {name:24s} interior CV {cv:.3f}   SNR {snr:.1f}")
    return {"cv": cv, "snr": snr}


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    vol_f = np.abs(loadmat(os.path.join(folder, "faraz", "faraz_recon.mat"),
                           squeeze_me=True)["img_gp"])
    d = ar.load_steve_npy(folder)
    meta = json.load(open(os.path.join(folder, "meta.json")))

    # --- 1. alpha from the actual trajectory (his formula) ---
    MS, IS = meta["MS"], meta["IS"]
    kdk = d["traj"] - MS / 2                       # prop. to delta-k (grid units)
    kmax_meas_invmm = np.linalg.norm(kdk, axis=1).max() / (MS / IS) / FOV_MM
    alpha = ((80 / 2 - 1) / FOV_MM) / kmax_meas_invmm
    print(f"alpha (his resizing formula on this trajectory) = {alpha:.4f}  "
          f"-> magnification x{1/alpha:.3f}")

    # --- recons ---
    print("our CG recon ...")
    vol_cg = np.abs(ar.recon(d["traj"], d["acq"], method="cg", cg_iters=15))
    print("steve-equiv recon ...")
    vol_st = np.abs(steve_recon(d["traj"], d["acq"], npts=meta["npts"], MS=MS,
                                IS=IS, smoothing=meta["gplb"], axes="xyz"))

    # --- 2. proof: alpha-emulated CG must match his volume at s = 1.00 ---
    traj_a = kdk * alpha + MS / 2
    print("alpha-emulated CG recon (reproduces his geometry) ...")
    vol_cga = np.abs(ar.recon(traj_a, d["acq"], method="cg", cg_iters=15))
    vf_xyz = np.flip(np.transpose(vol_f, (0, 2, 1)), axis=1)   # his (x,z,y) -> xyz

    def corr_at_scale(target, s):
        f = center_zoom(vf_xyz, 100 / 80 * s, out_shape=(100, 100, 100))
        return float(np.corrcoef(f.ravel(), target.ravel())[0, 1])

    scales = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10]
    sweep = {f"{s:.2f}": corr_at_scale(vol_cga, s) for s in scales}
    s_best = max(sweep, key=sweep.get)
    print(f"scale sweep vs alpha-emulated CG: best at s = {s_best} "
          f"(corr {sweep[s_best]:.4f})  [1.00 = bug fully explains geometry]")

    # --- 3. corrected Faraz volume (true geometry) ---
    vol_f_fix = center_zoom(vol_f, alpha, out_shape=vol_f.shape)

    # --- 5. metrics ---
    print("\nobject extents:")
    ext = {
        "faraz_original": extent_mm(vol_f, "faraz original"),
        "faraz_corrected": extent_mm(vol_f_fix, "faraz zoom-corrected"),
        "ours_cg": extent_mm(vol_cg, "ours CG"),
        "steve_equiv": extent_mm(vol_st, "steve-equiv"),
    }
    print("\nquality (scale-free):")
    qual = {
        "faraz_original": quality(vol_f, "faraz original"),
        "faraz_corrected": quality(vol_f_fix, "faraz zoom-corrected"),
        "ours_cg": quality(vol_cg, "ours CG"),
        "steve_equiv": quality(vol_st, "steve-equiv"),
    }

    # --- 4. montage: 4 rows ---
    idx_f = [j - 1 for j in SLICES_1B]
    idx_o = fraction_matched_slices(100, 80)
    print("\norienting our volumes to his layout for display:")
    cg_d = orient_to_match(vol_cg, vol_f_fix)
    st_d = orient_to_match(vol_st, vol_f_fix)
    panels = [
        (montage(vol_f, idx_f), f"Faraz ORIGINAL (zoom bug, x{1/alpha:.2f} magnified)"),
        (montage(vol_f_fix, idx_f), "Faraz CORRECTED (true geometry)"),
        (montage(cg_d, idx_o), "Ours (FINUFFT CG)"),
        (montage(st_d, idx_o), "Steve-equiv (numpy kernel)"),
    ]
    fig, axes = plt.subplots(4, 1, figsize=(16, 45))
    for ax, (img, t) in zip(axes, panels):
        ax.imshow(img, cmap="gray")
        ax.set_title(t, fontsize=14)
        ax.axis("off")
    fig.tight_layout()
    out = os.path.join(folder, "montage_zoomfix.png")
    fig.savefig(out, dpi=100)
    print(f"\nwrote {out}")

    with open(os.path.join(folder, "zoom_check_metrics.json"), "w") as f:
        json.dump({"alpha": alpha, "magnification": 1 / alpha,
                   "scale_sweep_vs_alpha_emulated_cg": sweep,
                   "extent_mm": ext, "acr_true_mm": ACR_TRUE_MM,
                   "quality": qual}, f, indent=2)
    print(f"wrote {os.path.join(folder, 'zoom_check_metrics.json')}")


if __name__ == "__main__":
    main()

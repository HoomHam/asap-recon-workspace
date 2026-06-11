"""Arbiter experiment: FINUFFT baseline vs Steve's saved output.

Usage:
    .venv/bin/python compare_baseline.py <folder-with-npy> [savedbin0.npy]

<folder> must contain trajx/y/z.npy + acq.npy (dumped by results.dyn_recon,
results.py:258-262). If a savedbin*.npy (Steve's IS^3 output) is given or
found in the folder, prints agreement metrics and writes side-by-side slices.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import asap_recon as ar


def norm(img):
    a = np.abs(img)
    return a / a.max() if a.max() > 0 else a


def best_aligned_corr(a, b):
    """Correlation of |a| vs |b| over axis flips (FFT-direction/orientation
    differences show up as flips — comparison doc, shared property #8/11)."""
    best = (-1.0, None)
    for fx in (False, True):
        for fy in (False, True):
            for fz in (False, True):
                t = b
                for ax, f in enumerate((fx, fy, fz)):
                    if f:
                        t = np.flip(t, axis=ax)
                c = np.corrcoef(norm(a).ravel(), norm(t).ravel())[0, 1]
                if c > best[0]:
                    best = (c, (fx, fy, fz), t)
    return best


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    d = ar.load_steve_npy(folder)
    traj, acq = d["traj"], d["acq"]
    print(f"traj {traj.shape}  acq {acq.shape}  "
          f"({len(acq)/traj.shape[0]:.1f}x trajectory repeats)")
    print(f"traj grid range [{traj.min():.1f}, {traj.max():.1f}] (MS={ar.MS_DEFAULT})")

    print("\n[1/2] plain adjoint (density-biased — expect center-heavy blur) ...")
    img_adj = ar.recon(traj, acq, method="adjoint")
    print("[2/2] CG inverse (20 iters) ...")
    img_cg = ar.recon(traj, acq, method="cg", cg_iters=20)

    recons = [("adjoint", img_adj), ("CG", img_cg)]

    # Steve's output: GPU savedbin0.npy if present, else our numpy
    # reimplementation of his kernel (steve_kernel_numpy, axes-matched)
    steve_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(folder, "savedbin0.npy")
    meta_path = os.path.join(folder, "meta.json")
    steve = None
    steve_label = "steve"
    if os.path.exists(steve_path):
        steve = np.transpose(np.load(steve_path), (2, 1, 0))  # his (z,y,x) -> xyz
        print(f"\nSteve GPU output {steve.shape} from {steve_path}")
    elif os.path.exists(meta_path):
        import json
        from steve_kernel_numpy import steve_recon
        meta = json.load(open(meta_path))
        print(f"\nno savedbin0.npy — computing Steve-equivalent (numpy kernel, "
              f"gplb={meta['gplb']}) ...")
        steve = steve_recon(d["traj"], acq, npts=meta["npts"], MS=meta["MS"],
                            IS=meta["IS"], smoothing=meta["gplb"], axes="xyz")
        steve_label = "steve-equiv (numpy)"
    if steve is not None:
        for name, img in recons:
            c, flips, _ = best_aligned_corr(img, steve)
            print(f"  corr(|{name}|, |{steve_label}|) = {c:.4f}   (flips {flips})")
    else:
        print("\nno savedbin0.npy and no meta.json — showing our recons only")

    # figure: center slices, 3 orientations x recons (+ steve)
    panels = recons + ([(steve_label, steve)] if steve is not None else [])
    fig, axes = plt.subplots(3, len(panels), figsize=(4 * len(panels), 11))
    for j, (name, img) in enumerate(panels):
        a = norm(img)
        n2 = a.shape[0] // 2
        for i, (sl, lbl) in enumerate([(a[n2, :, :], "x-cut"),
                                       (a[:, n2, :], "y-cut"),
                                       (a[:, :, n2], "z-cut")]):
            ax = axes[i, j]
            ax.imshow(sl, cmap="gray")
            ax.set_title(f"{name} {lbl}", fontsize=9)
            ax.axis("off")
    out = os.path.join(folder, "baseline_comparison.png")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

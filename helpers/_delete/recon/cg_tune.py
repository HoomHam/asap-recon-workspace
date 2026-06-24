"""CG regularization sweep: lambda x iterations, jointly (early stopping is
itself a regularizer, so the two knobs trade off and must be swept together).

For each lambda (scaled relative to the estimated magnitude of diag(A^H A)),
runs CG once to max iters, snapshotting at checkpoints. Reports SNR, interior
CV, low-frequency CV, and object extent per (lambda, iters). Saves a
center-slice contact sheet for visual judgment and a JSON of all metrics.

Usage: .venv/bin/python cg_tune.py <recon_io folder>
Outputs: cg_tune_sheet.png, cg_tune_metrics.json
"""

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, binary_erosion, binary_dilation

import asap_recon as ar

LAM_RELS = [0.0, 1e-3, 1e-2, 3e-2, 1e-1]
SNAPS = [10, 15, 20, 30]
N = 100


def cg_with_snapshots(traj_rad, data, n, lam, snaps, eps=1e-7):
    """Plain CG on (A^H A + lam I) x = A^H s, yielding copies at snapshot iters."""
    def normal_op(x):
        return ar.adjoint(traj_rad, ar.forward(traj_rad, x, eps=eps), n, eps=eps) + lam * x
    b = ar.adjoint(traj_rad, data, n, eps=eps)
    x = np.zeros_like(b)
    r = b.copy()
    p = r.copy()
    rs = np.vdot(r, r).real
    out = {}
    for it in range(1, max(snaps) + 1):
        Ap = normal_op(p)
        x = x + (rs / np.vdot(p, Ap).real) * p
        r = r - (rs / np.vdot(p, Ap).real) * Ap
        rs_new = np.vdot(r, r).real
        p = r + (rs_new / rs) * p
        rs = rs_new
        if it in snaps:
            out[it] = x.copy()
    return out


def metrics(vol):
    v = np.abs(vol)
    thr = 0.4 * np.percentile(v, 99.5)
    obj = binary_erosion(v > thr, iterations=3)
    bg = ~binary_dilation(v > thr, iterations=5)
    if not obj.any() or not bg.any():  # over-regularized: object mask empty
        return {"snr": 0.0, "cv": float("nan"), "lowfreq_cv": float("nan"),
                "extent_mm": [0, 0, 0]}
    mu = v[obj].mean()
    lo = gaussian_filter(v, 2.5)
    ext = []
    mask = v > 0.25 * np.percentile(v, 99)
    for ax in range(3):
        proj = mask.any(axis=tuple(i for i in range(3) if i != ax))
        idx = np.where(proj)[0]
        ext.append(float((idx[-1] - idx[0] + 1) / v.shape[0] * 250)
                   if idx.size else 0.0)
    return {"snr": float(mu / v[bg].std()),
            "cv": float(v[obj].std() / mu),
            "lowfreq_cv": float(lo[obj].std() / mu),
            "extent_mm": [round(e) for e in ext]}


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    d = ar.load_steve_npy(folder)
    traj_rad = ar.grid_to_radians(ar.tile_traj(d["traj"].astype(float), len(d["acq"])))
    data = np.ascontiguousarray(d["acq"], dtype=np.complex128)

    # scale reference for lambda: ||A^H A x|| / ||x|| on a random image
    rng = np.random.default_rng(0)
    xt = rng.standard_normal((N, N, N)) + 1j * rng.standard_normal((N, N, N))
    scale = np.linalg.norm(ar.adjoint(traj_rad, ar.forward(traj_rad, xt), N)) / np.linalg.norm(xt)
    print(f"diag(A^H A) scale estimate: {scale:.3e}")

    results, slices = {}, {}
    for lr in LAM_RELS:
        lam = lr * scale
        print(f"lambda_rel = {lr:g} (lam = {lam:.3e}) ...")
        snaps = cg_with_snapshots(traj_rad, data, N, lam, SNAPS)
        for it, vol in snaps.items():
            key = f"lam{lr:g}_it{it}"
            results[key] = {"lam_rel": lr, "iters": it, **metrics(vol)}
            slices[key] = np.abs(vol)[:, :, N // 2]
            m = results[key]
            print(f"  it={it:2d}  SNR {m['snr']:5.1f}  CV {m['cv']:.3f}  "
                  f"lowfreqCV {m['lowfreq_cv']:.3f}  extent {m['extent_mm']}")

    # contact sheet: rows = lambda, cols = iters
    fig, axes = plt.subplots(len(LAM_RELS), len(SNAPS),
                             figsize=(4 * len(SNAPS), 4 * len(LAM_RELS)))
    for i, lr in enumerate(LAM_RELS):
        for j, it in enumerate(SNAPS):
            key = f"lam{lr:g}_it{it}"
            ax = axes[i, j]
            ax.imshow(slices[key], cmap="gray")
            m = results[key]
            ax.set_title(f"λrel={lr:g} it={it}\nSNR {m['snr']:.1f} CV {m['cv']:.3f}",
                         fontsize=9)
            ax.axis("off")
    fig.tight_layout()
    out = os.path.join(folder, "cg_tune_sheet.png")
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")

    with open(os.path.join(folder, "cg_tune_metrics.json"), "w") as f:
        json.dump({"scale": scale, "results": results}, f, indent=2)
    print(f"wrote {os.path.join(folder, 'cg_tune_metrics.json')}")

    # reference points from earlier runs
    print("\nreference: steve-equiv SNR 28.7 CV 0.154 | faraz SNR 22.8 CV 0.139 | "
          "old CG(15, lam=0) SNR 19.6 CV 0.160")


if __name__ == "__main__":
    main()

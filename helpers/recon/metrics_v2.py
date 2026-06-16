"""metrics_v2 — fixed-ROI recon quality metrics (supersedes cg_tune.metrics).

Why a NEW module (handoff 2026-06-15, roadmap step 1): the original
`cg_tune.metrics` (cg_tune.py:52) has two *measured* blind spots —

  * SNR inflated when sparsity priors zero the background. The old noise
    estimate uses an ADAPTIVE background mask (~dilate(object > 0.4*p99)). When
    a CS prior zeroes the air, sigma_bg collapses toward 0 and SNR explodes:
    wavelet t0.03/t0.1 scored SNR ~68 on visibly WORSE images.
  * extent_mm is threshold-fragile. The old mask is `v > 0.25*p99`; residual
    noise clears that threshold, so CG-20 reads 250/250/250 mm (the full FOV)
    when the ACR phantom truth is ~190/190/148 mm.

Fixes here:
  1. Noise from a FIXED set of 8 corner air ROIs, not an adaptive mask. Fixed
     location -> comparable across recons and not gameable by mask shift. When
     a prior truly zeroes the corners we *flag* it (`bg_collapsed=True`) and
     report snr=inf, instead of silently emitting an inflated finite number.
  2. extent from a HALF-MAXIMUM object threshold (0.5 * p99.5 peak). Half-max
     sits far above the noise floor, so noise can't pad the bounding box.
  3. `edge_sharp`: boundary gradient energy normalized by object mean. Blurring
     an image to win SNR lowers edge_sharp -> the guard the handoff asked for
     ("so SNR can't be won by blurring").

Returns the original keys (snr, cv, lowfreq_cv, extent_mm) plus edge_sharp,
noise, signal, bg_collapsed. Drop-in: `from metrics_v2 import metrics`.

Assumes the object is roughly centered (ACR phantom) so the 8 grid corners are
air. True for v3_fov250. FOV default 250 mm = IS(100) voxels * 2.5 mm.

DO NOT port this back into cg_tune.metrics: the old sweep JSONs
(cs_sweep_metrics.json, cg_tune_metrics.json) were computed with the old metric
and must stay reproducible. New code imports from here.

Usage (demo the fix on CG-20): .venv/bin/python metrics_v2.py <recon_io folder>
"""

import numpy as np
from scipy.ndimage import (gaussian_filter, sobel,
                           binary_dilation, binary_erosion)

FOV_MM = 250.0  # IS=100 spans one FOV = 250 mm (v3_fov250)


def corner_rois(shape, k=12):
    """Boolean mask of the 8 corner cubes (side k) — fixed air ROIs."""
    nx, ny, nz = shape
    m = np.zeros(shape, bool)
    for sx in (slice(0, k), slice(nx - k, nx)):
        for sy in (slice(0, k), slice(ny - k, ny)):
            for sz in (slice(0, k), slice(nz - k, nz)):
                m[sx, sy, sz] = True
    return m


def object_mask(v, frac=0.5):
    """Half-maximum object mask. peak = p99.5 (robust to single hot voxels)."""
    peak = np.percentile(v, 99.5)
    return v > frac * peak


def metrics(vol, fov_mm=FOV_MM, corner=12):
    """Quality metrics for a recon volume. See module docstring for rationale.

    Keys: snr, cv, lowfreq_cv, extent_mm[3], edge_sharp, noise, signal,
          bg_collapsed.
    """
    v = np.abs(vol)
    obj = object_mask(v)
    if not obj.any():  # nothing above half-max -> emptied by over-regularization
        return {"snr": 0.0, "cv": float("nan"), "lowfreq_cv": float("nan"),
                "extent_mm": [0, 0, 0], "edge_sharp": 0.0, "noise": 0.0,
                "signal": 0.0, "bg_collapsed": False}

    signal = float(v[obj].mean())
    bg = corner_rois(v.shape, corner)
    noise = float(v[bg].std())
    collapsed = bool(noise < 1e-6 * signal)
    snr = float("inf") if collapsed else float(signal / noise)

    lo = gaussian_filter(v, 2.5)

    # extent from the half-max bounding box (robust to noise floor)
    ext = []
    for ax in range(3):
        proj = obj.any(axis=tuple(i for i in range(3) if i != ax))
        idx = np.where(proj)[0]
        ext.append(float((idx[-1] - idx[0] + 1) / v.shape[ax] * fov_mm)
                   if idx.size else 0.0)

    # edge sharpness: mean |grad| on a 1-voxel boundary shell, normalized by
    # object mean. Blurring an image lowers this -> SNR can't be bought by blur.
    gm = np.sqrt(sum(sobel(v, axis=a) ** 2 for a in range(3)))
    shell = binary_dilation(obj, iterations=1) & ~binary_erosion(obj, iterations=1)
    edge = float(gm[shell].mean() / signal) if shell.any() else 0.0

    return {"snr": snr,
            "cv": float(v[obj].std() / signal),
            "lowfreq_cv": float(lo[obj].std() / signal),
            "extent_mm": [round(e) for e in ext],
            "edge_sharp": round(edge, 4),
            "noise": noise, "signal": signal,
            "bg_collapsed": collapsed}


def _demo():
    """Recompute CG-20 metrics with old vs new metric to show the fix."""
    import sys
    import asap_recon as ar
    from cg_tune import metrics as metrics_old

    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    d = ar.load_steve_npy(folder)
    print("reconstructing CG-20 (method of record) ...")
    img = ar.recon(d["traj"], np.ascontiguousarray(d["acq"], dtype=complex),
                   method="cg", cg_iters=20)
    old = metrics_old(img)
    new = metrics(img)
    print("\n            OLD (cg_tune)        NEW (metrics_v2)")
    print(f"  SNR        {old['snr']:8.1f}            {new['snr']:8.1f}")
    print(f"  extent_mm  {str(old['extent_mm']):>16}   {str(new['extent_mm']):>16}"
          "   (ACR truth ~190/190/148)")
    print(f"  lowfreq_cv {old['lowfreq_cv']:8.3f}            {new['lowfreq_cv']:8.3f}")
    print(f"  edge_sharp     {'(n/a)':>8}            {new['edge_sharp']:8.4f}")
    print(f"  bg_collapsed   {'(n/a)':>8}            {str(new['bg_collapsed']):>8}")


if __name__ == "__main__":
    _demo()

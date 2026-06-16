"""BART `pics` CS vs our CS — the independence test (handoff 2026-06-15 step 2).

Our CS layer (cs_recon.py) is built on sigpy, which is Lustig-lab (Frank Ong);
the MATLAB Lustig comparison shares that lineage, so it is a personal baseline,
not an independence test. BART (Uecker/Lustig, C) is the independent benchmark:
different codebase, native 3D non-Cartesian `pics`, its own ell1-wavelet/TV.

This is NEW code. It imports our helpers read-only and writes only into the
recon_io folder's `bart/` subdir. It does NOT touch cs_recon.py, cg_tune.py, or
asap_recon.py.

What it does
------------
  1. Build BART inputs from the same recon_io arrays:
       traj  -> bart_traj = traj_rad * N/(2*pi), dims (3, M, 1)  [BART units:
                 the N-point grid spans [-N/2, N/2]; ours is grid-index, the
                 conversion is the radian traj scaled by N/2pi]
       ksp   -> dims (1, M, 1)
       sens  -> ones(N,N,N,1)  (single coil)
  2. Sweep `bart pics` over regularizer x lambda:
       wavelet : -R W:7:0:lam     TV : -R T:7:0:lam     (flags 7 = x,y,z)
  3. Reconstruct our CG-20 + our wavelet finalists (t=0.003, 0.01) with the
     SAME primitives as cs_recon.
  4. Score every volume with metrics_v2 (the fixed metrics), orient BART output
     to our frame, write a montage + JSON into recon_io/bart/.

Usage:
  ../.venv/bin/python bart_compare.py ../../data/v3_fov250/recon_io \
      --bart /Users/hoomham/bin/bart-src/bart --iters 100
"""

import argparse
import json
import os
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sigpy as sp

import asap_recon as ar
import cs_recon as cs
from bartio import readcfl, writecfl
from metrics_v2 import metrics

N = cs.N  # 100
LAMBDAS = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]   # BART pics regularization strength
OUR_T = [0.003, 0.01]                       # our wavelet finalists (prior sweep)


# --------------------------------------------------------------- our pipeline

def our_setup(folder):
    """Mirror cs_recon.main(): normalized DCF-weighted operator A, y_n, and the
    CG-20 image used as both the method-of-record anchor and the threshold ref."""
    d = ar.load_steve_npy(folder)
    traj_rad = ar.grid_to_radians(
        ar.tile_traj(d["traj"].astype(float), len(d["acq"])))
    y = np.ascontiguousarray(d["acq"], dtype=np.complex128)

    A_raw = cs.FinufftForward(traj_rad, N)
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
    return {"traj_rad": traj_rad, "y": y, "A": A, "y_n": y_n,
            "img_cg": img_cg, "t_ref": t_ref}


# --------------------------------------------------------------- BART pipeline

def write_bart_inputs(out_dir, traj_rad, y):
    """traj_rad (M,3) -> BART traj (3,M,1); y (M,) -> ksp (1,M,1); ones sens."""
    bart_traj = (traj_rad.T * (N / (2 * np.pi))).astype(np.complex64)  # (3,M)
    writecfl(os.path.join(out_dir, "traj"), bart_traj[:, :, None])     # (3,M,1)
    writecfl(os.path.join(out_dir, "ksp"),
             y.astype(np.complex64)[None, :, None])                    # (1,M,1)
    writecfl(os.path.join(out_dir, "sens"),
             np.ones((N, N, N, 1), dtype=np.complex64))


def run_pics(bart, out_dir, reg, lam, iters):
    """`bart pics -t traj -R <reg>:7:0:lam ksp sens out`. reg in {W,T}.

    Wavelet (W) is forced onto ADMM (-m): pics' default FISTA path diverged to
    ~1e32 here (unnormalized ones-sensitivity -> bad Lipschitz step), while ADMM
    is stable and matches the TV path's solver. Measured 2026-06-15."""
    out = os.path.join(out_dir, f"bart_{reg}_l{lam:g}")
    cmd = [bart, "pics", "-i", str(iters)]
    if reg == "W":
        cmd.append("-m")  # ADMM — FISTA diverges on the wavelet prox here
    cmd += ["-t", os.path.join(out_dir, "traj"),
            "-R", f"{reg}:7:0:{lam}",
            os.path.join(out_dir, "ksp"),
            os.path.join(out_dir, "sens"), out]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    img = readcfl(out)
    return np.abs(np.squeeze(img))


def orient_to_ours(vol_bart, ref):
    """BART and finufft differ by axis/flip/conjugate conventions. Pick the
    flip/transpose of |vol_bart| that best matches |ref| (our CG-20) by
    normalized correlation. Magnitude-only, so conjugation drops out."""
    a = np.abs(ref)
    a = (a - a.mean()) / (a.std() + 1e-12)
    best, best_c = vol_bart, -np.inf
    cands = []
    for flips in [(), (0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]:
        v = vol_bart
        for ax in flips:
            v = np.flip(v, axis=ax)
        cands.append(v)
        cands.append(np.transpose(v, (1, 0, 2)))
    for v in cands:
        b = (v - v.mean()) / (v.std() + 1e-12)
        cc = float(np.mean(a * b))
        if cc > best_c:
            best_c, best = cc, v
    return best, best_c


# ---------------------------------------------------------------------- main

def norm01(vol):
    p = np.percentile(np.abs(vol), 99.5)
    return np.clip(np.abs(vol) / (p + 1e-12), 0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("--bart", default="/Users/hoomham/bin/bart-src/bart")
    ap.add_argument("--iters", type=int, default=100)
    ap.add_argument("--lambdas", type=float, nargs="+", default=LAMBDAS)
    args = ap.parse_args()

    out_dir = os.path.join(args.recon_io, "bart")
    os.makedirs(out_dir, exist_ok=True)

    print("== our pipeline setup (DCF-weighted normalized operator) ==")
    o = our_setup(args.recon_io)
    results, vols = {}, {}

    results["ours_cg20"] = metrics(o["img_cg"])
    vols["ours_cg20"] = o["img_cg"]
    print(f"  ours_cg20: {results['ours_cg20']}")

    for t in OUR_T:
        img = cs.wavelet_recon(o["A"], o["y_n"], t * o["t_ref"], args.iters)
        results[f"ours_wav_t{t:g}"] = {"t_rel": t, **metrics(img)}
        vols[f"ours_wav_t{t:g}"] = img
        print(f"  ours_wav_t{t:g}: {results[f'ours_wav_t{t:g}']}")

    print("== BART inputs ==")
    write_bart_inputs(out_dir, o["traj_rad"], o["y"])

    print("== BART pics sweep ==")
    for reg, name in [("W", "wav"), ("T", "tv")]:
        for lam in args.lambdas:
            key = f"bart_{name}_l{lam:g}"
            try:
                raw = run_pics(args.bart, out_dir, reg, lam, args.iters)
            except subprocess.CalledProcessError as e:
                print(f"  {key}: FAILED\n{e.stderr[-400:]}")
                continue
            v, cc = orient_to_ours(raw, o["img_cg"])
            results[key] = {"reg": name, "lam": lam, "orient_cc": round(cc, 3),
                            **metrics(v)}
            vols[key] = v
            m = results[key]
            print(f"  {key}: SNR {m['snr']:6.1f} lfCV {m['lowfreq_cv']:.3f} "
                  f"edge {m['edge_sharp']:.2f} ext {m['extent_mm']} "
                  f"cc {cc:.2f} bgcol {m['bg_collapsed']}")

    # montage: central axis-2 slice of each volume, shared [0,1] axis
    keys = list(vols.keys())
    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow))
    for ax in axes.ravel():
        ax.axis("off")
    for i, k in enumerate(keys):
        ax = axes.ravel()[i]
        ax.imshow(norm01(vols[k])[:, :, N // 2], cmap="gray", vmin=0, vmax=1)
        m = results[k]
        snr = m["snr"]
        snr_s = "inf" if snr == float("inf") else f"{snr:.1f}"
        ax.set_title(f"{k}\nSNR {snr_s} lfCV {m['lowfreq_cv']:.3f} "
                     f"edge {m.get('edge_sharp', 0):.1f}", fontsize=8)
    fig.tight_layout()
    png = os.path.join(out_dir, "bart_vs_ours_montage.png")
    fig.savefig(png, dpi=110)
    print(f"wrote {png}")

    with open(os.path.join(out_dir, "bart_compare_metrics.json"), "w") as f:
        json.dump({"iters": args.iters, "lambdas": args.lambdas,
                   "metric": "metrics_v2", "results": results}, f, indent=2)
    print(f"wrote {os.path.join(out_dir, 'bart_compare_metrics.json')}")
    print("\nNOTE metrics_v2: SNR uses fixed corner-air ROIs (bg_collapsed flags "
          "a prior-zeroed background -> SNR=inf is not a quality win); extent is "
          "half-max; edge_sharp guards against winning SNR by blurring.")


if __name__ == "__main__":
    main()

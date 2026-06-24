"""Compressed-sensing recon on the validated FINUFFT operator (handoff CS step 1-2).

min_x ||W^(1/2) (A x - y)||^2 + lam R(x)

A      : asap_recon.forward/adjoint wrapped as sigpy Linops (keeps the
         Steve-units, isign=-1 convention in one place — handoff preference
         over sigpy's own NUFFT).
W      : DCF preconditioning weights, w = 1/density (density via |A A^H 1|).
         Required for gradient-type solvers: unweighted A^H A has maxeig
         dominated by the oversampled spiral center (measured density spread
         1.8e6) so FISTA/PDHG amplitudes stall ~1000x low while CG cuts
         through. Weighted FISTA reaches the CG solution scale in ~20 iters
         (corr 0.945 vs unweighted CG-20). Standard CS-MRI practice (sigpy
         mri apps use dcf the same way). Bin weights later multiply into the
         same W (architecture D3).
R      : 'wavelet' — L1 of db4 coefficients, synthesis form, FISTA
         'tv'      — L1 of finite differences, PDHG
lam    : parameterized by the FISTA soft-threshold t = alpha*lam (alpha =
         1/maxeig(A^H A)), with t swept relative to the p99 coefficient
         magnitude of the CG-20 solution per regularizer. Scaling to
         max|A^H y| (first attempt 2026-06-12) put the threshold 4-5 orders
         below the coefficients — a measured no-op. Choice gated by Hooman's
         visual verdict on the contact sheet, not the scalar metrics
         (standing rule).

Static single bin: sample_weights = ones (omitted). Binned recon later =
same call with W_b folded into A and y (architecture D3 unchanged).

Usage: .venv/bin/python cs_recon.py <recon_io folder> [--max-iter N]
Outputs in <folder>: cs_sweep_sheet.png, cs_sweep_metrics.json
"""

import argparse
import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sigpy as sp

import asap_recon as ar
from cg_tune import metrics

N = 100
T_RELS = [0.003, 0.01, 0.03, 0.1]  # threshold / p99(|coeffs of x_cg|)
# (first grid [0.1..3] over-regularized: t=0.1 already shrank the object,
#  t=1 emptied it. t=0.1 kept as the known-too-strong anchor.)
EPS = 1e-7  # finufft tolerance, same as cg_recon


# ------------------------------------------------------- finufft as sigpy Linop

class FinufftForward(sp.linop.Linop):
    """A: (N,N,N) image -> (M,) nonuniform samples, times `scale`.

    scale = 1/sqrt(maxeig(A^H A)) normalizes the operator so
    maxeig(A_n^H A_n) = 1: FISTA's soft-threshold equals lam directly, and
    PDHG's primal/dual steps are sane (raw ||A|| ~ 8e4 left tau ~ 1e-10 —
    measured: TV output was byte-identical across 4 decades of lam).
    Normalizing y by the same factor keeps the least-squares minimizer
    unchanged.
    """

    def __init__(self, traj_rad, n, eps=EPS, scale=1.0):
        self.traj_rad = traj_rad
        self.eps = eps
        self.scale = scale
        super().__init__((traj_rad.shape[0],), (n, n, n))

    def _apply(self, input):
        return self.scale * ar.forward(self.traj_rad, input, eps=self.eps)

    def _adjoint_linop(self):
        return FinufftAdjoint(self.traj_rad, self.ishape[0], eps=self.eps,
                              scale=self.scale)


class FinufftAdjoint(sp.linop.Linop):
    """A^H: (M,) nonuniform samples -> (N,N,N) image, times `scale`."""

    def __init__(self, traj_rad, n, eps=EPS, scale=1.0):
        self.traj_rad = traj_rad
        self.eps = eps
        self.scale = scale
        super().__init__((n, n, n), (traj_rad.shape[0],))

    def _apply(self, input):
        return self.scale * ar.adjoint(self.traj_rad, input, self.oshape[0],
                                       eps=self.eps)

    def _adjoint_linop(self):
        return FinufftForward(self.traj_rad, self.oshape[0], eps=self.eps,
                              scale=self.scale)


# ------------------------------------------------------------------ CS recons

def wavelet_recon(A, y, lam, max_iter):
    """Synthesis-form L1-wavelet: min_a ||A W^H a - y||^2 + lam ||a||_1.
    Expects normalized A (maxeig = 1, so alpha = 1 and threshold = lam).
    W is orthogonal (db4): maxeig(W A^H A W^H) = maxeig(A^H A) = 1."""
    W = sp.linop.Wavelet(A.ishape, wave_name="db4")
    a = sp.app.LinearLeastSquares(A * W.H, y,
                                  proxg=sp.prox.L1Reg(W.oshape, lam),
                                  alpha=1.0,
                                  max_iter=max_iter, show_pbar=False).run()
    return W.H(a)


def tv_recon(A, y, lam, max_iter):
    """min_x ||A x - y||^2 + lam ||grad x||_1 via PDHG. Expects normalized A."""
    G = sp.linop.FiniteDifference(A.ishape)
    return sp.app.LinearLeastSquares(A, y,
                                     proxg=sp.prox.L1Reg(G.oshape, lam), G=G,
                                     max_iter=max_iter, show_pbar=False).run()


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", nargs="?", default=".")
    ap.add_argument("--max-iter", type=int, default=80)
    args = ap.parse_args()

    d = ar.load_steve_npy(args.folder)
    traj_rad = ar.grid_to_radians(
        ar.tile_traj(d["traj"].astype(float), len(d["acq"])))
    y = np.ascontiguousarray(d["acq"], dtype=np.complex128)

    A_raw = FinufftForward(traj_rad, N)
    M = len(y)

    # DCF preconditioning weights (see module docstring)
    dens = np.abs(A_raw(A_raw.H(np.ones(M, dtype=complex))))
    w = 1.0 / np.clip(dens, dens.max() * 1e-4, None)
    w /= w.mean()
    Aw_raw = sp.linop.Multiply((M,), np.sqrt(w)) * A_raw
    L = sp.app.MaxEig(Aw_raw.H * Aw_raw, dtype=np.complex128,
                      show_pbar=False).run()
    print(f"M = {M}   density spread {dens.max()/dens.min():.2e}   "
          f"maxeig(A_w^H A_w) = {L:.4e}")
    c = 1.0 / np.sqrt(L)
    A = sp.linop.Multiply((M,), c * np.sqrt(w)) * A_raw  # normalized: maxeig=1
    y_n = c * np.sqrt(w) * y

    print("baseline CG-20 (method of record) ...")
    img_cg = ar.recon(d["traj"], y, method="cg", cg_iters=20)

    results = {"cg20": metrics(img_cg)}
    slices = {"cg20": np.abs(img_cg)[:, :, N // 2]}
    print(f"  cg20: {results['cg20']}")

    # per-regularizer threshold reference: p99 coefficient magnitude of x_cg
    W = sp.linop.Wavelet(A.ishape, wave_name="db4")
    G = sp.linop.FiniteDifference(A.ishape)
    t_ref = {"wavelet": float(np.percentile(np.abs(W(img_cg)), 99)),
             "tv": float(np.percentile(np.abs(G(img_cg)), 99))}
    print(f"threshold refs: {t_ref}")

    for reg, fn in [("wavelet", wavelet_recon), ("tv", tv_recon)]:
        for tr in T_RELS:
            key = f"{reg}_t{tr:g}"
            lam = tr * t_ref[reg]  # normalized A: FISTA threshold = lam
            t0 = time.time()
            img = fn(A, y_n, lam, args.max_iter)
            m = metrics(img)
            results[key] = {"t_rel": tr, "lam": lam,
                            "seconds": round(time.time() - t0, 1), **m}
            slices[key] = np.abs(img)[:, :, N // 2]
            print(f"  {key}: SNR {m['snr']:5.1f}  CV {m['cv']:.3f}  "
                  f"lowfreqCV {m['lowfreq_cv']:.3f}  extent {m['extent_mm']}  "
                  f"({results[key]['seconds']}s)")

    # contact sheet: row 0 = CG baseline + reference bar, rows 1-2 = wavelet/TV
    ncols = len(T_RELS)
    fig, axes = plt.subplots(3, ncols, figsize=(4 * ncols, 12))
    for ax in axes.ravel():
        ax.axis("off")
    m = results["cg20"]
    axes[0, 0].imshow(slices["cg20"], cmap="gray")
    axes[0, 0].set_title(f"CG-20 baseline\nSNR {m['snr']:.1f} "
                         f"lfCV {m['lowfreq_cv']:.3f}", fontsize=9)
    axes[0, 1].text(0, 0.5,
                    "Bar (handoff):\n"
                    "SNR > 28.7 (Steve-equiv)\n"
                    "lowfreq-CV < 0.093 (Faraz)\n"
                    "extent 190/190/148 mm (ACR)\n"
                    "edges: judge visually",
                    fontsize=11, family="monospace", va="center")
    for i, reg in enumerate(["wavelet", "tv"], start=1):
        for j, tr in enumerate(T_RELS):
            key = f"{reg}_t{tr:g}"
            m = results[key]
            axes[i, j].imshow(slices[key], cmap="gray")
            axes[i, j].set_title(f"{reg} t={tr:g}·p99\nSNR {m['snr']:.1f} "
                                 f"lfCV {m['lowfreq_cv']:.3f}", fontsize=9)
    fig.tight_layout()
    out = os.path.join(args.folder, "cs_sweep_sheet.png")
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")

    with open(os.path.join(args.folder, "cs_sweep_metrics.json"), "w") as f:
        json.dump({"maxeig": L, "t_ref": t_ref, "max_iter": args.max_iter,
                   "results": results}, f, indent=2)
    print(f"wrote {os.path.join(args.folder, 'cs_sweep_metrics.json')}")

    print("\nreference: steve-equiv SNR 28.7 lowfreqCV 0.110 | "
          "faraz lowfreqCV 0.093 | CG-20 prior runs SNR ~19.6")


if __name__ == "__main__":
    main()

"""Wavelet two-way (ours vs BART) + the DCF-grid proof. handoff step-2 close.

NEW assembler. Two figures into recon_io/:

  wavelet_twoway_montage.png / .json
    ours db4-wavelet (FISTA, DCF-weighted) vs BART pics -R W (ADMM), several
    lambda, all scored metrics_v2, BART oriented to ours, shared [0,1] axis.
    Lustig is absent on purpose — its Wavelet operator is FWT2_PO (2D), unusable
    on the 3D volume (BART_Comparison.md).

  texture_compare.png
    THREE panels: ours+DCF | ours WITHOUT DCF | BART pics no-reg. BY EYE
    (the reliable instrument here — standing rule): BART shows a REGULAR
    cross-hatch grid; ours shows only irregular mottle, DCF or not. So the grid
    is real and BART-specific, and NOT explained by DCF (ours-no-DCF is still
    grid-free). Cause = BART nufft gridding internals, UNRESOLVED (ruled out:
    traj over-range, Toeplitz, kernel width/oversampling, fftmod, DCF). The
    grid-HF / peakiness scalars printed below CANNOT separate grid from mottle
    here (phantom spectrum overlaps the band) — do not trust them for this; they
    fooled the read twice. See [[eye_vs_metric]].

    ../.venv/bin/python wavelet_twoway.py ../../data/v3_fov250/recon_io \
        --bart /Users/hoomham/bin/bart-src/bart
"""

import argparse
import glob
import json
import os
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import asap_recon as ar
import cs_recon as cs
from bartio import readcfl
from bart_compare import our_setup, orient_to_ours, norm01
from metrics_v2 import metrics

N = cs.N


def grid_hf(vol):
    """Fraction of central-slice energy in the spiral-grid band (r 6-22)."""
    v = np.abs(np.squeeze(vol)); sl = v[:, :, N // 2]; sl = sl / sl.max()
    F = np.abs(np.fft.fftshift(np.fft.fft2(sl - sl.mean()))) ** 2
    yy, xx = np.mgrid[0:N, 0:N]; r = np.hypot(yy - N // 2, xx - N // 2)
    return float(F[(r > 6) & (r < 22)].sum() / F[r < N // 2 - 2].sum())


def dcf_proof(recon_io, bart, o):
    """ours+DCF | ours no-DCF | BART no-reg — the grid = unweighted-LS proof."""
    d = ar.load_steve_npy(recon_io)
    y = np.ascontiguousarray(d["acq"], complex)

    # ours WITH DCF (the wavelet champion's operator -> just CG on weighted op)
    A = o["A"]; ours_dcf = np.abs(cs.tv_recon(A, o["y_n"], 0.0, 40))  # lam0 = weighted CG
    # ours WITHOUT DCF (plain CG, method of record)
    ours_nodcf = np.abs(ar.recon(d["traj"], y, method="cg", cg_iters=30))
    # BART pics, no regularization (pure CG, ones sens, no DCF)
    out = os.path.join(recon_io, "bart", "pics_noreg")
    subprocess.run([bart, "pics", "-i", "30",
                    "-t", os.path.join(recon_io, "bart", "traj"),
                    os.path.join(recon_io, "bart", "ksp"),
                    os.path.join(recon_io, "bart", "sens"), out],
                   check=True, capture_output=True, text=True)
    bart_nr, _ = orient_to_ours(np.abs(np.squeeze(readcfl(out))), o["img_cg"])

    def peakiness(v):
        vv = np.abs(np.squeeze(v)); sl = vv[:, :, N // 2]; sl = sl / sl.max()
        F = np.abs(np.fft.fftshift(np.fft.fft2(sl - sl.mean())))
        yy, xx = np.mgrid[0:N, 0:N]; r = np.hypot(yy - N // 2, xx - N // 2)
        b = F[(r > 6) & (r < 22)]
        return b.max() / (np.median(b) + 1e-9)

    panels = [("ours + DCF\n(weighted LS, D6)", ours_dcf),
              ("ours NO DCF\n(unweighted CG)", ours_nodcf),
              ("BART pics no-reg\n(unweighted)", bart_nr)]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.6))
    for ax, (t, v) in zip(axes, panels):
        sl = np.abs(v)[:, :, N // 2]
        ax.imshow(sl / np.percentile(sl, 99.5), cmap="gray", vmin=0, vmax=1)
        ax.set_title(f"{t}\ngrid-HF {grid_hf(v):.3f}  peakiness {peakiness(v):.1f}",
                     fontsize=11)
        ax.axis("off")
    fig.suptitle("BY EYE: BART (right) has a REGULAR grid; ours (DCF or not) only "
                 "irregular mottle -> grid is real & BART-specific, NOT from DCF. "
                 "HF/peakiness below can't separate grid from mottle - don't trust "
                 "them here.", fontsize=10)
    fig.tight_layout()
    png = os.path.join(recon_io, "texture_compare.png")
    fig.savefig(png, dpi=120)
    print(f"wrote {png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("--bart", default="/Users/hoomham/bin/bart-src/bart")
    ap.add_argument("--t-rel", type=float, nargs="+", default=[0.003, 0.01, 0.03])
    ap.add_argument("--iters", type=int, default=100)
    a = ap.parse_args()

    o = our_setup(a.recon_io)
    results, vols = {}, {}

    # ours wavelet (db4 FISTA, DCF-weighted)
    import sigpy as sp
    W = sp.linop.Wavelet(o["A"].ishape, wave_name="db4")
    t_ref = float(np.percentile(np.abs(W(o["img_cg"])), 99))
    print("== ours wavelet ==")
    for t in a.t_rel:
        img = cs.wavelet_recon(o["A"], o["y_n"], t * t_ref, a.iters)
        k = f"ours_wav_t{t:g}"
        results[k] = {"t_rel": t, "grid_hf": round(grid_hf(img), 3), **metrics(img)}
        vols[k] = np.abs(img)
        print(f"  {k}: lfCV {results[k]['lowfreq_cv']:.3f} edge {results[k]['edge_sharp']:.2f} "
              f"gridHF {results[k]['grid_hf']:.3f}")

    # BART wavelet (ADMM)
    print("== BART wavelet ==")
    for f in sorted(glob.glob(os.path.join(a.recon_io, "bart", "bart_W_l*.cfl"))):
        raw = np.abs(np.squeeze(readcfl(f[:-4])))
        v, cc = orient_to_ours(raw, o["img_cg"])
        k = os.path.basename(f)[:-4].replace("bart_W", "bart_wav")
        results[k] = {"orient_cc": round(cc, 3), "grid_hf": round(grid_hf(v), 3),
                      **metrics(v)}
        vols[k] = v
        print(f"  {k}: lfCV {results[k]['lowfreq_cv']:.3f} edge {results[k]['edge_sharp']:.2f} "
              f"gridHF {results[k]['grid_hf']:.3f} cc {cc:.2f}")

    # montage
    keys = list(vols.keys())
    ncol = 3
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow))
    for ax in np.atleast_1d(axes).ravel():
        ax.axis("off")
    for i, k in enumerate(keys):
        ax = np.atleast_1d(axes).ravel()[i]
        ax.imshow(norm01(vols[k])[:, :, N // 2], cmap="gray", vmin=0, vmax=1)
        m = results[k]
        ax.set_title(f"{k}\nlfCV {m['lowfreq_cv']:.3f} edge {m['edge_sharp']:.1f} "
                     f"gridHF {m['grid_hf']:.2f}", fontsize=8)
    fig.tight_layout()
    png = os.path.join(a.recon_io, "wavelet_twoway_montage.png")
    fig.savefig(png, dpi=110)
    print(f"wrote {png}")
    with open(os.path.join(a.recon_io, "wavelet_twoway_metrics.json"), "w") as f:
        json.dump({"metric": "metrics_v2", "results": results}, f, indent=2)
    print(f"wrote {os.path.join(a.recon_io, 'wavelet_twoway_metrics.json')}")

    print("== DCF grid proof ==")
    dcf_proof(a.recon_io, a.bart, o)


if __name__ == "__main__":
    main()

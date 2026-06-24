"""TV three-way: ours vs BART vs Lustig, ALL scored with metrics_v2, same data.

NEW assembler (handoff step 2 close). TV is the one regularizer all three
pipelines do cleanly (Lustig's wavelet operator is 2D-only, so wavelet stays
ours-vs-BART — see bart_compare.py / BART_Comparison.md). This puts the three TV
recons on ONE ruler (metrics_v2) and one montage.

Sources (run the producers first):
  ours   : recomputed here — cs.tv_recon at a few t_rel (PDHG on our finufft op)
  BART   : recon_io/bart/bart_T_l*.{cfl,hdr}        (from bart_compare.py)
  Lustig : recon_io/lustig/lustig_tv_sweep.mat       (from run_lustig_sweep.py)

    ../.venv/bin/python tv_threeway.py ../../data/v3_fov250/recon_io \
        --t-rel 0.003 0.01 0.03

Outputs in recon_io/: tv_threeway_montage.png, tv_threeway_metrics.json
"""

import argparse
import glob
import json
import os

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import itertools

import cs_recon as cs
from bartio import readcfl
from bart_compare import our_setup, orient_to_ours, norm01
from metrics_v2 import metrics

N = cs.N


def best_orient_full(vol, ref):
    """Full 48-orientation search (6 axis perms x 8 flips) — needed for Lustig,
    whose per-slice rot90 mixes axes beyond orient_to_ours' flip+single-transpose
    set. Magnitude correlation, so conjugation is irrelevant."""
    a = np.abs(ref); a = (a - a.mean()) / (a.std() + 1e-12)
    best, best_c = vol, -np.inf
    for perm in itertools.permutations(range(3)):
        vp = np.transpose(np.abs(vol), perm)
        for fl in itertools.product([False, True], repeat=3):
            v = vp
            for ax, do in enumerate(fl):
                if do:
                    v = np.flip(v, axis=ax)
            b = (v - v.mean()) / (v.std() + 1e-12)
            cc = float(np.mean(a * b))
            if cc > best_c:
                best_c, best = cc, v
    return best, best_c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("--t-rel", type=float, nargs="+", default=[0.003, 0.01, 0.03])
    ap.add_argument("--iters", type=int, default=100)
    a = ap.parse_args()

    results, vols = {}, {}

    # ours TV (PDHG on our operator), scored with metrics_v2
    print("== ours TV ==")
    o = our_setup(a.recon_io)
    G = __import__("sigpy").linop.FiniteDifference(o["A"].ishape)
    t_ref = float(np.percentile(np.abs(G(o["img_cg"])), 99))
    for t in a.t_rel:
        img = cs.tv_recon(o["A"], o["y_n"], t * t_ref, a.iters)
        k = f"ours_tv_t{t:g}"
        results[k] = {"t_rel": t, **metrics(img)}
        vols[k] = np.abs(img)
        print(f"  {k}: {results[k]['snr']:.1f} lfCV {results[k]['lowfreq_cv']:.3f} "
              f"edge {results[k]['edge_sharp']:.2f}")

    # BART TV (already reconstructed; orient to ours, score with metrics_v2)
    print("== BART TV ==")
    for f in sorted(glob.glob(os.path.join(a.recon_io, "bart", "bart_T_l*.cfl"))):
        tag = os.path.basename(f)[:-4]
        raw = np.abs(np.squeeze(readcfl(f[:-4])))
        v, cc = orient_to_ours(raw, o["img_cg"])
        k = tag.replace("bart_T", "bart_tv")
        results[k] = {"orient_cc": round(cc, 3), **metrics(v)}
        vols[k] = v
        print(f"  {k}: {results[k]['snr']:.1f} lfCV {results[k]['lowfreq_cv']:.3f} "
              f"edge {results[k]['edge_sharp']:.2f} cc {cc:.2f}")

    # Lustig TV sweep
    print("== Lustig TV ==")
    mat = os.path.join(a.recon_io, "lustig", "lustig_tv_sweep.mat")
    if os.path.exists(mat):
        g = sio.loadmat(mat)
        gf = np.asarray(g["gas_final"], dtype=float)
        tvW = np.asarray(g["tvW"], dtype=float).ravel()
        for i, lam in enumerate(tvW):
            k = f"lustig_tv_l{lam:g}"
            # Lustig saved with its own per-slice rot90 -> orient to our frame,
            # else [:, :, N//2] is a different plane than ours/BART. metrics are
            # orientation-invariant (global), but the montage plane must match.
            v, cc = best_orient_full(gf[i], o["img_cg"])
            results[k] = {"TVWeight": float(lam), "orient_cc": round(cc, 3),
                          **metrics(v)}
            vols[k] = v
            print(f"  {k}: {results[k]['snr']:.1f} lfCV {results[k]['lowfreq_cv']:.3f} "
                  f"edge {results[k]['edge_sharp']:.2f} cc {cc:.2f}")
    else:
        print(f"  SKIP — {mat} not found (run run_lustig_sweep.py first)")

    # montage: central axis-2 slice, shared [0,1] axis, grouped by pipeline
    keys = list(vols.keys())
    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow))
    for ax in np.atleast_1d(axes).ravel():
        ax.axis("off")
    for i, k in enumerate(keys):
        ax = np.atleast_1d(axes).ravel()[i]
        ax.imshow(norm01(vols[k])[:, :, N // 2], cmap="gray", vmin=0, vmax=1)
        m = results[k]
        snr = "inf" if m["snr"] == float("inf") else f"{m['snr']:.1f}"
        ax.set_title(f"{k}\nSNR {snr} lfCV {m['lowfreq_cv']:.3f} "
                     f"edge {m['edge_sharp']:.1f}", fontsize=8)
    fig.tight_layout()
    png = os.path.join(a.recon_io, "tv_threeway_montage.png")
    fig.savefig(png, dpi=110)
    print(f"wrote {png}")

    with open(os.path.join(a.recon_io, "tv_threeway_metrics.json"), "w") as f:
        json.dump({"metric": "metrics_v2", "results": results}, f, indent=2)
    print(f"wrote {os.path.join(a.recon_io, 'tv_threeway_metrics.json')}")


if __name__ == "__main__":
    main()

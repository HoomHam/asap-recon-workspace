"""Six-row z-registered fine-z montage: Faraz, Steve, our CG, ours wav, BART, Lustig.

Extends z_register_compare.py to all pipelines. Every volume is (a) oriented to
our-CG (perm+flip), (b) resampled to 100^3 if needed (Faraz is 80^3), (c) z-affine
registered (scale+shift) to our-CG so the ACR resolution insert lands on the same
z across rows (phantoms differ ~7% in size — see [[slice-matching-zaffine]]).

    ../.venv/bin/python zreg_sixway_montage.py ../../data/v3_fov250/recon_io \
        --bart /Users/hoomham/bin/bart-src/bart --z0 44 --z1 61

Outputs: recon_io/zreg_sixway_montage.png
"""

import argparse
import itertools
import json
import os

import numpy as np
import scipy.io as sio
from scipy.ndimage import zoom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sigpy as sp

import asap_recon as ar
import cs_recon as cs
from steve_kernel_numpy import steve_recon
from bartio import readcfl
from bart_compare import our_setup
from z_register_compare import fit_z_affine, resample_z

N = cs.N


def orient_any(vol, ref):
    """Perm+flip of |vol| best matching |ref| by correlation (handles any shape
    via a zoomed copy for scoring); returns oriented full-res |vol|."""
    a = np.abs(vol)
    r = np.abs(ref)
    a_s = zoom(a, np.array(r.shape) / np.array(a.shape), order=1)
    best = (-2.0, None)
    for perm in itertools.permutations(range(3)):
        for fl in itertools.product((False, True), repeat=3):
            t = np.transpose(a_s, perm)
            for ax, f in enumerate(fl):
                if f:
                    t = np.flip(t, axis=ax)
            c = np.corrcoef(t.ravel(), r.ravel())[0, 1]
            if c > best[0]:
                best = (c, (perm, fl))
    _, (perm, fl) = best
    out = np.transpose(np.abs(vol), perm)
    for ax, f in enumerate(fl):
        if f:
            out = np.flip(out, axis=ax)
    return out


def to_ref_frame(vol, ref, z0, z1, label):
    """Orient -> resample to 100^3 -> z-affine register to ref."""
    v = orient_any(vol, ref)
    if v.shape != (N, N, N):
        v = zoom(v, np.array([N, N, N]) / np.array(v.shape), order=1)
    s, d, c = fit_z_affine(v, np.abs(ref), z0, z1)
    print(f"  {label:12s} z-affine scale {s:.3f} shift {d:+.2f}  corr {c:.3f}")
    return resample_z(v, s, d)


def n01(s):
    p = np.percentile(s, 99.5)
    return np.clip(s / (p + 1e-12), 0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("--bart", default="/Users/hoomham/bin/bart-src/bart")
    ap.add_argument("--our-t", type=float, default=0.003)
    ap.add_argument("--bart-l", type=float, default=0.0003)
    ap.add_argument("--lustig-tv", type=float, default=0.001)
    ap.add_argument("--z0", type=int, default=44)
    ap.add_argument("--z1", type=int, default=61)
    ap.add_argument("--iters", type=int, default=100)
    a = ap.parse_args()
    D = a.recon_io

    o = our_setup(D)
    ref = np.abs(o["img_cg"])             # our CG = reference frame
    d = ar.load_steve_npy(D)
    meta = json.load(open(os.path.join(D, "meta.json")))

    print("recon + register each pipeline ...")
    faraz = np.abs(sio.loadmat(os.path.join(D, "faraz", "faraz_recon.mat"),
                               squeeze_me=True)["img_gp"])
    steve = steve_recon(d["traj"], d["acq"], npts=meta["npts"], MS=meta["MS"],
                        IS=meta["IS"], smoothing=meta["gplb"], axes="xyz")
    W = sp.linop.Wavelet(o["A"].ishape, wave_name="db4")
    tref = float(np.percentile(np.abs(W(o["img_cg"])), 99))
    ourwav = np.abs(cs.wavelet_recon(o["A"], o["y_n"], a.our_t * tref, a.iters))
    bart = np.squeeze(readcfl(os.path.join(D, "bart", f"bart_W_l{a.bart_l:g}")))
    g = sio.loadmat(os.path.join(D, "lustig", "lustig_tv_sweep.mat"))
    gf = np.asarray(g["gas_final"], float)
    tvW = np.asarray(g["tvW"], float).ravel()
    lus = gf[int(np.argmin(np.abs(tvW - a.lustig_tv)))]

    rows = [
        ("Faraz (MATLAB)", to_ref_frame(faraz, ref, a.z0, a.z1, "Faraz")),
        ("Steve-equiv",    to_ref_frame(steve, ref, a.z0, a.z1, "Steve")),
        ("our CG-20",      to_ref_frame(ref,   ref, a.z0, a.z1, "ourCG")),
        (f"our wav t{a.our_t:g}", to_ref_frame(ourwav, ref, a.z0, a.z1, "ourWav")),
        (f"BART wav l{a.bart_l:g}", to_ref_frame(bart, ref, a.z0, a.z1, "BART")),
        (f"Lustig TV w{a.lustig_tv:g}", to_ref_frame(lus, ref, a.z0, a.z1, "Lustig")),
    ]

    zs = list(range(a.z0, a.z1, 2))
    fig, axes = plt.subplots(len(rows), len(zs), figsize=(2.0 * len(zs), 2.0 * len(rows)))
    for ri, (lab, v) in enumerate(rows):
        for ci, z in enumerate(zs):
            ax = axes[ri, ci]
            ax.imshow(n01(np.abs(v[:, :, z])), cmap="gray", vmin=0, vmax=1)
            ax.axis("off")
            if ri == 0:
                ax.set_title(f"z={z}", fontsize=9)
            if ci == 0:
                ax.text(-0.22, 0.5, lab, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=10)
    fig.suptitle("z-affine registered to our CG — ACR grid insert aligns across all "
                 "6 pipelines (Faraz / Steve / our CG / our wav / BART / Lustig)",
                 fontsize=12)
    fig.tight_layout()
    png = os.path.join(D, "zreg_sixway_montage.png")
    fig.savefig(png, dpi=115)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()

"""Z-affine register ours/BART/Lustig then matched-slice compare.

Hooman (2026-06-15): the grid IS visible in all three, but the phantoms are
slightly different SIZES, so a flip-only orientation aligns the bulk yet leaves
the thin grid-insert slice at different z in each -> at a fixed z one pipeline is
on the grid, the others on a uniform plane.

Fix: after the 48-orientation (flips/perms), register each volume to ours along z
with a 1-D AFFINE (scale s + shift d), found by maximizing the mean slicewise
correlation. A scale (not just a shift) is needed because the size differs, so a
single offset can't align both ends of the phantom. Then a fine-z montage: if the
grid appears in the SAME column across all three rows, the slices are matched.

    ../.venv/bin/python z_register_compare.py ../../data/v3_fov250/recon_io \
        --bart /Users/hoomham/bin/bart-src/bart \
        --our-t 0.003 --bart-l 0.0003 --lustig-tv 0.001 --z0 44 --z1 60

Outputs: recon_io/zreg_finez_montage.png  (+ prints fitted s,d per pipeline)
"""

import argparse
import itertools
import os

import numpy as np
import scipy.io as sio
from scipy.ndimage import map_coordinates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sigpy as sp

import asap_recon as ar
import cs_recon as cs
from bartio import readcfl
from bart_compare import our_setup

N = cs.N


def orient_flips(vol, ref):
    a = np.abs(ref); a = (a - a.mean()) / (a.std() + 1e-12)
    best, bc = np.abs(vol), -9.0
    for p in itertools.permutations(range(3)):
        vp = np.transpose(np.abs(vol), p)
        for fl in itertools.product([0, 1], repeat=3):
            w = vp
            for ax, dd in enumerate(fl):
                if dd:
                    w = np.flip(w, ax)
            b = (w - w.mean()) / (w.std() + 1e-12)
            c = float((a * b).mean())
            if c > bc:
                bc, best = c, w
    return best


def resample_z(vol, s, d):
    """Resample vol along z at z_src = s*(z - N/2) + N/2 + d for each output z."""
    zout = np.arange(N)
    zsrc = s * (zout - N / 2) + N / 2 + d
    out = np.empty_like(vol)
    for z in range(N):
        # bilinear in-plane is identity; only z is resampled -> interp slices
        zl = zsrc[z]
        out[:, :, z] = map_coordinates(
            vol, np.array([
                np.repeat(np.arange(N), N),
                np.tile(np.arange(N), N),
                np.full(N * N, zl)]),
            order=1, mode="nearest").reshape(N, N)
    return out


def fit_z_affine(vol, ref, z0, z1, smin=0.80, smax=1.30):
    """Grid-search scale s, shift d maximizing mean slicewise corr over [z0,z1].
    Wider scale range (0.80-1.30) covers Faraz, whose phantom is larger in z."""
    best, bs, bd = -9.0, 1.0, 0.0
    for s in np.linspace(smin, smax, 34):
        for d in np.linspace(-6, 6, 25):
            zsrc = s * (np.arange(z0, z1) - N / 2) + N / 2 + d
            cc, n = 0.0, 0
            for zi, zs in zip(range(z0, z1), zsrc):
                if 1 <= zs <= N - 2:
                    sl = (1 - (zs % 1)) * vol[:, :, int(zs)] + (zs % 1) * vol[:, :, int(zs) + 1]
                    a = ref[:, :, zi]
                    if a.std() > 0 and sl.std() > 0:
                        cc += np.corrcoef(a.ravel(), sl.ravel())[0, 1]; n += 1
            if n and cc / n > best:
                best, bs, bd = cc / n, s, d
    return bs, bd, best


def n01(s):
    p = np.percentile(s, 99.5); return np.clip(s / (p + 1e-12), 0, 1)


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
    ref = np.abs(o["img_cg"])
    W = sp.linop.Wavelet(o["A"].ishape, wave_name="db4")
    tref = float(np.percentile(np.abs(W(o["img_cg"])), 99))
    ours = np.abs(cs.wavelet_recon(o["A"], o["y_n"], a.our_t * tref, a.iters))

    bart = orient_flips(np.squeeze(readcfl(os.path.join(D, "bart", f"bart_W_l{a.bart_l:g}"))), ref)
    g = sio.loadmat(os.path.join(D, "lustig", "lustig_tv_sweep.mat"))
    gf = np.asarray(g["gas_final"], float); tvW = np.asarray(g["tvW"], float).ravel()
    li = int(np.argmin(np.abs(tvW - a.lustig_tv)))
    lus = orient_flips(gf[li], ref)

    # register BART, Lustig to OURS (ours is the reference frame)
    sb, db, cb = fit_z_affine(bart, ours, a.z0, a.z1)
    sl_, dl, cl = fit_z_affine(lus, ours, a.z0, a.z1)
    print(f"BART  z-affine: scale {sb:.3f} shift {db:+.2f}  corr {cb:.3f}")
    print(f"Lustig z-affine: scale {sl_:.3f} shift {dl:+.2f}  corr {cl:.3f}")
    bart_r = resample_z(bart, sb, db)
    lus_r = resample_z(lus, sl_, dl)

    zs = list(range(a.z0, a.z1, 2))
    rows = [(f"ours wav t{a.our_t:g}", ours),
            (f"BART wav l{a.bart_l:g} (z-reg)", bart_r),
            (f"Lustig TV w{tvW[li]:g} (z-reg)", lus_r)]
    fig, axes = plt.subplots(3, len(zs), figsize=(2.1 * len(zs), 6.8))
    for ri, (lab, v) in enumerate(rows):
        for ci, z in enumerate(zs):
            ax = axes[ri, ci]
            ax.imshow(n01(np.abs(v[:, :, z])), cmap="gray", vmin=0, vmax=1)
            ax.axis("off")
            if ri == 0:
                ax.set_title(f"z={z}", fontsize=9)
            if ci == 0:
                ax.text(-0.2, 0.5, lab, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=10)
    fig.suptitle("After z-affine registration to ours — grid insert should align "
                 "in the SAME column across all rows", fontsize=11)
    fig.tight_layout()
    png = os.path.join(D, "zreg_finez_montage.png")
    fig.savefig(png, dpi=115)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()

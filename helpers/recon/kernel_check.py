"""Show the 1D S-I 'kernel' (profile) the diaphragm detector works on, for a few
windows across one breath, so we can SEE what each position metric finds:

  - centroid  : signal centre-of-mass (sits in the MIDDLE of the lung)
  - edge      : half-max crossings = the lung boundaries; the diaphragm is the
                INFERIOR one (the boundary that moves with breathing)

Left column: coronal nav image (S-I vertical) with centroid (orange) + both
half-max edges (red/green). Right column: the 1D profile with the same markers.

Usage: ../.venv/bin/python kernel_check.py <dump_dir> [--win-ilv 20] [--nav-n 80]
       [--start auto|<win>] [--n 6]
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import asap_recon as ar

SI = 2  # S-I axis (empirically; see surrogates.SI_AXIS)


def profile(img):
    other = tuple(a for a in range(3) if a != SI)
    return np.sum(np.abs(img), axis=other)


def centroid(p0):
    s = p0.sum()
    return float((np.arange(len(p0)) * p0).sum() / s) if s > 0 else None


def half_max_edges(p0):
    """First/last indices where the bg-subtracted profile exceeds half its max."""
    if p0.max() <= 0:
        return None, None
    above = np.where(p0 > 0.5 * p0.max())[0]
    return (float(above[0]), float(above[-1])) if above.size else (None, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--win-ilv", type=int, default=20)
    ap.add_argument("--nav-n", type=int, default=80)
    ap.add_argument("--start", default="auto")
    ap.add_argument("--n", type=int, default=6)
    args = ap.parse_args()

    meta = json.load(open(os.path.join(args.dump_dir, "meta.json")))
    npts, nilv, MS = meta["npts"], meta["ntotalilvs"], meta["MS"]
    N = args.nav_n
    traj = np.stack([np.load(os.path.join(args.dump_dir, f"traj{a}.npy")) for a in "xyz"], 1)
    y = np.ascontiguousarray(np.load(os.path.join(args.dump_dir, "acq_dyn.npy"))[0],
                             dtype=np.complex128)
    tt = ar.tile_traj(traj.astype(float), npts * nilv)
    n_win = nilv // args.win_ilv
    start = n_win // 3 if args.start == "auto" else int(args.start)
    wins = list(range(start, min(start + args.n, n_win)))

    fig, axes = plt.subplots(len(wins), 2, figsize=(8, 2.3 * len(wins)))
    for r, w in enumerate(wins):
        sl = slice(w * args.win_ilv * npts, (w + 1) * args.win_ilv * npts)
        img = np.abs(ar.recon(tt[sl], y[sl], method="cg", MS=MS, IS=N, cg_iters=20))
        p = profile(img)
        p0 = p - np.median(p); p0[p0 < 0] = 0
        cen = centroid(p0)
        e0, e1 = half_max_edges(p0)

        # coronal image: project over axis0 -> (axis1, axis2); show S-I (axis2) vertical
        cor = np.abs(img).max(axis=0)                 # (axis1, axis2)
        ax_im = axes[r, 0]
        ax_im.imshow(cor.T, cmap="gray", origin="lower", aspect="auto",
                     vmax=np.percentile(cor, 99.5))
        if cen is not None:
            ax_im.axhline(cen, color="orange", lw=1.3, label="centroid")
        for e, col, lab in [(e0, "red", "edge lo"), (e1, "lime", "edge hi")]:
            if e is not None:
                ax_im.axhline(e, color=col, lw=1.2, ls="--", label=lab)
        ax_im.set_ylabel(f"win {w}\nS-I (axis2)", fontsize=8)
        if r == 0:
            ax_im.legend(fontsize=6, loc="upper right")
            ax_im.set_title("coronal nav (axis1 x S-I)", fontsize=9)

        ax_p = axes[r, 1]
        ax_p.plot(np.arange(N), p0, "k-", lw=1)
        ax_p.axhline(0.5 * p0.max(), color="gray", ls=":", lw=.8)
        if cen is not None:
            ax_p.axvline(cen, color="orange", lw=1.3)
        for e, col in [(e0, "red"), (e1, "lime")]:
            if e is not None:
                ax_p.axvline(e, color=col, lw=1.2, ls="--")
        if r == 0:
            ax_p.set_title("1D S-I profile (the 'kernel')", fontsize=9)
        ax_p.set_xlabel("S-I index" if r == len(wins) - 1 else "")

    fig.suptitle("What the diaphragm detector sees: centroid=middle vs half-max edges",
                 fontsize=10)
    fig.tight_layout()
    out = os.path.join(args.dump_dir, "kernel_check.png")
    fig.savefig(out, dpi=120); plt.close(fig)
    print(f"wrote {out}  (windows {wins})")


if __name__ == "__main__":
    main()

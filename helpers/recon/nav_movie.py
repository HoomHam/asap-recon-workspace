"""Record the per-window CS-nav undersampled lung images as a movie (QA for the
diaphragm surrogate).

Reconstructs each consecutive-interleave window (same windows the diaphragm
surrogate uses), projects a coronal view (S-I vertical), overlays the detected
diaphragm position, and writes an MP4 + GIF + a montage grid. Watch this to see
whether the nav images are real moving lungs or noise -- if the diaphragm curve
looks wrong, this shows why.

Usage:
    ../.venv/bin/python nav_movie.py <dump_dir> [--win-ilv 40] [--nav-n 80]
        [--method cg|wavelet] [--metric centroid|edge] [--stride 1] [--fps 8]
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import asap_recon as ar
import surrogates as S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--win-ilv", type=int, default=20)
    ap.add_argument("--nav-n", type=int, default=80)
    ap.add_argument("--method", choices=["cg", "wavelet"], default="cg")
    ap.add_argument("--metric", choices=["centroid", "edge"], default="centroid")
    ap.add_argument("--view", choices=["coronal", "sagittal"], default="sagittal",
                    help="DISPLAY projection. sagittal (over axis0)=more lung/SNR; "
                         "coronal (over axis1)=nicer to eye. Detection is the same "
                         "full S-I profile either way.")
    ap.add_argument("--stride", type=int, default=1, help="use every Nth window")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(args.dump_dir, f"nav_movie_{args.view}")
    os.makedirs(out, exist_ok=True)

    meta = json.load(open(os.path.join(args.dump_dir, "meta.json")))
    npts, nilv, MS = meta["npts"], meta["ntotalilvs"], meta["MS"]
    N = args.nav_n
    traj = np.stack([np.load(os.path.join(args.dump_dir, f"traj{a}.npy")) for a in "xyz"], 1)
    y = np.ascontiguousarray(np.load(os.path.join(args.dump_dir, "acq_dyn.npy"))[0],
                             dtype=np.complex128)
    ilvtime = np.load(os.path.join(args.dump_dir, "ilvtime.npy"))
    traj_tiled = ar.tile_traj(traj.astype(float), npts * nilv)
    n_win = nilv // args.win_ilv

    frames, cens, los, his, tlabels, wlabels = [], [], [], [], [], []
    for w in range(0, n_win, args.stride):
        i0, i1 = w * args.win_ilv, (w + 1) * args.win_ilv
        sl = slice(i0 * npts, i1 * npts)
        yw = y[sl]
        if np.count_nonzero(yw) < 0.5 * len(yw):
            continue
        img = ar.recon(traj_tiled[sl], yw, method="cg", MS=MS, IS=N, cg_iters=20)
        # DETECTION (_halfmax_edges / centroid) uses the full S-I profile (sum over
        # both in-plane axes), so it keeps all signal regardless of the DISPLAY view.
        # DISPLAY projection: coronal=sum over axis1 -> (axis0, S-I); sagittal=sum
        # over axis0 -> (axis1, S-I). Both keep S-I=axis2 as the 2nd dim.
        proj = (np.abs(img).sum(axis=0) if args.view == "sagittal"
                else np.abs(img).sum(axis=1))         # (other-axis, S-I=axis2)
        cen = S._centroid_position(img)
        lo, hi = S._halfmax_edges(img)                # lung boundaries on S-I; hi = diaphragm
        frames.append(proj)
        cens.append(cen); los.append(lo); his.append(hi)
        tlabels.append(float(np.mean(ilvtime[i0:i1])))
        wlabels.append(w)
        if w % 25 == 0:
            print(f"  nav window {w}/{n_win}  centroid={cen} edges=({lo},{hi})")

    frames = np.array(frames)
    vmax = np.percentile(frames, 99.5)
    print(f"{len(frames)} frames, nav N={N}, vmax={vmax:.3e}")

    def _h(ax, val, **kw):
        return ax.axhline(val if val is not None else N / 2, **kw)

    # --- animation: coronal lung (S-I vertical, inferior at bottom). hi edge =
    #     diaphragm (solid red), lo edge = apex (cyan dashed), centroid (orange). ---
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(frames[0].T, cmap="gray", vmin=0, vmax=vmax, origin="lower", aspect="auto")
    ax.invert_yaxis()                                 # display-only: inferior (hi) at bottom
    l_hi = _h(ax, his[0], color="cyan", lw=1.4, ls="--", label="diaphragm (hi edge)")
    ax.legend(fontsize=6, loc="upper right")
    txt = ax.set_title("")
    ax.set_xlabel("axis1" if args.view == "sagittal" else "axis0")
    ax.set_ylabel("S-I (axis2), inferior down")

    def update(k):
        im.set_data(frames[k].T)
        if his[k] is not None:
            l_hi.set_ydata([his[k], his[k]])
        hi = f"{his[k]:.1f}" if his[k] is not None else "n/a"
        txt.set_text(f"win {wlabels[k]}  t={tlabels[k]:.1f}s  diaph={hi}")
        return im, l_hi, txt

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
    mp4 = os.path.join(out, "nav_movie.mp4")
    gif = os.path.join(out, "nav_movie.gif")
    anim.save(mp4, writer=animation.FFMpegWriter(fps=args.fps))
    anim.save(gif, writer=animation.PillowWriter(fps=args.fps))
    plt.close(fig)
    print(f"wrote {mp4}\nwrote {gif}")

    # --- montage of a subset ---
    npanel = min(24, len(frames))
    idx = np.linspace(0, len(frames) - 1, npanel).astype(int)
    cols = 6; rows = int(np.ceil(npanel / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(2 * cols, 2 * rows))
    for a in axes.ravel():
        a.axis("off")
    for n, k in enumerate(idx):
        a = axes.ravel()[n]
        a.imshow(frames[k].T, cmap="gray", vmin=0, vmax=vmax, origin="lower", aspect="auto")
        a.invert_yaxis()
        if his[k] is not None:
            a.axhline(his[k], color="cyan", lw=1.1, ls="--")  # diaphragm (only line)
        a.set_title(f"t={tlabels[k]:.0f}s", fontsize=7)
    fig.suptitle(f"nav {args.view}: cyan dashed = diaphragm (hi edge) "
                 f"[{args.win_ilv} ilv/win, N={N}]")
    fig.tight_layout()
    mont = os.path.join(out, "nav_montage.png")
    fig.savefig(mont, dpi=110); plt.close(fig)
    print(f"wrote {mont}")


if __name__ == "__main__":
    main()

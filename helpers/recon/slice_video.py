"""All-slice cine video: grid of every SI slice animated across respiratory bins.

Each frame = one bin; slices arranged in a grid (10 per row along SI axis).
Reads cine_joint.npy (or cine_baseline.npy) from cine_<surrogate>/.

Usage:
    ../.venv/bin/python slice_video.py <dump_dir> [--surrogate diaphragm]
        [--use baseline|joint] [--axis 0|1|2] [--slices-per-row 10] [--fps 4]
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# SI_AXIS=2: third spatial axis of the (N,N,N) image volume.
# cine shape is (B, N, N, N); spatial axis 2 = cine dim 3.
SI_AXIS = 2


def slice_grid(vol, axis, slices_per_row):
    """Return a (rows*H, cols*W) grid image from vol (N,N,N), sliced along axis."""
    n = vol.shape[axis]
    n_rows = int(np.ceil(n / slices_per_row))
    # per-slice shape after taking along axis
    other = [vol.shape[a] for a in range(3) if a != axis]
    H, W = other[0], other[1]
    grid = np.zeros((n_rows * H, slices_per_row * W))
    for s in range(n):
        r, c = s // slices_per_row, s % slices_per_row
        sl = np.take(vol, s, axis=axis)   # shape (H, W)
        grid[r * H:(r + 1) * H, c * W:(c + 1) * W] = sl
    return grid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--surrogate", default="diaphragm",
                    choices=["signal", "pneumo", "diaphragm"])
    ap.add_argument("--use", default="joint", choices=["baseline", "joint"])
    ap.add_argument("--axis", type=int, default=SI_AXIS,
                    help="spatial axis to slice along (0/1/2); default=2 (SI)")
    ap.add_argument("--slices-per-row", type=int, default=10)
    ap.add_argument("--fps", type=int, default=4)
    args = ap.parse_args()

    cine_dir = os.path.join(args.dump_dir, f"cine_{args.surrogate}")
    cine_path = os.path.join(cine_dir, f"cine_{args.use}.npy")
    if not os.path.exists(cine_path):
        raise SystemExit(f"not found: {cine_path}")

    cine = np.abs(np.load(cine_path))          # (B, N, N, N)
    B = cine.shape[0]
    N = cine.shape[1 + args.axis]
    vmax = np.percentile(cine, 99.5)
    print(f"cine shape {cine.shape}  B={B}  slicing axis {args.axis}  N={N}  vmax={vmax:.3e}")

    # figure sizing: each slice cell is displayed at 1 inch / 100 dpi = 1px per voxel
    cell = cine.shape[1 + ((args.axis + 1) % 3)]   # one of the non-slice dims (square anyway)
    n_rows = int(np.ceil(N / args.slices_per_row))
    fig_w = args.slices_per_row * cell / 100
    fig_h = n_rows * cell / 100 + 0.35   # +0.35 for title

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    frame0 = slice_grid(cine[0], axis=args.axis, slices_per_row=args.slices_per_row)
    im = ax.imshow(frame0, cmap="gray", vmin=0, vmax=vmax, origin="lower", aspect="equal")
    ttl = ax.set_title(
        f"bin 0/{B-1} — {args.surrogate} surrogate, axis {args.axis} ({args.use})",
        fontsize=8)
    fig.tight_layout(pad=0.3)

    def update(b):
        im.set_data(slice_grid(cine[b], axis=args.axis, slices_per_row=args.slices_per_row))
        ttl.set_text(f"bin {b}/{B-1} — {args.surrogate} surrogate, axis {args.axis} ({args.use})")
        return [im, ttl]

    ani = animation.FuncAnimation(fig, update, frames=B, blit=False, interval=1000 // args.fps)

    stem = f"all_slices_{args.use}_axis{args.axis}"
    mp4 = os.path.join(cine_dir, f"{stem}.mp4")
    ani.save(mp4, writer="ffmpeg", fps=args.fps, dpi=100)
    print(f"wrote {mp4}")

    gif = os.path.join(cine_dir, f"{stem}.gif")
    ani.save(gif, writer="pillow", fps=args.fps)
    print(f"wrote {gif}")

    plt.close(fig)


if __name__ == "__main__":
    main()

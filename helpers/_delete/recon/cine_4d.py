"""Driver: dynamic dump -> respiratory binning -> 4D CS cine + montage.

Runs the whole Stage 1-4 path for one surrogate:
    load dump -> get per-interleave surrogate -> soft bins -> Stage-3 baseline
    (per-bin wavelet) and/or Stage-4 joint (wavelet + temporal TV) -> save + montage.

Usage:
    ../.venv/bin/python cine_4d.py <dump_dir> --surrogate signal|pneumo|diaphragm \
        [--bins 16] [--sigma-bins 0.75] [--stage both] [--max-iter 60] \
        [--lam-s-rel 0.01] [--lam-t-rel 0.05] [--out <dir>]
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import binning
import cs_recon_4d as c4


def load_surrogate(dump_dir, which, npts, ntotalilvs):
    """Per-interleave volume in [0,1] for the chosen surrogate."""
    if which == "signal":
        return np.load(os.path.join(dump_dir, "ilvvol_signal.npy"))
    if which == "pneumo":
        p = os.path.join(dump_dir, "ilvvol_pneumo.npy")
        if not os.path.exists(p):
            raise SystemExit("no ilvvol_pneumo.npy in dump (rerun dump with --pneumotach)")
        return np.load(p)
    if which == "diaphragm":
        p = os.path.join(dump_dir, "ilvvol_diaphragm.npy")
        if os.path.exists(p):
            return np.load(p)
        import surrogates
        return surrogates.diaphragm_surrogate(dump_dir, metric="edge", win_ilv=20, smooth_win=5)
    raise SystemExit(f"unknown surrogate {which!r}")


def montage(cine, path, title, slice_axis=1):
    """Row of B panels: one mid-slice per bin, magnitude, shared p99.5 scale."""
    B = cine.shape[0]
    mag = np.abs(cine)
    vmax = np.percentile(mag, 99.5)
    sl = cine.shape[slice_axis] // 2
    fig, axes = plt.subplots(1, B, figsize=(1.6 * B, 2.2))
    for b in range(B):
        img = np.take(mag[b], sl, axis=slice_axis)
        axes[b].imshow(img.T, cmap="gray", vmin=0, vmax=vmax, origin="lower")
        axes[b].set_title(f"bin {b}", fontsize=7)
        axes[b].axis("off")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--surrogate", choices=["signal", "pneumo", "diaphragm"], default="signal")
    ap.add_argument("--bins", type=int, default=None, help="default = meta nbins")
    ap.add_argument("--sigma-bins", type=float, default=0.75)
    ap.add_argument("--stage", choices=["baseline", "joint", "both"], default="both")
    ap.add_argument("--max-iter", type=int, default=60)
    ap.add_argument("--lam-s-rel", type=float, default=0.01)
    ap.add_argument("--lam-t-rel", type=float, default=0.05)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    meta = json.load(open(os.path.join(args.dump_dir, "meta.json")))
    npts, nilv, MS, IS = meta["npts"], meta["ntotalilvs"], meta["MS"], meta["IS"]
    B = args.bins or meta["nbins"]
    out = args.out or os.path.join(args.dump_dir, f"cine_{args.surrogate}")
    os.makedirs(out, exist_ok=True)

    traj = np.stack([np.load(os.path.join(args.dump_dir, f"traj{a}.npy")) for a in "xyz"], 1)
    y = np.ascontiguousarray(np.load(os.path.join(args.dump_dir, "acq_dyn.npy"))[0],
                             dtype=np.complex128)
    exclude = np.load(os.path.join(args.dump_dir, "exclude_ilv.npy"))

    ilvvol = load_surrogate(args.dump_dir, args.surrogate, npts, nilv)
    phi, M_ilv = binning.membership_from_surrogate(ilvvol, n_bins=B,
                                                   sigma_bins=args.sigma_bins,
                                                   exclude=exclude)
    occ = M_ilv.sum(1)
    print(f"surrogate={args.surrogate}  B={B}  bin occupancy (eff. interleaves): "
          f"{np.array2string(occ, precision=0)}")
    M_samp = binning.tile_to_samples(M_ilv, npts)
    np.save(os.path.join(out, "membership_ilv.npy"), M_ilv)
    np.save(os.path.join(out, "phase.npy"), phi)

    baseline = info = None
    if args.stage in ("baseline", "both"):
        print("=== Stage 3: per-bin wavelet baseline ===")
        baseline, info = c4.recon_4d_baseline(traj, y, M_samp, N=IS, MS=MS,
                                              lam_s_rel=args.lam_s_rel,
                                              max_iter=args.max_iter)
        np.save(os.path.join(out, "cine_baseline.npy"), baseline.astype(np.complex64))
        montage(baseline, os.path.join(out, "montage_baseline.png"),
                f"{args.surrogate} baseline (wavelet, no temporal)")

    if args.stage in ("joint", "both"):
        print("=== Stage 4: joint wavelet + temporal TV ===")
        w_dcf, c = (info[0], info[1]) if info else (None, None)
        joint, baseline = c4.recon_4d_joint(traj, y, M_samp, N=IS, MS=MS,
                                            lam_s_rel=args.lam_s_rel,
                                            lam_t_rel=args.lam_t_rel,
                                            max_iter=args.max_iter, baseline=baseline,
                                            w_dcf=w_dcf, c=c)
        np.save(os.path.join(out, "cine_joint.npy"), joint.astype(np.complex64))
        montage(joint, os.path.join(out, "montage_joint.png"),
                f"{args.surrogate} joint (wavelet + temporal TV, "
                f"lam_t_rel={args.lam_t_rel})")

    print(f"done -> {out}")


if __name__ == "__main__":
    main()

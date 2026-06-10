"""Produce Steve's npy dumps WITHOUT GPU — pure CPU, uses his own loaders.

Replicates exactly what results.dyn_recon saves (results.py:258-262):
    trajx.npy / trajy.npy / trajz.npy   grid-unit trajectory (GPDYN)
    acq.npy                             channel-0 data, F-order flat, x mask
    bins.npy                            ilvbin * nbins (SIGNAL binning)
plus meta.json (npts, nuniquesmp, ntotalilvs, nch, MS, IS, TR) needed by
steve_kernel_numpy.py.

Usage:
    .venv/bin/python dump_inputs.py <meas.dat> <gp_traj.npy> [dp_traj.npy] [outdir]

Runs Steve's raw.load() (noise normalization, rephasing, spike filter,
exclusion ranges) so the dumped acq is bit-faithful to what his recon grids.
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")  # raw.py plots spectra fits; never block headless
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)  # Steve's modules live at repo root

from gtypes import gvar, imgtype, graddir, bintype  # noqa: E402
import raw as steve_raw                              # noqa: E402


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    datfile = sys.argv[1]
    gptraj = sys.argv[2]
    dptraj = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3].endswith(".npy") else "/nonexistent"
    outdir = sys.argv[-1] if not sys.argv[-1].endswith((".npy", ".dat")) else "."
    os.makedirs(outdir, exist_ok=True)

    g = gvar()

    trajec = steve_raw.traj()
    trajec.load(gptraj, dptraj, nusimg=1)
    if trajec.npts == 0:
        sys.exit(f"trajectory load failed for {gptraj}")
    print(f"traj: npts/ilv={trajec.npts}  nuniquesmp={trajec.nuniquesmp}")

    r = steve_raw.raw()
    r.load(g, trajec, [datfile], [os.path.getsize(datfile)], "siemens", [], [])
    print(f"raw: npts={r.npts} nch={r.nch} ntotalilvs={r.ntotalilvs} TR={r.TR}")

    trajec.rescale_to_MS(g.MS, g.IS)

    # mask from exclusion ranges — verbatim from results.dyn_recon
    mask = np.ones(r.npts * r.ntotalilvs)
    for iilv in range(r.ntotalilvs):
        for el in r.excluderanges:
            if iilv * r.TR >= el.start and iilv * r.TR <= el.stop:
                mask[(iilv * r.npts):((iilv + 1) * r.npts)] = 0

    np.save(os.path.join(outdir, "trajx"), trajec.gettraj(imgtype.GPDYN, graddir.X))
    np.save(os.path.join(outdir, "trajy"), trajec.gettraj(imgtype.GPDYN, graddir.Y))
    np.save(os.path.join(outdir, "trajz"), trajec.gettraj(imgtype.GPDYN, graddir.Z))
    np.save(os.path.join(outdir, "acq"),
            np.reshape(r.getimg(imgtype.GPDYN)[:, 0, :],
                       r.npts * r.ntotalilvs, order="F") * mask)
    ilvbin = r.ilvbin[bintype.SIGNAL]
    if len(ilvbin) > 0:
        b = np.asarray(ilvbin, dtype=float)
        b[np.isnan(b)] = -1.0 / g.nbins
        np.save(os.path.join(outdir, "bins"), b * g.nbins)

    with open(os.path.join(outdir, "meta.json"), "w") as f:
        json.dump({"npts": int(r.npts), "nuniquesmp": int(trajec.nuniquesmp),
                   "ntotalilvs": int(r.ntotalilvs), "nch": int(r.nch),
                   "MS": int(g.MS), "IS": int(g.IS), "TR": float(r.TR),
                   "gplb": float(g.gplb), "datfile": os.path.basename(datfile)},
                  f, indent=2)
    print(f"wrote trajx/y/z.npy, acq.npy, bins.npy, meta.json -> {outdir}")
    print("NOTE: acq is channel 0 only (matching results.py:261).")


if __name__ == "__main__":
    main()

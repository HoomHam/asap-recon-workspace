"""recon_io npy  ->  ACR_test.mat  (exact replica of spiral3d_frames_mat_hoom.ipynb)

Collapses the old 3-stage dance (MATLAB ACR_data -> Colab notebook -> ACR_test)
into one call, reading directly from a v3_fov250-style recon_io folder.

EXACTNESS — every step mirrors the Colab notebook so the CS output reproduces
the old analysis:
  * trajectory: recenter grid-index -> k centered at 0 (subtract MS/2), then
    the notebook's max-radius normalization  k = k / |k|.max() * pi   ([-pi,pi])
  * kdata: real + 1j*imag, raw (MATLAB script does the /max scaling)
  * DCF: torchkbnufft pipe, tkbn.calc_density_compensation_function with
    im_size=(100,100,100) — the FIXED size (the (256,64,256) cardiac-leftover
    bug is not present in this notebook version)
  * saved keys/shapes match the existing ACR_test.mat exactly:
        ktrajs (1,3,N) f8 [-pi,pi] | kdatas (1,N) c128 | kcomps (1,N) c128

The notebook builds ACR_data from trajectory in 1/mm; recon_io stores it in
grid-index units (g = k*MS/IS + MS/2). The notebook's max-normalization is
scale-invariant, so only the CENTER offset matters — handled by subtracting
MS/2. Geometry then matches both the old Lustig run and our finufft pipeline.
"""

import argparse
import json
import os

import numpy as np
import scipy.io as sio
import torch
import torchkbnufft as tkbn

IM_SIZE = (100, 100, 100)  # notebook im_size (k_size=slices=100)


def build(recon_io, out_mat):
    meta = json.load(open(os.path.join(recon_io, "meta.json")))
    MS = float(meta["MS"])

    p = lambda f: os.path.join(recon_io, f)
    acq = np.load(p("acq.npy"))
    tx = np.load(p("trajx.npy")).astype(np.float64)
    ty = np.load(p("trajy.npy")).astype(np.float64)
    tz = np.load(p("trajz.npy")).astype(np.float64)
    assert len(acq) == len(tx) == len(ty) == len(tz)

    # grid-index -> k centered at 0 (matches asap_recon.grid_to_radians offset)
    dkx = tx - MS / 2.0
    dky = ty - MS / 2.0
    dkz = tz - MS / 2.0

    # notebook normalization: max radius -> pi
    dk_mag = np.sqrt(dkx**2 + dky**2 + dkz**2)
    s = np.pi / dk_mag.max()
    ktraj = np.stack([dkx * s, dky * s, dkz * s], axis=0)  # (3, N) in [-pi,pi]

    device = torch.device("cpu")
    ktraj_t = torch.tensor(ktraj, dtype=torch.float64).to(device)
    dcomp = tkbn.calc_density_compensation_function(
        ktraj=ktraj_t, im_size=IM_SIZE)
    dcomp = dcomp.cpu().squeeze().numpy()

    kdata = np.ascontiguousarray(acq, dtype=np.complex128)

    # match ACR_test.mat layout: cell-of-frames -> (1, ...) leading dim
    sio.savemat(out_mat, {
        "ktrajs": ktraj[None, :, :],            # (1,3,N)
        "kdatas": kdata[None, :],               # (1,N)
        "kcomps": dcomp.astype(np.complex128)[None, :],  # (1,N)
    })
    print(f"wrote {out_mat}: ktrajs (1,3,{len(acq)}) kdatas (1,{len(acq)}) "
          f"kcomps (1,{len(acq)})  | krange [{ktraj.min():.3f},{ktraj.max():.3f}]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("out_mat")
    a = ap.parse_args()
    build(a.recon_io, a.out_mat)

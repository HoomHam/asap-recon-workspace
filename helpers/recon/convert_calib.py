"""Convert a calibration .mat (cal struct, k in 1/mm) to Steve-format .npy.

Steve's traj.load does np.load(file) * FOV, so the .npy must hold k in 1/mm
(cycles/mm); x FOV[mm] then yields delta-k (cycles/FOV) as his code expects.
The cal struct already stores 1/mm -> straight passthrough with sanity checks.

Usage: .venv/bin/python convert_calib.py <calib.mat> <out_gp_traj.npy>
"""

import sys
import numpy as np
from scipy.io import loadmat


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    m = loadmat(sys.argv[1], squeeze_me=True, struct_as_record=False)
    cal = m["cal"]
    assert getattr(cal.units, "k") == "1/mm", f"unexpected k units: {cal.units.k}"
    k = np.stack([cal.kx, cal.ky, cal.kz], axis=1).astype(np.float64)
    fov = float(cal.FOV_mm)
    npts = int(cal.nsamples)
    assert k.shape[0] % npts == 0, "k length not a multiple of nsamples"
    kdk = np.linalg.norm(k, axis=1).max() * fov
    print(f"{k.shape[0]} pts = {k.shape[0]//npts} ilvs x {npts}; "
          f"FOV {fov:.0f} mm; |k|max x FOV = {kdk:.1f} delta-k "
          f"(imgsize/2 = {int(cal.imgsize)/2:.0f})")
    np.save(sys.argv[2], k)
    print(f"wrote {sys.argv[2]}  — run dump_inputs.py with --fov {fov:.0f}")


if __name__ == "__main__":
    main()

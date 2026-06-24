"""Export cine_joint.npy to NIfTI (.nii.gz) and MATLAB (.mat).

NIfTI: shape (N,N,N,B) = (100,100,100,16), float32 magnitude, 2.4mm isotropic
       (FOV=240mm / IS=100 vox). Affine = identity scaled by voxel size.
MATLAB: dict with keys 'cine' (B,N,N,N), 'voxel_mm', 'surrogate', 'patient'.

Usage:
    ../.venv/bin/python export_4d.py <cine_joint.npy> <out_dir>
        [--voxel-mm 2.4] [--surrogate diaphragm] [--patient 025JC]
"""

import argparse
import os
import numpy as np
import nibabel as nib
from scipy.io import savemat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cine_npy")
    ap.add_argument("out_dir")
    ap.add_argument("--voxel-mm", type=float, default=2.4)
    ap.add_argument("--surrogate", default="diaphragm")
    ap.add_argument("--patient", default="025JC")
    args = ap.parse_args()

    cine = np.abs(np.load(args.cine_npy)).astype(np.float32)  # (B,N,N,N)
    B, N1, N2, N3 = cine.shape
    print(f"cine shape {cine.shape}  voxel {args.voxel_mm} mm")

    os.makedirs(args.out_dir, exist_ok=True)
    stem = f"cine_4d_{args.patient}_{args.surrogate}"

    # NIfTI: (N,N,N,B)
    nii_data = np.transpose(cine, (1, 2, 3, 0))   # (N,N,N,B)
    affine = np.diag([args.voxel_mm, args.voxel_mm, args.voxel_mm, 1.0])
    img = nib.Nifti1Image(nii_data, affine)
    img.header.set_xyzt_units("mm", "sec")
    nii_path = os.path.join(args.out_dir, f"{stem}.nii.gz")
    nib.save(img, nii_path)
    print(f"wrote {nii_path}")

    # MATLAB
    mat_path = os.path.join(args.out_dir, f"{stem}.mat")
    savemat(mat_path, {
        "cine": cine,                    # (B,N,N,N) float32
        "voxel_mm": args.voxel_mm,
        "surrogate": args.surrogate,
        "patient": args.patient,
        "B": B,
        "N": np.array([N1, N2, N3]),
    }, do_compression=True)
    print(f"wrote {mat_path}")


if __name__ == "__main__":
    main()

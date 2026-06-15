"""One-shot Lustig CS on a recon_io folder. Replaces MATLAB->Colab->MATLAB.

    python run_lustig.py <recon_io> [--out-dir DIR] [--skip-dcf]

Steps:
  1. build_acrtest.py : recon_io npy -> ACR_test.mat (exact notebook DCF)
     (uses this folder's torch venv: .venv_lustig)
  2. MATLAB -batch run_cs : exact spiral3d_cs_3D_hoom.m -> gas(15,100,100,100)
  3. print per-iter + final cg_tune metrics

Outputs land in <recon_io>/lustig/ by default:
    ACR_test.mat, lustig_cs.mat (gas)
"""

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODES_ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "codes"))
DCF_PY = os.path.join(HERE, ".venv_lustig", "bin", "python")
RECON_VENV_PY = os.path.normpath(
    os.path.join(HERE, "..", "recon", "..", ".venv", "bin", "python"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--skip-dcf", action="store_true",
                    help="reuse existing ACR_test.mat (skip torch/DCF step)")
    ap.add_argument("--matlab", default="matlab")
    a = ap.parse_args()

    recon_io = os.path.abspath(a.recon_io)
    out_dir = a.out_dir or os.path.join(recon_io, "lustig")
    os.makedirs(out_dir, exist_ok=True)
    acrtest = os.path.join(out_dir, "ACR_test.mat")
    out_mat = os.path.join(out_dir, "lustig_cs.mat")

    # 1. DCF / ACR_test
    if not a.skip_dcf:
        py = DCF_PY if os.path.exists(DCF_PY) else sys.executable
        print(f"[1/3] build ACR_test ({py})")
        subprocess.run([py, os.path.join(HERE, "build_acrtest.py"),
                        recon_io, acrtest], check=True)
    elif not os.path.exists(acrtest):
        sys.exit(f"--skip-dcf but {acrtest} missing")

    # 2. MATLAB CS
    print(f"[2/3] MATLAB CS -> {out_mat}")
    cmd = (f"run_cs('{acrtest}','{out_mat}','{CODES_ROOT}')")
    subprocess.run([a.matlab, "-batch", cmd], check=True, cwd=HERE)

    # 3. metrics (use recon venv which has cg_tune)
    print("[3/3] metrics")
    metric_py = RECON_VENV_PY if os.path.exists(RECON_VENV_PY) else sys.executable
    snippet = (
        "import numpy as np,sys;"
        "sys.path.insert(0,r'%s');" % os.path.join(HERE, "..", "recon") +
        "from cg_tune import metrics;"
        "f=r'%s';" % out_mat +
        "import scipy.io as sio\n"
        "try:\n g=sio.loadmat(f)['gas']\n"
        "except NotImplementedError:\n"
        " import h5py;g=np.array(h5py.File(f,'r')['gas']).transpose(3,2,1,0)\n"
        "print('iters',g.shape)\n"
        "[print(i+1,metrics(np.asarray(g[i],float))) for i in range(g.shape[0])]")
    subprocess.run([metric_py, "-c", snippet], check=False)
    print(f"done. gas at {out_mat}")


if __name__ == "__main__":
    main()

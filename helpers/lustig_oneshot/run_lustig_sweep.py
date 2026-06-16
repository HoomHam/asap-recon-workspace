"""Lustig TV lambda-sweep on a recon_io folder, scored with metrics_v2.

NEW driver — does NOT modify run_lustig.py / run_cs.m. Makes the Lustig baseline
match BART's lambda sweep (handoff step 2): one fnlCg run per TVWeight, pure TV
(xfmWeight=0), final iterate scored with the FIXED metric (metrics_v2), so all
three pipelines (ours / BART / Lustig) finally sit on one ruler.

Reuses the existing ACR_test.mat (the torch/DCF step) unless --rebuild-dcf.

    .venv_lustig/bin/python run_lustig_sweep.py <recon_io> \
        --tvweights 1e-3 3e-3 1e-2 3e-2 1e-1

Outputs in <recon_io>/lustig/:
    lustig_tv_sweep.mat   gas_final(numW,100,100,100), tvW
    lustig_tv_sweep_metrics.json
"""

import argparse
import json
import os
import subprocess
import sys

import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
CODES_ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "codes"))
RECON_DIR = os.path.normpath(os.path.join(HERE, "..", "recon"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recon_io")
    ap.add_argument("--tvweights", type=float, nargs="+",
                    default=[1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
    ap.add_argument("--matlab", default="matlab")
    ap.add_argument("--rebuild-dcf", action="store_true",
                    help="rerun build_acrtest.py (torch); default reuses ACR_test.mat")
    a = ap.parse_args()

    recon_io = os.path.abspath(a.recon_io)
    out_dir = os.path.join(recon_io, "lustig")
    os.makedirs(out_dir, exist_ok=True)
    acrtest = os.path.join(out_dir, "ACR_test.mat")
    out_mat = os.path.join(out_dir, "lustig_tv_sweep.mat")

    if a.rebuild_dcf or not os.path.exists(acrtest):
        print(f"[1/3] build ACR_test")
        subprocess.run([sys.executable, os.path.join(HERE, "build_acrtest.py"),
                        recon_io, acrtest], check=True)
    else:
        print(f"[1/3] reuse {acrtest}")

    tvw = "[" + " ".join(repr(float(x)) for x in a.tvweights) + "]"
    cmd = f"run_cs_sweep('{acrtest}','{out_mat}','{CODES_ROOT}',{tvw})"
    print(f"[2/3] MATLAB TV sweep {a.tvweights} -> {out_mat}")
    subprocess.run([a.matlab, "-batch", cmd], check=True, cwd=HERE)

    print("[3/3] metrics_v2")
    g = sio.loadmat(out_mat)
    vols = np.asarray(g["gas_final"], dtype=float)        # (numW,100,100,100)
    tvW = np.asarray(g["tvW"], dtype=float).ravel()
    sys.path.insert(0, RECON_DIR)
    from metrics_v2 import metrics

    results = {}
    for i, lam in enumerate(tvW):
        m = metrics(vols[i])
        results[f"lustig_tv_l{lam:g}"] = {"TVWeight": float(lam), **m}
        snr = m["snr"]
        snr_s = "inf" if snr == float("inf") else f"{snr:.1f}"
        print(f"  TVWeight {lam:<7g} SNR {snr_s:>6} lfCV {m['lowfreq_cv']:.3f} "
              f"edge {m['edge_sharp']:.2f} ext {m['extent_mm']} "
              f"bgcol {m['bg_collapsed']}")

    with open(os.path.join(out_dir, "lustig_tv_sweep_metrics.json"), "w") as f:
        json.dump({"metric": "metrics_v2", "tvweights": list(map(float, tvW)),
                   "results": results}, f, indent=2)
    print(f"wrote {os.path.join(out_dir, 'lustig_tv_sweep_metrics.json')}")


if __name__ == "__main__":
    main()

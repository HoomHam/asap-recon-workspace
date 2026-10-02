"""Scratch: gradient-delay sweep on the measured 1H water trajectory (006KL), end-exp lines, autofocus metrics."""
import sys, os, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _proton_local as L
N, DX = 128, 3.0
raw, k, tr = L.load()
lines, b = L.endexp_lines(raw.shape[2])
print(f'end-exp lines: {len(lines)} of {raw.shape[2]}', flush=True)
t0 = time.time(); dcf = L.pipe_dcf(k, N, DX); print(f'DCF {time.time()-t0:.0f} s', flush=True)
x0 = L.recon(raw, k, lines, dcf, N, DX); mask = x0 > 0.15 * np.percentile(x0, 99)
res = {}
for tau in np.arange(-4, 4.01, 0.5):
    x = L.recon(raw, L.delay_traj(k, tau), lines, dcf, N, DX)
    e, g = L.sharpness(x, mask); res[float(tau)] = (e, g)
    print(f'tau {tau:+.1f} samples ({tau*5:+.1f} us): entropy {e:.4f}  grad {g:.5f}', flush=True)
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
json.dump(res, open(os.path.join(out, 'delay_iso.json'), 'w'), indent=1)

"""Scratch: off-resonance autofocus for a pre-v2 gas twix (shift 0). Demodulate exp(-i2pi f t_samp),
t_samp = n*dwell (5 us), recon (virtual coil) for f in a sweep, score = L4/L2 peakiness.
Usage: _prev2_freq_scan.py <gas.npz> <tag> [start_gas_scan] [n_gas_scans]
"""
import sys, pathlib
import numpy as np, finufft
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from _prev2_static_recon import pipe_dcf, FOV, N, TRAJ, OUT
npz, tag = sys.argv[1], sys.argv[2]
z = np.load(npz); g = z['g']
a0 = int(sys.argv[3]) if len(sys.argv) > 3 else 0
n = int(sys.argv[4]) if len(sys.argv) > 4 else g.shape[0] - a0
traj = np.load(TRAJ).reshape(-1, 3); ns = g.shape[2]; narms = traj.shape[0] // ns
idx = np.arange(a0, a0 + n); arm = idx % narms; gw = g[idx]
u, s, vh = np.linalg.svd(gw[:, :, 1:40].transpose(1, 0, 2).reshape(gw.shape[1], -1), full_matrices=False)
vc = np.einsum('c,ncs->ns', np.conj(u[:, 0]), gw)
kd0 = np.stack([vc[arm == i].mean(0) if (arm == i).any() else np.zeros(ns) for i in range(narms)])
kr = (traj * FOV) * 2 * np.pi / N
keep = np.ones((narms, ns), bool); keep[:, :2] = False; keep = keep.ravel()
kr = np.ascontiguousarray(kr[keep]); w = pipe_dcf(kr)
plan = finufft.Plan(1, (N, N, N), eps=1e-4); plan.setpts(kr[:, 0].copy(), kr[:, 1].copy(), kr[:, 2].copy())
tsamp = np.arange(ns) * 5e-6
res, ims = [], {}
for f in np.arange(-400, 401, 25):
    y = (kd0 * np.exp(-2j * np.pi * f * tsamp)[None]).ravel()[keep] * w
    im = np.abs(plan.execute(y.astype(np.complex128)))
    res.append((f, float((im ** 4).sum() / (im ** 2).sum() ** 2 * im.size))); ims[f] = im
r = np.array(res); best = r[np.argmax(r[:, 1]), 0]
print(f'{tag} window {a0}+{n}: best f {best:+.0f} Hz peaky {r[:,1].max():.2f} | f=0 peaky {r[r[:,0]==0,1][0]:.2f}')
print('  sweep:', [(int(a), round(b, 2)) for a, b in r[::2]])
np.save(OUT / f'freqscan_{tag}_{a0}.npy', np.stack([ims[0], ims[best]]))

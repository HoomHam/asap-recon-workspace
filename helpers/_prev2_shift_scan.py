"""Scratch: find the arm phase (shift) of a pre-v2 gas twix vs a candidate trajectory by brute force.
Input: gas npz from _prev2_static_recon.read_scans (g: gas scans x ch x 512). One SVD virtual coil,
NUFFT adjoint with Pipe DCF for every cyclic shift of the arm-averaged data; score = normalized
gradient energy (sharpness) and L4/L2 peakiness. Usage: _prev2_shift_scan.py <gas.npz> <tag> [traj.npy] [narms]
"""
import sys, pathlib, time
import numpy as np, finufft
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from _prev2_static_recon import pipe_dcf, FOV, N, TRAJ, OUT
npz, tag = sys.argv[1], sys.argv[2]
traj = np.load(sys.argv[3] if len(sys.argv) > 3 else TRAJ).reshape(-1, 3)
z = np.load(npz); g = z['g']
ns = g.shape[2]; narms = traj.shape[0] // ns
arm = np.arange(g.shape[0]) % narms
# virtual coil: first SVD component over coils
u, s, vh = np.linalg.svd(g[:, :, 1:40].transpose(1, 0, 2).reshape(g.shape[1], -1), full_matrices=False)
vc = np.einsum('c,ncs->ns', np.conj(u[:, 0]), g)
kd0 = np.stack([vc[arm == i].mean(0) for i in range(narms)])          # narms x ns
kr = (traj * FOV) * 2 * np.pi / N
keep = np.ones((narms, ns), bool); keep[:, :2] = False; keep = keep.ravel()
kr = np.ascontiguousarray(kr[keep]); w = pipe_dcf(kr)
plan = finufft.Plan(1, (N, N, N), eps=1e-3); plan.setpts(kr[:, 0].copy(), kr[:, 1].copy(), kr[:, 2].copy())
res = []
t0 = time.time()
for sft in range(narms):
    y = np.roll(kd0, sft, axis=0).ravel()[keep] * w
    im = np.abs(plan.execute(y.astype(np.complex128)))
    gm = sum(x ** 2 for x in np.gradient(im))
    res.append((sft, float(gm.sum() / (im ** 2).sum()), float((im ** 4).sum() / (im ** 2).sum() ** 2 * im.size)))
r = np.array(res); np.save(OUT / f'shiftscan_{tag}.npy', r)
o1 = np.argsort(r[:, 1])[::-1][:8]; o2 = np.argsort(r[:, 2])[::-1][:8]
print(f'{tag}: {narms} shifts in {time.time()-t0:.0f}s; median sharp {np.median(r[:,1]):.4f} peaky {np.median(r[:,2]):.2f}')
print(' top sharpness:', [(int(r[i,0]), round(r[i,1],4)) for i in o1])
print(' top peakiness:', [(int(r[i,0]), round(r[i,2],2)) for i in o2])
print(' shift0:', round(r[0,1],4), round(r[0,2],2))

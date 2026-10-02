"""Scratch (h1 challenge): first look. Coil whitening numbers, static gridded image, orientation vs Tyger p,
and how well local soft-bin gridding reproduces Tyger p."""
import os, sys, json, itertools, time
import numpy as np
import scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H

os.makedirs(H.EXT, exist_ok=True); os.makedirs(H.OUT, exist_ok=True)
D = H.load(); raw, k, vol = D['raw'], D['k'], D['vol']
L = raw.shape[1]
v, info = H.prewhiten(raw)
C = info['noise_cov']
print('noise cov diag', np.diag(C).real, 'corr', abs(C[0, 1]) / np.sqrt(C[0, 0].real * C[1, 1].real))
print('whitened signal eigenvalues (virtual ch 0, 1):', info['signal_eig'], 'ratio', info['signal_eig'][1] / info['signal_eig'][0])
# after whitening noise var per channel = 1: signal/noise per sample along the readout, each virtual channel
ss = np.arange(L) >= H.NILV
for s in (5, 20, 50, 100, 200, 300, 400, 500):
    print(f'  sample {s:3d}: rms |v0| {np.sqrt((np.abs(v[0, ss, s])**2).mean()):8.2f}   rms |v1| {np.sqrt((np.abs(v[1, ss, s])**2).mean()):6.2f}  (noise = 1)')
t0 = time.time(); dcf = H.pipe_dcf(k); print('dcf %.1f s' % (time.time() - t0)); np.save(os.path.join(H.EXT, 'dcf_n100.npy'), dcf)
lines = np.where(ss)[0]
op = H.Nufft(k, lines)
w = dcf[lines % H.NILV][:, H.KILL:]
t0 = time.time(); x0 = op.adj(v[0, lines, H.KILL:] * w); x1 = op.adj(v[1, lines, H.KILL:] * w); print('adjoint x2 %.1f s' % (time.time() - t0))
np.save(os.path.join(H.EXT, 'static_grid_v01.npy'), np.stack([x0, x1]).astype(np.complex64))
P = sio.loadmat(H.TYGER_P)['gas_phase'].astype(np.float32)            # (16, SI, AP, LR)
ref = P.mean(0); a = np.abs(x0)
best = (-1, None, None)
for perm in itertools.permutations(range(3)):
    ap = np.transpose(a, perm)
    for flip in itertools.product((0, 1), repeat=3):
        y = ap
        for ax, f in enumerate(flip):
            if f: y = np.flip(y, axis=ax)
        c = np.corrcoef(y.ravel(), ref.ravel())[0, 1]
        if c > best[0]: best = (c, perm, flip)
print('orientation: corr %.4f perm %s flip %s' % best)
np.savez(os.path.join(H.EXT, 'orient.npz'), perm=np.array(best[1]), flip=np.array(best[2]))
# sub-voxel shift check between local and Tyger grids
y = H.to_tyger(a)
from skimage.registration import phase_cross_correlation
sh, err, _ = phase_cross_correlation(ref, y, upsample_factor=20); print('shift local->tyger (vox)', sh)
# reproduce Tyger p: soft bins exp(-d^2/0.44) on Steve's cyclic coordinate, SOS of the two RAW channels vs virtual ch 0
c = H.cyc_coord(vol)
for s2 in (0.44, 2.0):
    W = H.soft_w(c, 16, s2)
    cs = []
    for b in (0, 4, 8, 12):
        wl = (W[lines, b][:, None] * w)
        xb = np.abs(H.to_tyger(op.adj(v[0, lines, H.KILL:] * wl)))
        cs.append(np.corrcoef(xb.ravel(), P[b].ravel())[0, 1])
    print(f's2 {s2}: corr local bin vs Tyger p bins 0,4,8,12:', np.round(cs, 4))
json.dump(dict(noise_cov_abs=np.abs(C).tolist(), signal_eig=info['signal_eig'].tolist(), orient_corr=float(best[0]),
               perm=list(map(int, best[1])), flip=list(map(int, best[2]))), open(os.path.join(H.OUT, 'orient.json'), 'w'), indent=1)

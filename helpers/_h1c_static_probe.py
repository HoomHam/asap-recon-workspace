"""Scratch probe: can one channel's pass-averaged data be fitted to the noise level by a 700 mm static image?"""
import os, sys, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H
NB = 200
KIT = np.load('/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/h1challenge/shared/judging_kit_v1.npz')
D = H.load(); k = D['k']; v, _ = H.prewhiten(D['raw']); L = v.shape[1]
train = np.setdiff1d(np.arange(H.NILV, L), KIT['heldout_lines']); il = train % H.NILV
cnt = np.bincount(il, minlength=H.NILV).astype(float)
ybar = np.zeros((2, H.NILV, H.NPTS - H.KILL), np.complex128)
for c in range(2): np.add.at(ybar[c], il, v[c, train, H.KILL:])
ybar /= cnt[None, :, None]
np.save(os.path.join(H.EXT, 'ybar_train.npy'), ybar.astype(np.complex64))
op = H.Nufft(k, np.arange(H.NILV), n=NB)
dcf = np.load(os.path.join(H.EXT, 'dcf_n200.npy'))
bands = ((0, 8), (8, 48), (48, 198), (198, 510))
def report(x, c, tag):
    r = op.fwd(x) - ybar[c]
    print(tag, ' '.join(f'[{a+2}-{b+2}] {(np.abs(r[:, a:b])**2).mean()*13:9.2f} (sig {(np.abs(ybar[c][:, a:b])**2).mean()*13:.3g})' for a, b in bands), flush=True)
for p in (1.0, 0.5):
    Wt = dcf ** p; Wt = Wt / Wt.sum()
    for c in (0, 1):
        rhs = op.adj(Wt * ybar[c]); sc = np.vdot(rhs, op.adj(Wt * op.fwd(rhs))).real / np.vdot(rhs, rhs).real
        nrm = lambda x: op.adj(Wt * op.fwd(x)) + 1e-4 * sc * x
        t0 = time.time()
        cb = lambda i, x, r: report(x, c, f'p={p} ch{c} it {i+1:3d} rr {r:.4f}') if (i + 1) in (10, 30, 80) else None
        x = H.cg(nrm, rhs, it=80, cb=cb)
        print('  %.0f s' % (time.time() - t0))
        if p == 1.0: np.save(os.path.join(H.EXT, f'static700_ls_ch{c}.npy'), x.astype(np.complex64))

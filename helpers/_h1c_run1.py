"""Scratch: first 4D least-squares runs. Held-out (line-wise) error against CG iteration and temporal weight."""
import os, sys, time, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
lam_t = float(sys.argv[1]); lam_s = float(sys.argv[2]); it = int(sys.argv[3]); tag = sys.argv[4]
d = DY.Dyn()
train = np.setdiff1d(d.ss, d.heldout)
log = []
def cb(i, x):
    if (i + 1) in (5, 10, 15, 20, 30, 40):
        U = d.fwd_all(x)
        ho = d.chi2(x, d.heldout, U); tr = d.chi2(x, train[::8], U)
        log.append(dict(it=i + 1, heldout=ho.tolist(), train=tr.tolist()))
        print(f'  it {i+1}: held-out ch0 {np.round(ho[0], 2)} ch1 {np.round(ho[1], 2)} | train ch0 {np.round(tr[0], 2)} ch1 {np.round(tr[1], 2)}', flush=True)
        np.save(os.path.join(H.EXT, f'x_{tag}_it{i+1}_crop.npy'), DY.crop(x).astype(np.complex64))
t0 = time.time()
x = d.recon(train, lam_s=lam_s, lam_t=lam_t, it=it, cb=cb)
print('total %.0f s' % (time.time() - t0))
json.dump(log, open(os.path.join(H.OUT, f'run1_{tag}.json'), 'w'), indent=1)

"""Scratch (h1 challenge): the ASAP arm pipeline, identical on the real scan and on the phantom.

1. weakly coupled 4D least squares (16 nodes, receiver filter, body support, principal channel)
2. breathing displacement field from it (demons end-expiration -> end-inspiration, rank 1, per-node amplitude)
3. optional: model-based self-navigation (each line's breathing coordinate corrected from the data), repeat 1-2
4. 4D least squares with the temporal coupling taken ALONG the motion (lam_t), for half A, half B, scored training set
5. optional linear post-filters: radial k-space weighting (= readout filter exp(-t/tau)), MP-PCA across the phases
6. kit v2 scores; on the phantom also the error against the truth

usage: _h1c_pipeline.py <real|v1|v1b> <tag> [--lt 0.5] [--ls 1e-3] [--it 15] [--it1 12] [--wpow 1.0] [--nav 0] [--apod 0] [--mppca 0] [--prior] [--sym] [--motion file.npz] [--display]
  --display: one extra reconstruction from ALL steady-state lines (the cine to show; not scored on held-out data)
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY, _h1c_score as SC, _h1c_motion as MO, _h1c_nav as NV, _h1c_post as PO, _h1c_phantom_score as PS
from _h1c_mppca import mppca

src, tag = sys.argv[1], sys.argv[2]
def opt(name, default, typ=float):
    return typ(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default
lt, it, it1, wpow, nnav, tau, mp, ls = opt('--lt', 0.5), opt('--it', 15, int), opt('--it1', 12, int), opt('--wpow', 1.0), opt('--nav', 0, int), opt('--apod', 0.0), opt('--mppca', 0, int), opt('--ls', 1e-3)
PD = os.path.join(H.EXT, 'phantom'); OD = os.path.join(H.EXT, 'arms'); os.makedirs(OD, exist_ok=True)
if src == 'real':
    d = DY.Dyn()
else:
    z = np.load(os.path.join(PD, f'cache_phantom_{src}.npz'))
    d = DY.Dyn(data=H.prewhiten(z['raw'])[0], vol=z['vol'].astype(np.float64))
d.wpow = wpow
kit = d.kit; tr = kit['scored_train_lines']; hold = np.concatenate([kit['heldout_lines'], d.ho_ilv_lines])
t0 = time.time(); log = dict(src=src, tag=tag, prior_sigma_vox=opt('--psig', 2.5), sym='--sym' in sys.argv, prior='--prior' in sys.argv, lam_s=ls, lt=lt, it=it, it1=it1, wpow=wpow, nav=nnav, apod_ms=tau, mppca=mp)

mfile = sys.argv[sys.argv.index('--motion') + 1] if '--motion' in sys.argv else None
if mfile is not None:                                       # reuse first pass, motion field, coordinates and prior of an earlier run
    mz = np.load(mfile)
    d.c = mz['c']; dfield, amp = mz['d'], mz['a']; prior_all = mz['prior']
else:
    x1 = d.recon(tr, lam_s=1e-3, lam_t=0.05, it=it1, verbose=False)
    for r in range(nnav):
        c_new, delta, info = NV.refine(d, x1, tr, exclude=hold)
        print(f'  self-navigation round {r+1}: rms coordinate change {info["delta_rms"]:.4f} cycles', flush=True)
        d.c = c_new
        x1 = d.recon(tr, lam_s=1e-3, lam_t=0.05, it=it1, verbose=False)
    dfield, amp = MO.estimate(x1, d.M)
    from scipy import ndimage as ndi                        # prior power per node = smoothed first-pass magnitude (9 mm), clipped, squared
    pw = []
    for b_ in range(d.nb):
        a_ = ndi.gaussian_filter(np.abs(x1[b_]), opt('--psig', 2.5)); ref_ = np.percentile(a_[d.M], 90)   # prior map smoothing, voxels
        q_ = 1 / np.clip(a_ / ref_, 0.1, 1.0) ** 2; pw.append((q_ / q_[d.M].mean()).astype(np.float32))
    prior_all = np.stack(pw)
    np.savez(os.path.join(OD, f'{src}_{tag}_motion.npz'), d=dfield.astype(np.float32), a=amp, c=d.c, prior=prior_all)
    del x1
mats = MO.warp_mats(dfield, amp, d.M, sym='--sym' in sys.argv)
log['node_amplitude'] = np.asarray(amp).tolist(); log['disp_p99_vox'] = float(np.percentile(np.sqrt((dfield ** 2).sum(0))[d.M], 99))
prior = prior_all if '--prior' in sys.argv else None        # Tikhonov weighted by 1 / prior power (one fixed reweighting, never iterated)


def post(x):
    if tau > 0:
        x = PO.apod(x, tau)
    return x


def rec(lines):
    return post(d.recon(lines, lam_s=ls, lam_t=lt, it=it, verbose=False, moco=mats, prior=prior))


def crop(x):
    c = DY.crop(x)
    if mp:
        c = mppca(c, p=mp, stride=2)[0]
    return c


xa = crop(rec(np.intersect1d(tr, kit['half_a_lines']))); xb = crop(rec(np.intersect1d(tr, kit['half_b_lines'])))
m = SC.metrics(xa, xb, kit)
x = d.recon(tr, lam_s=ls, lam_t=lt, it=it, verbose=False, moco=mats, prior=prior)
U = d.fwd_all(x)                                             # held-out scores from the model itself (before post-filters)
m['heldout_lines'] = d.chi2(x, kit['heldout_lines'], U).tolist(); m['heldout_interleaves'] = d.chi2(x, d.ho_ilv_lines, U).tolist(); m['train'] = d.chi2(x, tr[::6], U).tolist()
xc = crop(post(x))
m.update(tag=f'{src}_{tag}', seconds=time.time() - t0); log.update(m)
np.save(os.path.join(OD, f'{src}_{tag}_scored_crop.npy'), xc.astype(np.complex64))
np.save(os.path.join(OD, f'{src}_{tag}_halfA_crop.npy'), xa.astype(np.complex64)); np.save(os.path.join(OD, f'{src}_{tag}_halfB_crop.npy'), xb.astype(np.complex64))
SC.show(m)
if src != 'real':
    pm = PS.score(np.abs(xc), truth_file=os.path.join(PD, f'truth_{src}.npz')); PS.show(f'{src}_{tag}', pm); log['phantom'] = pm
if '--display' in sys.argv:
    xd = crop(post(d.recon(d.ss, lam_s=ls, lam_t=lt, it=it, verbose=False, moco=mats, prior=prior)))
    np.save(os.path.join(OD, f'{src}_{tag}_display_crop.npy'), xd.astype(np.complex64))
os.makedirs(os.path.join(H.OUT, 'scores'), exist_ok=True)
json.dump(log, open(os.path.join(H.OUT, 'scores', f'{src}_{tag}.json'), 'w'), indent=1)
print(f'total {time.time()-t0:.0f} s')

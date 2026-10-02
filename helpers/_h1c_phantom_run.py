"""Scratch (h1 challenge): run the ASAP 4D least-squares pipeline on the phantom and score it against the truth.
usage: _h1c_phantom_run.py <tag> <lam_t> <iterations> [variant=v1]"""
import os, sys, json, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY, _h1c_phantom_score as PS
tag, lam_t, it = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]); var = sys.argv[4] if len(sys.argv) > 4 else 'v1'
PD = os.path.join(H.EXT, 'phantom')
z = np.load(os.path.join(PD, f'cache_phantom_{var}.npz'))
v, info = H.prewhiten(z['raw'])
print('phantom thermal covariance estimate (diag):', np.diag(info['noise_cov']).real, ' (real data: 1.12e-11, 9.4e-12)')
d = DY.Dyn(data=v, vol=z['vol'].astype(np.float64))
d.wpow = float(os.environ.get('WPOW', 1.0))
t0 = time.time(); x = d.recon(d.kit['scored_train_lines'], lam_s=1e-3, lam_t=lam_t, it=it, verbose=False)
xc = DY.crop(x)
os.makedirs(os.path.join(PD, 'asap'), exist_ok=True); np.save(os.path.join(PD, 'asap', f'{tag}_{var}_crop.npy'), xc.astype(np.complex64))
m = PS.score(np.abs(xc), truth_file=os.path.join(PD, f'truth_{var}.npz')); PS.show(f'{tag} ({var})', m)
U = d.fwd_all(x)
print('   held-out lines', np.round(d.chi2(x, d.kit['heldout_lines'], U)[0], 2), ' held-out interleaves', np.round(d.chi2(x, d.ho_ilv_lines, U)[0], 2), ' train', np.round(d.chi2(x, d.kit['scored_train_lines'][::6], U)[0], 2), f' ({time.time()-t0:.0f} s)')
os.makedirs(os.path.join(H.OUT, 'scores'), exist_ok=True); json.dump(m, open(os.path.join(H.OUT, 'scores', f'phantom_{var}_{tag}.json'), 'w'), indent=1)

"""Scratch: which knob costs vessel contrast? Phantom v1, scored set only, truth scores + held-out interleaves.
Motion field from an earlier pipeline run; first pass redone once to build prior maps at two smoothing widths."""
import os, sys, time, numpy as np
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY, _h1c_motion as MO, _h1c_phantom_score as PS
src = 'v1'; PD = os.path.join(H.EXT, 'phantom')
z = np.load(os.path.join(PD, f'cache_phantom_{src}.npz')); d = DY.Dyn(data=H.prewhiten(z['raw'])[0], vol=z['vol'].astype(np.float64))
tr = d.kit['scored_train_lines']
mz = np.load(os.path.join(H.EXT, 'arms', 'v1_sharp_sym_motion.npz')); mats = MO.warp_mats(mz['d'], mz['a'], d.M)
x1 = d.recon(tr, lam_s=1e-3, lam_t=0.05, it=12, verbose=False)
def prior(sig, floor=0.1):
    pw = []
    for b in range(d.nb):
        a_ = ndi.gaussian_filter(np.abs(x1[b]), sig) if sig > 0 else np.abs(x1[b]); ref_ = np.percentile(a_[d.M], 90)
        q_ = 1 / np.clip(a_ / ref_, floor, 1.0) ** 2; pw.append((q_ / q_[d.M].mean()).astype(np.float32))
    return np.stack(pw)
P = {2.5: prior(2.5), 1.0: prior(1.0)}
del x1
def run(tag, lt=0.5, ls=0.01, pr=2.5, wpow=0.75, it=20, moco=True):
    d.wpow = wpow; t0 = time.time()
    x = d.recon(tr, lam_s=ls, lam_t=lt, it=it, verbose=False, moco=mats if moco else None, prior=None if pr is None else P[pr])
    m = PS.score(np.abs(DY.crop(x))); hi = d.chi2(x, d.ho_ilv_lines)[0]
    print(f"{tag:46s} err box {m['err_box']:.3f} err lung {m['err_lung']:.3f} bias {m['lung_bias']:+.3f} vessel {m['vessel_contrast']:.3f} exc {m['dome_excursion_mm']:.2f} | held-out ilv {hi[0]:5.0f} {hi[1]:5.1f} {hi[2]:5.2f} ({time.time()-t0:.0f} s)", flush=True)
run('A reference: lt .5, prior 9 mm, lam .01')
run('B prior 3.5 mm smoothing, lam .01', pr=1.0)
run('C prior 3.5 mm, lam .003', pr=1.0, ls=0.003)
run('D uniform Tikhonov .003', pr=None, ls=0.003)
run('E uniform Tikhonov .0005, 40 iterations', pr=None, ls=0.0005, it=40)
run('F prior 3.5 mm, lam .01, coupling .2', pr=1.0, lt=0.2)
run('G prior 3.5 mm, lam .01, 40 iterations', pr=1.0, it=40)
run('H prior 3.5 mm, lam .01, weights p = 1', pr=1.0, wpow=1.0)

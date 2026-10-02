"""Scratch: motion-compensated temporal coupling arms (kit v2 scores). usage: _h1c_moco_run.py <motion.npz> <lam_t,lam_t,...> [it]"""
import os, sys, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY, _h1c_score as SC, _h1c_motion as MO
mz = np.load(sys.argv[1]); lts = [float(v) for v in sys.argv[2].split(',')]; it = int(sys.argv[3]) if len(sys.argv) > 3 else 15
d = DY.Dyn()
t0 = time.time(); mats, idx = MO.warp_mats(mz['d'], mz['a'], d.M); print('warp matrices %.0f s' % (time.time() - t0), flush=True)
for lt in lts:
    SC.show(SC.run_arm(f'moco_lt{lt:g}_it{it}', lt, 1e-3, it, d=d, moco=(mats, idx)))

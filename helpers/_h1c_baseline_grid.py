"""Scratch: kit items 1-5 for a run-808-equivalent baseline (plain gridding on the 350 mm grid, Steve's soft bins
exp(-d^2/0.44) on the same cyclic coordinate, principal channel, complex), so the old cine has numbers on the same scale."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_score as SC
kit = np.load(SC.KIT2); D = H.load(); k = D['k']; v = H.prewhiten(D['raw'])[0][0]; c = H.cyc_coord(D['vol'])
dcf = np.load(os.path.join(H.EXT, 'dcf_n100.npy'))
tr = kit['scored_train_lines']
def grid(lines):
    op = H.Nufft(k, lines, n=100); W = H.soft_w(c[lines], 16, 0.44); w = dcf[lines % H.NILV][:, H.KILL:]
    return np.stack([H.to_tyger(op.adj(v[lines, H.KILL:] * w * W[:, b][:, None])) for b in range(16)])
a = grid(np.intersect1d(tr, kit['half_a_lines'])); b = grid(np.intersect1d(tr, kit['half_b_lines']))
m = SC.metrics(a, b, kit)
print(f"gridding 350 mm (run-808 equivalent): noise/tis {m['noise_lung_over_tissue']:.4f}  lung/tis {m['lung_over_tissue']:.3f}  tisSNR {m['tissue_snr']:.1f}  bg/tis {m['artefact_bg_over_tissue']:.3f}  dome exc {m['dome_excursion_mm']:.2f} mm  width {m['dome_width_mm']:.2f} mm")
os.makedirs(os.path.join(H.EXT, 'deliver', 'real'), exist_ok=True)
np.save(os.path.join(H.EXT, 'deliver', 'real', 'baseline_gridding350_halfA.npy'), a.astype(np.complex64)); np.save(os.path.join(H.EXT, 'deliver', 'real', 'baseline_gridding350_halfB.npy'), b.astype(np.complex64))
json.dump(m, open(os.path.join(H.OUT, 'scores', 'real_baseline_gridding350.json'), 'w'), indent=1)

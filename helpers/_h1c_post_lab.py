"""Scratch: linear post-filters on a finished arm (half A / half B / scored crops): readout-filter equivalent and
MP-PCA across phases; kit items 1-5 for each. usage: _h1c_post_lab.py <arm prefix in arms/>"""
import os, sys, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_score as SC, _h1c_post as PO
from _h1c_mppca import mppca
A = os.path.join(H.EXT, 'arms'); pre = sys.argv[1]
a, b, x = [np.load(os.path.join(A, f'{pre}_{s}_crop.npy')) for s in ('halfA', 'halfB', 'scored')]
kit = np.load(SC.KIT2)
def show(tag, m):
    print(f"{tag:34s} noise/tis {m['noise_lung_over_tissue']:.4f}  lung/tis {m['lung_over_tissue']:.3f}  tisSNR {m['tissue_snr']:6.1f}  bg/tis {m['artefact_bg_over_tissue']:.3f}  dome exc {m['dome_excursion_mm']:5.2f} mm width {m['dome_width_mm']:5.2f} mm", flush=True)
show('as reconstructed', SC.metrics(a, b, kit))
for tau in (2.0, 1.5, 1.0):
    fa, fb = PO.apod(a, tau), PO.apod(b, tau); show(f'readout filter exp(-t/{tau} ms)', SC.metrics(fa, fb, kit))
    if tau == 1.5: np.save(os.path.join(A, f'{pre}_apod1.5_scored_crop.npy'), PO.apod(x, tau))
t0 = time.time(); ma, sa, ra = mppca(a, 5, 2); mb, sb, rb = mppca(b, 5, 2); print('MP-PCA 2 volumes %.0f s; mean rank kept %.2f of 16' % (time.time() - t0, ra[ra > 0].mean()))
show('MP-PCA 5x5x5', SC.metrics(ma, mb, kit)); np.save(os.path.join(A, f'{pre}_mppca_scored_crop.npy'), mppca(x, 5, 2)[0])
fa, fb = PO.apod(ma, 1.5), PO.apod(mb, 1.5); show('MP-PCA + readout filter 1.5 ms', SC.metrics(fa, fb, kit))

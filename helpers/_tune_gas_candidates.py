"""Scratch: real-data bin-0 gas images of the tuning candidates (prod / prod+deapod / MS200+deapod / MS160+deapod),
common absolute window per subject (prod lung p99), coronal through the lung centroid + an axial slice.
usage: python _tune_gas_candidates.py KEY1 KEY2 -> outputs/tune_gas/candidates_bin0.png"""
import os, sys
import numpy as np
from scipy import ndimage as ndi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'recon'))
from steve_kernel_numpy import steve_rolloff

O = os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon/workspace/outputs/tune_gas')
cands = [('prod', 'prod', 240, None), ('prod + deapod', 'prod', 240, 0.2),
         ('MS200 + deapod', 'MS200', 200, 0.2), ('MS160 + deapod', 'MS160', 160, 0.2)]
keys = sys.argv[1:]
fig, ax = plt.subplots(2 * len(keys), len(cands), figsize=(3.4 * len(cands), 3.3 * 2 * len(keys)))
for i, k in enumerate(keys):
    z = np.load(os.path.join('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tune_gas', k, 'montage_bin0.npz'))
    lung = z['lung']
    cz, cy, cx = (int(round(c)) for c in ndi.center_of_mass(lung))
    vmax = np.percentile(z['prod'][lung], 99)
    for j, (lab, name, ms, kd2) in enumerate(cands):
        img = z[name] / (steve_rolloff(ms, 100, kd2) if kd2 else 1.0)
        m = img[lung].mean() / z['prod'][lung].mean()
        ax[2 * i, j].imshow(img[:, cy, :], cmap='gray', vmin=-0.1 * vmax, vmax=vmax)
        ax[2 * i, j].set_title(f'{k[-5:]} {lab}\nlung mean ×{m:.3f}', fontsize=9)
        ax[2 * i + 1, j].imshow(img[cz, :, :], cmap='gray', vmin=-0.1 * vmax, vmax=vmax)
for a in ax.ravel():
    a.axis('off')
fig.suptitle('real data, bin 0: coronal (rows 1, 3) and axial (rows 2, 4) through the lung centroid; one window per subject (prod lung p99)')
fig.tight_layout()
fig.savefig(os.path.join(O, 'candidates_bin0.png'), dpi=110)
print(os.path.join(O, 'candidates_bin0.png'))

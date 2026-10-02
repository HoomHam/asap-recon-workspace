"""Scratch: evaluate readout windows on bias-corrected lung-parenchyma SNR (reads apod_images.npz).
Lung = low-signal voxels inside the per-coronal-slice filled chest hull; background = outside dilated hull.
Noise floor from background: E[x^2]=2*Nc*sigma^2 (2-coil SOS); S_lung = sqrt(<x^2>_lung - <x^2>_bg)."""
import sys, os, json, numpy as np
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
Z = np.load('/Volumes/HoomHamExt/AIkill_Dynamic_1H/2024-03-12_006KL/apod_images.npz'); names = list(Z.keys())
NC = 2
x0 = Z[names[0]].mean(0); xs = ndi.gaussian_filter(x0, 2.0); p99 = np.percentile(xs, 99); N = x0.shape[0]
hull = np.zeros_like(xs, bool)
for ap in range(N):                                           # axis1 = AP; fill holes in each coronal slice
    hull[:, ap, :] = ndi.binary_fill_holes(xs[:, ap, :] > 0.33 * p99)
hull = ndi.binary_opening(hull, iterations=2)
lung = hull & (xs < 0.28 * p99)
lab, n = ndi.label(lung); sz = ndi.sum(lung, lab, range(1, n + 1)); big = np.argsort(sz)[::-1][:2] + 1
lung = ndi.binary_erosion(np.isin(lab, big), iterations=2)
yy = np.indices(xs.shape) - N / 2; inside = np.sqrt((yy ** 2).sum(0)) < 0.45 * N
bg = inside & (xs < 0.2 * p99) & ~ndi.binary_dilation(hull, iterations=3)
tissue = hull & (xs > 0.5 * p99)
print(f'masks: lung {lung.sum()} vox ({lung.sum()*27/1000:.0f} mL), background {bg.sum()}, tissue {tissue.sum()}')
res = {}
for nm in names:
    xa, xb = Z[nm]; x = 0.5 * (xa + xb)
    d = 0.5 * (xa - xb); sigma = d[tissue].std()              # high-SNR tissue: |x| noise ~ per-component sigma (aliasing cancels)
    floor2 = 2 * NC * sigma ** 2
    s_lung = np.sqrt(max((x[lung] ** 2).mean() - floor2, 0)); s_tis = np.sqrt(max((x[tissue] ** 2).mean() - floor2, 0))
    sig_diff = float((x[bg] ** 2).mean() / (2 * NC)) ** 0.5 if bg.any() else float('nan')   # cross-check from background
    res[nm] = dict(lung_SNR=float(s_lung / sigma), tissue_SNR=float(s_tis / sigma), sigma_bg=float(sigma), sigma_passdiff=float(sig_diff),
                   lung_over_floor=float(np.sqrt((x[lung] ** 2).mean() / floor2)))
    print(f'{nm:22s} lung SNR(bias-corr) {s_lung/sigma:5.2f}  tissue SNR {s_tis/sigma:6.1f}  lung rms/floor {res[nm]["lung_over_floor"]:.2f}  sigma bg {sigma:.3g} vs passdiff {sig_diff:.3g}')
json.dump(res, open(os.path.join(sys.argv[1], 'apod_eval.json'), 'w'), indent=1)
ap = int(np.round(np.argwhere(lung)[:, 1].mean()))
fig, ax = plt.subplots(1, len(names) + 1, figsize=(3.8 * (len(names) + 1), 4.4))
ax[0].imshow(x0[:, ap, :].T, cmap='gray', origin='lower'); ax[0].contour(lung[:, ap, :].T, [.5], colors='c', linewidths=.6)
ax[0].contour(bg[:, ap, :].T, [.5], colors='y', linewidths=.6); ax[0].set_title('masks: lung (c), bg (y)', fontsize=9); ax[0].axis('off')
for a, nm in zip(ax[1:], names):
    x = Z[nm].mean(0); a.imshow(x[:, ap, :].T, cmap='gray', origin='lower', vmin=0, vmax=np.percentile(x[hull], 99))
    a.set_title(f'{nm.replace("_"," ")}\nlung SNR {res[nm]["lung_SNR"]:.2f}', fontsize=9); a.axis('off')
fig.tight_layout(); fig.savefig(os.path.join(sys.argv[1], 'apod_eval.png'), dpi=80); print('ok')

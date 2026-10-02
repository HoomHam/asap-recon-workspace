"""Scratch: readout apodization / truncation vs lung-parenchyma SNR (006KL 1H, end-exp lines, 3 mm grid).
Noise from even- vs odd-pass difference (same k-locations -> aliasing cancels)."""
import sys, os, json, numpy as np
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _proton_local as L
N, DX, DW = 128, 3.0, 5e-6
raw, k, tr = L.load(); lines, _ = L.endexp_lines(raw.shape[2])
A, B = lines[(lines // 832) % 2 == 0], lines[(lines // 832) % 2 == 1]
dcf = L.pipe_dcf(k, N, DX)
t = np.arange(512) * DW                                  # time after TE (s)
kr = np.linalg.norm(k - k[:, :1], axis=-1).mean(0)       # mean |k| vs sample (offset removed)
def tukey(n0): w = np.ones(512); m = np.arange(n0, 512); w[n0:] = np.cos(0.5 * np.pi * (m - n0) / (511 - n0)) ** 2; return w
W = {'full (2.56 ms)': np.ones(512),
     'trunc 1.5 ms': (t <= 1.5e-3).astype(float),
     'trunc 1.0 ms': (t <= 1.0e-3).astype(float),
     'matched exp(-t/2ms)': np.exp(-t / 2e-3),
     'tukey from 1.0 ms': tukey(200)}
img = {}
for name, w in W.items():
    xa = L.recon(raw, k, A, dcf, N, DX, win=w); xb = L.recon(raw, k, B, dcf, N, DX, win=w)
    img[name] = (xa, xb)
    print('done', name, flush=True)
np.savez('/Volumes/HoomHamExt/AIkill_Dynamic_1H/2024-03-12_006KL/apod_images.npz', **{n.replace(' ','_'): np.stack(v).astype(np.float32) for n, v in img.items()})
if os.environ.get('RECON_ONLY'): sys.exit(0)
# masks from the full-readout image
x0 = 0.5 * (img['full (2.56 ms)'][0] + img['full (2.56 ms)'][1]); xs = ndi.gaussian_filter(x0, 1.5)
body = ndi.binary_fill_holes(xs > 0.2 * np.percentile(xs, 99))
body = ndi.binary_erosion(body, iterations=4)
dark = body & (xs < 0.45 * np.median(xs[body]))
lab, n = ndi.label(dark); sizes = ndi.sum(dark, lab, range(1, n + 1)); keep = np.argsort(sizes)[::-1][:2] + 1
lung = ndi.binary_erosion(np.isin(lab, keep), iterations=2)
tissue = body & ~ndi.binary_dilation(np.isin(lab, keep), iterations=3) & (xs > np.median(xs[body]))
print(f'lung mask {lung.sum()} vox, tissue mask {tissue.sum()} vox', flush=True)
res = {}
for name, (xa, xb) in img.items():
    x = 0.5 * (xa + xb); sig = (xa - xb) / 2                 # noise of the mean image
    sn = sig[body].std()
    w = W[name]; keff = np.sqrt((w ** 2 * kr ** 2).sum() / (w ** 2).sum())
    res[name] = dict(lung_SNR=float(x[lung].mean() / sn), tissue_SNR=float(x[tissue].mean() / sn),
                     lung_tissue_ratio=float(x[lung].mean() / x[tissue].mean()), k_rms=float(keff),
                     kmax_used=float(kr[w > 0.01].max()))
    print(f'{name:22s} lung SNR {res[name]["lung_SNR"]:6.1f}  tissue SNR {res[name]["tissue_SNR"]:6.1f}  '
          f'lung/tissue {res[name]["lung_tissue_ratio"]:.3f}  k_rms {keff:.4f}  kmax {res[name]["kmax_used"]:.4f} ({1/(2*res[name]["kmax_used"]):.1f} mm)', flush=True)
out = sys.argv[1]; json.dump(res, open(os.path.join(out, 'apod.json'), 'w'), indent=1)
ap = int(np.round(np.argwhere(lung)[:, 1].mean()))
fig, ax = plt.subplots(1, len(img), figsize=(4.2 * len(img), 4.6))
for a, (name, (xa, xb)) in zip(ax, img.items()):
    x = 0.5 * (xa + xb); sl = x[:, ap, :].T
    a.imshow(sl, cmap='gray', origin='lower', vmin=0, vmax=np.percentile(x[body], 99))
    a.contour(lung[:, ap, :].T, [0.5], colors='c', linewidths=.5)
    a.set_title(f'{name}\nlung SNR {res[name]["lung_SNR"]:.1f}', fontsize=9); a.axis('off')
fig.tight_layout(); fig.savefig(os.path.join(out, 'apod_compare.png'), dpi=80)
np.save('/Volumes/HoomHamExt/AIkill_Dynamic_1H/2024-03-12_006KL/apod_full_vs_matched.npy', np.stack([0.5*sum(img['full (2.56 ms)']), 0.5*sum(img['matched exp(-t/2ms)'])]).astype(np.float32))

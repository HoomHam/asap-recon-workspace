"""Scratch: sharper 1H test recon at a chosen k-scale s and voxel DX (low-pass to grid Nyquist)."""
import sys, numpy as np, mapvbvd, finufft, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os
TRAJ = os.environ.get('TRAJ', 'data/xe/human/traj/fa_spiral_dyn_fancy_v3_20240130_gp.npy')
dat, out, s, DX, N = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
k = np.load(TRAJ).reshape(832, 512, 3) * s
tw = mapvbvd.mapVBVD(dat, quiet=True); tw = tw[-1] if isinstance(tw, list) else tw
tw.image.flagRemoveOS = False; tw.image.squeeze = True
raw = tw.image.unsorted(); raw = raw if raw.ndim == 3 else raw[:, None, :]
nl = raw.shape[2]; kk = k[np.arange(nl) % 832]                  # (nl, 512, 3)
d = raw.transpose(2, 0, 1)                                       # (nl, 512, ch)
kr = np.linalg.norm(kk, axis=-1); keep = kr <= 1 / (2 * DX)      # grid Nyquist
if os.environ.get('DCF', 'pipe') == 'pipe':
    import sigpy.mri as smr                                      # Pipe-Menon on the 832 unique interleaves, tiled
    ku = k.reshape(-1, 3) * N * DX                               # cycles/FOV
    wu = np.abs(smr.pipe_menon_dcf(ku, img_shape=(N, N, N), max_iter=30, show_pbar=False)).reshape(832, 512)
    w = wu[np.arange(nl) % 832][keep]
else:
    w = kr[keep] ** 2
kx = 2 * np.pi * DX * kk[keep]
sos = 0
for c in range(d.shape[2]):
    im = finufft.nufft3d1(kx[:, 0], kx[:, 1], kx[:, 2], (d[..., c][keep] * w).astype(np.complex128), (N, N, N), eps=1e-4)
    sos = sos + np.abs(im) ** 2
im = np.sqrt(sos); np.save(out.replace('.png', '.npy'), im.astype(np.float32))
fig, ax = plt.subplots(1, 3, figsize=(15, 5)); ext = [-N * DX / 2, N * DX / 2] * 2
for a, (sl, nm) in zip(ax, ((im[N // 2], 'X-mid'), (im[:, N // 2], 'Y-mid'), (im[:, :, N // 2], 'Z-mid'))):
    a.imshow(sl.T, cmap='gray', origin='lower', extent=ext, vmax=np.percentile(im, 99.9)); a.set_title(nm); a.set_xlabel('mm')
fig.suptitle(f'{dat.split("/")[-3]}  traj={os.path.basename(TRAJ)[:40]}  s={s}  DX={DX} mm  kept {keep.mean()*100:.0f}%'); fig.tight_layout(); fig.savefig(out, dpi=80)
print('wrote', out, 'kept', keep.mean())

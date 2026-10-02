"""Scratch (2026-10-01): undersampled-frame movie of a 1H dynamic (26 ilv/frame for v3, 20 for v2).

Frame i = lines [i*ILV, (i+1)*ILV) with traj interleaves (line % NILV), adjoint NUFFT,
Pipe-Menon DCF computed per Thomson set (frames cycle through NILV/ILV distinct sets),
2-coil SOS, low-res grid. Saves frames .npy + mp4/gif of coronal slices and the k0 trace.

usage: _proton_usimg_video.py <dat> <traj.npy> <outdir> [ILV=26] [N=64] [DX=5.5]
"""
import os, sys
import numpy as np
import mapvbvd, finufft
import sigpy.mri as smr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
sys.path.insert(0, os.path.dirname(__file__))
from _proton_resp_pneumo import k0_resp

dat, trajf, out = sys.argv[1:4]
ILV = int(sys.argv[4]) if len(sys.argv) > 4 else 26
N = int(sys.argv[5]) if len(sys.argv) > 5 else 64
DX = float(sys.argv[6]) if len(sys.argv) > 6 else 5.5
os.makedirs(out, exist_ok=True)

k = np.load(trajf)
NPTS = 512
NILV = k.shape[0] // NPTS
k = k.reshape(NILV, NPTS, 3)
NSET = NILV // ILV
kmax = 1 / (2 * DX)

tw = mapvbvd.mapVBVD(dat, quiet=True)
tw = tw[-1] if isinstance(tw, list) else tw
tw.image.flagRemoveOS = False
tw.image.squeeze = True
tr = float(tw.hdr.MeasYaps[('alTR', '0')]) * 1e-6
raw = tw.image.unsorted()
raw = raw if raw.ndim == 3 else raw[:, None, :]                 # (512, ch, lines)
nfr = raw.shape[2] // ILV
print(f'{nfr} frames x {ILV} ilv ({ILV*tr:.3f} s/frame), {NSET} Thomson sets, grid {N}^3 x {DX} mm', flush=True)

# per-set geometry + DCF
sets = []
for s in range(NSET):
    ks = k[s * ILV:(s + 1) * ILV]                                # (ILV, 512, 3)
    keep = np.linalg.norm(ks, axis=-1) <= kmax
    kk = ks[keep]
    w = np.abs(smr.pipe_menon_dcf(kk * N * DX, img_shape=(N, N, N), max_iter=30, show_pbar=False))
    sets.append((keep, 2 * np.pi * DX * kk, w))

frames = np.zeros((nfr, N, N, N), np.float32)
for i in range(nfr):
    s = i % NSET                                                 # frame i plays interleaves i*ILV.. (mod NILV)
    keep, kx, w = sets[s]
    d = raw[:, :, i * ILV:(i + 1) * ILV].transpose(2, 0, 1)      # (ILV, 512, ch)
    sos = 0
    for c in range(d.shape[2]):
        im = finufft.nufft3d1(kx[:, 0], kx[:, 1], kx[:, 2], (d[..., c][keep] * w).astype(np.complex128), (N, N, N), eps=1e-3)
        sos = sos + np.abs(im) ** 2
    frames[i] = np.sqrt(sos)
fdir = os.environ.get('FRAMES_DIR', out)                       # frames are ~0.5 GB -> keep on Ext
os.makedirs(fdir, exist_ok=True)
np.save(os.path.join(fdir, f'usimg{ILV}_frames.npy'), frames)

# coronal = (axis0 horizontal, axis2 vertical) at fixed axis1 index; pick 3 AP slices around lung centre
mean = frames.mean(0)
body = mean > np.percentile(mean, 70)
ap_c = int(np.round(np.argwhere(body)[:, 1].mean()))
aps = [ap_c - 3, ap_c, ap_c + 3]
t, r, _ = k0_resp(dat)
r = -r                                                           # volume-like (plateaus = end-exp)
tf = (np.arange(nfr) + 0.5) * ILV * tr
vmax = np.percentile(frames[:, :, aps, :], 99.5)

fig = plt.figure(figsize=(13, 6.2))
gs = fig.add_gridspec(2, 3, height_ratios=[4, 1.3])
axs = [fig.add_subplot(gs[0, j]) for j in range(3)]
axt = fig.add_subplot(gs[1, :])
ims = [a.imshow(frames[0][:, ap, :].T, cmap='gray', origin='lower', vmin=0, vmax=vmax) for a, ap in zip(axs, aps)]
for a, ap in zip(axs, aps):
    a.set_title(f'coronal, AP slice {ap} ({(ap - N // 2) * DX:+.0f} mm)'); a.axis('off')
axt.plot(t, r, lw=.6, color='C0'); axt.set_xlim(0, t[-1]); axt.set_ylabel('k0 resp'); axt.set_xlabel('s')
cur = axt.axvline(0, color='r')
ttl = fig.suptitle('')
fig.tight_layout()


def upd(i):
    for im, ap in zip(ims, aps):
        im.set_data(frames[i][:, ap, :].T)
    cur.set_xdata([tf[i], tf[i]])
    ttl.set_text(f'{dat.split("/")[-3]} 1H, {ILV}-interleave frames: frame {i}/{nfr}  t = {tf[i]:.1f} s')
    return ims + [cur, ttl]


anim = FuncAnimation(fig, upd, frames=nfr, blit=False)
mp4 = os.path.join(out, f'usimg{ILV}_coronal.mp4')
anim.save(mp4, writer=FFMpegWriter(fps=6, bitrate=3000))
print('wrote', mp4, flush=True)

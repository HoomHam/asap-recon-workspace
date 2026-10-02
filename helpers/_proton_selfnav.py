"""Scratch (2026-10-01): image-based diaphragm self-navigation for 1H dynamics.

The 1H ASAP readout winds from the start (no radial segment -> no 1D projections), but
the low-k ball is densely sampled: ILV interleaves reach Nyquist out to |k| ~ KCUT.
Per ILV-interleave frame: keep only |k| <= KCUT, Pipe DCF per Thomson set, coarse 3D
adjoint NUFFT (N^3 x DX), 2-coil SOS. Dome = top temporal-std voxels; SI shift of the
ROI-averaged SI profile vs the mean profile (upsampled xcorr) -> dz(t).
Compares dz with the k0 surrogate; saves dz trace (.npz), figure, coronal mp4.

usage: _proton_selfnav.py <dat> <traj.npy> <outdir> [ILV=26] [KCUT=0.019] [N=32] [DX=11]
"""
import os, sys
import numpy as np
import mapvbvd, finufft
import sigpy.mri as smr
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
sys.path.insert(0, os.path.dirname(__file__))
from _proton_resp_pneumo import k0_resp

dat, trajf, out = sys.argv[1:4]
ILV = int(sys.argv[4]) if len(sys.argv) > 4 else 26
KCUT = float(sys.argv[5]) if len(sys.argv) > 5 else 0.019
N = int(sys.argv[6]) if len(sys.argv) > 6 else 32
DX = float(sys.argv[7]) if len(sys.argv) > 7 else 11.0
SI = 2                                              # traj/image axis 2 = superior-inferior
os.makedirs(out, exist_ok=True)
tag = f'ilv{ILV}_k{KCUT:g}'

k = np.load(trajf)
NILV = k.shape[0] // 512
k = k.reshape(NILV, 512, 3)
NSET = NILV // ILV

tw = mapvbvd.mapVBVD(dat, quiet=True)
tw = tw[-1] if isinstance(tw, list) else tw
tw.image.flagRemoveOS = False
tw.image.squeeze = True
tr = float(tw.hdr.MeasYaps[('alTR', '0')]) * 1e-6
raw = tw.image.unsorted()
raw = raw if raw.ndim == 3 else raw[:, None, :]
nfr = raw.shape[2] // ILV

sets = []
for s in range(NSET):
    ks = k[s * ILV:(s + 1) * ILV]
    keep = np.linalg.norm(ks, axis=-1) <= KCUT
    keep[:, :2] = False                             # samples 0+1 are ADC transients
    kk = ks[keep]
    w = np.abs(smr.pipe_menon_dcf(kk * N * DX, img_shape=(N, N, N), max_iter=40, show_pbar=False))
    sets.append((keep, 2 * np.pi * DX * kk, w))
print(f'{nfr} frames x {ILV} ilv ({ILV*tr*1e3:.0f} ms), |k|<={KCUT} ({1/(2*KCUT):.0f} mm res), '
      f'{sets[0][0].sum()} samples/frame, grid {N}^3 x {DX} mm', flush=True)

fr = np.zeros((nfr, N, N, N), np.float32)
for i in range(nfr):
    keep, kx, w = sets[i % NSET]
    d = raw[:, :, i * ILV:(i + 1) * ILV].transpose(2, 0, 1)
    sos = 0
    for c in range(d.shape[2]):
        im = finufft.nufft3d1(kx[:, 0], kx[:, 1], kx[:, 2], (d[..., c][keep] * w).astype(np.complex128), (N, N, N), eps=1e-4)
        sos = sos + np.abs(im) ** 2
    fr[i] = np.sqrt(sos)
skip = int(np.ceil(NILV / ILV))                     # first pass = approach to steady state
F = fr[skip:]
tf = (np.arange(nfr) + 0.5) * ILV * tr
tF = tf[skip:]
fdir = os.environ.get('FRAMES_DIR', out)                           # keep frame stacks off git/laptop
os.makedirs(fdir, exist_ok=True)
np.save(os.path.join(fdir, f'selfnav_{tag}_frames.npy'), fr)

# dome ROI: voxels with the largest temporal std (after removing per-set mean -> kills set-to-set PSF wobble)
Fc = F.copy()
for s in range(NSET):
    idx = (np.arange(skip, nfr) % NSET) == s
    Fc[idx] -= F[idx].mean(0)
sd = Fc.std(0)
mean = F.mean(0)
thr = np.percentile(sd, 99)
roi = sd >= thr
c = np.argwhere(roi).mean(0).round().astype(int)
half = 3
sl = [slice(max(ci - half, 0), ci + half + 1) for ci in c]
sl[SI] = slice(None)
other = tuple(a for a in range(3) if a != SI)
prof = np.array([f[tuple(sl)].mean(axis=other) for f in F])          # (frames, N)
UP = 20
x = np.arange(N)
xf = np.linspace(0, N - 1, N * UP)
P = np.array([np.interp(xf, x, p) for p in prof])
ref = P.mean(0)
rw = gaussian_filter1d(np.gradient(ref), UP)                         # weight the edge
dz = []
for p in P:
    g = gaussian_filter1d(np.gradient(p), UP)
    cc = np.correlate(g, rw, 'full')
    m = np.argmax(cc[len(rw) - 1 - 4 * UP:len(rw) + 4 * UP]) - 4 * UP
    dz.append(m / UP * DX)
dz = np.array(dz)
dz -= np.median(dz)

t, r, _ = k0_resp(dat)
v = np.interp(tF, t, -r)                                             # k0 volume surrogate at frame times
v = (v - v.mean()) / v.std()
dzn = (dz - dz.mean()) / dz.std()
lags = np.arange(-4, 5)
cor = [np.corrcoef(np.roll(dzn, L), v)[0, 1] for L in lags]
bl = lags[int(np.argmax(np.abs(cor)))]
print(f'dome ROI centre {c.tolist()}, dz p2p (5-95%) {np.percentile(dz,95)-np.percentile(dz,5):.1f} mm; '
      f'corr(dz, k0 volume) = {cor[4]:+.2f} at lag 0, best {cor[bl+4]:+.2f} at lag {bl} frames', flush=True)
np.savez(os.path.join(out, f'selfnav_{tag}.npz'), t=tF, dz=dz, k0vol=v, roi_centre=c, prof=prof)

fig, ax = plt.subplots(3, 1, figsize=(14, 9), gridspec_kw=dict(height_ratios=[2, 1, 1]))
pr = (prof - prof.mean(1, keepdims=True)) / prof.std(1, keepdims=True)
ax[0].imshow(pr.T, aspect='auto', origin='lower', cmap='gray', extent=[tF[0], tF[-1], 0, N * DX])
ax[0].set_ylabel('SI (mm)'); ax[0].set_title(f'ROI SI profile vs time ({ILV} ilv = {ILV*tr*1e3:.0f} ms/frame, |k|<={KCUT})')
ax[1].plot(tF, dz, lw=.8); ax[1].set_ylabel('dome dz (mm)')
ax[2].plot(tF, dzn, lw=.8, label='dome dz (z)'); ax[2].plot(tF, v, lw=.8, label='k0 volume (z)')
ax[2].legend(loc='upper right'); ax[2].set_xlabel('s'); ax[2].set_title(f'corr = {cor[4]:+.2f}')
for a in ax[1:]:
    a.set_xlim(tF[0], tF[-1])
fig.tight_layout()
fig.savefig(os.path.join(out, f'selfnav_{tag}.png'), dpi=75)

# coronal movie through the ROI centre
ap = c[1]
vmax = np.percentile(F[:, :, ap, :], 99.5)
fig2, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw=dict(width_ratios=[1, 1.6]))
im = a1.imshow(F[0][:, ap, :].T, cmap='gray', origin='lower', vmin=0, vmax=vmax); a1.axis('off')
a1.plot([c[0]], [c[2]], 'r+', ms=12)
a2.plot(tF, dz, lw=.6); cur = a2.axvline(tF[0], color='r'); a2.set_xlabel('s'); a2.set_ylabel('dome dz (mm)')
ttl = fig2.suptitle('')
fig2.tight_layout()


def upd(i):
    im.set_data(F[i][:, ap, :].T)
    cur.set_xdata([tF[i], tF[i]])
    ttl.set_text(f'low-k navigator {ILV} ilv ({ILV*tr*1e3:.0f} ms)  t={tF[i]:.1f} s')
    return [im, cur, ttl]


FuncAnimation(fig2, upd, frames=len(F)).save(os.path.join(out, f'selfnav_{tag}_coronal.mp4'), writer=FFMpegWriter(fps=6, bitrate=2000))
print('wrote', out, flush=True)

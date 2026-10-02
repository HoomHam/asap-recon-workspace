"""Scratch: does the Xe v3_20240130 trajectory fit the 1H dynamic scans as-is,
or does k need the gamma ratio (42.577/11.777 = 3.615)?

Time-averaged adjoint NUFFT of the whole dynamic run, two k scalings:
  A  traj as-is                       (sequence rescales gradients by gamma)
  B  traj x 3.615, |k| <= kmax(A)     (fixed gradient waveform; low-passed to same res)
The right one shows an anatomically sized torso; the wrong one is 3.6x off.
"""
import sys
import numpy as np
import mapvbvd
import finufft
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TRAJ = '/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/data/xe/human/traj/fa_spiral_dyn_fancy_v3_20240130_gp.npy'
GRATIO = 42.577 / 11.777
N, DX = 128, 4.0                     # grid 128 x 4 mm = 512 mm FOV


def npts_from_traj(k):
    a = np.abs(np.fft.fft(np.sum(k**2, 1)))
    idx = next(i for i in range(1, len(a)) if a[i] > a[0] / 10)
    return int(len(k) / idx + 0.1)


def recon(dat, out_png, title, skip=0, stride=1, scales=(1.0, GRATIO)):
    k = np.load(TRAJ)                                   # cycles/mm
    npts = npts_from_traj(k)
    nilv_traj = k.shape[0] // npts
    k = k.reshape(nilv_traj, npts, 3)

    tw = mapvbvd.mapVBVD(dat)
    tw = tw[-1] if isinstance(tw, list) else tw
    tw.image.flagRemoveOS = False
    tw.image.squeeze = True
    raw = tw.image.unsorted().astype('complex64')       # (samples, ch, lines)
    if raw.ndim == 2:
        raw = raw[:, None, :]
    ns, nch, nl = raw.shape
    print(f'{dat.split("/")[-1]}: raw {raw.shape}, traj npts {npts} x {nilv_traj} ilv', flush=True)
    m = min(ns, npts)
    raw = raw[:, :, skip::stride]                       # drop spectra lines, pick gas lines
    nl = raw.shape[2]
    lines = np.arange(nl)
    kk = k[lines % nilv_traj, :m, :]                    # (nl, m, 3)
    d = raw[:m].transpose(2, 0, 1)                      # (nl, m, ch)

    imgs = {}
    for s in scales:
        tag = f's={s:.3g}'
        ks = kk * s
        kr = np.linalg.norm(ks, axis=-1)
        keep = kr <= np.linalg.norm(kk, axis=-1).max()        # s>1: low-pass to the unscaled kmax
        kx = (2 * np.pi * DX * ks[keep]).astype(np.float64)
        w = np.minimum(kr[keep], np.percentile(kr[keep], 99)) ** 2   # 3D center-out DCF ~ |k|^2
        sos = np.zeros((N, N, N))
        for c in range(nch):
            c_ = (d[..., c][keep] * w).astype(np.complex128)
            im = finufft.nufft3d1(kx[:, 0], kx[:, 1], kx[:, 2], c_, (N, N, N), eps=1e-4)
            sos += np.abs(im) ** 2
        imgs[tag] = np.sqrt(sos)
        print(f'  {tag}: kept {keep.mean()*100:.0f}% samples', flush=True)

    fig, ax = plt.subplots(len(imgs), 3, figsize=(12, 4 * len(imgs)), squeeze=False)
    for r, (tag, im) in enumerate(imgs.items()):
        for c, (sl, name) in enumerate(((im[N // 2], 'X-mid'), (im[:, N // 2], 'Y-mid'), (im[:, :, N // 2], 'Z-mid'))):
            ax[r, c].imshow(np.abs(sl).T, cmap='gray', origin='lower', extent=[-N*DX/2, N*DX/2]*2,
                            vmax=np.percentile(im, 99.9))
            ax[r, c].set_title(f'{tag}  {name}')
            ax[r, c].set_xlabel('mm')
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=90)
    print('wrote', out_png, flush=True)


if __name__ == '__main__':
    a = sys.argv
    recon(a[1], a[2], a[3], skip=int(a[4]) if len(a) > 4 else 0, stride=int(a[5]) if len(a) > 5 else 1,
          scales=tuple(float(x) for x in a[6].split(',')) if len(a) > 6 else (1.0, GRATIO))

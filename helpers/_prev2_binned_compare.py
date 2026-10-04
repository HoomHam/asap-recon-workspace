"""Scratch: breathing-binned local recon of a 20220812 acqOrder-7 twix (entry-4 = v2 gp trajectory, gas-counter
arms) vs Faraz's own 16-bin recon of the same file — fair (binned vs binned) image-quality comparison.

Phase: virtual-coil gas k0 -> smoothed, detrended -> Hilbert phase; bin 0 = end-expiration (k0 minimum), as
Faraz's bins (min signal at bin 0). Each recon bin takes scans within +-1/16 cycle of the bin centre (2/16 wide).
Orientation: brute-force 48 axis permutations/flips of our grid onto Faraz's (correlation of time means).

Usage: helpers/.venv/bin/python helpers/_prev2_binned_compare.py <twix> <faraz_img_dyn.mat> <tag> <smooth_scans> <detrend_scans>
Writes: outputs/prev2_calib_2026-10-03/binned_<tag>.png
"""
import itertools
import pathlib
import sys

import numpy as np
import scipy.io as sio
from scipy.ndimage import uniform_filter1d, zoom
from scipy.signal import hilbert

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from _prev2_static_recon import read_scans, recon, TRAJ, OUT  # noqa: E402


def best_orient(ours, ref):
    a = zoom(ours, 0.5, order=1)
    b = zoom(ref, 0.5, order=1)
    b = (b - b.mean()) / b.std()
    best = None
    for perm in itertools.permutations(range(3)):
        t = np.transpose(a, perm)
        for fl in itertools.product((False, True), repeat=3):
            u = t[tuple(slice(None, None, -1) if f else slice(None) for f in fl)]
            u = (u - u.mean()) / u.std()
            c = float((u * b).mean())
            if best is None or c > best[0]:
                best = (c, perm, fl)
    return best


def apply_orient(v, perm, fl):
    t = np.transpose(v, perm)
    return t[tuple(slice(None, None, -1) if f else slice(None) for f in fl)]


def main():
    twix, fz_mat, tag, sm, dtw = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    d, ts, _, _ = read_scans(twix)
    g = d[0::3]
    del d
    traj = np.load(TRAJ).reshape(-1, 3)
    narms = traj.shape[0] // g.shape[2]
    arm = np.arange(g.shape[0]) % narms
    u, s, vh = np.linalg.svd(g[:, :, 1:40].transpose(1, 0, 2).reshape(g.shape[1], -1), full_matrices=False)
    k0 = np.abs(np.einsum('c,ncs->ns', np.conj(u[:, 0]), g[:, :, 1:6]).mean(1))
    ks = uniform_filter1d(k0, sm)
    x = ks - uniform_filter1d(ks, dtw)          # detrend window > 1 breath
    ph = np.angle(hilbert(x))                       # 0 at signal max (EI), +-pi at min (EE)
    cyc = ((ph + np.pi) / (2 * np.pi)) % 1.0        # 0 = EE ... 0.5 = EI
    fz = sio.loadmat(fz_mat, variable_names=['img_gp'])['img_gp']
    fsig = fz.reshape(-1, fz.shape[-1]).sum(0)
    fEE, fEI = int(np.argmin(fsig)), int(np.argmax(fsig))
    res = {}
    for name, c0 in [('EE', 0.0), ('EI', 0.5)]:
        dist = np.abs(((cyc - c0 + 0.5) % 1.0) - 0.5)
        sel = dist <= 1 / 16
        kd = np.zeros((narms, g.shape[1], g.shape[2]), np.complex128)
        for i in range(narms):
            m = sel & (arm == i)
            if m.any():
                kd[i] = g[m].mean(0)
        res[name] = recon(kd, traj)
        print(f'{tag} {name}: {int(sel.sum())} scans, arms {len(np.unique(arm[sel]))}/{narms}', flush=True)
    c, perm, fl = best_orient(res['EE'] + res['EI'], fz.mean(-1))
    print(f'{tag}: orientation corr {c:.3f} perm {perm} flips {fl}; Faraz EE bin {fEE}, EI bin {fEI}')
    rows = [(f'ours EE', apply_orient(res['EE'], perm, fl)), (f'Faraz bin {fEE} (EE)', fz[..., fEE]),
            (f'ours EI', apply_orient(res['EI'], perm, fl)), (f'Faraz bin {fEI} (EI)', fz[..., fEI])]
    ks_ = list(range(56, 23, -8))                   # Faraz axis 2 = A-P, high = anterior
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(len(rows), len(ks_), figsize=(1.9 * len(ks_), 2.0 * len(rows)), facecolor='k')
    for r, (lab, v) in enumerate(rows):
        vm = np.percentile(v, 99.7)
        for j, k in enumerate(ks_):
            ax[r, j].imshow(v[6:77, 2:66, k].T[::-1], cmap='gray', vmin=0, vmax=vm)
            ax[r, j].axis('off')
            ax[r, j].set_title(f'{lab}  sl{k}' if j == 0 else f'sl{k}', color='w', fontsize=7, loc='left')
    fig.suptitle(f'{tag}: breathing-binned, entry-4 traj (ours, 2/16-wide bins) vs Faraz 16-bin; coronal A→P',
                 color='w', fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / f'binned_{tag}.png', dpi=90, facecolor='k')
    print('wrote', OUT / f'binned_{tag}.png')


if __name__ == '__main__':
    main()

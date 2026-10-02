"""Scratch (h1 challenge): joint figure and cine, both sides, same slices, each row its own window.

Rows (real): run 808, ASAP sharp, ASAP SNR, XeCS sharp, XeCS SNR (display cines, all steady-state lines).
Rows (phantom): truth, then each side's scored reconstruction.
Reads <side>/deliver/<src>/<side>_<arm>_display.npy (or _scored.npy on the phantom); missing rows are skipped.
usage: _h1c_joint.py <real|v1|v1b> [--win 0.45]
"""
import os, sys, subprocess
import numpy as np, scipy.io as sio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H

src = sys.argv[1]; win = float(sys.argv[sys.argv.index('--win') + 1]) if '--win' in sys.argv else 0.45
DA = os.path.join(H.EXT, 'deliver', src); DX_ = os.path.join('/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/h1challenge/deliver', src)
part = 'display' if src == 'real' else 'scored'
cand = []
if src == 'real':
    cand.append(('run 808 (gridding, 350 mm, as delivered 2026-10-01)', H.TYGER_P))
else:
    cand.append(('TRUTH (limited to the acquired k-sphere)', os.path.join(H.EXT, 'phantom', f'truth_{src}.npz')))
for side, d in (('ASAP', DA), ('XeCS', DX_)):
    for arm, nm in (('sharp', 'sharp'), ('snr', 'SNR')):
        f = os.path.join(d, f'{side.lower()}_{arm}_{part}.npy')
        if not os.path.exists(f) and part == 'scored':
            f = os.path.join(d, f'{side.lower()}_{arm}_display.npy')
        cand.append((f'{side} {nm}', f))
rows = [(n, f) for n, f in cand if os.path.exists(f)]
print('rows:', [n for n, _ in rows], '| missing:', [n for n, f in cand if not os.path.exists(f)])


def load(f):
    if f.endswith('.mat'): return np.abs(sio.loadmat(f)['gas_phase']).astype(np.float32)
    if f.endswith('.npz'): return np.abs(np.load(f)['truth_ksphere']).astype(np.float32)
    return np.abs(np.load(f)).astype(np.float32)


vols = [load(f) for _, f in rows]
views = (((slice(20, 95), 48, slice(5, 95)), 'coronal AP 48'), ((slice(20, 95), 58, slice(5, 95)), 'coronal AP 58'), ((slice(20, 95), slice(25, 90), 30), 'sagittal LR 30'), ((45, slice(20, 90), slice(5, 95)), 'axial SI 45'))
fig, ax = plt.subplots(len(rows), 2 * len(views), figsize=(4.6 * 2 * len(views), 4.4 * len(rows)), squeeze=False)
for r, ((nm, _), v) in enumerate(zip(rows, vols)):
    vm = win * np.percentile(v[0], 99.5)
    for c, (sl, t) in enumerate(views):
        for j, ph in enumerate((0, 8)):
            a_ = ax[r, 2 * c + j]; a_.imshow(v[ph][sl], cmap='gray', vmin=0, vmax=vm, interpolation='nearest'); a_.axis('off')
            a_.set_title(f'{nm}\n{t}, phase {ph} ({"end expiration" if ph == 0 else "end inspiration"})', fontsize=9)
fig.tight_layout(); out = os.path.join(H.OUT, f'joint_{src}.png'); fig.savefig(out, dpi=50); print('wrote', out)
cine = os.path.join(H.OUT, f'joint_{src}.mp4')
subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), '_h1c_cine.py'), cine] + [f'{n}={f}' for n, f in rows] + ['--win', str(win)], check=True)

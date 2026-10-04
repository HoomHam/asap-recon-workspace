"""Scratch: 2022-10-26_000HH slow vs fast two-row coronal video from OUR Tyger recons (fork ee3c91f/6b01071,
DIAPHRAGM, Faraz calib entry 4 = v2 gp, pipeline/convert_mixed_adc.py), same layout as _hh_coronal_video.py
(Faraz recon version): 10 coronal slices per row, anterior (left) -> posterior (right), rows matched in A-P.

Orientation: each Tyger gas_phase (16,100,100,100) is mapped onto Faraz's img_gp layout with the 48-way
permutation/flip that best correlates with Faraz's own recon of the same run (helpers/_prev2_binned_compare);
then Faraz's checked convention applies: axis0 L-R (high = patient left), axis1 S-I (high = superior -> flipud),
axis2 A-P (high = anterior).
Slice matching: lung A-P extent per row (time-mean mask) -> 10 positions at the same fractional depths.
Phase: each row rolled so its min-signal bin (end-exp) is frame 0.

Usage: helpers/.venv/bin/python helpers/_hh_coronal_video_tyger.py
Writes: outputs/hh_fastslow_2026-10-03/HH_2022-10-26_coronal_slow_vs_fast_TYGER.{mp4,gif} + _TYGER_still_EI.png
"""
import pathlib
import subprocess
import sys

import numpy as np
import scipy.io as sio
from scipy.ndimage import gaussian_filter, zoom

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from _prev2_binned_compare import best_orient, apply_orient  # noqa: E402

T = pathlib.Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_2022')
F = pathlib.Path('/Volumes/HoomHamExt/Work/Images/MRI/Human/2022-10-26_000HH')
OUT = pathlib.Path(__file__).resolve().parents[1] / 'outputs' / 'hh_fastslow_2026-10-03'
RUNS = [('SLOW / deep breathing  ~6.2 bpm (9.7 s)  (MID00018)', '2022-10-26_000HHslow', 'img_dyn_16ph_lowbpm.mat'),
        ('FAST breathing  ~27 bpm (2.2 s)  (MID00017, coil 1 dropped)', '2022-10-26_000HHfast_drop1', 'img_dyn_16ph_highbpm.mat')]
# fast: 7 coils beats 8 (median gas SNR 9.1 vs 8.4); slow: 8 coils kept (16.0 vs 15.9 with coil 1 dropped)
NSL, FPS, CYCLES = 10, 8, 4


def load(key, fz_name):
    t = np.abs(sio.loadmat(T / key / 'd' / 'recon.mat', variable_names=['gas_phase'])['gas_phase']).astype(np.float32)
    fz = sio.loadmat(F / fz_name, variable_names=['img_gp'])['img_gp'].mean(-1)
    c, perm, fl = best_orient(zoom(t.mean(0), 0.8, order=1), fz)
    v = np.stack([apply_orient(t[b], perm, fl) for b in range(t.shape[0])], -1)     # (x, y, z, bin), Faraz layout
    v = np.roll(v, -int(np.argmin(v.reshape(-1, v.shape[-1]).sum(0))), axis=-1)
    print(f'{key}: orient corr {c:.3f} perm {perm} flips {fl}', flush=True)
    return v


def ap_slices(v):
    m = gaussian_filter(v.mean(-1), 1.5)
    bg = np.median(m)
    mask = m > bg + 0.35 * (np.percentile(m, 99.5) - bg)
    prof = np.nonzero(mask.sum((0, 1)) > 20)[0]
    lo, hi = prof.min(), prof.max()
    frac = np.linspace(0.92, 0.08, NSL)                      # anterior (high index) -> posterior
    return [int(round(lo + f * (hi - lo))) for f in frac], (lo, hi), mask


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = []
    for lab, key, fzn in RUNS:
        v = load(key, fzn)
        ks, ext, mask = ap_slices(v)
        ys = np.nonzero(mask.sum((0, 2)) > 15)[0]
        xs = np.nonzero(mask.sum((1, 2)) > 15)[0]
        rows.append(dict(lab=lab, v=v, ks=ks, ext=ext, box=(xs.min(), xs.max(), ys.min(), ys.max())))
        print(f'{key}: A-P extent {ext}, slices {ks}')
    x0 = max(min(r['box'][0] for r in rows) - 4, 0)
    x1 = min(max(r['box'][1] for r in rows) + 5, rows[0]['v'].shape[0])
    y0 = max(min(r['box'][2] for r in rows) - 4, 0)
    y1 = min(max(r['box'][3] for r in rows) + 5, rows[0]['v'].shape[1])
    for r in rows:
        sel = np.stack([r['v'][..., k, :] for k in r['ks']])
        r['vmax'] = np.percentile(sel, 99.8)

    def panel(v, k, b):
        return v[x0:x1, y0:y1, k, b].T[::-1]              # rows = S-I (superior top), cols = L-R (radiological)

    fig, ax = plt.subplots(2, NSL, figsize=(NSL * 1.55, 2 * 1.9 + 0.7), facecolor='k')
    fig.subplots_adjust(left=0.005, right=0.995, top=0.86, bottom=0.02, wspace=0.03, hspace=0.22)
    ims = []
    for i, r in enumerate(rows):
        for j, k in enumerate(r['ks']):
            a = ax[i, j]
            ims.append(a.imshow(panel(r['v'], k, 0), cmap='gray', vmin=0, vmax=r['vmax'], interpolation='bilinear'))
            a.set_xticks([]); a.set_yticks([])
            for s in a.spines.values():
                s.set_visible(False)
            a.set_title(f'sl {k}', color='0.75', fontsize=7, pad=2)
        ax[i, 0].text(0.0, 1.17, r['lab'], transform=ax[i, 0].transAxes, color='w', fontsize=9, ha='left', va='bottom')
    fig.text(0.005, 0.975, '000HH 2022-10-26  gas, OUR Tyger recon (fork ee3c91f, DIAPHRAGM, calib entry 4 = v2; fast: coil 1 dropped), '
             '16 bins/cycle (phase-normalised)   coronal: ANTERIOR (left) → POSTERIOR (right); R on viewer left; '
             'rows matched by fractional A-P depth', color='0.85', fontsize=7, va='top')
    btxt = fig.text(0.995, 0.975, '', color='y', fontsize=9, ha='right', va='top')

    def draw(b):
        i = 0
        for r in rows:
            for k in r['ks']:
                ims[i].set_data(panel(r['v'], k, b))
                i += 1
        btxt.set_text(f'bin {b + 1:2d}/16' + ('  end-exp' if b == 0 else ''))

    tmp = OUT / '_frames_tyger'
    tmp.mkdir(parents=True, exist_ok=True)
    for b in range(16):
        draw(b)
        fig.savefig(tmp / f'f{b:02d}.png', dpi=130, facecolor='k')
    draw(int(np.argmax(rows[0]['v'].reshape(-1, 16).sum(0))))
    fig.savefig(OUT / '_TYGER_still_EI.png', dpi=130, facecolor='k')
    plt.close(fig)
    lst = tmp / 'list.txt'
    lst.write_text(''.join(f"file 'f{b:02d}.png'\nduration {1 / FPS}\n" for _ in range(CYCLES) for b in range(16))
                   + "file 'f00.png'\n")
    stem = OUT / 'HH_2022-10-26_coronal_slow_vs_fast_TYGER'
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-f', 'concat', '-safe', '0', '-i', str(lst),
                    '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-pix_fmt', 'yuv420p', '-c:v', 'libx264',
                    '-r', str(FPS), f'{stem}.mp4'], check=True)
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(FPS), '-i', str(tmp / 'f%02d.png'),
                    '-vf', 'split[a][b];[a]palettegen[p];[b][p]paletteuse', '-loop', '0', f'{stem}.gif'], check=True)
    import shutil
    shutil.rmtree(tmp)
    print('wrote', f'{stem}.mp4', f'{stem}.gif', OUT / '_TYGER_still_EI.png')


if __name__ == '__main__':
    main()

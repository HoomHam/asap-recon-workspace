"""Scratch: 2022-10-26_000HH Faraz dynamic gas recon, slow (deep, MID00018) vs fast (high-bpm, MID00017).

Two-row coronal video, 10 slices each, anterior (left) -> posterior (right), slices matched in A-P.

Faraz img_gp layout (80,80,80,16), checked by anatomy 2026-10-03:
  axis 0 = L-R  (high = patient LEFT -> raw columns already radiological, no flip)
  axis 1 = S-I  (high = SUPERIOR: trachea/carina there, diaphragm edge moves at low index -> flipud)
  axis 2 = A-P  (high = ANTERIOR: V-shaped lungs around heart; low = two columns beside spine)
  axis 3 = 16 respiratory bins, one cycle, starting near end-expiration.
Matching: lung A-P extent low 22-63 / high 21-62, centroid offset 0.75 -> fast slice = slow slice - 1.
Phase: both rows rolled so their min-signal bin (end-exp) is frame 0 (slow: 0, fast: 1).

Usage: helpers/.venv/bin/python helpers/_hh_coronal_video.py
Writes: outputs/hh_fastslow_2026-10-03/HH_2022-10-26_coronal_slow_vs_fast.{mp4,gif} + _still_EI.png
"""
import pathlib
import subprocess

import numpy as np
import scipy.io as sio

SRC = pathlib.Path('/Volumes/HoomHamExt/Work/Images/MRI/Human/2022-10-26_000HH')
OUT = pathlib.Path(__file__).resolve().parents[1] / 'outputs' / 'hh_fastslow_2026-10-03'
SLOW_K = list(range(60, 23, -4))          # anterior -> posterior: 60,56,...,24
SHIFT = -1                                # fast slice = slow slice + SHIFT
ROWS = (2, 66)                            # crop on axis 1 (S-I)
COLS = (6, 77)                            # crop on axis 0 (L-R)
FPS = 8
CYCLES = 4


def load(name):
    g = sio.loadmat(SRC / f'img_dyn_16ph_{name}.mat', variable_names=['img_gp'])['img_gp'].astype(np.float32)
    sig = g.reshape(-1, 16).sum(0)
    return np.roll(g, -int(np.argmin(sig)), axis=3)


def panel(vol, k, b):
    sl = vol[COLS[0]:COLS[1], ROWS[0]:ROWS[1], k, b].T    # rows = S-I, cols = L-R
    return sl[::-1]                                       # superior on top


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [('SLOW / deep breathing  ~6.2 bpm (9.7 s)  (MID00018)', load('lowbpm'), SLOW_K),
            ('FAST breathing  ~27 bpm (2.2 s)  (MID00017)', load('highbpm'), [k + SHIFT for k in SLOW_K])]
    vmax = [np.percentile(np.stack([v[..., k, :] for k in ks]), 99.8) for _, v, ks in rows]
    nk = len(SLOW_K)
    fig, ax = plt.subplots(2, nk, figsize=(nk * 1.55, 2 * 1.75 + 0.7), facecolor='k')
    fig.subplots_adjust(left=0.005, right=0.995, top=0.86, bottom=0.02, wspace=0.03, hspace=0.22)
    ims = []
    for r, (lab, v, ks) in enumerate(rows):
        for j, k in enumerate(ks):
            a = ax[r, j]
            ims.append(a.imshow(panel(v, k, 0), cmap='gray', vmin=0, vmax=vmax[r], interpolation='bilinear'))
            a.set_xticks([]); a.set_yticks([])
            for s in a.spines.values():
                s.set_visible(False)
            a.set_title(f'sl {k}', color='0.75', fontsize=7, pad=2)
        ax[r, 0].text(0.0, 1.17, lab, transform=ax[r, 0].transAxes, color='w', fontsize=9, ha='left', va='bottom')
    fig.text(0.005, 0.975, '000HH 2022-10-26  gas, Faraz recon, 16 bins/cycle (phase-normalised, not real time)   '
             'coronal: ANTERIOR (left) → POSTERIOR (right); R on viewer left; fast slice = slow − 1',
             color='0.85', fontsize=7.5, va='top')
    btxt = fig.text(0.995, 0.975, '', color='y', fontsize=9, ha='right', va='top')

    def draw(b):
        i = 0
        for _, v, ks in rows:
            for k in ks:
                ims[i].set_data(panel(v, k, b))
                i += 1
        btxt.set_text(f'bin {b + 1:2d}/16' + ('  end-exp' if b == 0 else ''))

    tmp = OUT / '_frames'
    tmp.mkdir(exist_ok=True)
    for b in range(16):
        draw(b)
        fig.savefig(tmp / f'f{b:02d}.png', dpi=130, facecolor='k')
    sig = rows[0][1].reshape(-1, 16).sum(0)
    draw(int(np.argmax(sig)))
    fig.savefig(OUT / '_still_EI.png', dpi=130, facecolor='k')
    plt.close(fig)
    lst = OUT / '_frames' / 'list.txt'
    lst.write_text(''.join(f"file 'f{b:02d}.png'\nduration {1 / FPS}\n" for _ in range(CYCLES) for b in range(16))
                   + "file 'f00.png'\n")
    stem = OUT / 'HH_2022-10-26_coronal_slow_vs_fast'
    vf = 'scale=trunc(iw/2)*2:trunc(ih/2)*2'
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-f', 'concat', '-safe', '0', '-i', str(lst),
                    '-vf', vf, '-pix_fmt', 'yuv420p', '-c:v', 'libx264', '-r', str(FPS), f'{stem}.mp4'], check=True)
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(FPS), '-i', str(tmp / 'f%02d.png'),
                    '-vf', 'split[a][b];[a]palettegen[p];[b][p]paletteuse', '-loop', '0', f'{stem}.gif'], check=True)
    print('wrote', stem.with_suffix('.mp4'), stem.with_suffix('.gif'), OUT / '_still_EI.png')


if __name__ == '__main__':
    main()

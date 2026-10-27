"""Scratch (2026-10-27, Hooman): four-row video per session — s (SIGNAL, b44d_sig) / d (DIAPHRAGM, b44d) /
p (PNEUMOTACH, clock-mapped, image 40c248b) / mix (MIX, sigma 0.5 PoU, image 36f67fb), top to bottom, same coronal
slice, one frame per respiratory bin, shared grey scale. Conventions as post_process.make_diaphragm_gif.
Per-row annotation: bin, median real time after the k0 EE (s: <key>_bins_tyger.csv; d: <key>_bins_diaphragm_tyger.csv;
p / mix: recon.mat ilvbin + ee_times), hard-bin interleave count, background sigma of the frame, cyan = 15 % diaphragm
estimate.
usage: python helpers/_fourway_video.py <key> ...  -> outputs/mix_video_2026-10-27/<key>_fourway.{mp4,gif}
"""
import csv
import os
import subprocess
import sys

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon'
BT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'bin_timing_2026-10-26')
OUT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'mix_video_2026-10-27')
NB = 16


def diaphragm_pos(g):
    pos = np.zeros(g.shape[0])
    for b in range(g.shape[0]):
        prof = g[b].mean(axis=(1, 2)); above = np.where(prof > prof.max() * 0.15)[0]
        pos[b] = float(above[-1]) if len(above) else g.shape[1] / 2.0
    return pos


def times_from_mat(m):
    bins = m['ilvbin'].ravel().astype(float); t = m['ilvtime'].ravel().astype(float); ee = m['ee_times'].ravel().astype(float)
    k = np.searchsorted(ee, t, side='right') - 1; inside = (k >= 0) & (k < len(ee) - 1)
    tee = np.full(len(t), np.nan); tee[inside] = t[inside] - ee[k[inside]]
    hard = np.where(bins >= 0, np.floor(bins * NB + 0.5).astype(int) % NB, -1)
    return (np.array([np.nanmedian(tee[hard == b]) if (hard == b).any() else np.nan for b in range(NB)]),
            np.array([int((hard == b).sum()) for b in range(NB)]))


def times_from_csv(path, tcol, ncol, binning=None):
    rows = [r for r in csv.DictReader(open(path)) if binning is None or r.get('binning', binning) == binning]
    return np.array([float(r[tcol]) for r in rows]), np.array([int(float(r[ncol])) for r in rows])


def main(key):
    os.makedirs(OUT, exist_ok=True)
    stacks = [
        ('s  SIGNAL (phase fraction, σ 1)', f'{EXT}/AIkill_Dynamic_b44d_sig/{key}/s/recon.mat'),
        ('d  DIAPHRAGM (nav rank, σ 1)', f'{EXT}/AIkill_Dynamic_b44d/{key}/d/recon.mat'),
        ('p  PNEUMOTACH (volume rank, clock-mapped, σ 1)', f'{EXT}/AIkill_Dynamic_mix/{key}/p/recon.mat'),
        ('mix  MIX (insp abs clock + exhale nav rank, σ 0.5 PoU)', f'{EXT}/AIkill_Dynamic_mix/{key}/m/recon.mat'),
    ]
    mats = [sio.loadmat(p) for _, p in stacks]
    gs = [m['gas_phase'].astype(np.float64) for m in mats]
    nb, nz, ny, nx = gs[1].shape
    y0 = int(np.argmax(gs[1].mean(0).sum(axis=(0, 2))))
    vmax = np.percentile(np.concatenate([g[:, :, y0, :].ravel() for g in gs]), 99.7)
    dps = [diaphragm_pos(g) for g in gs]
    sigs = [g[:, :8, :8, :8].std(axis=(1, 2, 3)) for g in gs]
    tm = [None] * 4; nn = [None] * 4
    tm[0], nn[0] = times_from_csv(f'{BT}/{key}_bins_tyger.csv', 't_after_EE_med_s', 'n_hard_kept')
    tm[1], nn[1] = times_from_csv(f'{BT}/{key}_bins_diaphragm_tyger.csv', 't_after_EE_med_s', 'n_hard_kept')
    tm[2], nn[2] = times_from_mat(mats[2]); tm[3], nn[3] = times_from_mat(mats[3])
    extra = ['', '', f'clock chk a {float(mats[2]["pneumo_clock_chk_a"].ravel()[0]):.4f} b {float(mats[2]["pneumo_clock_chk_b"].ravel()[0]):+.2f} s',
             f'τ_I {float(mats[3]["mix_tau_I"].ravel()[0]):.2f} s n_insp {int(mats[3]["mix_n_insp"].ravel()[0])}']
    ninsp = int(mats[3]['mix_n_insp'].ravel()[0])
    frames = []
    for b in range(nb):
        fig, axs = plt.subplots(4, 1, figsize=(4.4, 16.4), dpi=100)
        for i, (ax, g, (lab, _)) in enumerate(zip(axs, gs, stacks)):
            ax.imshow(np.clip(g[b, :, y0, :] / vmax, 0, 1), cmap='gray', origin='upper', aspect='equal', vmin=0, vmax=1)
            ax.axhline(dps[i][b], color='cyan', linewidth=1.0, linestyle='--')
            ax.axis('off')
            limb = ''
            if i == 3:
                limb = '  insp' if b < ninsp else '  exhale'
            elif i in (1, 2):
                limb = '  inhale' if b < 8 else '  exhale'
            ax.set_title(f'{lab}{limb}\nbin {b:2d}/{nb - 1} · t ≈ {tm[i][b]:.2f} s after EE · n = {nn[i][b]} · σ_bg {sigs[i][b]:.0f}'
                         + (f' · {extra[i]}' if extra[i] else ''), fontsize=7.5)
        fig.suptitle(f'{key} · coronal y = {y0} · shared grey scale · cyan = 15 % diaphragm estimate', fontsize=8, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()))
        plt.close(fig)
    base = f'{OUT}/{key}_fourway'
    frames[0].save(base + '.gif', save_all=True, append_images=frames[1:], duration=400, loop=0)
    tmp = f'{OUT}/_frames4_{key}'
    os.makedirs(tmp, exist_ok=True)
    for i, f in enumerate(frames):
        f.save(f'{tmp}/f{i:02d}.png')
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '2.5', '-i', f'{tmp}/f%02d.png',
                    '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', base + '.mp4'], check=True)
    print(f'{key}: y0 {y0} -> {base}.mp4 / .gif')


if __name__ == '__main__':
    for k in sys.argv[1:]:
        main(k)

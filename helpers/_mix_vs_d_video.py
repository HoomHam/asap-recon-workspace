"""Scratch (2026-10-27, Hooman): video per session — top row MIX stack, bottom row DIAPHRAGM stack, same coronal slice,
one frame per respiratory bin (0..15), shared grey scale. Conventions as post_process.make_diaphragm_gif: gas_phase
(nbins, Z, Y, X), coronal = gas[b, :, y0, :], origin upper (diaphragm = high z = bottom), no flip.
y0 = the coronal slice with the most lung signal in the bin-mean of the DIAPHRAGM stack (same y0 for both rows).
Annotations: bin index, median real time after the k0 EE (MIX from recon.mat ilvbin/ee_times; DIAPHRAGM from
<key>_bins_diaphragm_tyger.csv), interleave count per hard bin, estimated diaphragm z (post_process rule, cyan line),
background sigma of the frame.
usage: python helpers/_mix_vs_d_video.py <key> ...  -> outputs/mix_video_2026-10-27/<key>_mix_vs_diaphragm.{mp4,gif}
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


def mix_times(m):
    bins = m['ilvbin'].ravel().astype(float); t = m['ilvtime'].ravel().astype(float); ee = m['ee_times'].ravel().astype(float)
    k = np.searchsorted(ee, t, side='right') - 1; inside = (k >= 0) & (k < len(ee) - 1)
    tee = np.full(len(t), np.nan); tee[inside] = t[inside] - ee[k[inside]]
    hard = np.where(bins >= 0, np.floor(bins * NB + 0.5).astype(int) % NB, -1)
    tm = np.array([np.nanmedian(tee[hard == b]) if (hard == b).any() else np.nan for b in range(NB)])
    n = np.array([int((hard == b).sum()) for b in range(NB)])
    return tm, n


def main(key):
    os.makedirs(OUT, exist_ok=True)
    m = sio.loadmat(f'{EXT}/AIkill_Dynamic_mix/{key}/m/recon.mat')
    d = sio.loadmat(f'{EXT}/AIkill_Dynamic_b44d/{key}/d/recon.mat', variable_names=['gas_phase'])
    gm, gd = m['gas_phase'].astype(np.float64), d['gas_phase'].astype(np.float64)
    nb, nz, ny, nx = gm.shape
    y0 = int(np.argmax(gd.mean(0).sum(axis=(0, 2))))
    vmax = np.percentile(np.concatenate([gm[:, :, y0, :].ravel(), gd[:, :, y0, :].ravel()]), 99.7)
    dpm, dpd = diaphragm_pos(gm), diaphragm_pos(gd)
    tmix, nmix = mix_times(m)
    tau, ninsp = float(m['mix_tau_I'].ravel()[0]), int(m['mix_n_insp'].ravel()[0])
    rows = list(csv.DictReader(open(f'{BT}/{key}_bins_diaphragm_tyger.csv')))
    td = np.array([float(r['t_after_EE_med_s']) for r in rows]); nd = np.array([int(r['n_hard_kept']) for r in rows])
    zd = np.array([float(r['nav_z_vox_wmean']) for r in rows])
    sig_m = gm[:, :8, :8, :8].std(axis=(1, 2, 3)); sig_d = gd[:, :8, :8, :8].std(axis=(1, 2, 3))
    frames = []
    for b in range(nb):
        fig, axs = plt.subplots(2, 1, figsize=(4.2, 8.6), dpi=110)
        for ax, g, lab, dp in ((axs[0], gm, 'MIX', dpm), (axs[1], gd, 'DIAPHRAGM', dpd)):
            ax.imshow(np.clip(g[b, :, y0, :] / vmax, 0, 1), cmap='gray', origin='upper', aspect='equal', vmin=0, vmax=1)
            ax.axhline(dp[b], color='cyan', linewidth=1.0, linestyle='--')
            ax.axis('off')
        limb = 'insp (abs clock)' if b < ninsp else 'exhale (nav rank)'
        axs[0].set_title(f'{key}  MIX  bin {b:2d}/{nb - 1}  {limb}\n'
                         f't ≈ {tmix[b]:.2f} s after EE · n = {nmix[b]} · σ_bg {sig_m[b]:.0f}  (τ_I {tau:.2f} s, n_insp {ninsp}, σ 0.5 PoU)',
                         fontsize=8)
        axs[1].set_title(f'DIAPHRAGM  bin {b:2d}/{nb - 1}  {"inhale" if b < 8 else "exhale"} rank\n'
                         f't ≈ {td[b]:.2f} s after EE · n = {nd[b]} · nav z {zd[b]:.1f} vox · σ_bg {sig_d[b]:.0f}  (σ 1, production)',
                         fontsize=8)
        fig.text(0.5, 0.005, f'coronal slice y = {y0} · shared grey scale (99.7 % of both stacks) · cyan = 15 % diaphragm estimate',
                 ha='center', fontsize=6.5)
        fig.tight_layout(rect=(0, 0.02, 1, 1))
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()))
        plt.close(fig)
    base = f'{OUT}/{key}_mix_vs_diaphragm'
    frames[0].save(base + '.gif', save_all=True, append_images=frames[1:], duration=400, loop=0)
    tmp = f'{OUT}/_frames_{key}'
    os.makedirs(tmp, exist_ok=True)
    for i, f in enumerate(frames):
        f.save(f'{tmp}/f{i:02d}.png')
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '2.5', '-i', f'{tmp}/f%02d.png',
                    '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', base + '.mp4'], check=True)
    print(f'{key}: y0 {y0}, vmax {vmax:.0f}, tau_I {tau:.2f} n_insp {ninsp} -> {base}.mp4 / .gif')


if __name__ == '__main__':
    for k in sys.argv[1:]:
        main(k)

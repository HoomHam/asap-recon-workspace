#!/usr/bin/env python3
"""
Post-process an ASAP output directory after Tyger recon:
  - saves recon.mat        (gas_phase, dissolved_phase_*, diaphragm_pos)
  - saves signal_pneumo.npz (fid_signal, pneumo_time, pneumo_vol — from input.mrd)
  - creates fig/
      axial.gif            all lung Z slices, 10-wide montage, 1 frame/bin (flipud)
      coronal.gif          all lung Y slices, 10-wide montage, 1 frame/bin (rev order)
      sagittal.gif         all lung X slices, 10-wide montage, 1 frame/bin (rot90 ccw)
      diaphragm.gif        coronal mid-slice per bin + estimated diaphragm line
      navigator.gif        low-res navigator coronal projection + tracked z line
      resp_traces.png      FID + pneumotach + navigator/image-derived diaphragm
    Slice montages share a fixed color axis (vmin=EE bin, vmax=EI bin).
    navigator.gif + the navigator traces appear only for DIAPHRAGM-binned recons.
    (gas_montage.png disabled per request)

gas_phase in .mat: shape (nbins, Z, Y, X) → in MATLAB: gas_phase(bin, z, y, x)
diaphragm_pos: (nbins,) z-index of estimated inferior lung edge per bin.

Usage:
    python post_process.py <output_dir> [--input-mrd <input.mrd>]
    python post_process.py <output.mrd> [--input-mrd <input.mrd>]
"""
import argparse
import io
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import mrd
from PIL import Image


# ---------------------------------------------------------------------------
# MRD readers

NAV_KEYS = ('nav_coronal', 'nav_diaphragm_z', 'nav_time',
            'nav_volume', 'nav_ilvtime', 'nav_volmeastime')


def read_output_mrd(mrd_path):
    """Returns (gas, dissolved, nav) where nav is a dict of navigator arrays
    (empty if the recon was not DIAPHRAGM-binned)."""
    gas = dissolved = None
    nav = {}
    with mrd.BinaryMrdReader(str(mrd_path)) as reader:
        reader.read_header()
        for item in reader.read_data():
            if isinstance(item, mrd.StreamItem.NdArrayFloat):
                if item.value.meta.get('gas_phase_image'):
                    gas = item.value.data
                else:
                    for k in NAV_KEYS:
                        if item.value.meta.get(k):
                            nav[k] = item.value.data
                            break
            elif isinstance(item, mrd.StreamItem.NdArrayComplexFloat):
                if item.value.meta.get('dissolved_phase_image'):
                    dissolved = item.value.data
    return gas, dissolved, nav


def read_input_mrd(mrd_path):
    """Extract FID signal and pneumotach from input.mrd.
    Returns (fid_signal, pneumo_time, pneumo_vol) — any may be None.
    fid_signal: (N_interleaves,) normalized [0,1]
    pneumo_time: (M,) seconds
    pneumo_vol:  (M,) integrated volume, normalized [0,1]
    """
    dyn_acq = None
    pneumo = None
    with mrd.BinaryMrdReader(str(mrd_path)) as reader:
        reader.read_header()
        for item in reader.read_data():
            if isinstance(item, mrd.StreamItem.NdArrayComplexFloat):
                if item.value.meta.get('dynamic_acquisition'):
                    dyn_acq = item.value.data  # (channels, samples, lines)
            elif isinstance(item, (mrd.StreamItem.NdArrayFloat, mrd.StreamItem.NdArrayDouble)):
                if item.value.meta.get('pneumotach'):
                    pneumo = item.value.data   # (2, N): row0=time, row1=pressure

    fid_signal = pneumo_time = pneumo_vol = None

    if dyn_acq is not None:
        # FID signal = sum of first 8 k-space point magnitudes across channels,
        # per interleave — matches raw.py:429
        sig = np.sum(np.abs(dyn_acq[:, :8, :]), axis=(0, 1)).astype(float)
        mn, mx = sig.min(), sig.max()
        fid_signal = (sig - mn) / (mx - mn + 1e-9)

    if pneumo is not None:
        t = pneumo[0].astype(float)
        P = pneumo[1].astype(float)
        vol = np.cumsum(P)   # integrate pressure → volume (unnormalized)
        mn, mx = vol.min(), vol.max()
        pneumo_vol = (vol - mn) / (mx - mn + 1e-9)
        pneumo_time = t

    return fid_signal, pneumo_time, pneumo_vol


# ---------------------------------------------------------------------------
# Diaphragm estimator

def estimate_diaphragm_pos(gas):
    """Inferior-most z with signal > 15 % of max, per bin.
    gas shape: (nbins, Z, Y, X). Returns (nbins,) float."""
    nbins = gas.shape[0]
    pos = np.zeros(nbins)
    for b in range(nbins):
        profile = gas[b].mean(axis=(1, 2))   # (Z,)
        thresh = profile.max() * 0.15
        above = np.where(profile > thresh)[0]
        pos[b] = float(above[-1]) if len(above) else gas.shape[1] / 2.0
    return pos


# ---------------------------------------------------------------------------
# GIF helpers

def _norm01(arr):
    mn, mx = arr.min(), arr.max()
    return (arr - mn) / (mx - mn + 1e-9)


def _fig_to_pil(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=100)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).copy().convert('RGB')


def _save_gif(frames, out_path, duration=80):
    frames[0].save(str(out_path), save_all=True, append_images=frames[1:],
                   loop=0, duration=duration, format='GIF')


def _ee_ei_bins(gas):
    """End-expiration / end-inspiration bin indices by total gas signal.
    EI (inspiration) = most HP gas in lung = brightest; EE = dimmest."""
    bin_sig = gas.reshape(gas.shape[0], -1).sum(axis=1)
    return int(np.argmin(bin_sig)), int(np.argmax(bin_sig))   # EE, EI


def make_allslices_video(gas, axis, out_path, ncols=10, thresh=0.10, pad=2,
                         duration=200):
    """Animated montage: every lung-bearing slice on `axis`, tiled ncols-wide,
    one frame per respiratory bin. Slice set + FOV bbox are fixed across bins
    (computed from the bin-mean), so the video shows ventilation change at a
    constant color axis.

    Color axis fixed: vmin from the EE (end-expiration) bin, vmax from the EI
    (end-inspiration) bin. Last row dropped if it holds < 4 tiles.
    Orientation: axial flipud, coronal slice-order reversed, sagittal rot90 ccw.
    axis: 1=Z(axial), 2=Y(coronal), 3=X(sagittal). gas: (bins,Z,Y,X).
    """
    nbins = gas.shape[0]
    ax0 = axis - 1
    mean_vol = gas.mean(axis=0)
    mask = mean_vol > thresh * (mean_vol.max() + 1e-9)
    msk = np.moveaxis(mask, ax0, 0)                   # (n_slice, A, B)

    idx = np.where(msk.any(axis=(1, 2)))[0]
    if len(idx) == 0:
        idx = np.arange(msk.shape[0])

    # shared in-plane lung bbox (union over kept slices), padded
    plane = msk[idx].any(axis=0)
    rows = np.where(plane.any(axis=1))[0]
    cols = np.where(plane.any(axis=0))[0]
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, msk.shape[1])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, msk.shape[2])

    if axis == 2:                                     # coronal: reverse order
        idx = idx[::-1]

    # drop the final row if it would hold fewer than 4 tiles
    rem = len(idx) % ncols
    if 0 < rem < 4:
        idx = idx[:len(idx) - rem]
    n = len(idx)
    nrows = int(np.ceil(n / ncols))

    def orient(t):
        if axis == 1:
            return np.flipud(t)
        if axis == 3:
            return np.rot90(t, 1)                     # ccw
        return t

    # fixed color axis from EE / EI bins
    ee, ei = _ee_ei_bins(gas)
    vmin = float(gas[ee].min())
    vmax = float(gas[ei].max())

    th, tw = orient(np.empty((r1 - r0, c1 - c0))).shape
    g = 1
    axis_labels = {1: ('Axial', 'z'), 2: ('Coronal', 'y'), 3: ('Sagittal', 'x')}
    title, ax_name = axis_labels[axis]
    fig_w = ncols * 0.9
    fig_h = nrows * 0.9 * (th / tw) + 0.4

    frames = []
    for b in range(nbins):
        volb = np.moveaxis(gas[b], ax0, 0)
        grid = np.full((nrows * th + (nrows - 1) * g,
                        ncols * tw + (ncols - 1) * g), vmin, dtype=float)
        for k, i in enumerate(idx):
            t = orient(volb[i, r0:r1, c0:c1])
            rr, cc = divmod(k, ncols)
            y0 = rr * (th + g)
            x0 = cc * (tw + g)
            grid[y0:y0 + th, x0:x0 + tw] = t
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=110)
        ax.imshow(grid, cmap='gray', origin='upper', aspect='equal',
                  vmin=vmin, vmax=vmax)
        tag = '  [EE]' if b == ee else ('  [EI]' if b == ei else '')
        ax.set_title(f'{title} — {n} lung slices, {nrows}×{ncols} — bin {b}{tag}',
                     fontsize=9)
        ax.axis('off')
        fig.tight_layout()
        frames.append(_fig_to_pil(fig))
    _save_gif(frames, out_path, duration)
    print(f'[post_process] saved {out_path}  ({n} slices, {nrows}x{ncols}, '
          f'{nbins} bins, EE={ee} EI={ei}, vmin={vmin:.3g} vmax={vmax:.3g})')


def make_diaphragm_gif(gas, diaphragm_pos, out_path, duration=150):
    """Coronal mid-slice per respiratory bin with estimated diaphragm line."""
    nbins, nz, ny, nx = gas.shape
    mid_y = ny // 2
    vmax = gas.max() + 1e-9
    frames = []
    for b in range(nbins):
        slc = gas[b, :, mid_y, :] / vmax
        fig, ax = plt.subplots(figsize=(3, 3), dpi=100)
        ax.imshow(slc, cmap='gray', origin='upper', aspect='equal')
        ax.axhline(diaphragm_pos[b], color='cyan', linewidth=1.5, linestyle='--')
        ax.set_title(f'bin {b}  diaphragm z={diaphragm_pos[b]:.0f}', fontsize=8)
        ax.axis('off')
        frames.append(_fig_to_pil(fig))
    _save_gif(frames, out_path, duration)
    print(f'[post_process] saved {out_path}')


def make_navigator_video(nav, out_path, duration=120):
    """Video of the low-res navigator coronal projection with the tracked
    diaphragm line overlaid, one frame per undersampled interleave group.
    nav['nav_coronal']: (nframes, z, x); nav['nav_diaphragm_z']: (nframes,)."""
    imgs = np.asarray(nav.get('nav_coronal'))
    if imgs.ndim != 3 or imgs.shape[0] == 0:
        print('[post_process] no nav_coronal frames — skipping navigator video')
        return
    lines = np.asarray(nav.get('nav_diaphragm_z', np.full(imgs.shape[0], np.nan)))
    times = np.asarray(nav.get('nav_time', np.arange(imgs.shape[0])))
    # nav_coronal is exported already oriented apex-up / diaphragm-down (tyger_recon
    # flips the navigator GPDYN before the edge-finder), and nav_diaphragm_z is the
    # diaphragm (high-z) edge in that frame — so no flip or coordinate transform here.
    # fixed contrast across frames (robust percentiles)
    vmin, vmax = np.percentile(imgs, 1), np.percentile(imgs, 99.5)
    frames = []
    for i in range(imgs.shape[0]):
        fig, ax = plt.subplots(figsize=(3.2, 3.2), dpi=110)
        ax.imshow(imgs[i], cmap='gray', origin='upper', aspect='auto',
                  vmin=vmin, vmax=vmax)
        z = lines[i]
        if np.isfinite(z):
            ax.axhline(z, color='cyan', linewidth=1.5)
        t = times[i] if i < len(times) else i
        ax.set_title(f'navigator {i}  t={t:.2f}s  z={z:.1f}' if np.isfinite(z)
                     else f'navigator {i}  t={t:.2f}s  (no fit)', fontsize=8)
        ax.axis('off')
        frames.append(_fig_to_pil(fig))
    _save_gif(frames, out_path, duration)
    print(f'[post_process] saved {out_path}  ({imgs.shape[0]} nav frames)')


def _tnorm(t):
    return (t - t[0]) / (t[-1] - t[0] + 1e-9)


def make_resp_traces(diaphragm_pos, fid_signal, pneumo_time, pneumo_vol,
                     nbins, out_path, nav=None):
    """Plot FID signal, pneumotach, and diaphragm position on one figure.
    If `nav` holds the real navigator waveform, plot that (and the raw tracked
    measurements); otherwise fall back to the image-derived per-bin estimate."""
    fig, ax = plt.subplots(figsize=(10, 4), dpi=120)

    if fid_signal is not None:
        x_sig = np.linspace(0, 1, len(fid_signal))
        ax.plot(x_sig, fid_signal, color='steelblue', lw=0.6, alpha=0.85,
                label='FID signal (k=0 envelope)')

    if pneumo_vol is not None and pneumo_time is not None:
        ax.plot(_tnorm(pneumo_time), pneumo_vol, color='tomato', lw=0.8, alpha=0.85,
                label='Pneumotach volume')

    nav_vol = None if nav is None else np.asarray(nav.get('nav_volume', []))
    nav_ilvt = None if nav is None else np.asarray(nav.get('nav_ilvtime', []))
    nav_z = None if nav is None else np.asarray(nav.get('nav_diaphragm_z', []))
    nav_t = None if nav is None else np.asarray(nav.get('nav_time', []))
    have_nav = nav_vol is not None and nav_vol.size and nav_ilvt is not None and nav_ilvt.size

    if have_nav:
        # navigator respiratory waveform (interpolated, normalized) vs time
        m = np.isfinite(nav_vol)
        ax.plot(_tnorm(nav_ilvt)[m], nav_vol[m], color='green', lw=1.3,
                label='Navigator volume (tracked)')
        # raw tracked diaphragm-z measurements
        if nav_z is not None and nav_z.size and nav_t is not None and nav_t.size:
            zf = np.isfinite(nav_z)
            if zf.any():
                z_norm = _norm01(nav_z[zf])
                ax.scatter(_tnorm(nav_t)[zf], z_norm, color='limegreen', s=12,
                           zorder=5, alpha=0.7, label='Diaphragm z (navigator meas.)')
    # image-derived per-bin estimate (faint when navigator present)
    a = 0.3 if have_nav else 1.0
    diaphragm_norm = _norm01(diaphragm_pos)
    bin_centers = (np.arange(nbins) + 0.5) / nbins
    ax.scatter(bin_centers, diaphragm_norm, color='cyan', s=50, zorder=4, alpha=a,
               label='Diaphragm pos (image-derived, per bin)')
    ax.step(bin_centers, diaphragm_norm, color='cyan', lw=1.2, alpha=a * 0.6, where='mid')

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('Normalized time (acquisition fraction)', fontsize=10)
    ax.set_ylabel('Normalized value', fontsize=10)
    ax.set_title('Respiratory Traces', fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=120)
    plt.close(fig)
    print(f'[post_process] saved {out_path}')


# ---------------------------------------------------------------------------
# Main

def run(out_dir, input_mrd_path=None):
    out_dir = Path(out_dir)
    mrd_path = out_dir / 'output.mrd'
    if not mrd_path.is_file():
        sys.exit(f'[post_process] ERROR: {mrd_path} not found')

    print(f'[post_process] reading {mrd_path} ...')
    gas, dissolved, nav = read_output_mrd(mrd_path)
    if gas is None:
        sys.exit('[post_process] ERROR: no gas_phase_image in output.mrd')
    print(f'[post_process] gas shape {gas.shape} dtype {gas.dtype}')
    if nav:
        print(f'[post_process] navigator arrays: '
              + ', '.join(f'{k}{np.asarray(v).shape}' for k, v in nav.items()))

    # signal / pneumotach from input.mrd (if available)
    fid_signal = pneumo_time = pneumo_vol = None
    npz_path = out_dir / 'signal_pneumo.npz'
    if npz_path.is_file():
        d = np.load(npz_path, allow_pickle=True)
        fid_signal  = d['fid_signal']  if 'fid_signal'  in d else None
        pneumo_time = d['pneumo_time'] if 'pneumo_time' in d else None
        pneumo_vol  = d['pneumo_vol']  if 'pneumo_vol'  in d else None
        print('[post_process] loaded signal_pneumo.npz')
    elif input_mrd_path and Path(input_mrd_path).is_file():
        print(f'[post_process] extracting signal/pneumo from {input_mrd_path} ...')
        fid_signal, pneumo_time, pneumo_vol = read_input_mrd(input_mrd_path)
        np.savez(str(npz_path),
                 fid_signal=fid_signal if fid_signal is not None else np.array([]),
                 pneumo_time=pneumo_time if pneumo_time is not None else np.array([]),
                 pneumo_vol=pneumo_vol if pneumo_vol is not None else np.array([]))
        print(f'[post_process] saved {npz_path}')
    else:
        print('[post_process] WARNING: no input.mrd — signal/pneumotach unavailable for graph')

    # --- .mat ---------------------------------------------------------------
    mat_vars = {'gas_phase': gas.astype(np.float32)}
    if dissolved is not None:
        mat_vars['dissolved_phase_magnitude'] = np.abs(dissolved).astype(np.float32)
        mat_vars['dissolved_phase_real'] = dissolved.real.astype(np.float32)
        mat_vars['dissolved_phase_imag'] = dissolved.imag.astype(np.float32)
    diaphragm_pos = estimate_diaphragm_pos(gas)
    mat_vars['diaphragm_pos'] = diaphragm_pos
    if fid_signal is not None:
        mat_vars['fid_signal'] = fid_signal.astype(np.float32)
    if pneumo_vol is not None:
        mat_vars['pneumo_time'] = pneumo_time.astype(np.float32)
        mat_vars['pneumo_vol']  = pneumo_vol.astype(np.float32)
    for k, v in (nav or {}).items():     # real navigator arrays, if present
        mat_vars[k] = np.asarray(v).astype(np.float32)
    mat_path = out_dir / 'recon.mat'
    scipy.io.savemat(str(mat_path), mat_vars)
    print(f'[post_process] saved {mat_path}')

    # --- fig/ ---------------------------------------------------------------
    fig_dir = out_dir / 'fig'
    fig_dir.mkdir(exist_ok=True)

    # All-lung-slice montage videos (one frame per respiratory bin), 10-wide
    make_allslices_video(gas, axis=1, out_path=fig_dir / 'axial.gif')
    make_allslices_video(gas, axis=2, out_path=fig_dir / 'coronal.gif')
    make_allslices_video(gas, axis=3, out_path=fig_dir / 'sagittal.gif')

    # Static all-bin gas montage — disabled per request (was gas_montage.png)
    # gp_png = out_dir / 'output_gp.png'
    # if gp_png.is_file():
    #     shutil.copy2(gp_png, fig_dir / 'gas_montage.png')
    #     print(f'[post_process] copied gas_montage.png to {fig_dir}')

    # Diaphragm image GIF (coronal, per bin, with diaphragm line)
    make_diaphragm_gif(gas, diaphragm_pos, fig_dir / 'diaphragm.gif')

    # Navigator coronal-projection video with the tracked diaphragm line
    if nav and np.asarray(nav.get('nav_coronal', [])).ndim == 3:
        make_navigator_video(nav, fig_dir / 'navigator.gif')

    # Respiratory traces graph (uses real navigator waveform when available)
    make_resp_traces(diaphragm_pos, fid_signal, pneumo_time, pneumo_vol,
                     nbins=gas.shape[0], out_path=fig_dir / 'resp_traces.png',
                     nav=nav)

    print(f'[post_process] DONE -> {fig_dir}')
    return fig_dir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output', help='output dir or output.mrd path')
    ap.add_argument('--input-mrd', default=None, metavar='FILE',
                    help='input.mrd from runs/ dir (for signal + pneumotach)')
    args = ap.parse_args()
    p = Path(args.output)
    run(p if p.is_dir() else p.parent, input_mrd_path=args.input_mrd)


if __name__ == '__main__':
    main()

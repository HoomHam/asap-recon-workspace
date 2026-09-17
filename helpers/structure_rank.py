#!/usr/bin/env python3
"""
Rank the dynamic Xe-129 gas recons (AIkill_Dynamic) by intra-lung "patchiness":
how much bright/dark edge structure sits inside the lung on 2D coronal slices,
at the max-signal bin, normalized so lung size does not matter.

Per session (binning folder `d`, gas_phase (bins, Z, Y, X)):
  bin      gas SNR-max bin from outputs/snr_2026-09-13/snr_table.csv (falls back to
           the max-lung-volume bin computed here). Hooman: "max signal, usually b6-9".
  noise    sigma = std of the 3D corner cubes (10^3) minus the 3-voxel-dilated lung.
  mask     3D lung = gas > 5 sigma, blobs >= 500 voxels (same as snr_calc.py).
  per coronal slice y (image I = gas[bin, :, y, :], apex up):
           mask2D = 3D mask slice, closed with a disk r=8 + holes filled, so dark defects
                  (below 5 sigma, up to 16 px wide) stay INSIDE the lung and their edges count.
                  Closing is done per side of the midline so it never bridges the mediastinum.
                  (A looser 2 sigma extent was tried: it adds a partial-volume band around
                  healthy lungs whose ramp reads as a rim edge, so it was dropped.)
           P    = perimeter(bright 5-sigma mask) / perimeter(mask2D): fragmentation. 1 = one
                  solid lung, >1 = patches with edges. solidity = bright area / extent area.
           w    = rim weight: 0 within 3 px of the lung boundary, 1 beyond 8 px (the signal
                  ramp at the rim is PSF/partial volume, not a patch edge); skip if sum w < MIN_PIX
           R    = I / localmean(I, sigma=10)  (mask-normalized, floored at 0.1*lung mean):
                  divides out coil shading / gravity ramp; only relative patch edges remain
           R_s  = mask-normalized gaussian(R, SIG_SMOOTH=3): patch scale, kills noise AND
                  vessel-scale detail; pixels outside the lung never leak in
           G    = w-mean |grad R_s|             dimensionless edge density (per px)
           Gn   = same on the slice's background corners with I/mean_lung
                  (what pure noise contributes at this smoothing)
           CVr  = std_{m2} R_s                  relative-image CV, trend removed
           G_fine = same at sigma 1.5 (diagnostic; healthy high-SNR lungs score high here
                  from vessels, so it is NOT the headline)
           E    = Canny edge pixels inside m2 / |m2|      edge-length density
                  (Canny on I_s/mean, sigma=1, hysteresis 0.05/0.15 of the range)
           CV   = std_{m2} I_s / mean_{m2} I_s            noise-corrected: sqrt(max(CV^2 - 1/SNR_s^2, 0))
  session  every per-slice value is averaged over slices weighted by |m2|, i.e. a mean
           over all lung pixels -> a big lung does not win by size. Also per-slice max.
           G_struct = G - Gn is the headline number.
  show     coronal slice through the CARINA (Hooman 2026-09-17: bifurcation governs): on the
           max-over-bins volume, Y template = trachea column + two diagonal arms anchored at the
           lung apex / midline; slice maximising column-mean x arms-mean (see carina_slice).
           '*' = trachea column < 5 sigma (weak). Also the max-lung-area slice.
           (The earlier trachea-column-only picker is kept as trachea_slice, unused.)

Outputs (workspace/outputs/structure_rank/):
  structure_table.csv                one row per session
  panel_bifurcation_slice.png        all sessions, carina coronal slice at max-SNR bin, sorted
  panel_maxarea_slice.png            same, max-lung-area slice
  panel_bifurcation_slice.mp4/.gif   same panel, one frame per respiratory bin (all 16 bins),
  panel_maxarea_slice.mp4/.gif       fixed grey scale per session; frames in *_frames/
  scatter_snr_vs_G.png               gas SNR (x) vs G_struct (y), labelled
  tiles.npz                          per-session (bins, z, x) crops for --figs-only

Usage:
  python structure_rank.py                 # all sessions
  python structure_rank.py --only 2023-03-27_02BB,2023-10-10_002BB
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy import ndimage
from skimage import feature

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DYN = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
WS = Path(__file__).resolve().parents[1]
OUT = WS / 'outputs' / 'structure_rank'
SNR_CSV = WS / 'outputs' / 'snr_2026-09-13' / 'snr_table.csv'

CORNER = 10
THR = 5.0
MIN_BLOB = 500
EXCL_DILATE = 3
EDGE_D0, EDGE_D1 = 3.0, 8.0   # rim weight ramp (px inside the lung boundary)
CLOSE_R = 8
MIN_BLOB2D = 100
_yy, _xx = np.ogrid[-CLOSE_R:CLOSE_R + 1, -CLOSE_R:CLOSE_R + 1]
DISK = (_yy ** 2 + _xx ** 2) <= CLOSE_R ** 2
MIN_PIX = 150            # eroded-mask pixels needed for a slice to count
SIG_SMOOTH = 3.0         # headline: COPD patch scale (~5-15 px); vessels/noise gone
SIG_FINE = 1.5           # diagnostic: keeps vessel-scale detail
SIG_TREND = 10.0         # local-mean scale divided out (coil shading, gravity ramp)
CANNY_SIG = 1.0
CANNY_LO, CANNY_HI = 0.05, 0.15
TRACHEA_HALF_W = 4
TRACHEA_ROWS = (-2, 25)
AREA_FRAC = 0.4
AREA_FRAC_BIF = 0.3
FFMPEG = '/opt/homebrew/bin/ffmpeg'


# ----------------------------------------------------------------------------- masks
def corner_box(shape, k=CORNER):
    box = np.zeros(shape, dtype=bool)
    for a in (slice(0, k), slice(shape[0] - k, shape[0])):
        for b in (slice(0, k), slice(shape[1] - k, shape[1])):
            for c in (slice(0, k), slice(shape[2] - k, shape[2])):
                box[a, b, c] = True
    return box


def lung_mask(img, sigma):
    m = img > THR * sigma
    lab, n = ndimage.label(m)
    if n == 0:
        return m
    sizes = ndimage.sum(m, lab, range(1, n + 1))
    return np.isin(lab, 1 + np.flatnonzero(sizes >= MIN_BLOB))


def corner_box2d(shape, k=CORNER):
    box = np.zeros(shape, dtype=bool)
    for a in (slice(0, k), slice(shape[0] - k, shape[0])):
        for b in (slice(0, k), slice(shape[1] - k, shape[1])):
            box[a, b] = True
    return box


# ----------------------------------------------------------------------------- per-slice metrics
def lung_extent_2d(bright, x_mid):
    """Lung extent including dark defects between bright patches: the 5-sigma mask closed
    with a disk r=CLOSE_R (bridges defects <= 2r wide) + holes filled, on each side of the
    midline separately (never bridges the mediastinum), 2D blobs < MIN_BLOB2D dropped.
    Tight on a healthy lung (no partial-volume band -> no false rim edge); defects that
    open onto the lung boundary are still lost (no proton mask here)."""
    ext = bright
    out = np.zeros_like(ext)
    for sl in (slice(0, x_mid), slice(x_mid, ext.shape[1])):
        part = ext[:, sl]
        if part.any():
            out[:, sl] = ndimage.binary_fill_holes(ndimage.binary_closing(part, structure=DISK))
    lab, n = ndimage.label(out)
    if n:
        sizes = ndimage.sum(out, lab, range(1, n + 1))
        out = np.isin(lab, 1 + np.flatnonzero(sizes >= MIN_BLOB2D))
    return out


def perimeter(m):
    """4-connected boundary length of a binary 2D mask (pixels with a non-mask neighbour)."""
    er = ndimage.binary_erosion(m, iterations=1, border_value=0)
    return float((m & ~er).sum())


def local_mean(I, m, sig):
    """Mask-normalized Gaussian local mean (no leakage from outside the lung)."""
    num = ndimage.gaussian_filter(I * m, sig)
    den = ndimage.gaussian_filter(m.astype(float), sig)
    return num / np.maximum(den, 1e-3)


def grad_mag(I, sig, m=None):
    """Gradient magnitude of the smoothed image (central-difference scale). With a mask,
    the smoothing is mask-normalized so pixels outside never leak in (no false rim edge)."""
    I_s = local_mean(I, m, sig) if m is not None else ndimage.gaussian_filter(I, sig)
    return I_s, np.hypot(ndimage.sobel(I_s, axis=0) / 8.0, ndimage.sobel(I_s, axis=1) / 8.0)


def edge_weight(m2d):
    """Per-pixel weight rising from 0 at EDGE_D0 px inside the lung boundary to 1 at EDGE_D1:
    the real signal ramp at the lung rim (PSF / partial volume) is not a patch edge."""
    d = ndimage.distance_transform_edt(m2d)
    return np.clip((d - EDGE_D0) / (EDGE_D1 - EDGE_D0), 0.0, 1.0)


def slice_metrics(I, m2d, sigma, x_mid):
    """I: (Z, X) coronal gas slice; m2d: 2D lung mask on that slice; sigma: noise std;
    x_mid: midline column. Returns dict or None if the eroded mask is too small."""
    bright = m2d                                         # 5-sigma mask slice (bright patches)
    m2d = lung_extent_2d(bright, x_mid)                 # extent incl. bridged defects
    if not m2d.any():
        return None
    bright = bright & m2d
    w = edge_weight(m2d)                                 # 0 at the rim -> 1 in the core
    wsum = float(w.sum())
    if wsum < MIN_PIX:
        return None
    mean_lung = float(I[m2d].mean())
    if mean_lung <= 0:
        return None
    bg = corner_box2d(I.shape) & ~ndimage.binary_dilation(m2d, iterations=EXCL_DILATE)
    m2 = w > 0.5

    def wmean(a):
        return float((a * w).sum() / wsum)

    # relative image: divide out the slow trend (coil shading, gravity ramp) so only
    # mid-scale bright/dark patch edges remain. Floor at 10% of lung mean so a defect
    # region does not blow up the ratio.
    lm = local_mean(I, m2d, SIG_TREND)
    R = np.where(m2d, I / np.maximum(lm, 0.1 * mean_lung), 0.0)
    # noise-only reference: same pipeline on the background corners, in lung-mean units
    Rn = np.where(bg, I / mean_lung, 0.0)

    R_s, gmag = grad_mag(R, SIG_SMOOTH, m2d)             # patch scale (headline)
    G = wmean(gmag)
    _, gn = grad_mag(Rn, SIG_SMOOTH, bg)
    Gn = float(gn[bg].mean()) if bg.sum() > 20 else np.nan
    R_f, gfine = grad_mag(R, SIG_FINE, m2d)              # vessel/fine scale (diagnostic)
    Gf = wmean(gfine)
    _, gfn = grad_mag(Rn, SIG_FINE, bg)
    Gfn = float(gfn[bg].mean()) if bg.sum() > 20 else np.nan

    edges = feature.canny(np.where(m2d, R_f, 1.0), sigma=CANNY_SIG,
                          low_threshold=CANNY_LO, high_threshold=CANNY_HI)
    E = wmean(edges.astype(float))

    I_s = local_mean(I, m2d, SIG_SMOOTH)
    cv = float(I_s[m2].std() / I_s[m2].mean())
    snr_s = I_s[m2].mean() / (sigma / (2 * np.sqrt(np.pi) * SIG_SMOOTH))   # smoothed-noise sigma, 2D gaussian
    cv_c = float(np.sqrt(max(cv ** 2 - 1.0 / snr_s ** 2, 0.0)))
    cv_r = float(np.sqrt(wmean((R_s - wmean(R_s)) ** 2)))   # relative-image CV, trend removed
    # fragmentation: how much boundary the bright (5 sigma) patches have relative to the
    # smooth lung outline; 1 = one solid lung, >1 = patches. Bright mask lightly closed
    # (r=1) so single-pixel noise raggedness does not count.
    b_s = ndimage.binary_closing(bright, structure=np.ones((3, 3)))
    p_ext = perimeter(m2d)
    P = perimeter(b_s) / p_ext if p_ext else np.nan
    solidity = float(b_s.sum() / m2d.sum())             # bright fraction of the extent (1-VDP-like)
    return dict(npx=wsum, G=G, Gn=Gn, Gf=Gf, Gfn=Gfn, E=E, CV=cv, CVc=cv_c, CVr=cv_r,
                P=P, solidity=solidity, ext=float(m2d.sum()))


# ----------------------------------------------------------------------------- carina slice
def carina_slice(gas, mask3d, sigma):
    """Coronal slice through the tracheal bifurcation (Hooman 2026-09-17: the carina governs).
    Detection volume V = max over bins (airways are brightest at their own phase), sigma-1
    smoothed. Y template anchored at the lung apex / midline: trachea column (rows apex-2..
    apex+18, |dx|<=3) and two diagonal arms (rows apex+10..apex+30, x offset 4+0.45*dz, +-2).
    score_y = mean V in column * mean V in arms, over slices with >= 30 % of the max lung
    area. Weak flag when the column mean is < 5 sigma (trachea not really seen).
    Returns (y, y_maxarea, column_signal_in_sigma, weak)."""
    V = ndimage.gaussian_filter(gas.max(axis=0), 1.0)
    nz, ny, nx = V.shape
    area = mask3d.sum(axis=(0, 2))
    y_max = int(np.argmax(area))
    zz, _, xx = np.nonzero(mask3d)
    if len(zz) == 0:
        return y_max, y_max, 0.0, True
    za, xm = int(zz.min()), int(round(xx.mean()))
    col = np.zeros((nz, nx), bool)
    col[max(za - 2, 0):min(za + 18, nz), max(xm - 3, 0):min(xm + 4, nx)] = True
    arms = np.zeros((nz, nx), bool)
    for dz in range(10, 30):
        z = za + dz
        if z >= nz:
            break
        off = 4 + int(dz * 0.45)
        for sgn in (-1, 1):
            x = xm + sgn * off
            arms[z, max(x - 2, 0):min(x + 3, nx)] = True
    ok = area >= AREA_FRAC_BIF * area.max()
    score = np.full(ny, -np.inf)
    colsig = np.zeros(ny)
    for y in np.flatnonzero(ok):
        s = V[:, y, :]
        colsig[y] = s[col].mean()
        score[y] = colsig[y] * s[arms].mean()
    y = int(np.argmax(score))
    weak = bool(colsig[y] < THR * sigma)
    return y, y_max, float(colsig[y] / sigma), weak


# ----------------------------------------------------------------------------- trachea slice (old)
def trachea_slice(vol, mask3d, sigma):
    """vol, mask3d: (Z, Y, X). Returns (y_trachea, y_maxarea, trachea_score, used_fallback)."""
    area = mask3d.sum(axis=(0, 2))                      # per y
    y_max = int(np.argmax(area))
    zz, _, xx = np.nonzero(mask3d)
    if len(zz) == 0:
        return y_max, y_max, 0.0, True
    z_apex = int(zz.min())
    x_mid = int(round(xx.mean()))
    z0 = max(z_apex + TRACHEA_ROWS[0], 0)
    z1 = min(z_apex + TRACHEA_ROWS[1], vol.shape[0])
    x0, x1 = max(x_mid - TRACHEA_HALF_W, 0), min(x_mid + TRACHEA_HALF_W + 1, vol.shape[2])
    band = vol[z0:z1, :, x0:x1]                          # (dz, Y, dx)
    score = band.mean(axis=(0, 2)) / sigma               # per y, in sigma units
    ok = area >= AREA_FRAC * area.max()
    score = np.where(ok, score, -np.inf)
    y_t = int(np.argmax(score))
    fallback = not np.isfinite(score[y_t]) or score[y_t] < THR
    if fallback:
        y_t = y_max
    return y_t, y_max, float(score[y_t]) if np.isfinite(score[y_t]) else 0.0, fallback


# ----------------------------------------------------------------------------- session
def load_snr_bins():
    bins = {}
    if not SNR_CSV.is_file():
        return bins
    with open(SNR_CSV) as f:
        for r in csv.DictReader(f):
            if r.get('error'):
                continue
            try:
                bins[r['key']] = dict(gas_max_bin=int(r['gas_max_bin']), insp_bin=int(r['insp_bin']),
                                      gas_snr_max=float(r['gas_snr_max']), gas_snr_insp=float(r['gas_snr_insp']))
            except (KeyError, ValueError):
                pass
    return bins


def session_structure(key, snr_bins):
    mat = DYN / key / 'd' / 'recon.mat'
    if not mat.is_file():
        return dict(key=key, error='no d/recon.mat'), None
    gas = sio.loadmat(mat, variable_names=['gas_phase'])['gas_phase'].astype(np.float64)
    nb = gas.shape[0]
    corners = corner_box(gas.shape[1:])
    gmean = gas.mean(0)
    mmean = lung_mask(gmean, gmean[corners].std())
    bg3 = corners & ~ndimage.binary_dilation(mmean, iterations=EXCL_DILATE)

    info = snr_bins.get(key)
    if info is not None:
        b = info['gas_max_bin']
        bin_src = 'snr_table gas_max_bin'
    else:
        vols = [lung_mask(gas[i], gas[i][bg3].std()).sum() for i in range(nb)]
        b = int(np.argmax(vols))
        bin_src = 'max lung volume (no snr_table row)'
    vol = gas[b]
    sigma = float(vol[bg3].std())
    m3 = lung_mask(vol, sigma)
    snr = float(vol[m3].mean() / sigma) if m3.any() else np.nan

    per = []
    _, _, xx = np.nonzero(m3)
    x_mid = int(round(xx.mean())) if len(xx) else vol.shape[2] // 2
    for y in range(vol.shape[1]):
        m2d = m3[:, y, :]
        if m2d.sum() < MIN_PIX:
            continue
        r = slice_metrics(vol[:, y, :], m2d, sigma, x_mid)
        if r is not None:
            r['y'] = y
            per.append(r)
    if not per:
        return dict(key=key, bin=b, error='no usable coronal slices'), None

    w = np.array([r['npx'] for r in per], dtype=float)
    def wmean(k):
        v = np.array([r[k] for r in per], dtype=float)
        ok = np.isfinite(v)
        return float(np.sum(v[ok] * w[ok]) / np.sum(w[ok])) if ok.any() else np.nan
    def vmax(k):
        v = np.array([r[k] for r in per], dtype=float)
        return float(np.nanmax(v))

    y_b, y_max, bscore, weak = carina_slice(gas, m3, sigma)
    row = dict(key=key, bin=b, bin_src=bin_src, gas_snr_bin=snr, sigma=sigma,
               lung_vox=int(m3.sum()), n_slices=len(per),
               G=wmean('G'), Gn=wmean('Gn'), G_struct=wmean('G') - wmean('Gn'),
               G_max_slice=vmax('G'), G_fine=wmean('Gf') - wmean('Gfn'),
               E=wmean('E'), E_max_slice=vmax('E'),
               CV=wmean('CV'), CVc=wmean('CVc'), CVr=wmean('CVr'),
               P=wmean('P'), P_max_slice=vmax('P'), solidity=wmean('solidity'),
               y_bif=y_b, y_maxarea=y_max, bif_col_sigma=bscore, bif_weak=int(weak), error='')
    # slices to show — ALL bins — cropped to the 3D lung bbox in (z, x), padded
    zz, _, xx = np.nonzero(m3)
    p = 4
    z0, z1 = max(zz.min() - p, 0), min(zz.max() + p + 1, vol.shape[0])
    x0, x1 = max(xx.min() - p, 0), min(xx.max() + p + 1, vol.shape[2])
    tiles = dict(bif=gas[:, z0:z1, y_b, x0:x1].astype(np.float32),        # (bins, z, x)
                 maxarea=gas[:, z0:z1, y_max, x0:x1].astype(np.float32),
                 vmax=float(np.percentile(vol[m3], 99.5)))
    return row, tiles


# ----------------------------------------------------------------------------- figures
def load_results(out):
    """Rows (numeric) from structure_table.csv + tiles from tiles.npz, for --figs-only."""
    rows = []
    with open(out / 'structure_table.csv') as f:
        for r in csv.DictReader(f):
            if r.get('error') or 'REJECTED' in r['key']:
                continue
            for k, v in r.items():
                try:
                    r[k] = int(v) if k in ('bin', 'y_bif', 'y_maxarea', 'bif_weak') else float(v)
                except (TypeError, ValueError):
                    pass
            rows.append(r)
    z = np.load(out / 'tiles.npz')
    tiles = {k: dict(bif=z[k + '__bif'], maxarea=z[k + '__maxarea'], vmax=float(z[k + '__vmax']))
             for k in {n.rsplit('__', 1)[0] for n in z.files}}
    return rows, tiles


LABEL = {'bif': 'carina (bifurcation) coronal slice (* = trachea weak, < 5 sigma)',
         'maxarea': 'max-lung-area coronal slice'}


def _draw_panel(rows, tiles, which, frame, ncols, sort, subtitle, head=''):
    """One panel figure. frame = bin index, or None for each session's own max-SNR bin.
    rows may carry 'label' (tile title prefix) set by the cohort ordering."""
    n = len(rows)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 1.75, nrows * 1.95), dpi=120, squeeze=False)
    for ax in axes.ravel():
        ax.axis('off')
    for i, r in enumerate(rows):
        ax = axes.ravel()[i]
        t = tiles[r['key']]
        b = r['bin'] if frame is None else min(frame, t[which].shape[0] - 1)
        ax.imshow(t[which][b], cmap='gray', vmin=0, vmax=t['vmax'], origin='upper', aspect='equal',
                  interpolation='nearest')
        y = r['y_bif'] if which == 'bif' else r['y_maxarea']
        tag = '*' if (which == 'bif' and r['bif_weak']) else ''
        name = r.get('label') or f"#{i + 1} {r['key']}"
        ax.set_title(f"{name}\nG={r['G_struct']:.3f} P={r['P']:.2f} "
                     f"SNR={r['gas_snr_bin']:.0f} b{r['bin']} y{y}{tag}", fontsize=5.5, pad=2)
    title = (f'{head}Gas recon, {subtitle}, {LABEL[which]} — sorted by {sort}\n'
             f'G = noise-subtracted gradient of I/localmean at sigma 3, P = perimeter(bright)/perimeter(lung); '
             f'area-weighted means over all coronal slices (lung size cancels)')
    fig.suptitle(title, fontsize=8 if ncols >= 8 else 6.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96 if nrows > 1 else 0.90))
    return fig


def _ordered(rows, tiles, sort):
    rows = [r for r in rows if not r.get('error') and r['key'] in tiles]
    if any('label' in r for r in rows):
        return rows                                   # cohort ordering already applied
    return sorted(rows, key=lambda r: -r[sort])


def panel(rows, tiles, which, path, ncols=11, sort='G_struct', head=''):
    """Static panel at each session's max-SNR bin."""
    rows = _ordered(rows, tiles, sort)
    fig = _draw_panel(rows, tiles, which, None, ncols, sort, 'max-SNR bin', head)
    fig.savefig(path)
    plt.close(fig)
    print(f'panel -> {path}')


def panel_video(rows, tiles, which, path_base, ncols=11, sort='G_struct', fps=4, head=''):
    """Same panel, one frame per respiratory bin (fixed grey scale per session) -> .mp4 + .gif
    via ffmpeg. Frame PNGs go to <path_base>_frames/ (kept, so a single frame can be reused)."""
    import subprocess
    rows = _ordered(rows, tiles, sort)
    nb = max(tiles[r['key']][which].shape[0] for r in rows)
    fdir = Path(str(path_base) + '_frames')
    fdir.mkdir(parents=True, exist_ok=True)
    for b in range(nb):
        fig = _draw_panel(rows, tiles, which, b, ncols, sort, f'bin {b}/{nb - 1}', head)
        fig.savefig(fdir / f'frame_{b:02d}.png')
        plt.close(fig)
        print(f'  frame {b + 1}/{nb}', flush=True)
    pat = str(fdir / 'frame_%02d.png')
    mp4, gif = str(path_base) + '.mp4', str(path_base) + '.gif'
    # even dimensions for yuv420p
    subprocess.run([FFMPEG, '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', pat,
                    '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
                    '-crf', '20', mp4], check=True)
    subprocess.run([FFMPEG, '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', pat,
                    '-vf', 'scale=1400:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];'
                           '[b][p]paletteuse=dither=none', gif], check=True)
    print(f'video -> {mp4}  (+ {gif}, {nb} frames @ {fps} fps)')


def scatter(rows, path):
    rows = [r for r in rows if not r.get('error')]
    x = np.array([r['gas_snr_bin'] for r in rows])
    y = np.array([r['G_struct'] for r in rows])
    fig, ax = plt.subplots(figsize=(9, 7), dpi=120)
    ax.scatter(x, y, s=18, color='#3b6ea5')
    for r in rows:
        ax.annotate(r['key'][5:], (r['gas_snr_bin'], r['G_struct']), fontsize=5, alpha=0.8,
                    xytext=(2, 2), textcoords='offset points')
    ax.set_xlabel('gas SNR at max-SNR bin')
    ax.set_ylabel('G_struct  (mean |grad R_s|, R = I / local mean, sigma 3, noise-subtracted)')
    ax.set_title('Patchiness vs SNR — want upper right', fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f'scatter -> {path}')


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', help='comma-separated session keys (date_id)')
    ap.add_argument('--out', default=str(OUT))
    ap.add_argument('--sort', default='G_struct', help='column to sort the panels by')
    ap.add_argument('--figs-only', action='store_true',
                    help='rebuild panels/scatter from structure_table.csv + tiles.npz (no recompute)')
    ap.add_argument('--no-video', action='store_true', help='skip the per-bin panel videos')
    ap.add_argument('--cohort-map', default=str(WS / 'notes' / 'calspec_package_2026-09-16' / 'cohort_map.csv'),
                    help='sid,cohort csv -> per-cohort panels in outputs/structure_rank/cohorts/ ("" = off)')
    ap.add_argument('--cohorts-only', action='store_true', help='with --figs-only: rebuild only the cohort panels')
    ap.add_argument('--fps', type=int, default=4, help='panel video frame rate (bins per second)')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.figs_only:
        rows, tiles = load_results(out)
        make_figs(rows, tiles, out, args)
        return 0

    keys = sorted(p.name for p in DYN.iterdir()
                  if p.is_dir() and not p.name.startswith(('.', '_')) and 'REJECTED' not in p.name)
    if args.only:
        keys = [k for k in keys if k in set(args.only.split(','))]
    snr_bins = load_snr_bins()

    fields = ['key', 'bin', 'bin_src', 'gas_snr_bin', 'G_struct', 'G', 'Gn', 'G_max_slice', 'G_fine',
              'E', 'E_max_slice',
              'CV', 'CVc', 'CVr', 'P', 'P_max_slice', 'solidity', 'y_bif', 'y_maxarea', 'bif_col_sigma', 'bif_weak',
              'lung_vox', 'n_slices', 'sigma', 'error']
    rows, tiles = [], {}
    with open(out / 'structure_table.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for i, k in enumerate(keys, 1):
            t0 = time.time()
            try:
                r, t = session_structure(k, snr_bins)
            except Exception as e:
                r, t = dict(key=k, error=f'{type(e).__name__}: {e}'), None
            rows.append(r)
            if t is not None:
                tiles[k] = t
            w.writerow({kk: (f'{v:.4f}' if isinstance(v, float) else v) for kk, v in r.items()})
            f.flush()
            if r.get('error'):
                print(f'[{i}/{len(keys)}] {k}: ERROR {r["error"]}', flush=True)
            else:
                print(f'[{i}/{len(keys)}] {k}: b{r["bin"]} SNR {r["gas_snr_bin"]:.1f}  G {r["G"]:.3f} '
                      f'Gn {r["Gn"]:.3f} -> G_struct {r["G_struct"]:.3f}  E {r["E"]:.3f}  CVc {r["CVc"]:.3f}  '
                      f'carina y{r["y_bif"]}{"*" if r["bif_weak"] else ""} '
                      f'({r["bif_col_sigma"]:.0f}σ)  {time.time() - t0:.0f}s', flush=True)
    np.savez_compressed(out / 'tiles.npz',
                        **{k + '__bif': t['bif'] for k, t in tiles.items()},
                        **{k + '__maxarea': t['maxarea'] for k, t in tiles.items()},
                        **{k + '__vmax': np.array(t['vmax']) for k, t in tiles.items()})
    make_figs(rows, tiles, out, args)
    good = sorted([r for r in rows if not r.get('error')], key=lambda r: -r['G_struct'])
    print('\nTop 15 by G_struct:')
    for r in good[:15]:
        print(f"  {r['key']:22s} G_struct {r['G_struct']:.3f}  E {r['E']:.3f}  CVc {r['CVc']:.3f}  SNR {r['gas_snr_bin']:.0f}")


COHORT_ORDER = ['HC', 'HC_OLD', 'EBV', 'LTX', 'RT']
COHORT_NAME = {'HC': 'healthy (young)', 'HC_OLD': 'healthy (>45)', 'EBV': 'EBV (pre/post valve)',
               'LTX': 'lung transplant (longitudinal)', 'RT': 'radiotherapy (pre/post)'}


def load_cohort_map(path):
    with open(path) as f:
        return {r['sid']: r['cohort'] for r in csv.DictReader(f)}


def cohort_rows(rows, tiles, cmap, sort):
    """Split rows by cohort. Inside a cohort: one block per subject (sessions adjacent, date
    ascending, a `_merged` session right after its original), blocks ordered by the subject's
    max `sort` value descending. Adds r['cohort'] and r['label']."""
    import collections
    out = collections.OrderedDict((c, []) for c in COHORT_ORDER)
    by = collections.defaultdict(list)
    for r in rows:
        if r.get('error') or r['key'] not in tiles:
            continue
        sid = r['key'].split('_')[1]
        r['cohort'] = cmap.get(sid, '?')
        by[(r['cohort'], sid)].append(r)
    blocks = collections.defaultdict(list)
    for (c, sid), rs in by.items():
        rs.sort(key=lambda r: (r['key'][:10], r['key'].endswith('_merged')))   # date, merged after
        for k, r in enumerate(rs, 1):
            date = r['key'][:10]
            r['label'] = f"{sid} {k}/{len(rs)}  {date}{' MERGED' if r['key'].endswith('_merged') else ''}"
        blocks[c].append(rs)
    for c, bl in blocks.items():
        bl.sort(key=lambda rs: -max(r[sort] for r in rs))
        out.setdefault(c, [])
        out[c] = [r for rs in bl for r in rs]
    return out


def make_figs(rows, tiles, out, args):
    if not args.cohorts_only:
        panel(rows, tiles, 'bif', out / 'panel_bifurcation_slice.png', sort=args.sort)
        panel(rows, tiles, 'maxarea', out / 'panel_maxarea_slice.png', sort=args.sort)
        scatter(rows, out / 'scatter_snr_vs_G.png')
        if not args.no_video:
            panel_video(rows, tiles, 'bif', out / 'panel_bifurcation_slice', sort=args.sort, fps=args.fps)
            panel_video(rows, tiles, 'maxarea', out / 'panel_maxarea_slice', sort=args.sort, fps=args.fps)
    if args.cohort_map and Path(args.cohort_map).is_file():
        cmap = load_cohort_map(args.cohort_map)
        cdir = out / 'cohorts'
        cdir.mkdir(exist_ok=True)
        for c, crows in cohort_rows([dict(r) for r in rows], tiles, cmap, args.sort).items():
            if not crows:
                continue
            head = f'{c} = {COHORT_NAME.get(c, c)}: {len(crows)} sessions, same subject adjacent (date order) · '
            for which, nm in (('bif', 'bifurcation'), ('maxarea', 'maxarea')):
                panel(crows, tiles, which, cdir / f'panel_{c}_{nm}.png', sort=args.sort, head=head)
                if not args.no_video:
                    panel_video(crows, tiles, which, cdir / f'panel_{c}_{nm}', sort=args.sort,
                                fps=args.fps, head=head)


if __name__ == '__main__':
    sys.exit(main())

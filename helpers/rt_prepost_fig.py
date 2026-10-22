#!/usr/bin/env python
"""RT pre/post coronal comparison figures.

For each RT-study subject: 10 coronal lung slices per visit, one row per
visit (top = pre-RT, bottom = post-RT). Gas-phase image from the SIGNAL-binned
recon (s/recon.mat), end-inspiration bin (brightest bin = fullest lungs).

Slice picking: lung-bearing Y slices (bin-mean > 10% max), 10 evenly spaced,
anterior->posterior order matching fig/coronal.gif. Crop = union lung bbox
across the subject's visits (pad 2). Each visit normalized to its own EI max
(HP gas intensity is polarization-dependent, not comparable across visits).

Outputs -> workspace/outputs/rt_prepost/
  rt_<ID>.png            one per subject
  rt_all_pairs.png       all pre/post pairs stacked
  rt_all_singles.png     single-visit RT studies (one row each)
"""
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic')
OUT = Path(__file__).resolve().parents[1] / 'outputs' / 'rt_prepost'
OUT.mkdir(parents=True, exist_ok=True)

NSLICES = 10
THRESH = 0.10
PAD = 2
TARGET_GRAY = 0.5      # display gray for the mean lung voxel, every row

# manual per-visit column shift, fine-tune on top of the automatic profile
# alignment: +n shifts that row's content n tiles LEFT
COL_SHIFT = {}


def _profile_shift(prof_ref, prof, max_lag=20):
    """Integer y-shift aligning `prof` to `prof_ref` (normalized xcorr)."""
    a = (prof_ref - prof_ref.mean()) / (prof_ref.std() + 1e-9)
    b = (prof - prof.mean()) / (prof.std() + 1e-9)
    lags = range(-max_lag, max_lag + 1)
    best, best_s = -np.inf, 0
    for s in lags:
        if s >= 0:
            c = float(np.dot(a[s:], b[:len(b) - s])) / (len(b) - s)
        else:
            c = float(np.dot(a[:s], b[-s:])) / (len(b) + s)
        if c > best:
            best, best_s = c, s
    return best_s


def _shift_pick(pick, s, ny):
    if s == 0:
        return pick
    step = (pick[-1] - pick[0]) / (len(pick) - 1)
    new = list(pick[s:]) + [pick[-1] + step * (k + 1) for k in range(s)]
    return np.clip(np.round(new).astype(int), 0, ny - 1)

# (subject label, [(date, folder-id, snr), ...] chronological)
SUBJECTS = [
    ('001BB', [('2023-10-10', '002BB', 33.9)]),          # pre (2023-03-27) not on drive
    ('002ZS', [('2023-11-06', '002ZS', 30.7), ('2024-05-13', '002ZS', 14.1)]),
    ('003PM', [('2023-12-06', '003PM', 58.7), ('2024-06-18', '003PM', 48.4)]),
    ('004DS', [('2024-01-16', '004DS', 22.8), ('2024-07-29', '004DS', 30.1)]),
    ('005JJ', [('2024-05-13', '005JJ', 16.1)]),
    ('006MM', [('2024-08-15', '006MM', 38.9)]),
    # 007IT 2024-08-28: recon failed container-side (curve_fit), no image
    ('008TP', [('2024-09-10', '008TP', 57.8)]),
    ('009ML', [('2024-09-26', '009ML', 49.3)]),
]


def load_ei(date, fid):
    m = sio.loadmat(DATA / f'{date}_{fid}' / 's' / 'recon.mat',
                    variable_names=['gas_phase'])
    gas = m['gas_phase']                       # (bins, Z, Y, X)
    ei = int(np.argmax(gas.reshape(gas.shape[0], -1).sum(axis=1)))
    return gas.mean(axis=0), gas[ei]           # bin-mean (for masks), EI vol


def lung_y_slices(mean_vol):
    """Lung-bearing Y slices by masked-voxel COUNT per slice (a bare .any()
    keeps near-empty edge slices where one noise voxel clears threshold)."""
    mask = mean_vol > THRESH * (mean_vol.max() + 1e-9)
    count = np.moveaxis(mask, 1, 0).sum(axis=(1, 2))   # (Y,)
    keep = count >= max(30, 0.10 * count.max())
    idx = np.where(keep)[0]
    if len(idx) == 0:
        idx = np.arange(mask.shape[1])
    return idx, mask


def bbox_zx(mask, idx):
    plane = np.moveaxis(mask, 1, 0)[idx].any(axis=0)   # (Z, X)
    rows = np.where(plane.any(axis=1))[0]
    cols = np.where(plane.any(axis=0))[0]
    return rows[0], rows[-1] + 1, cols[0], cols[-1] + 1


def subject_rows(subj, visits):
    """-> list of (row_label, [10 cropped slices]) for one subject."""
    loaded = [load_ei(d, f) for d, f, _ in visits]
    # union crop across visits
    z0 = x0 = 10 ** 9
    z1 = x1 = 0
    masks_idx = []
    for mean_vol, _ in loaded:
        idx, mask = lung_y_slices(mean_vol)
        masks_idx.append(idx)
        a, b, c, d = bbox_zx(mask, idx)
        z0, z1, x0, x1 = min(z0, a), max(z1, b), min(x0, c), max(x1, d)

    ny = loaded[0][1].shape[1]
    if len(visits) == 2:
        # same voxel spacing in both rows; post row offset by y-profile
        # cross-correlation. Base picks span the INTERSECTION of the two
        # lung extents (post extent mapped into pre coords by the shift),
        # so every column holds real lung in both rows.
        prof = [mv.mean(axis=(0, 2)) for mv, _ in loaded]
        s = _profile_shift(prof[0], prof[1])
        i0, i1 = masks_idx
        lo = max(i0[0], i1[0] + s)
        hi = min(i0[-1], i1[-1] + s)
        if hi - lo < NSLICES:                  # degenerate overlap fallback
            lo, hi = i0[0], i0[-1]
        base = np.round(np.linspace(lo, hi, NSLICES)).astype(int)[::-1]
        print(f'  {subj}: post y-shift {-s:+d} vox, extent pre {i0[0]}-{i0[-1]} '
              f'post {i1[0]}-{i1[-1]} -> base {lo}-{hi}')
        picks = [base, np.clip(base - s, 0, ny - 1)]
    else:
        idx = masks_idx[0]
        picks = [idx[np.round(np.linspace(0, len(idx) - 1, NSLICES))
                     .astype(int)][::-1]]
    picks = [_shift_pick(p, COL_SHIFT.get((subj, d), 0), ny)
             for p, (d, _, _) in zip(picks, visits)]
    z0, x0 = max(z0 - PAD, 0), max(x0 - PAD, 0)
    z1 = min(z1 + PAD, loaded[0][1].shape[0])
    x1 = min(x1 + PAD, loaded[0][1].shape[2])

    rows = []
    tags = (['Pre-RT', 'Post-RT'] if len(visits) == 2 else [''])
    for (date, _, snr), (mean_vol, ei_vol), pick, tag in zip(visits, loaded,
                                                             picks, tags):
        # equalized color axis: scale so the mean lung voxel lands at the same
        # gray in every row — rows then match in average brightness across
        # visits/subjects; airway hotspots clip to white
        lung = ei_vol[mean_vol > THRESH * (mean_vol.max() + 1e-9)]
        v = np.clip(ei_vol * (TARGET_GRAY / (lung.mean() + 1e-9)), 0, 1)
        tiles = [v[z0:z1, y, x0:x1] for y in pick]
        label = f'{tag}\n{date}\nSNR {snr:.1f}'.lstrip('\n')
        rows.append((label, tiles))
    return rows


def draw(rows, title, out_path):
    nrows = len(rows)
    th, tw = rows[0][1][0].shape
    fig_w = NSLICES * 1.25 + 1.0
    fig_h = nrows * 1.25 * (th / tw) + 0.55
    fig, axes = plt.subplots(nrows, NSLICES, figsize=(fig_w, fig_h),
                             facecolor='black', squeeze=False)
    for r, (label, tiles) in enumerate(rows):
        for c in range(NSLICES):
            ax = axes[r, c]
            ax.imshow(tiles[c], cmap='gray', vmin=0, vmax=1,
                      origin='upper', aspect='equal')
            ax.axis('off')
        axes[r, 0].text(-0.08, 0.5, label, transform=axes[r, 0].transAxes,
                        ha='right', va='center', color='white', fontsize=8)
    fig.suptitle(title, color='white', fontsize=11, y=0.995)
    fig.subplots_adjust(left=0.075, right=0.998, top=1 - 0.42 / fig_h,
                        bottom=0.005, wspace=0.02, hspace=0.02)
    fig.savefig(out_path, dpi=200, facecolor='black')
    plt.close(fig)
    print(f'saved {out_path}')


def stack(figs, title, out_path):
    """Stack multiple subjects into one figure: figs = [(subj, rows), ...]."""
    all_rows, seps = [], []
    for subj, rows in figs:
        for label, tiles in rows:
            all_rows.append((f'{subj}\n{label}' if '\n' not in label or True
                             else label, tiles))
        seps.append(len(all_rows))
    nrows = len(all_rows)
    th, tw = all_rows[0][1][0].shape
    fig_w = NSLICES * 1.25 + 1.1
    fig_h = nrows * 1.25 * (th / tw) + 0.55
    fig, axes = plt.subplots(nrows, NSLICES, figsize=(fig_w, fig_h),
                             facecolor='black', squeeze=False)
    for r, (label, tiles) in enumerate(all_rows):
        for c in range(NSLICES):
            ax = axes[r, c]
            ax.imshow(tiles[c], cmap='gray', vmin=0, vmax=1,
                      origin='upper', aspect='equal')
            ax.axis('off')
        axes[r, 0].text(-0.08, 0.5, label, transform=axes[r, 0].transAxes,
                        ha='right', va='center', color='white', fontsize=7)
    for s in seps[:-1]:
        y = 1 - (0.42 / fig_h) - (1 - 0.42 / fig_h - 0.005) * s / nrows
        fig.add_artist(plt.Line2D([0.01, 0.99], [y, y], color='0.35',
                                  lw=0.8, transform=fig.transFigure))
    fig.suptitle(title, color='white', fontsize=11, y=0.995)
    fig.subplots_adjust(left=0.085, right=0.998, top=1 - 0.42 / fig_h,
                        bottom=0.005, wspace=0.02, hspace=0.02)
    fig.savefig(out_path, dpi=200, facecolor='black')
    plt.close(fig)
    print(f'saved {out_path}')


def main():
    pairs, singles = [], []
    for subj, visits in SUBJECTS:
        rows = subject_rows(subj, visits)
        title = (f'{subj} — RT pre/post, coronal (EI bin, SIGNAL recon)'
                 if len(rows) == 2 else
                 f'{subj} — RT single visit, coronal (EI bin, SIGNAL recon)')
        draw(rows, title, OUT / f'rt_{subj}.png')
        (pairs if len(rows) == 2 else singles).append((subj, rows))
    if pairs:
        stack(pairs, 'RT studies — pre (top) / post (bottom) per subject, '
              '10 coronal slices, EI bin', OUT / 'rt_all_pairs.png')
    if singles:
        stack(singles, 'RT studies — single visit, 10 coronal slices, EI bin',
              OUT / 'rt_all_singles.png')


if __name__ == '__main__':
    main()

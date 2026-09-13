#!/usr/bin/env python3
"""
Scratch check (2026-09-12): does the fork's DIAPHRAGM z-flip (commit db80f16) only
undo an orientation difference introduced by Kento's MRD input path, or does it
change Steve's navigator behaviour on the same data?

Part A — raw data: Steve's original reader (mapVBVD inside raw.load, exact calls)
vs the dynamic array Kento's converter wrote to input.mrd (after the same killpts
trim + transpose that tyger_recon.py / raw.load_from_arr apply). Identical arrays
=> identical z-order downstream (dyn_usimg_recon is the same function in both).

Part B — edge-finder: exported navigator projections (fork orientation = z flipped
back) -> dropoff profile -> Steve's edge-finder verbatim, run (1) as the fork runs
it and (2) on the unflipped profile as Steve's GUI does. Diaphragm moves with
breathing; apex barely moves -> excursion tells which edge each one tracked.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mapvbvd
import mrd
import numpy as np
import scipy.io as sio

SUBJ = '2024-03-06_030DN'
DAT = Path('/Volumes/HoomHamExt/_5t_images_roundtrip/Images/2024-03-06/030DN/'
           'meas_MID00543_FID08467_fa_spiral_dyn_fancy_v3_20240130.dat')
OUT = Path('/Volumes/HoomHamExt/AIkill_Dynamic') / SUBJ / 'd'
FIG_DIR = Path(__file__).resolve().parents[1] / 'outputs' / 'zorder_check_2026-09-12'
KILLPTS = 2


# ---------------- Part A: raw data equivalence ----------------
def steve_read(dat):
    """Steve's raw.load (siemens branch), dynamic file, verbatim calls."""
    tw = mapvbvd.mapVBVD(str(dat))
    try:
        tw.image.flagRemoveOS = False
    except Exception:
        tw = tw[1]
    tw.image.flagRemoveOS = False
    tw.image.squeeze = True
    a = tw.image.unsorted()
    if a.ndim == 2:                       # Steve's unsqueeze()
        a = a[:, None, :]
    return a.astype('complex64')[KILLPTS:, :, :]   # (samples, ch, lines)


def kento_read(input_mrd):
    """dynamic_acquisition from input.mrd -> tyger_recon killpts trim -> load_from_arr transpose."""
    dyn = None
    with mrd.BinaryMrdReader(str(input_mrd)) as r:
        r.read_header()
        for it in r.read_data():
            if isinstance(it, mrd.StreamItem.NdArrayComplexFloat) and \
                    it.value.meta.get('dynamic_acquisition'):
                dyn = it.value.data
    dyn = dyn[:, KILLPTS:, :]                          # (ch, samples, lines)
    return np.ascontiguousarray(dyn.transpose(1, 0, 2))  # (samples, ch, lines)


print(f'== Part A: raw data, {DAT.name}')
s = steve_read(DAT)
k = kento_read(OUT / 'input.mrd')
print(f'  Steve  (mapVBVD) shape {s.shape} {s.dtype}')
print(f'  Kento  (MRD)     shape {k.shape} {k.dtype}')
if s.shape == k.shape:
    eq = np.array_equal(s, k)
    print(f'  identical: {eq}   max|diff| = {np.max(np.abs(s - k)):.3g}')
else:
    n = min(s.shape[2], k.shape[2])
    d = np.max(np.abs(s[:, :, :n] - k[:, :, :n]))
    print(f'  SHAPE MISMATCH; first {n} lines max|diff| = {d:.3g}')
    eq = False
# acquisition-order fingerprint used by raw.py to split gas/dissolved
for name, arr in (('Steve', s), ('Kento', k)):
    pat = np.abs(np.fft.fft(arr[0, 0, :]))
    n = len(pat)
    third = max(pat[n // 3 - 5:n // 3 + 5]) > pat[0] / 5
    half = max(pat[n // 2 - 5:n // 2 + 5]) > pat[0] / 3
    print(f'  {name}: pattern {"gas-dis-dis" if third else "gas-dissolved" if half else "gas-only"}')


# ---------------- Part B: edge-finder with and without flip ----------------
def steve_edge(dropoff):
    """Steve's main.py calcLVcb edge-finder, verbatim logic. Returns z or nan."""
    IS = len(dropoff)
    lo, hi = np.min(dropoff), np.max(dropoff)
    for p1 in range(IS - 1, 0, -1):
        if dropoff[p1] > lo + 0.25 * (hi - lo):
            break
    for p2 in range(p1, 0, -1):
        if dropoff[p2] > lo + 0.75 * (hi - lo):
            break
    if p1 - p2 < 2:
        return np.nan, p1, p2
    p = np.polyfit(np.arange(p2, p1), dropoff[p2:p1], 2)
    det = np.sqrt(p[1] ** 2 - 4 * p[0] * (p[2] - (hi + lo) / 2))
    m1 = (-p[1] + det) / (2 * p[0])
    m2 = (-p[1] - det) / (2 * p[0])
    z = m1 if p2 < m1 < p1 else (m2 if p2 < m2 < p1 else np.nan)
    return z, p1, p2


print(f'\n== Part B: navigator edge-finder, {SUBJ} DIAPHRAGM run')
m = sio.loadmat(OUT / 'recon.mat', variable_names=['nav_coronal', 'nav_diaphragm_z', 'nav_time'])
nav = m['nav_coronal']                     # (frames, z, x), fork orientation
t = m['nav_time'].ravel()
stored = m['nav_diaphragm_z'].ravel()
IS = nav.shape[1]
drop_fork = nav.sum(axis=2)                # == np.sum(proj, (1, 2)) in both versions

z_fork = np.array([steve_edge(d)[0] for d in drop_fork])
z_steve_raw = np.array([steve_edge(d[::-1])[0] for d in drop_fork])   # unflipped, Steve's GUI
z_steve = (IS - 1) - z_steve_raw                                      # mapped to fork coords

ok = ~np.isnan(stored)
print(f'  recomputed fork z matches exported nav_diaphragm_z: '
      f'{np.allclose(z_fork[ok], stored[ok], equal_nan=True)} ({ok.sum()}/{len(ok)} frames with a fit)')


def summary(name, z):
    good = z[~np.isnan(z)]
    if len(good) == 0:
        print(f'  {name}: no fits')
        return
    print(f'  {name}: fits {len(good)}/{len(z)}, median z {np.median(good):.1f}, '
          f'excursion p5-p95 {np.percentile(good, 95) - np.percentile(good, 5):.2f} px')


summary('fork (flipped)   ', z_fork)
summary('Steve (unflipped)', z_steve)

# figure
FIG_DIR.mkdir(parents=True, exist_ok=True)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
mean_img = nav.mean(axis=0)
ax[0].imshow(mean_img, cmap='gray', aspect='auto', origin='upper')
for z, c, lab in ((z_fork, 'tab:orange', 'fork edge (flipped)'),
                  (z_steve, 'tab:cyan', 'Steve edge (unflipped)')):
    if np.any(~np.isnan(z)):
        ax[0].axhline(np.nanmedian(z), color=c, lw=2, label=lab)
ax[0].set_title('mean navigator projection (fork orientation)')
ax[0].set_xlabel('x'); ax[0].set_ylabel('z (row)')
ax[0].legend(loc='lower right', fontsize=8)

i = int(np.nanargmax(drop_fork.sum(1)))
ax[1].plot(drop_fork[i], np.arange(IS), color='k')
ax[1].invert_yaxis()
for z, c in ((z_fork[i], 'tab:orange'), (z_steve[i], 'tab:cyan')):
    if not np.isnan(z):
        ax[1].axhline(z, color=c, lw=2)
ax[1].set_title(f'dropoff profile, frame {i}')
ax[1].set_xlabel('summed signal'); ax[1].set_ylabel('z (row)')

ax[2].plot(t, z_fork, color='tab:orange', label='fork (flipped)')
ax[2].plot(t, z_steve, color='tab:cyan', label='Steve (unflipped)')
ax[2].invert_yaxis()
ax[2].set_title('tracked edge vs time')
ax[2].set_xlabel('time (s)'); ax[2].set_ylabel('z (row, fork orientation)')
ax[2].legend(fontsize=8)
fig.suptitle(f'{SUBJ}: DIAPHRAGM edge-finder with vs without db80f16 z-flip'
             f'   |   raw data Steve==Kento: {eq}')
fig.tight_layout()
out_png = FIG_DIR / f'{SUBJ}_zorder_edgefinder.png'
fig.savefig(out_png, dpi=110)
print(f'\n  figure -> {out_png}')

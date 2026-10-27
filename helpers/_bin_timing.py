"""Scratch (2026-10-26, PCA Step-1 cross-audit): bin timing of Steve's DIAPHRAGM and SIGNAL binning.

DIAPHRAGM: nav_volume / nav_ilvtime from the b44d recon.mat -> raw.bin() (verbatim copy) -> rank bins.
SIGNAL:    k0 signal per line from the old-tree signal_pneumo.npz -> Steve's two-pass EE-minimum finder
           (raw.py load_from_arr, N = 50) -> time-phase bins.
Both: hard count per bin, soft-kernel effective weight per bin (exp(-d^2/2), cyclic), breath period,
median seconds-after-EE per bin, spread. Prints a table per study.
"""
import sys
import numpy as np
import scipy.io as sio

EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon'
NBINS = 16
TR = 0.0172


def steve_bin(val):
    """raw.bin() verbatim (rank bins, inhale ascending then exhale reflected; nan -> -1; result in [0,1))."""
    bins = val.copy()
    v = np.array([x for x in val if not np.isnan(x)])
    if len(v) > 2:
        breathdir = np.zeros(len(v))
        if v[1] < v[0]:
            breathdir[0] = 1
        for iilv in range(1, len(v)):
            if v[iilv] < v[iilv - 1]:
                breathdir[iilv] = 1
        iv = np.sort(v[breathdir == 1])
        ev = np.sort(v[breathdir == 0])
        maxiv = 0.0 if len(iv) == 0 else max(iv) * 1.00001
        maxev = 0.0 if len(ev) == 0 else max(ev)
        v[breathdir == 1] = maxiv + maxev - v[breathdir == 1]
        vsort = np.sort(v)
        vidx = 0
        for iilv in range(0, len(bins)):
            if np.isnan(bins[iilv]):
                continue
            bins[iilv] = np.searchsorted(vsort, v[vidx])
            vidx += 1
        bins[np.isnan(bins)] = -1
        bins /= len(v)
    return bins


def steve_ee_minima(sig, t, N=50):
    """raw.py two-pass strict local minimum within +-N samples."""
    ee = []
    for j in range(N, len(t) - N - 1):
        ok = True
        for k in range(j - N, j + N + 1):
            if k != j and sig[k] <= sig[j]:
                ok = False
                break
        if ok:
            ee.append(j)
    return np.array(ee, dtype=int)


def kernel_weff(binctr):
    """effective weight per bin: sum over interleaves of exp(-d^2/bindist0sq), cyclic (recon.py)."""
    w = np.zeros(NBINS)
    bc = binctr[binctr >= 0] * NBINS
    for ib in range(NBINS):
        d = np.minimum.reduce([np.abs(ib - bc), np.abs(ib + NBINS - bc), np.abs(ib - NBINS - bc)])
        w[ib] = np.sum(np.exp(-d * d / 2.0))
    return w


def per_bin_time(binfrac, t, ee_t):
    """for each hard bin: median and IQR of (t - t_EE of own breath) / period of own breath, and in s."""
    ib = np.floor(binfrac * NBINS + 0.5).astype(int) % NBINS   # nearest integer bin = kernel centre
    rows = []
    k = np.searchsorted(ee_t, t, side='right') - 1
    valid = (binfrac >= 0) & (k >= 0) & (k < len(ee_t) - 1)
    per = np.diff(ee_t)
    phase = np.full(len(t), np.nan)
    phase[valid] = (t[valid] - ee_t[k[valid]]) / per[k[valid]]
    secs = np.full(len(t), np.nan)
    secs[valid] = t[valid] - ee_t[k[valid]]
    for b in range(NBINS):
        m = valid & (ib == b)
        if m.sum() == 0:
            rows.append((b, 0, np.nan, np.nan, np.nan, np.nan))
            continue
        p = phase[m]
        s = secs[m]
        rows.append((b, int(m.sum()), np.median(s), np.percentile(s, 25), np.percentile(s, 75),
                     np.percentile(p, 75) - np.percentile(p, 25)))
    return rows


def report(name, binfrac, t, ee_t, used):
    per = np.diff(ee_t)
    w = kernel_weff(binfrac)
    print(f'\n== {name}: used {used} of {len(t)} gas interleaves, {len(ee_t)} EE -> {len(per)} breaths, '
          f'period median {np.median(per):.2f} s (IQR {np.percentile(per,25):.2f}-{np.percentile(per,75):.2f}), '
          f'hard count/bin ~{used/NBINS:.0f}, kernel w_eff/bin mean {w.mean():.0f} (= {w.mean()/(used/NBINS):.2f} x hard)')
    print('bin  n_hard  w_eff  t_med(s)  t_q25  t_q75  phaseIQR')
    for (b, n, tm, q1, q3, piqr) in per_bin_time(binfrac, t, ee_t):
        print(f'{b:3d} {n:7d} {w[b]:6.0f} {tm:8.2f} {q1:6.2f} {q3:6.2f} {piqr:8.3f}')


def diaphragm(study):
    m = sio.loadmat(f'{EXT}/AIkill_Dynamic_b44d/{study}/d/recon.mat',
                    variable_names=['nav_volume', 'nav_ilvtime', 'nav_time', 'nav_diaphragm_z'])
    v = m['nav_volume'].ravel().astype(float)
    t = m['nav_ilvtime'].ravel().astype(float)
    used = int(np.sum(~np.isnan(v)))
    binfrac = steve_bin(v)
    # EE of the navigator itself: strict minima of nav_volume within +-50 interleaves (same rule as SIGNAL)
    vv = np.where(np.isnan(v), np.inf, v)
    ee = steve_ee_minima(vv, t, 50)
    ee = ee[np.isfinite(vv[ee])]
    report(f'{study} DIAPHRAGM (b44d)', binfrac, t, t[ee], used)
    return t, binfrac


def signal(study, numspec, ilvperTR=2):
    z = np.load(f'{EXT}/AIkill_Dynamic/{study}/d/signal_pneumo.npz')
    s = z['fid_signal'].astype(float)
    gas = s[numspec::ilvperTR]
    t = np.arange(len(gas)) * TR
    ee = steve_ee_minima(gas, t, 50)
    ee_t = t[ee]
    binfrac = np.full(len(t), -1.0)
    for j in range(len(t)):
        if len(ee_t) and ee_t[0] <= t[j] <= ee_t[-1]:
            k = np.searchsorted(ee_t, t[j], side='right') - 1
            k = min(k, len(ee_t) - 2)
            binfrac[j] = (t[j] - ee_t[k]) / (ee_t[k + 1] - ee_t[k])
    used = int(np.sum(binfrac >= 0))
    report(f'{study} SIGNAL (k0, N=50)', binfrac, t, ee_t, used)
    if 'pneumo_time' in z.files and z['pneumo_time'].size:
        pt = z['pneumo_time']
        print(f'   pneumotach span {pt.min():.1f}-{pt.max():.1f} s (gas clock ends {t[-1]:.1f} s)')
    return t, binfrac, ee_t


if __name__ == '__main__':
    for study in sys.argv[1:]:
        try:
            diaphragm(study)
        except Exception as e:
            print(f'{study} DIAPHRAGM: {e}')
        try:
            signal(study, numspec=780)
        except Exception as e:
            print(f'{study} SIGNAL: {e}')

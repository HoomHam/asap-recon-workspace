"""Scratch (2026-10-26, PCA ask): DIAPHRAGM-bin equivalent of <key>_bins_tyger.csv with Tyger's EXACT membership.

tyger_recon._diaphragm_navigator exports nav_volume = g_raw.ilvvol[DIAPHRAGM] AFTER raw.rescale, i.e. exactly the
array it passes to raw.bin(). So calling the ROOT raw.bin (not a re-implementation) on the exported nav_volume gives
ilvbin[DIAPHRAGM] as gridded (float32 export rounding can only swap exact ties by 1/N bin). Excluded interleaves
(auto low-signal sets) are NaN in nav_volume -> bin -1, so exclusions are honoured by construction.
Time column: DIAPHRAGM has no EE list; times are relative to Tyger's SIGNAL EE picks (<key>_ee_tyger.csv, same raw,
binning-independent), kernel-weighted mean/sd per bin. Inhale bins 0-7 (volume ascending), exhale 8-15 (descending).

Outputs (workspace/outputs/bin_timing_2026-10-26/):
  <key>_ilvbin_diaphragm_tyger.csv  ilv, t_s, ilvbin_frac, nav_volume, nav_z_vox, breathdir(1=volume falling)
  <key>_bins_diaphragm_tyger.csv    per bin: n_hard, w_eff, t_after_EE_wmean_s, t_after_EE_wsd_s, t_after_EE_med_s,
                                    t_before_nextEE_wmean_s, nav_z_vox_wmean, nav_z_vox_wsd, nav_volume_wmean,
                                    k0_gas_raw_wmean, k0_diss_raw_wmean, img_sum_all, img_sum_all_over_k0_norm
usage: cd <root> && python3 workspace/helpers/_bins_diaphragm_tyger.py <key> ...
"""
import csv
import os
import sys

import numpy as np
import scipy.io as sio

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(__file__))
from raw import raw as steve_raw            # noqa: E402  root (fork) raw.bin, verbatim
from _bin_timing import steve_bin           # noqa: E402  my earlier copy, for the comparison only

EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon'
OUT = os.path.join(ROOT, 'workspace', 'outputs', 'bin_timing_2026-10-26')
NB = 16


def main(key):
    m = sio.loadmat(f'{EXT}/AIkill_Dynamic_b44d/{key}/d/recon.mat',
                    variable_names=['nav_volume', 'nav_ilvtime', 'nav_diaphragm_z', 'nav_time', 'gas_phase'])
    navv = m['nav_volume'].ravel().astype(np.float64)
    t = m['nav_ilvtime'].ravel().astype(np.float64)
    fz, ft = m['nav_diaphragm_z'].ravel().astype(np.float64), m['nav_time'].ravel().astype(np.float64)
    g = m['gas_phase'].astype(np.float64)
    n = len(t)
    bins = steve_raw().bin(navv.copy())                      # root raw.bin on the exact array Tyger binned
    mine = steve_bin(navv.copy())
    good = np.isfinite(fz)
    navz = np.interp(t, ft[good], fz[good], left=np.nan, right=np.nan); navz[np.isnan(navv)] = np.nan
    v = navv.copy(); ok = np.isfinite(v)
    bd = np.zeros(n); vv = v[ok]; d = np.zeros(len(vv)); d[0] = 1 if vv[1] < vv[0] else 0
    d[1:] = (vv[1:] < vv[:-1]).astype(float); bd[ok] = d
    raw_rows = list(csv.DictReader(open(f'{OUT}/{key}_ilv_raw.csv')))
    assert len(raw_rows) == n, (len(raw_rows), n)
    k0g = np.array([float(r['k0_gas_sum8']) for r in raw_rows]); k0d = np.array([float(r['k0_diss_sum8']) for r in raw_rows])
    ee = np.loadtxt(f'{OUT}/{key}_ee_tyger.csv', skiprows=1)
    k = np.searchsorted(ee, t, side='right') - 1
    inside = (k >= 0) & (k < len(ee) - 1)
    tee = np.full(n, np.nan); tee[inside] = t[inside] - ee[k[inside]]
    nxt = np.full(n, np.nan); nxt[inside] = ee[k[inside] + 1] - t[inside]
    # auto-excluded interleaves (raw.py low-signal sets) are zeroed before gridding and skipped by the kernel, but
    # nav_volume is linearly INTERPOLATED across a mid-scan gap (np.interp inside the measured range), so they do get
    # a bin: they must carry zero weight here, exactly as in the image. Flags from _tyger_bins_local.py (same raw).
    excl = np.array([int(r['excluded']) for r in csv.DictReader(open(f'{OUT}/{key}_ilvbin_tyger.csv'))])
    assert len(excl) == n
    binned_excl = int(((bins >= 0) & (excl == 1)).sum())
    okb = (bins >= 0) & (excl == 0)
    c = bins * NB
    w = np.zeros((n, NB))
    for ib in range(NB):
        dd = np.minimum.reduce([np.abs(ib - c), np.abs(ib + NB - c), np.abs(ib - NB - c)])
        w[:, ib] = np.where(okb, np.exp(-dd * dd / 2.0), 0.0)
    hard = np.where(okb, np.floor(c + 0.5).astype(int) % NB, -1)

    def wstat(x):
        mm = np.isfinite(x); ww = w * mm[:, None]; s = ww.sum(0) + 1e-12
        mean = (ww * np.nan_to_num(x)[:, None]).sum(0) / s
        sd = np.sqrt(np.maximum((ww * (np.nan_to_num(x)[:, None] - mean) ** 2).sum(0) / s, 0))
        return mean, sd
    tm, tsd = wstat(tee); nm, _ = wstat(nxt); zm, zsd = wstat(navz); vm, _ = wstat(navv)
    kg, _ = wstat(k0g); kd, _ = wstat(k0d)
    img_all = g.sum((1, 2, 3)); ratio = img_all / kg; ratio /= ratio.mean()
    np.savetxt(f'{OUT}/{key}_ilvbin_diaphragm_tyger.csv', np.column_stack([np.arange(n), t, bins, navv, navz, bd]),
               delimiter=',', header='ilv,t_s,ilvbin_frac,nav_volume,nav_z_vox,breathdir', comments='',
               fmt='%d,%.4f,%.6f,%.5f,%.3f,%d')
    hard_all = np.where(bins >= 0, np.floor(c + 0.5).astype(int) % NB, -1)
    with open(f'{OUT}/{key}_bins_diaphragm_tyger.csv', 'w') as f:
        f.write('bin,n_hard_kept,w_eff_kept,n_hard_excluded,t_after_EE_wmean_s,t_after_EE_wsd_s,t_after_EE_med_s,'
                't_before_nextEE_wmean_s,nav_z_vox_wmean,nav_z_vox_wsd,nav_volume_wmean,k0_gas_raw_wmean,'
                'k0_diss_raw_wmean,img_sum_all,img_sum_all_over_k0_norm\n')
        for b in range(NB):
            hb = hard == b; hx = (hard_all == b) & (excl == 1)
            med = np.nanmedian(tee[hb]) if hb.any() else np.nan
            f.write(f'{b},{hb.sum()},{w[:, b].sum():.1f},{hx.sum()},{tm[b]:.3f},{tsd[b]:.3f},{med:.3f},{nm[b]:.3f},'
                    f'{zm[b]:.3f},{zsd[b]:.3f},{vm[b]:.4f},{kg[b]:.6g},{kd[b]:.6g},{img_all[b]:.6g},{ratio[b]:.4f}\n')
    both = (bins >= 0) & (mine >= 0)
    dif = np.abs(bins[both] - mine[both]) * NB
    print(f'== {key}: nav-binned {int((bins >= 0).sum())}/{n} (NaN nav: {int((bins < 0).sum())}), of which auto-excluded '
          f'{binned_excl} -> kept {int(okb.sum())}; root raw.bin vs my steve_bin: max |dbin| {dif.max():.4f}')
    print('    n_hard excluded per bin:', [int(((hard_all == b) & (excl == 1)).sum()) for b in range(NB)])
    print('    img_sum_all/k0 norm:', np.round(ratio, 3).tolist(), 'range', round(ratio.min(), 3), round(ratio.max(), 3))
    print('    t_after_EE wmean (s):', np.round(tm, 2).tolist())
    print('    t_after_EE wsd   (s):', np.round(tsd, 2).tolist())
    print('    nav_z wmean (vox):', np.round(zm, 2).tolist())
    print('    n_hard:', [int((hard == b).sum()) for b in range(NB)])


if __name__ == '__main__':
    for kk in sys.argv[1:]:
        try:
            main(kk)
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f'{kk}: FAILED {e!r}')

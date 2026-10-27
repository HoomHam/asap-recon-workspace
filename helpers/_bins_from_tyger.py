"""Scratch (2026-10-26): per-bin tables from Tyger's EXACT SIGNAL membership (+ auto exclusions), replacing the
re-derived <key>_bins.csv for the b44d_sig s/ images.

Inputs: <key>_ilvbin_tyger.csv (from _tyger_bins_local.py), <key>_ee_tyger.csv, <key>_ilv_raw.csv (raw |k0| per gas
interleave, _bin_timing_rawk0.py), <key>_ilv.csv (nav_volume, nav_z_vox per interleave, _bin_timing_export.py),
AIkill_Dynamic_b44d_sig/<key>/s/recon.mat (gas_phase) for the flatness test.

Output <key>_bins_tyger.csv per bin: n_hard_kept, w_eff_kept, n_hard_excluded, t_after_EE_wmean_s, t_after_EE_wsd_s,
t_after_EE_med_s, phase_wmean, nav_z_vox_wmean, nav_volume_wmean, k0_gas_raw_wmean, k0_diss_raw_wmean,
img_sum_all, img_sum_all_over_k0 (normalised to its mean).
Kernel weights = recon.py cudarecon: exp(-d^2/2), cyclic, 0 for excluded or unbinned interleaves.
"""
import csv
import os
import sys

import numpy as np
import scipy.io as sio

EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon'
OUT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'bin_timing_2026-10-26')
NB = 16


def rd(p):
    return list(csv.DictReader(open(p)))


def col(rows, name):
    return np.array([float(r[name]) for r in rows])


def pick(rows, *cands):
    for c in cands:
        if c in rows[0]:
            return col(rows, c)
    raise KeyError(f'{cands} not in {list(rows[0])}')


def main(key):
    ty = rd(f'{OUT}/{key}_ilvbin_tyger.csv')
    t = col(ty, 't_s'); bins = col(ty, 'ilvbin_frac'); excl = col(ty, 'excluded').astype(int)
    ee = np.loadtxt(f'{OUT}/{key}_ee_tyger.csv', skiprows=1)
    raw = rd(f'{OUT}/{key}_ilv_raw.csv')
    nav = rd(f'{OUT}/{key}_ilv.csv')
    n = len(t)
    assert len(raw) == n and len(nav) == n, (n, len(raw), len(nav))
    k0g = pick(raw, 'k0_gas_sum8', 'k0_gas'); k0d = pick(raw, 'k0_diss_sum8', 'k0_diss')
    navv = col(nav, 'nav_volume'); navz = col(nav, 'nav_z_vox')
    # my previous memberships, for the comparison
    sig_fid = col(nav, 'sigbin_frac'); sig_raw = pick(raw, 'sigbin_frac')
    ok = (bins >= 0) & (excl == 0)
    c = bins * NB
    w = np.zeros((n, NB))
    for ib in range(NB):
        d = np.minimum.reduce([np.abs(ib - c), np.abs(ib + NB - c), np.abs(ib - NB - c)])
        w[:, ib] = np.where(ok, np.exp(-d * d / 2.0), 0.0)
    hard = np.where(bins >= 0, np.floor(c + 0.5).astype(int) % NB, -1)
    k = np.searchsorted(ee, t, side='right') - 1
    inside = (k >= 0) & (k < len(ee) - 1)
    tee = np.full(n, np.nan); tee[inside] = t[inside] - ee[k[inside]]
    # real time per bin: the kernel of bin b also weights interleaves on the far side of the EE wrap (phase just
    # below 1 for bin 0/1, just above 0 for bin 15). For those, use time relative to the NEXT EE (negative), so
    # each bin's weighted mean is the time of its own phase centre relative to the EE nearest to that centre.
    # Anchor of bin ib = the EE preceding the bin's own phase centre. A member whose phase sits > 8 bins AHEAD of
    # the centre (c - ib > 8, e.g. c = 15.7 for bin 1) is in the previous breath just before the anchor EE: its
    # time relative to the anchor = -(time to its next EE). A member > 8 bins BEHIND (c - ib < -8, e.g. c = 0.3 for
    # bin 15) is in the next breath: anchor-relative time = its tee + its previous period. Else: tee.
    nxt = np.full(n, np.nan); nxt[inside] = ee[k[inside] + 1] - t[inside]
    prevper = np.full(n, np.nan); kk = inside & (k >= 1); prevper[kk] = ee[k[kk]] - ee[k[kk] - 1]
    tee_fold = np.zeros((n, NB))
    for ib in range(NB):
        rawd = c - ib
        tee_fold[:, ib] = np.where(rawd > NB / 2, -nxt, np.where(rawd < -NB / 2, tee + prevper, tee))

    def wstat(x, ww):
        m = np.isfinite(x); ww = ww * m[:, None]; s = ww.sum(0) + 1e-12
        mean = (ww * np.nan_to_num(x)[:, None]).sum(0) / s
        sd = np.sqrt(np.maximum((ww * (np.nan_to_num(x)[:, None] - mean) ** 2).sum(0) / s, 0))
        return mean, sd

    def wstat2(x2, ww):
        m = np.isfinite(x2); ww = ww * m; s = ww.sum(0) + 1e-12
        return (ww * np.nan_to_num(x2)).sum(0) / s
    g = sio.loadmat(f'{EXT}/AIkill_Dynamic_b44d_sig/{key}/s/recon.mat', variable_names=['gas_phase'])['gas_phase'].astype(float)
    img_all = g.sum((1, 2, 3))
    tm, tsd = wstat(tee, w); tfm = wstat2(tee_fold, w)
    phm, _ = wstat(bins, w); zm, _ = wstat(navz, w); vm, _ = wstat(navv, w)
    kg, _ = wstat(k0g, w); kd, _ = wstat(k0d, w)
    ratio = img_all / kg; ratio /= ratio.mean()
    with open(f'{OUT}/{key}_bins_tyger.csv', 'w') as f:
        f.write('bin,n_hard_kept,w_eff_kept,n_hard_excluded,t_after_EE_wmean_s,t_after_EE_wsd_s,t_after_anchorEE_wmean_s,'
                't_after_EE_med_s,phase_wmean,nav_z_vox_wmean,nav_volume_wmean,k0_gas_raw_wmean,k0_diss_raw_wmean,'
                'img_sum_all,img_sum_all_over_k0_norm\n')
        for b in range(NB):
            hk = (hard == b) & ok; hx = (hard == b) & (excl == 1)
            med = np.nanmedian(tee[hk]) if hk.any() else np.nan
            f.write(f'{b},{hk.sum()},{w[:, b].sum():.1f},{hx.sum()},{tm[b]:.3f},{tsd[b]:.3f},{tfm[b]:.3f},{med:.3f},'
                    f'{phm[b]:.4f},{zm[b]:.3f},{vm[b]:.4f},{kg[b]:.6g},{kd[b]:.6g},{img_all[b]:.6g},{ratio[b]:.4f}\n')
    # membership comparison with my earlier derivations (binned interleaves only)
    def cmp(mine, name):
        both = (bins >= 0) & (mine >= 0)
        dd = np.abs(bins[both] - mine[both]); dd = np.minimum(dd, 1 - dd) * NB
        print(f'    vs {name}: binned tyger {int((bins >= 0).sum())} mine {int((mine >= 0).sum())} both {int(both.sum())}; '
              f'|dbin| > 0.5: {int((dd > 0.5).sum())} ({100 * (dd > 0.5).mean():.1f} %), max {dd.max():.2f}')
    print(f'== {key}: EE {len(ee)}, excluded ilv {excl.sum()} (of which binned {int(((bins >= 0) & (excl == 1)).sum())})')
    cmp(sig_fid, 'fid-affine bins (bins.csv)'); cmp(sig_raw, 'raw-k0 bins (bins_raw.csv)')
    print('    img_sum_all/k0 (tyger membership, kept only) norm:', np.round(ratio, 3).tolist(), 'range', round(ratio.min(), 3), round(ratio.max(), 3))
    print('    t_after_EE wmean (s):', np.round(tm, 2).tolist())
    print('    t_after_anchorEE wmean (s):', np.round(tfm, 2).tolist())
    print('    n_hard kept:', [int(((hard == b) & ok).sum()) for b in range(NB)])
    # per 832-interleave set: k0 average vs scan mean (the auto-exclude statistic, raw.py:437-440)
    nset = int(np.ceil(n / 832))
    avg = np.array([k0g[s * 832:(s + 1) * 832].sum() for s in range(nset)])
    print('    per-set k0 sum / mean (auto-exclude drops < 0.5):', np.round(avg / avg.mean(), 2).tolist())


if __name__ == '__main__':
    for k in sys.argv[1:]:
        try:
            main(k)
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f'{k}: FAILED {e!r}')

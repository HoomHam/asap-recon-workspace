"""Scratch (2026-10-26, PCA Step 1): RAW |k0| per gas / dissolved line straight from d/input.mrd (no affine
scan normalisation), plus per-breath A / B / C with C = navigator maximum. Same layout as _bin_timing_export.py,
files suffixed _raw. Needs the `mrd` package -> run with conda base python:
    /opt/homebrew/Caskroom/miniforge/base/bin/python3 helpers/_bin_timing_rawk0.py <study> ...

k0 definitions (per line, summed over channels, samples AFTER killpts = 2, i.e. Steve's raw.py:461 after
tyger_recon's trim; raw ADC units, no noise normalisation):
  k0_sum8   = sum |s[2:10]|          (Steve's SIGNAL navigator quantity; Rician-biased at low signal)
  k0_coh8   = |sum s[2:10]| / 8      (coherent mean of the same 8 samples; unbiased at low signal)
  k0_s2     = |s[2]|                 (first kept sample)
Per breath (signal EE -> next EE, EE from k0_sum8 with Steve's +-50 rule): trough A, peak B (k0 max), C = navigator
maximum inside the breath (nav z interpolated per interleave from the b44d d/recon.mat frames), k0 at C, nav z at
A / B / C, and the dissolved k0 at A / B / C. Pneumotach: per-cycle tidal amplitude trend (its own clock) in the
summary line, for the dose-drift question.
"""
import os
import sys
import numpy as np
import scipy.io as sio
import mrd

sys.path.insert(0, os.path.dirname(__file__))
from _bin_timing import steve_bin, steve_ee_minima, NBINS, EXT
from _bin_timing_export import kernel_w, load_nav, OUT

KILLPTS = 2


def read_dyn(study):
    p = f'{EXT}/AIkill_Dynamic/{study}/d/input.mrd'
    pneumo = None
    with mrd.BinaryMrdReader(p) as r:
        hdr = r.read_header()
        for it in r.read_data():
            if isinstance(it, mrd.StreamItem.NdArrayComplexFloat) and it.value.meta.get('dynamic_acquisition'):
                dyn = np.asarray(it.value.data)          # (ch, samples, lines)
            elif isinstance(it, (mrd.StreamItem.NdArrayFloat, mrd.StreamItem.NdArrayDouble)) and it.value.meta.get('pneumotach'):
                pneumo = np.asarray(it.value.data, dtype=float)
    ul = {q.name: q.value for q in hdr.user_parameters.user_parameter_long} if hdr.user_parameters else {}
    ud = {q.name: q.value for q in hdr.user_parameters.user_parameter_double} if hdr.user_parameters else {}
    return dyn, pneumo, ul, ud


def main(study):
    dyn, pneumo, ul, ud = read_dyn(study)
    kp = int(ul.get('killpts', KILLPTS))
    TR = float(ud.get('TR', 0.0172))
    s8 = dyn[:, kp:kp + 8, :]
    k0_sum8 = np.sum(np.abs(s8), axis=(0, 1)).astype(float)
    k0_coh8 = np.sum(np.abs(np.sum(s8, axis=1)), axis=0).astype(float) / 8.0
    k0_s2 = np.sum(np.abs(dyn[:, kp, :]), axis=0).astype(float)
    # noise floor of the same estimators from the last 8 readout samples (high k, ~noise): per-line values
    tail = dyn[:, -8:, :]
    nf_sum8 = np.sum(np.abs(tail), axis=(0, 1)).astype(float)
    nf_coh8 = np.sum(np.abs(np.sum(tail, axis=1)), axis=0).astype(float) / 8.0
    tree, nav = load_nav(study)
    ngas = len(nav['nav_volume'])
    nlines = dyn.shape[2]
    numspec = 2 * ((nlines - 2 * ngas) // 2)
    # nav arrays can be longer than this input.mrd's gas lines (nav taken from a _merged / other-MID recon of the
    # same session): fall back to the protocol's spectral block (780 v3 / 500 v2) and truncate the nav arrays
    if numspec < 0 or numspec > 2000:
        numspec = 500 if abs(TR - 0.0149) < 5e-4 else 780
        ngas_new = (nlines - numspec) // 2
        print(f'{study}: nav length {ngas} > gas lines {ngas_new} (numspec {numspec} assumed); nav truncated')
        ngas = ngas_new
        nav = {k: v[:ngas] if k in ('nav_volume', 'nav_ilvtime') else v for k, v in nav.items()}
    g = slice(numspec, numspec + 2 * ngas, 2)
    d = slice(numspec + 1, numspec + 2 * ngas, 2)
    gas, dis = k0_sum8[g], k0_sum8[d]
    gasc, disc = k0_coh8[g], k0_coh8[d]
    gas2, dis2 = k0_s2[g], k0_s2[d]
    t = nav['nav_ilvtime'] if len(nav['nav_ilvtime']) == ngas else np.arange(ngas) * TR
    ee = steve_ee_minima(gas, t, 50)
    from scipy.signal import savgol_filter
    gsm = savgol_filter(gas, 11, 2)     # trough value without the min-of-noise pick: local quadratic over 11 ilv (0.19 s)
    ee_t = t[ee]
    sig = np.full(ngas, -1.0)
    kb = np.full(ngas, -1)
    for j in range(ngas):
        if len(ee_t) and ee_t[0] <= t[j] <= ee_t[-1]:
            k = min(np.searchsorted(ee_t, t[j], side='right') - 1, len(ee_t) - 2)
            sig[j] = (t[j] - ee_t[k]) / (ee_t[k + 1] - ee_t[k])
            kb[j] = k
    dia = steve_bin(nav['nav_volume'])
    fz, ft = nav['nav_diaphragm_z'], nav['nav_time']
    good = np.isfinite(fz)
    navz = np.interp(t, ft[good], fz[good], left=np.nan, right=np.nan)
    navz[np.isnan(nav['nav_volume'])] = np.nan
    t_ee = np.where(kb >= 0, t - ee_t[np.clip(kb, 0, None)], np.nan)
    os.makedirs(OUT, exist_ok=True)
    hdr = 'ilv,t_s,k0_gas_sum8,k0_diss_sum8,k0_gas_coh8,k0_diss_coh8,k0_gas_s2,k0_diss_s2,sigbin_frac,diabin_frac,nav_volume,nav_z_vox,breath,t_since_EE_s'
    np.savetxt(f'{OUT}/{study}_ilv_raw.csv',
               np.column_stack([np.arange(ngas), t, gas, dis, gasc, disc, gas2, dis2, sig, dia, nav['nav_volume'], navz, kb, t_ee]),
               delimiter=',', header=hdr, comments='',
               fmt='%d,%.4f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.5f,%.5f,%.5f,%.3f,%d,%.3f')
    # per-bin kernel-weighted raw k0
    with open(f'{OUT}/{study}_bins_raw.csv', 'w') as f:
        f.write('binning,bin,n_hard,w_eff,t_after_EE_med_s,nav_z_vox_wmean,k0_gas_sum8_wmean,k0_diss_sum8_wmean,k0_gas_coh8_wmean,k0_diss_coh8_wmean\n')
        for name, bf in (('SIGNAL', sig), ('DIAPHRAGM', dia)):
            w = kernel_w(bf)
            ib = np.where(bf >= 0, np.floor(bf * NBINS + 0.5).astype(int) % NBINS, -1)
            okz = np.isfinite(navz)
            for b in range(NBINS):
                wb = w[:, b]
                m = ib == b
                wm = lambda x, msk=None: np.sum(wb * (msk if msk is not None else 1) * np.nan_to_num(x)) / max(np.sum(wb * (msk if msk is not None else 1)), 1e-9)
                f.write('%s,%d,%d,%.1f,%.3f,%.3f,%.6g,%.6g,%.6g,%.6g\n' % (
                    name, b, int(m.sum()), wb.sum(), np.nanmedian(t_ee[m]) if m.any() else np.nan,
                    wm(navz, okz), wm(gas), wm(dis), wm(gasc), wm(disc)))
    # per-breath A / B / C
    with open(f'{OUT}/{study}_breaths_raw.csv', 'w') as f:
        f.write('breath,t_EE_s,period_s,k0_A,k0_B,t_B_after_EE_s,k0_C,t_C_after_EE_s,k0_nextEE,B_over_A,C_over_B,nextEE_over_C,'
                'navz_A,navz_B,navz_C,navz_min_in_breath,k0diss_A,k0diss_B,k0diss_C,k0diss_max,t_dissmax_after_EE_s,'
                'k0coh_A,k0coh_B,k0coh_C,k0coh_nextEE,B_over_A_coh,nextEE_over_C_coh,k0coh_A5,k0coh_nextEE5,nextEE_over_C_coh5,'
                't_C95_start_s,t_C95_mid_s,t_C95_end_s,k0_C95_start,k0_C95_mid,k0_C95_end,C95mid_over_B,k0_A_sm,k0_nextEE_sm,B_over_A_sm,nextEE_over_C_sm,k0_B_sm,k0_C_sm,B_over_A_smsm,C_over_B_smsm,nextEE_over_C_smsm\n')
        for k in range(len(ee) - 1):
            a, b = ee[k], ee[k + 1]
            seg = gas[a:b + 1]
            ip = a + int(np.argmax(seg))
            zseg = navz[a:b + 1]
            if np.all(np.isnan(zseg)):
                ic = -1
            else:
                ic = a + int(np.nanargmax(zseg))
            dseg = dis[a:b + 1]
            idm = a + int(np.argmax(dseg))
            # 95 % plateau of the navigator inside the breath: first / centre / last interleave with z >= min + 0.95*(max-min)
            if ic >= 0:
                zr = np.nanmin(zseg), np.nanmax(zseg)
                on = np.where(zseg >= zr[0] + 0.95 * (zr[1] - zr[0]))[0]
                i0, i1 = a + int(on[0]), a + int(on[-1]); im = (i0 + i1) // 2
            else:
                i0 = i1 = im = -1
            c5 = lambda i: np.mean(gasc[max(i - 2, 0):i + 3])
            vals = (k, t[a], t[b] - t[a], gas[a], gas[ip], t[ip] - t[a],
                    gas[ic] if ic >= 0 else np.nan, t[ic] - t[a] if ic >= 0 else np.nan, gas[b],
                    gas[ip] / gas[a], (gas[ic] / gas[ip]) if ic >= 0 else np.nan, (gas[b] / gas[ic]) if ic >= 0 else np.nan,
                    navz[a], navz[ip], navz[ic] if ic >= 0 else np.nan, np.nanmin(zseg) if ic >= 0 else np.nan,
                    dis[a], dis[ip], dis[ic] if ic >= 0 else np.nan, dis[idm], t[idm] - t[a],
                    gasc[a], gasc[ip], gasc[ic] if ic >= 0 else np.nan, gasc[b], gasc[ip] / gasc[a],
                    (gasc[b] / gasc[ic]) if ic >= 0 else np.nan, c5(a), c5(b), (c5(b) / gasc[ic]) if ic >= 0 else np.nan,
                    t[i0] - t[a] if ic >= 0 else np.nan, t[im] - t[a] if ic >= 0 else np.nan, t[i1] - t[a] if ic >= 0 else np.nan,
                    gas[i0] if ic >= 0 else np.nan, gas[im] if ic >= 0 else np.nan, gas[i1] if ic >= 0 else np.nan,
                    (gas[im] / gas[ip]) if ic >= 0 else np.nan,
                    gsm[a], gsm[b], gas[ip] / gsm[a], (gsm[b] / gas[ic]) if ic >= 0 else np.nan,
                    gsm[ip], gsm[ic] if ic >= 0 else np.nan, gsm[ip] / gsm[a], (gsm[ic] / gsm[ip]) if ic >= 0 else np.nan, (gsm[b] / gsm[ic]) if ic >= 0 else np.nan)
            f.write(','.join('%d' % v if i == 0 else '%.6g' % v for i, v in enumerate(vals)) + '\n')
    per = np.diff(ee_t)
    r = np.array([gas[ee[k] + int(np.argmax(gas[ee[k]:ee[k + 1] + 1]))] / gas[ee[k]] for k in range(len(ee) - 1)])
    nfs, nfc = np.median(nf_sum8[g]), np.median(nf_coh8[g])
    trA = gas[ee[:-1]]; trAc = gasc[ee[:-1]]
    msg0 = (f'   k0 noise floor (tail samples): sum8 {nfs:.3g} -> trough/floor {np.median(trA)/nfs:.1f}; coh8 {nfc:.3g} -> '
            f'{np.median(trAc)/nfc:.1f}; sum8/coh8 at A {np.median(trA/trAc):.3f}, at B '
            f'{np.median([gas[ee[k] + int(np.argmax(gas[ee[k]:ee[k + 1] + 1]))] / gasc[ee[k] + int(np.argmax(gas[ee[k]:ee[k + 1] + 1]))] for k in range(len(ee) - 1)]):.3f}')
    print(msg0)
    msg = (f'{study}: nav {tree}; {nlines} lines = {numspec} spec + 2x{ngas}; TR {TR*1e3:.1f} ms; {len(ee)} EE, period '
           f'{np.median(per):.2f} s; RAW k0 B/A median {np.median(r):.2f} (IQR {np.percentile(r,25):.2f}-{np.percentile(r,75):.2f}); '
           f'diss/gas sum8 at A {np.median(dis[ee]/gas[ee]):.3f}, coh8 {np.median(disc[ee]/gasc[ee]):.3f}')
    if pneumo is not None and pneumo.shape[0] == 2:
        pt, pp = pneumo[0], pneumo[1]
        try:
            from scipy.signal import savgol_filter, find_peaks
            P = savgol_filter(pp, 51, 2) if len(pp) > 51 else pp
            vol = np.cumsum(P)
            dt = np.median(np.diff(pt))
            # remove slow drift with a wide median, then per-cycle amplitude
            from scipy.ndimage import median_filter
            vol = vol - median_filter(vol, size=int(max(3, 20 / dt)) | 1)
            pk, _ = find_peaks(vol, distance=int(2.0 / dt))
            tr_, _ = find_peaks(-vol, distance=int(2.0 / dt))
            amp = np.array([vol[p] - vol[tr_[np.searchsorted(tr_, p) - 1]] for p in pk if np.searchsorted(tr_, p) > 0])
            tp = np.array([pt[p] for p in pk if np.searchsorted(tr_, p) > 0])
            if len(amp) > 5:
                sl = np.polyfit(tp, amp / np.median(amp), 1)[0] * 60 * 100
                msg += f'; pneumo tidal amplitude trend {sl:+.1f} %/min over {tp[0]:.0f}-{tp[-1]:.0f} s ({len(amp)} cycles, pneumo clock)'
        except Exception as e:
            msg += f'; pneumo trend failed ({e})'
    print(msg)


if __name__ == '__main__':
    for s in sys.argv[1:]:
        try:
            main(s)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f'{s}: FAILED {e}')

"""Scratch (2026-10-27): verify a MIX (m/) or PNEUMOTACH (p/) recon from the 508b329 image.

Checks per <key>/<sub>: recon.mat keys incl. the bin export; exported ilvbin vs an offline re-derivation
(tyger_recon._mix_binning on the fork loader + the d/ navigator volume) -> max |dbin| (expect 0.000 for m);
hard counts per bin; Sigma_FOV gas_phase(bin) / kernel-weighted raw k0 per bin (weights from the exported ilvbin with
the exported kernel settings) flat +-2 %; corner-noise sigma per bin; mix_tau_I / n_insp.
usage: cd <root> && python3 workspace/helpers/_mix_check.py <key> [m|p] ...
"""
import csv
import os
import sys

import numpy as np
import scipy.io as sio

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.dirname(__file__))
EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon'
BT = os.path.join(ROOT, 'workspace', 'outputs', 'bin_timing_2026-10-26')
NB = 16


def weights(bins, sigma, pou):
    ok = bins >= 0; c = bins * NB
    w = np.zeros((len(bins), NB))
    for ib in range(NB):
        d = np.minimum(np.minimum(np.abs(ib - c), np.abs(ib + NB - c)), np.abs(ib - NB - c))
        w[:, ib] = np.where(ok, np.exp(-d * d / (2 * sigma ** 2)), 0.0)
    if pou:
        w = w / np.maximum(w.sum(1, keepdims=True), 1e-12)
    return w


def main(key, sub):
    p = f'{EXT}/AIkill_Dynamic_mix/{key}/{sub}/recon.mat'
    m = sio.loadmat(p)
    keys = sorted(k for k in m if not k.startswith('__'))
    print(f'== {key}/{sub}: keys {keys}')
    g = m['gas_phase'].astype(np.float64)
    bins = m['ilvbin'].ravel().astype(np.float64)
    mode = int(m['bin_mode'].ravel()[0]) if 'bin_mode' in m else -1
    sigma = 0.5 if sub == 'm' else 1.0
    pou = 1 if sub == 'm' else 0
    print(f'   bin_mode {mode}; binned {int((bins >= 0).sum())}/{len(bins)}; ee_times n {m["ee_times"].size}; '
          + (f'tau_I {float(m["mix_tau_I"].ravel()[0]):.3f} s n_insp {int(m["mix_n_insp"].ravel()[0])} tB_med {float(m["mix_tB_med"].ravel()[0]):.2f} T_med {float(m["mix_T_med"].ravel()[0]):.2f}' if 'mix_tau_I' in m else ''))
    hard = np.where(bins >= 0, np.floor(bins * NB + 0.5).astype(int) % NB, -1)
    cnt = np.array([int((hard == b).sum()) for b in range(NB)])
    print(f'   hard counts {cnt.tolist()} cv {cnt.std() / cnt.mean():.3f} min/max {cnt.min() / cnt.max():.2f}')
    w = weights(bins, sigma, pou)
    weff = w.sum(0)
    print(f'   exported bin_weff == recomputed: {np.allclose(weff, m["bin_weff"].ravel(), rtol=1e-3)}')
    # flatness vs raw k0
    raw = list(csv.DictReader(open(f'{BT}/{key}_ilv_raw.csv')))
    k0 = np.array([float(r['k0_gas_sum8']) for r in raw])
    assert len(k0) == len(bins), (len(k0), len(bins))
    km = (w * k0[:, None]).sum(0) / np.maximum(weff, 1e-9)
    tot = g.sum((1, 2, 3)); ratio = tot / km; ratio /= ratio.mean()
    print(f'   Sigma_FOV gas / k0_wmean (norm): {np.round(ratio, 3).tolist()} range {ratio.min():.3f} {ratio.max():.3f}')
    sig = g[:, :8, :8, :8].std((1, 2, 3))
    print(f'   corner sigma per bin: {np.round(sig, 0).tolist()}')
    if sub == 'm':
        # offline re-derivation of the MIX bins (fork loader + d/ navigator volume)
        import _tyger_bins_local as L
        import tyger_recon as T
        from gtypes import bintype
        gg, g_raw, _ = L.load(f'{EXT}/AIkill_Dynamic_mix/{key}/m/input.mrd')
        g_raw.ilvvol[bintype.DIAPHRAGM] = m['nav_volume'].ravel().astype(float)
        bo, out = T._mix_binning(gg, g_raw)
        both = (bo >= 0) & (bins >= 0)
        print(f'   offline MIX re-derivation: binned {int((bo >= 0).sum())} vs exported {int((bins >= 0).sum())}; '
              f'max |dbin| {np.abs(bo[both] - bins[both]).max() * NB:.4f}; tau {float(out["mix_tau_I"][0]):.3f} n_insp {int(out["mix_n_insp"][0])}')
    # time per bin from ee_times
    t = m['ilvtime'].ravel().astype(float); ee = m['ee_times'].ravel().astype(float)
    k = np.searchsorted(ee, t, side='right') - 1; inside = (k >= 0) & (k < len(ee) - 1)
    tee = np.full(len(t), np.nan); tee[inside] = t[inside] - ee[k[inside]]
    tm = np.array([np.nanmedian(tee[hard == b]) if (hard == b).any() else np.nan for b in range(NB)])
    print(f'   median t after EE per bin (s): {np.round(tm, 2).tolist()}')
    navz = None
    if 'nav_diaphragm_z' in m:
        fz, ft = m['nav_diaphragm_z'].ravel().astype(float), m['nav_time'].ravel().astype(float)
        good = np.isfinite(fz)
        navz = np.interp(t, ft[good], fz[good], left=np.nan, right=np.nan)
        zsd = np.array([np.nanstd(navz[hard == b]) if (hard == b).sum() > 2 else np.nan for b in range(NB)])
        print(f'   within-bin nav z sd (vox): {np.round(zsd, 2).tolist()}')


if __name__ == '__main__':
    args = sys.argv[1:]
    keys = [a for a in args if '_' in a]; subs = [a for a in args if a in ('m', 'p')] or ['m']
    for k in keys:
        for s in subs:
            try:
                main(k, s)
            except Exception as e:
                import traceback; traceback.print_exc(); print(f'{k}/{s}: FAILED {e!r}')

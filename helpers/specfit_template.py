#!/usr/bin/env python3
"""
specfit_template.py -- template-basis LINEAR solve for the RBC/TP split inputs, and a stability study
against the nonlinear specfit (M3 lumped) on every cached cal block.

Why: the M3 two-membrane decomposition (16 free parameters) is multimodal on RBC-weak blocks: one early
sample or a different line subset flips the minimum (042DR: 61-91 deg, F68/F69). The split only needs one
complex coefficient per line (amplitude + phase of the RBC and of the LUMPED membrane phasor at the image
k0). So: freeze every line SHAPE at the cohort value (per protocol, from the stable specfit M3 fits) and
solve for the complex coefficients by linear least squares -- no initial guess, no local minima. Small
frequency offsets (subject shim) are absorbed by a coarse grid over the RBC and membrane shifts (variable
projection); every grid point is a global linear solve.

Model (t from ADC sample 0 = TE, samples 0-1 dropped as in specfit / killpts=2):
    y(t) = c_g T_g(t) + c_r T_r(t) + c_m T_m(t)
    T_g = e^{2pi i f_g t} e^{-pi w_g t}                       f_g from specfit.locate_gas, w_g = 26 Hz (cohort)
    T_r = e^{2pi i (f_g + (p_r + d_r) h) t} e^{-pi w_r t}     p_r, w_r = cohort RBC shift / FWHM (ppm -> Hz by h)
    T_m = e^{2pi i (f_g + (p_1 + d_m) h) t} e^{-pi w_1 t} + a21 e^{i psi} e^{2pi i (f_g + (p_2 + d_m) h) t} e^{-pi w_2 t}
          (cohort mem1/mem2 shifts, FWHMs, amplitude ratio a2/a1 and relative phase psi = -(phi1 - phi2))
    d_r, d_m on a grid of +-1.5 ppm (0.25 ppm steps); pick the minimum residual.
Outputs at t_k0 = 10 us (killpts * dt_img after sample 0):
    dphi_k0 = angle(c_r T_r(t_k0) / (c_m T_m(t_k0)));  ratio_lumped = |c_r T_r(t_k0)| / |c_m T_m(t_k0)|
    ratio_scalar = |c_r T_r(t_k0)| / (|c_m| (|L1(t_k0)| + a21 |L2(t_k0)|));  F_lump = ratio_lumped / ratio_scalar

Study (per block, both methods): drop-2 vs drop-3, odd vs even reps, first vs second half of the post-arrival
reps -> |delta dphi_k0|; template vs specfit agreement; where they disagree and why (SNR, scalar ratio).
Inputs : XeCS cache /Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/calspec/cache/<tag>.npz (cal, dwell, row)
         notes/calspec_package_2026-09-16/{resplit_inputs_2026-09-24.csv, specfit.py}
Outputs: outputs/resplit_2026-09-24/template_study.csv, fig/template_study.png
usage: python specfit_template.py [--no-specfit-halves]
"""
import argparse
import csv
import glob
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WS = HERE.parent
PKG = WS / 'notes' / 'calspec_package_2026-09-16'
OUT = WS / 'outputs' / 'resplit_2026-09-24'
CACHE = Path('/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/calspec/cache')
sys.path.insert(0, str(PKG))
import specfit  # noqa: E402

HZ = 17.61
T_K0 = 10e-6
W_GAS = 26.0                      # Hz, cohort median gas FWHM (subjects_pooled.csv)
GRID = np.arange(-1.5, 1.51, 0.25)   # ppm offsets for RBC and membrane
wrap = lambda d: (d + 180.0) % 360.0 - 180.0


def load_block(tag):
    z = np.load(CACHE / f'{tag}.npz', allow_pickle=True)
    if 'cal' not in z:
        return None
    row = json.loads(str(z['row']))
    cal = z['cal']                                    # (N, nlin, nrep)
    arr = int(row.get('arrival_rep', 0) or 0)
    post = cal[:, :, max(arr, 0):]                    # post-arrival reps
    return dict(cal=post, dwell=float(z['dwell']), nlin=cal.shape[1], nrep_post=post.shape[2], arr=arr)


def cohort_templates(rows):
    """Per-protocol cohort line shapes from the stable, valid, M3 specfit rows."""
    T = {}
    for seq in ('v3', 'v2'):
        sel = [r for r in rows if r['seq'].endswith(seq) and r['valid'] == 'yes' and r['stable'] == 'yes'
               and r['model_used'] == 'M3']
        if len(sel) < 5:
            sel = [r for r in rows if r['seq'].endswith(seq) and r['valid'] == 'yes' and r['model_used'] == 'M3']
        f = lambda k: float(np.median([float(r[k]) for r in sel]))
        T[seq] = dict(n=len(sel), p_r=f('rbc_ppm'), w_r=f('rbc_fwhm_ppm') * HZ, p_1=f('mem1_ppm'), w_1=f('mem1_fwhm_ppm') * HZ,
                      p_2=f('mem2_ppm'), w_2=f('mem2_fwhm_ppm') * HZ, a21=1.0 / f('a1_over_a2'),
                      psi=-np.radians(f('mem12_dphi_deg')))
    return T


def basis(t, f_g, tp, d_r, d_m):
    Tg = np.exp(2j * np.pi * f_g * t - np.pi * W_GAS * t)
    Tr = np.exp(2j * np.pi * (f_g + (tp['p_r'] + d_r) * HZ) * t - np.pi * tp['w_r'] * t)
    L1 = np.exp(2j * np.pi * (f_g + (tp['p_1'] + d_m) * HZ) * t - np.pi * tp['w_1'] * t)
    L2 = np.exp(2j * np.pi * (f_g + (tp['p_2'] + d_m) * HZ) * t - np.pi * tp['w_2'] * t)
    Tm = L1 + tp['a21'] * np.exp(1j * tp['psi']) * L2
    return Tg, Tr, Tm, L1, L2


ADAPTIVE = False                              # --adaptive: also grid the line widths and the mem2/mem1 amplitude ratio
GRID_A = np.arange(-2.0, 2.01, 0.25)          # ppm shifts (adaptive)
WSCALE = (0.7, 0.85, 1.0, 1.2, 1.4)           # width scale factors for RBC and (jointly) mem1/mem2
A21 = (0.8, 1.2, 1.6, 2.0, 2.5)               # mem2/mem1 amplitude ratio


def template_solve(y, t, tp, f_g=None, dwell=None):
    """Linear solve on the shift grid (and, with --adaptive, the width / a21 grid): variable projection,
    every grid point is a global least-squares solve for the three complex coefficients."""
    if f_g is None:
        f_g = specfit.locate_gas(y, dwell, HZ)
    best = None
    shifts = GRID_A if ADAPTIVE else GRID
    shapes = [(wr, wm, a) for wr in WSCALE for wm in WSCALE for a in A21] if ADAPTIVE else [(1.0, 1.0, tp['a21'])]
    for wr, wm, a21 in shapes:
        tq = dict(tp, w_r=tp['w_r'] * wr, w_1=tp['w_1'] * wm, w_2=tp['w_2'] * wm, a21=a21)
        for d_r in shifts:
            for d_m in shifts:
                Tg, Tr, Tm, L1, L2 = basis(t, f_g, tq, d_r, d_m)
                B = np.stack([Tg, Tr, Tm], axis=1)
                c, *_ = np.linalg.lstsq(B, y, rcond=None)
                res = float(np.linalg.norm(B @ c - y))
                if best is None or res < best['res']:
                    best = dict(res=res, c=c, d_r=d_r, d_m=d_m, f_g=f_g, tq=tq, wr=wr, wm=wm)
    c, d_r, d_m, tp = best['c'], best['d_r'], best['d_m'], best['tq']
    tk = np.array([T_K0])
    Tg, Tr, Tm, L1, L2 = basis(tk, f_g, tp, d_r, d_m)
    rbc, mem = c[1] * Tr[0], c[2] * Tm[0]
    scalar_mem = abs(c[2]) * (abs(L1[0]) + tp['a21'] * abs(L2[0]))
    return dict(dphi=float(wrap(np.degrees(np.angle(rbc / mem)))), ratio_lumped=float(abs(rbc) / abs(mem)),
                ratio_scalar=float(abs(rbc) / scalar_mem), F_lump=float(scalar_mem / abs(mem)),
                d_r=float(d_r), d_m=float(d_m), f_g=float(f_g), res_rel=float(best['res'] / np.linalg.norm(y)),
                gas_over_diss=float(abs(c[0]) / (abs(rbc) + abs(mem))), wr=float(best['wr']), wm=float(best['wm']), a21=float(tp['a21']))


def prep(cal3, drop):
    """(N, nlin, nrep) -> mean FID over lines and reps, trimmed; t from ADC sample 0."""
    y = cal3.reshape(cal3.shape[0], -1, order='F').mean(axis=1)
    idx = np.arange(y.size)
    return y[idx >= drop].astype(complex), idx[idx >= drop]


def run_specfit(cal3, dwell, drop=2):
    y, idx = prep(cal3, drop)
    s = specfit.fit_block(y, dwell, 620e-6, first_sample=int(idx[0]), t_k0_offset_s=T_K0, hz_per_ppm=HZ, drop_before=0)
    return dict(dphi=s['dphi_k0_deg'], ratio_lumped=s['ratio_lumped'], ratio_scalar=s['ratio_scalar'], model=s['model_used'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-specfit-halves', action='store_true')
    ap.add_argument('--adaptive', action='store_true', help='grid the line widths and mem2/mem1 ratio too')
    a = ap.parse_args()
    global ADAPTIVE
    ADAPTIVE = a.adaptive
    sfx = '_adaptive' if a.adaptive else ''
    rows = list(csv.DictReader(open(PKG / 'resplit_inputs_2026-09-24.csv')))
    tp = cohort_templates(rows)
    for seq, v in tp.items():
        print(f'template {seq} (n={v["n"]}): RBC {v["p_r"]:.2f} ppm / {v["w_r"] / HZ:.1f} ppm; mem1 {v["p_1"]:.2f}/{v["w_1"] / HZ:.1f}; '
              f'mem2 {v["p_2"]:.2f}/{v["w_2"] / HZ:.1f}; a2/a1 {v["a21"]:.2f}; psi {np.degrees(v["psi"]):.0f} deg')
    out = []
    for r in rows:
        blk = load_block(r['tag'])
        if blk is None:
            continue
        seq = 'v2' if r['seq'].endswith('v2') else 'v3'
        T = tp[seq]
        cal, dwell = blk['cal'], blk['dwell']
        rec = dict(session=r['session'], seq=seq, snr=float(r['snr_diss']), valid=r['valid'], stable=r['stable'],
                   model=r['model_used'], sf_dphi=float(r['dphi_k0_deg']), sf_ratio_l=float(r['ratio_lumped']),
                   sf_ratio_s=float(r['ratio_scalar']), sf_stab_drop=float(r['stab_dphi_deg']) if r['stab_dphi_deg'] else np.nan,
                   nrep=blk['nrep_post'])
        # template: full, drop 2 / drop 3, odd/even reps, first/second half
        y, idx = prep(cal, 2); t = idx * dwell
        full = template_solve(y, t, T, dwell=dwell)
        rec.update(tp_dphi=full['dphi'], tp_ratio_l=full['ratio_lumped'], tp_ratio_s=full['ratio_scalar'], tp_F=full['F_lump'],
                   tp_dr=full['d_r'], tp_dm=full['d_m'], tp_res=full['res_rel'], tp_wr=full['wr'], tp_wm=full['wm'], tp_a21=full['a21'])
        y3, idx3 = prep(cal, 3)
        d3 = template_solve(y3, idx3 * dwell, T, f_g=full['f_g'])
        rec['tp_stab_drop'] = abs(wrap(full['dphi'] - d3['dphi']))
        splits = {}
        nrep = cal.shape[2]
        if nrep >= 4:
            splits['oddeven'] = (cal[:, :, 0::2], cal[:, :, 1::2])
            splits['halves'] = (cal[:, :, :nrep // 2], cal[:, :, nrep // 2:])
        for name, (A, B) in splits.items():
            ya, ia = prep(A, 2); yb, ib = prep(B, 2)
            da = template_solve(ya, ia * dwell, T, f_g=full['f_g']); db = template_solve(yb, ib * dwell, T, f_g=full['f_g'])
            rec[f'tp_stab_{name}'] = abs(wrap(da['dphi'] - db['dphi']))
            if not a.no_specfit_halves:
                try:
                    sa, sb = run_specfit(A, dwell), run_specfit(B, dwell)
                    rec[f'sf_stab_{name}'] = abs(wrap(sa['dphi'] - sb['dphi']))
                except Exception as e:
                    rec[f'sf_stab_{name}'] = np.nan
        rec['tp_minus_sf'] = wrap(full['dphi'] - rec['sf_dphi'])
        out.append(rec)
        print(f"{rec['session']} {seq} snr {rec['snr']:5.1f} {rec['model']} stable={rec['stable']:3s} | specfit {rec['sf_dphi']:6.1f} "
              f"(drop {rec['sf_stab_drop']:4.1f} oe {rec.get('sf_stab_oddeven', np.nan):4.1f} hv {rec.get('sf_stab_halves', np.nan):4.1f}) | "
              f"template {full['dphi']:6.1f} (drop {rec['tp_stab_drop']:4.1f} oe {rec.get('tp_stab_oddeven', np.nan):4.1f} "
              f"hv {rec.get('tp_stab_halves', np.nan):4.1f}) ratio_l {full['ratio_lumped']:.3f}/{rec['sf_ratio_l']:.3f} "
              f"scalar {full['ratio_scalar']:.3f}/{rec['sf_ratio_s']:.3f} dr {full['d_r']:+.2f} dm {full['d_m']:+.2f} res {full['res_rel']:.3f}", flush=True)
    keys = list(out[0].keys())
    with open(OUT / f'template_study{sfx}.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in out:
            w.writerow({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items()})
    json.dump({k: {kk: float(vv) for kk, vv in v.items()} for k, v in tp.items()}, open(OUT / 'template_params.json', 'w'), indent=1)
    # summary
    st = lambda v: (np.nanmedian(v), np.nanpercentile(v, 90))
    g = lambda k, sel: np.array([r.get(k, np.nan) for r in out if sel(r)], float)
    good = lambda r: r['valid'] == 'yes' and r['stable'] == 'yes'
    bad = lambda r: not good(r)
    print('\n=== stable+valid rows (n=%d): template - specfit dphi median %+.1f MAD %.1f, |>10 deg| %d' % (
        sum(good(r) for r in out), np.nanmedian(g('tp_minus_sf', good)), np.nanmedian(np.abs(g('tp_minus_sf', good) - np.nanmedian(g('tp_minus_sf', good)))),
        int((np.abs(g('tp_minus_sf', good)) > 10).sum())))
    for k in ('stab_drop', 'stab_oddeven', 'stab_halves'):
        for m in ('sf', 'tp'):
            v = g(f'{m}_{k}', good); vb = g(f'{m}_{k}', bad)
            print(f'  {m}_{k:12s} stable rows median/p90 {st(v)[0]:.1f}/{st(v)[1]:.1f}   unstable-or-invalid rows median/p90 {st(vb)[0]:.1f}/{st(vb)[1]:.1f}')
    print('=== unstable/invalid rows:')
    for r in out:
        if bad(r):
            print(f"  {r['session']} snr {r['snr']:.0f} scalar {r['sf_ratio_s']:.3f}: specfit {r['sf_dphi']:.1f} (drop {r['sf_stab_drop']:.1f}, oe {r.get('sf_stab_oddeven', np.nan):.1f}) "
                  f"template {r['tp_dphi']:.1f} (drop {r['tp_stab_drop']:.1f}, oe {r.get('tp_stab_oddeven', np.nan):.1f}, hv {r.get('tp_stab_halves', np.nan):.1f})")
    # figure
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    col = ['tab:green' if good(r) else 'tab:red' for r in out]
    ax[0].scatter([r['sf_dphi'] for r in out], [r['tp_dphi'] for r in out], c=col, s=16); ax[0].plot([0, 120], [0, 120], 'k--', lw=.6)
    ax[0].set_xlabel('specfit M3-lumped dphi_k0 (deg)'); ax[0].set_ylabel('template linear dphi_k0 (deg)'); ax[0].set_title('angle: green = specfit stable+valid, red = flagged', fontsize=9)
    ax[1].scatter([r.get('sf_stab_oddeven', np.nan) for r in out], [r.get('tp_stab_oddeven', np.nan) for r in out], c=col, s=16); ax[1].plot([0, 40], [0, 40], 'k--', lw=.6)
    ax[1].set_xlabel('specfit |dphi odd - even reps| (deg)'); ax[1].set_ylabel('template |dphi odd - even| (deg)'); ax[1].set_title('stability under line-subset split', fontsize=9); ax[1].set_xlim(0, 40); ax[1].set_ylim(0, 40)
    ax[2].scatter([r['sf_ratio_s'] for r in out], [r.get('tp_stab_oddeven', np.nan) for r in out], c=col, s=16, label='template')
    ax[2].scatter([r['sf_ratio_s'] for r in out], [r.get('sf_stab_oddeven', np.nan) for r in out], facecolors='none', edgecolors='gray', s=16, label='specfit')
    ax[2].set_xlabel('RBC/mem scalar ratio (specfit)'); ax[2].set_ylabel('|dphi odd - even| (deg)'); ax[2].set_title('instability vs RBC strength', fontsize=9); ax[2].legend(fontsize=8); ax[2].set_ylim(0, 40)
    fig.suptitle('Template-basis linear solve vs specfit M3 (75 cal blocks, samples 0-1 dropped, t_k0 = 10 us)'); fig.tight_layout()
    p = OUT / 'fig' / f'template_study{sfx}.png'; fig.savefig(p, dpi=115); print('fig ->', p)


if __name__ == '__main__':
    main()

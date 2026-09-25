#!/usr/bin/env python
# ⚠ RETRACTED-DEPENDENCY [F74]: the imprint_index / static-pattern columns measured a weak statistic, not a smooth static phase field — the real cause is F79 (calcb b fine phase); grades A–D remain valid — see workspace/canon/facts.md
"""Where does the RBC/TM split fail, and how much can each session's separation be trusted?

Two questions from Hooman (2026-09-25, after the RBC/TM atlas):
  1. "I still see fake lung structures" — is the structure the static calcb phase reference
     (2steve/06: one b for all bins, phase measured only where |F_avg| > 10*noise, a global
     quadratic polynomial elsewhere)?  Test: the stored split is z = aRBC*exp(i*phi_k0) + aTP, so a
     lung voxel is PHYSICAL only if aRBC >= 0 and aTP >= 0 (phase inside the wedge [0, phi_k0]).
     Out-of-wedge voxels are pure error (noise or phase-reference error). If the error concentrates
     OUTSIDE calcb's measured-phase support (the polynomial region, at the moving rim), the static
     reference is imprinting; if it is uniform and matches the noise prediction, it is just SNR.
  2. "How much do we trust the separation?" — per voxel: phase noise sigma_phi = sigma_z/|z| (rad);
     the split resolves RBC vs TM only where sigma_phi << phi_k0. Per session: a grade from the
     angle conditioning, the fit stability, the resolvable-voxel fraction and the excess error.

Inputs : Ext AIkill_Dynamic/<key>/d/recon.mat (gas_phase), d/recon_resplit_xecs.mat (aRBC, aTP, ...),
         outputs/te90_2026-09-25/te90_table.csv (session list + fit numbers)
Outputs: outputs/split_trust_2026-09-25/split_trust.csv (one row per session),
         fig/<key>.png (EI + EE coronal: gas, |z|/sigma, phase-in-wedge, out-of-wedge + support contour),
         cohort_trust.png
Usage  : split_trust.py [--only KEY ...] [--no-fig]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import snr_calc as S  # noqa: E402  (corner_box, lung_mask, EXCL_DILATE — same masks as the re-split)

SRC = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
WS = HERE.parents[1]
OUT = WS / 'outputs' / 'split_trust_2026-09-25'
FIG = OUT / 'fig'
TABLE = WS / 'outputs' / 'te90_2026-09-25' / 'te90_table.csv'
SUPPORT_K = 10.0           # calcb: bmask = |b| > 10 * mean|b| in the corner (results.py:130)
RAYLEIGH = np.sqrt(np.pi / 2)   # mean|complex noise| / per-component sigma
RES_FRAC = 4.0             # resolvable voxel: sigma_phi = 1/SNR <= phi_k0/RES_FRAC  ->  SNR >= RES_FRAC/phi_k0
RIM_VOX = 3                # rim band width (voxels) around the support boundary
KERNEL_SIG2 = 0.1          # cudarecon: weight exp(-dsq/0.2) = Gaussian, sigma_k^2 = 0.1 cell^2 (2steve/03)


def apodization(shape):
    """Image-domain envelope of the gridding kernel, never divided out (2steve/03):
    A(r) = exp(-2 pi^2 sigma_k^2 r^2 / N^2). Noise and signal both carry it, so per-voxel
    sigma(x) = sigma_corner * A(x) / <A>_corner. Validated per session against the radial profile."""
    zz, yy, xx = np.indices(shape)
    c = (np.array(shape) - 1) / 2
    r2 = (zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2
    return np.exp(-2 * np.pi ** 2 * KERNEL_SIG2 * r2 / shape[0] ** 2), np.sqrt(r2)


def session_masks(gas):
    """Per-bin lung masks + background exactly as resplit_rbctp / snr_calc."""
    corners = S.corner_box(gas.shape[1:])
    gmean = gas.mean(0)
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gmean, gmean[corners].std()),
                                            iterations=S.EXCL_DILATE)
    masks = [S.lung_mask(gas[b], gas[b][bg].std()) for b in range(gas.shape[0])]
    return gmean, bg, masks


def calcb_support(gmean, bg):
    """Proxy for calcb's bmask: the cycle-average gas image above 10x the mean noise magnitude.
    gas_phase is real(F*b) ~ |F| inside the lung; corner noise per component sigma -> mean|noise| = 1.25 sigma."""
    sigma = gmean[bg].std()
    return gmean > SUPPORT_K * RAYLEIGH * sigma


ROOT = None   # --root: production tree with the container split in recon.mat (C49/C50)


def analyse(key, row, do_fig=True):
    p = (ROOT or SRC) / key
    gas = sio.loadmat(p / 'd' / 'recon.mat', variable_names=['gas_phase'])['gas_phase'].astype(np.float64)
    if ROOT:
        m = sio.loadmat(p / 'd' / 'recon.mat', variable_names=['dissolved_phase_real', 'dissolved_phase_imag',
                                                                 'rbc_tp_separated', 'rbc_tp_dphi_deg', 'rbc_tp_R'])
        if str(m['rbc_tp_separated'][0]).strip() != '1':
            raise RuntimeError('container left this session unsplit (gate)')
        aRBC, aTP = m['dissolved_phase_real'].astype(np.float64), m['dissolved_phase_imag'].astype(np.float64)
        phi = float(np.radians(float(str(m['rbc_tp_dphi_deg'][0]))))
        R_new = np.array([float(v) for v in str(m['rbc_tp_R'][0]).split()]) if 'rbc_tp_R' in m else np.full(gas.shape[0], np.nan)
        gsig = gas.reshape(gas.shape[0], -1).sum(axis=1); b_in = int(np.argmax(gsig))
    else:
        r = sio.loadmat(p / 'd' / 'recon_resplit_xecs.mat',
                        variable_names=['aRBC', 'aTP', 'dphi_new_deg', 'insp_bin', 'R_new', 'ph_new'])
        aRBC, aTP = r['aRBC'].astype(np.float64), r['aTP'].astype(np.float64)
        phi = float(np.radians(r['dphi_new_deg'].ravel()[0]))
        b_in = int(r['insp_bin'].ravel()[0])
        R_new = r['R_new'].ravel()
    gmean, bg, masks = session_masks(gas)
    support = calcb_support(gmean, bg)
    edge = support ^ ndimage.binary_erosion(support, iterations=RIM_VOX)
    rim = ndimage.binary_dilation(edge, iterations=RIM_VOX)      # band straddling the support boundary
    z = aRBC * np.exp(1j * phi) + aTP                            # stored complex dissolved image, TM axis = 0
    A, rad = apodization(gas.shape[1:])
    snr_res = RES_FRAC / phi
    out = {'key': key, 'phi_k0_deg': np.degrees(phi), 'gain': 1 / np.sin(phi), 'insp_bin': b_in,
           'snr_resolve': snr_res}
    # empirical check of the envelope: median |z| in radial rings outside the (dilated) lung, EI bin,
    # fitted as log|z| = a - r^2/L2 over r >= 40 (noise-dominated); theory L2 = N^2/(2 pi^2 sigma_k^2)
    far = ~ndimage.binary_dilation(masks[b_in], iterations=8)
    rr, mm = [], []
    for a in range(40, 80, 4):
        ring = (rad >= a) & (rad < a + 4) & far
        if ring.sum() > 200:
            rr.append((a + 2) ** 2); mm.append(np.log(np.median(np.abs(z[b_in][ring]))))
    slope = np.polyfit(rr, mm, 1)[0] if len(rr) > 3 else np.nan
    out['apod_L2_fit'] = -1 / slope if slope < 0 else np.nan
    out['apod_L2_theory'] = gas.shape[1] ** 2 / (2 * np.pi ** 2 * KERNEL_SIG2)
    L2 = out['apod_L2_fit'] if np.isfinite(out['apod_L2_fit']) else out['apod_L2_theory']
    A_fit = np.exp(-rad ** 2 / L2)                 # primary noise envelope (data-fitted)
    A_th = A                                        # secondary (kernel theory), steeper
    bins = {'EI': b_in, 'EE': 0}
    per = {}
    for tag, b in bins.items():
        mk = masks[b]
        if mk.sum() < 100:
            continue
        zb = z[b]
        sig0 = np.sqrt(0.5 * (zb.real[bg].var() + zb.imag[bg].var()))  # per-component noise, corners
        sig = sig0 * A_fit / A_fit[bg].mean()                            # envelope-scaled noise map (fit)
        sig_th = sig0 * A_th / A_th[bg].mean()                           # theory envelope (steeper -> larger sigma)
        snr = np.abs(zb) / sig
        oow = (aRBC[b] < 0) | (aTP[b] < 0)                                # outside the physical wedge
        ang = np.angle(zb)                                               # 0 = TM axis, phi = RBC axis
        # noise-only prediction of the out-of-wedge fraction: a voxel of true phase theta in [0, phi]
        # with phase noise 1/snr leaves the wedge with prob ~ Q((theta)/s) + Q((phi-theta)/s).
        # Use each voxel's own phase clipped into the wedge as theta.
        from scipy.stats import norm
        s = 1 / np.maximum(snr, 1e-3)
        th = np.clip(ang, 0, phi)
        pred = norm.sf(th / s) + norm.sf((phi - th) / s)
        s_th = sig_th / np.maximum(np.abs(zb), 1e-9)
        pred_th = norm.sf(th / s_th) + norm.sf((phi - th) / s_th)
        resm = snr >= snr_res
        regions = {'lung': mk, 'in': mk & support, 'out': mk & ~support, 'rim': mk & rim}
        d = {f'{tag}_lung_vox': int(mk.sum()),
             f'{tag}_snr_med': float(np.median(snr[mk])),
             f'{tag}_sigma0': float(sig0),
             f'{tag}_resolvable_frac': float((mk & resm).sum() / mk.sum()),
             f'{tag}_lung_out_support_frac': float(regions['out'].sum() / mk.sum())}
        for rn, rm in regions.items():
            if rm.sum() < 20:
                continue
            d[f'{tag}_oow_{rn}_pct'] = 100 * oow[rm].mean()
            d[f'{tag}_oow_{rn}_pred_pct'] = 100 * pred[rm].mean()
            d[f'{tag}_oow_{rn}_excess_pct'] = d[f'{tag}_oow_{rn}_pct'] - d[f'{tag}_oow_{rn}_pred_pct']
            d[f'{tag}_oow_{rn}_excess_th_pct'] = d[f'{tag}_oow_{rn}_pct'] - 100 * pred_th[rm].mean()
            d[f'{tag}_snr_{rn}_med'] = float(np.median(snr[rm]))
        # excess over noise among RESOLVABLE voxels = structured error, not SNR
        d[f'{tag}_oow_excess_resolvable_pct'] = 100 * float(
            (oow & resm)[mk].mean() - (pred * resm)[mk].mean())
        d[f'{tag}_phase_med_deg'] = float(np.degrees(np.median(ang[mk & resm]))) if (mk & resm).any() else np.nan
        # spatial imprint test: RBC fraction vs the static average-gas image, resolvable lung voxels
        sel = mk & resm
        if sel.sum() > 50:
            frac = np.clip(aRBC[b] / np.maximum(np.abs(zb), 1e-9), -1, 2)
            d[f'{tag}_corr_rbcfrac_gmean'] = float(np.corrcoef(frac[sel], gmean[sel])[0, 1])
            d[f'{tag}_corr_rbcfrac_gasbin'] = float(np.corrcoef(frac[sel], gas[b][sel])[0, 1])
        out.update(d)
        per[tag] = dict(b=b, mk=mk, snr=snr, ang=ang, oow=oow, zb=zb, sig=sig, resm=resm, pred=pred)
    out['R_bins_cv'] = float(np.nanstd(R_new) / np.nanmean(R_new)) if np.isfinite(R_new).sum() > 3 else np.nan
    # STATIC-STRUCTURE TEST: does the RBC-fraction pattern move with the lung (EI vs EE)?
    # global z-shift between the two gas bins by cross-correlation (crude, whole-lung), then compare the
    # correlation of the RBC-fraction maps at FIXED voxels vs after that shift. A static imprint scores
    # higher at fixed voxels; a pattern carried by the anatomy scores higher after the shift.
    if 'EI' in per and 'EE' in per:
        gi, ge = gas[per['EI']['b']], gas[per['EE']['b']]
        best, dz_best = -np.inf, 0
        for dz in range(-12, 13):
            ge_s = np.roll(ge, dz, axis=0)
            c = np.corrcoef(gi[masks[b_in]], ge_s[masks[b_in]])[0, 1]
            if c > best:
                best, dz_best = c, dz
        out['shift_dz_vox'] = dz_best
        out['gas_corr_fixed'] = float(np.corrcoef(gi[masks[b_in]], ge[masks[b_in]])[0, 1])
        out['gas_corr_shifted'] = float(best)
        fr = lambda t: np.clip(aRBC[per[t]['b']] / np.maximum(np.abs(per[t]['zb']), 1e-9), -1, 2)
        fi, fe = fr('EI'), fr('EE')
        sel = per['EI']['mk'] & per['EI']['resm']
        sel_s = sel & np.roll(per['EE']['mk'] & per['EE']['resm'], dz_best, axis=0)
        sel_f = sel & per['EE']['mk'] & per['EE']['resm']
        if sel_f.sum() > 200 and sel_s.sum() > 200:
            out['rbcfrac_corr_fixed'] = float(np.corrcoef(fi[sel_f], fe[sel_f])[0, 1])
            out['rbcfrac_corr_shifted'] = float(np.corrcoef(fi[sel_s], np.roll(fe, dz_best, axis=0)[sel_s])[0, 1])
            # same for the total dissolved magnitude (anatomy-carried by construction) as the yardstick
            mi, me = np.abs(per['EI']['zb']), np.abs(per['EE']['zb'])
            out['dismag_corr_fixed'] = float(np.corrcoef(mi[sel_f], me[sel_f])[0, 1])
            out['dismag_corr_shifted'] = float(np.corrcoef(mi[sel_s], np.roll(me, dz_best, axis=0)[sel_s])[0, 1])
            # imprint index: how much more the RBC fraction prefers FIXED coordinates than the magnitude does
            out['imprint_index'] = (out['rbcfrac_corr_fixed'] - out['rbcfrac_corr_shifted']) - \
                                   (out['dismag_corr_fixed'] - out['dismag_corr_shifted'])
            if do_fig:
                static_fig(key, out, per, gas, fi, fe, mi, me, dz_best)
    out['n_bins_split'] = int(np.isfinite(R_new).sum())
    if do_fig and per:
        make_fig(key, out, per, gas, gmean, support, phi)
    return out


def make_fig(key, out, per, gas, gmean, support, phi):
    ys = np.where(per['EI']['mk'].any(axis=(0, 2)))[0]
    y = int(ys[len(ys) // 2]) if len(ys) else gas.shape[2] // 2
    ncol, nrow = 5, len(per)
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 3.3 * nrow), facecolor='black')
    axs = np.atleast_2d(axs)
    sl = lambda v: v[:, y, :][::-1]       # coronal, apex up (post_process: y reversed -> here flip z for display)
    sup = sl(support)
    for i, (tag, d) in enumerate(per.items()):
        b, mk, snr, ang, oow, resm = d['b'], d['mk'], d['snr'], d['ang'], d['oow'], d['resm']
        g = sl(gas[b]); gv = np.percentile(g[sl(mk)], 99.5) if sl(mk).any() else g.max()
        panels = [
            (g, 'gray', 0, gv, f'gas bin {b} ({tag})'),
            (sl(snr), 'magma', 0, 12, '|z|/σ  (dissolved SNR)'),
            (np.where(sl(mk), np.degrees(sl(ang)), np.nan), 'twilight', -90, 90,
             f'phase of z in lung [deg]\nwedge 0 (TM) … {np.degrees(phi):.0f} (RBC)'),
            (np.where(sl(mk), 0.35 + 0.65 * sl(oow).astype(float), np.nan), 'Reds', 0, 1,
             f'out-of-wedge (RBC<0 or TM<0)\nobs {out[f"{tag}_oow_lung_pct"]:.1f}% · noise-pred {out[f"{tag}_oow_lung_pred_pct"]:.1f}%'),
            (np.where(sl(mk), 0.35 + 0.65 * sl(resm).astype(float), np.nan), 'Greens', 0, 1,
             f'resolvable (SNR≥{out["snr_resolve"]:.1f}: σφ≤φ_k0/{RES_FRAC:.0f})\n{100*out[f"{tag}_resolvable_frac"]:.0f}% of lung'),
        ]
        for j, (img, cm, lo, hi, title) in enumerate(panels):
            ax = axs[i, j]
            ax.set_facecolor('black')
            ax.imshow(img, cmap=cm, vmin=lo, vmax=hi, interpolation='nearest')
            ax.contour(sup, levels=[0.5], colors='cyan', linewidths=0.6)
            ax.set_title(title, fontsize=7, color='white')
            ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f'{key}  coronal y={y} · cyan = calcb static phase support (|F_avg| > 10·noise) · '
                 f'φ_k0={np.degrees(phi):.0f}°, gain {out["gain"]:.2f} · '
                 f'EI out-of-wedge obs/pred: in-support {out.get("EI_oow_in_pct", np.nan):.1f}/{out.get("EI_oow_in_pred_pct", np.nan):.1f}% · '
                 f'outside {out.get("EI_oow_out_pct", np.nan):.1f}/{out.get("EI_oow_out_pred_pct", np.nan):.1f}% · '
                 f'rim {out.get("EI_oow_rim_pct", np.nan):.1f}/{out.get("EI_oow_rim_pred_pct", np.nan):.1f}% · '
                 f'envelope L² fit/theory {out["apod_L2_fit"]:.0f}/{out["apod_L2_theory"]:.0f}',
                 fontsize=8, color='white')
    fig.tight_layout()
    fig.savefig(FIG / f'{key}.png', dpi=130, facecolor='black')
    plt.close(fig)


def static_fig(key, out, per, gas, fi, fe, mi, me, dz):
    """EI vs EE, coronal + sagittal: gas, |dissolved|, RBC fraction. The eye test for a static pattern."""
    mi_, me_ = per['EI']['mk'], per['EE']['mk']
    ys = np.where(mi_.any(axis=(0, 2)))[0]; y = int(ys[len(ys) // 2])
    xs = np.where(mi_.any(axis=(0, 1)))[0]; x = int(xs[len(xs) // 4])       # one lung, sagittal
    fig, axs = plt.subplots(4, 3, figsize=(10.5, 13), facecolor='black')
    for r, (tag, mk, f, m, b) in enumerate([('EI', mi_, fi, mi, per['EI']['b']), ('EE', me_, fe, me, per['EE']['b'])]):
        for c, (cut, name) in enumerate([(lambda v: v[:, y, :][::-1], f'coronal y={y}'),
                                         (lambda v: v[:, :, x][::-1], f'sagittal x={x}')]):
            row = 2 * c + r
            g = cut(gas[b]); gv = np.percentile(g[cut(mk)], 99.5)
            panels = [(g, 'gray', 0, gv, f'gas bin {b} ({tag}) {name}'),
                      (np.where(cut(mk), cut(m), np.nan), 'gray', 0, np.percentile(m[mk], 99.5), f'|dissolved| ({tag})'),
                      (np.where(cut(mk), cut(f), np.nan), 'RdYlBu_r', 0, 1, f'RBC fraction aRBC/|z| ({tag})')]
            for j, (img, cm, lo, hi, title) in enumerate(panels):
                ax = axs[row, j]; ax.set_facecolor('black')
                ax.imshow(img, cmap=cm, vmin=lo, vmax=hi, interpolation='nearest')
                ax.set_title(title, fontsize=8, color='white'); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f'{key} · lung shift EE→EI = {dz} vox (gas x-corr) · corr EI vs EE at fixed / shifted voxels: '
                 f'|dissolved| {out["dismag_corr_fixed"]:.2f} / {out["dismag_corr_shifted"]:.2f} · '
                 f'RBC fraction {out["rbcfrac_corr_fixed"]:.2f} / {out["rbcfrac_corr_shifted"]:.2f} · '
                 f'imprint index {out["imprint_index"]:+.2f}', fontsize=8, color='white')
    fig.tight_layout()
    fig.savefig(FIG / f'{key}_static.png', dpi=120, facecolor='black')
    plt.close(fig)


def grade(t):
    """A/B/C/D from four markers; thresholds stated here, nowhere else."""
    g = []
    g.append('A' if t.gain <= 1.35 else 'B' if t.gain <= 1.7 else 'C')                    # angle conditioning
    sd = t.dphi_rep_sd if np.isfinite(t.dphi_rep_sd) else 15
    g.append('A' if sd <= 4 else 'B' if sd <= 8 else 'C')                                  # fit stability (rep sd, deg)
    rf = t.EI_resolvable_frac
    g.append('A' if rf >= 0.6 else 'B' if rf >= 0.35 else 'C')                            # dissolved SNR
    ex = t.EI_oow_excess_resolvable_pct
    g.append('A' if ex <= 2 else 'B' if ex <= 5 else 'C')                                 # structured error beyond noise
    worst = max(g)
    return worst if worst != 'C' or g.count('C') == 1 else 'D', ''.join(g)


def cohort_fig(t):
    t = t.sort_values('key').reset_index(drop=True)
    x = np.arange(len(t))
    fig, axs = plt.subplots(5, 1, figsize=(16, 13.5), sharex=True)
    axs[0].bar(x - 0.2, t.EI_oow_in_excess_pct, 0.4, label='inside calcb support', color='tab:blue')
    axs[0].bar(x + 0.2, t.EI_oow_out_excess_pct, 0.4, label='outside support (polynomial region)', color='tab:red')
    axs[0].plot(x, t.EI_oow_rim_excess_pct, 'k_', ms=8, label='rim band (±3 vox of the boundary)')
    axs[0].axhline(0, color='grey', lw=0.8)
    axs[0].set_ylabel('out-of-wedge EXCESS over noise\n[% of region voxels], EI bin'); axs[0].legend(fontsize=8)
    axs[1].bar(x, t.EI_oow_excess_resolvable_pct, color='tab:orange')
    axs[1].axhline(2, color='grey', ls=':'); axs[1].axhline(5, color='grey', ls=':')
    axs[1].set_ylabel('excess error beyond noise,\nresolvable voxels [% of lung]')
    axs[2].bar(x, 100 * t.EI_resolvable_frac, color='tab:green')
    axs[2].axhline(60, color='grey', ls=':'); axs[2].axhline(35, color='grey', ls=':')
    axs[2].set_ylabel(f'resolvable lung fraction [%]\n(σφ ≤ φ_k0/{RES_FRAC:.0f})')
    cols = {'A': 'tab:green', 'B': 'tab:olive', 'C': 'tab:orange', 'D': 'tab:red'}
    axs[3].bar(x, t.gain, color=[cols[g] for g in t.grade])
    axs[3].set_ylabel('split noise gain 1/sin φ_k0\n(colour = grade)')
    for gname, c in cols.items():
        axs[3].bar([], [], color=c, label=f'{gname} (n={int((t.grade == gname).sum())})')
    axs[3].legend(fontsize=8, ncol=4)
    axs[4].bar(x - 0.2, t.rbcfrac_corr_fixed - t.rbcfrac_corr_shifted, 0.4, color='tab:red',
               label='RBC fraction: corr(EI,EE) fixed − shifted')
    axs[4].bar(x + 0.2, t.dismag_corr_fixed - t.dismag_corr_shifted, 0.4, color='tab:blue',
               label='|dissolved|: fixed − shifted (anatomy yardstick)')
    axs[4].axhline(0, color='grey', lw=0.8)
    axs[4].set_ylabel('static-pattern test\n(>0 = stays in scanner coords)'); axs[4].legend(fontsize=8)
    axs[4].set_xticks(x); axs[4].set_xticklabels(t.key, rotation=90, fontsize=5.5)
    fig.suptitle('RBC/TM split trust, 61 ON_RBC sessions — split_trust.py 2026-09-25', fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / 'cohort_trust.png', dpi=140)
    plt.close(fig)


def main():
    argv = sys.argv[1:]
    do_fig = '--no-fig' not in argv
    global ROOT, OUT, FIG
    if '--root' in argv:
        ROOT = Path(argv[argv.index('--root') + 1])
        tag = argv[argv.index('--tag') + 1] if '--tag' in argv else 'b44'
        OUT = OUT.parent / f'split_trust_{tag}'; FIG = OUT / 'fig'
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
    table = pd.read_csv(TABLE).set_index('key')
    keys = sorted(table[table.atlas_set.eq('yes')].index)
    if '--only' in argv:
        keys = [k for k in argv[argv.index('--only') + 1:] if not k.startswith('--')]
    rows = []
    for k in keys:
        try:
            o = analyse(k, table.loc[k], do_fig)
            o['dphi_rep_sd'] = table.loc[k, 'dphi_rep_sd']
            o['snr_diss_cal'] = table.loc[k, 'snr_diss']
            o['xecs_stable'] = table.loc[k, 'xecs_stable']
            rows.append(o)
            print(f'[trust] {k}: EI oow obs/pred lung {o.get("EI_oow_lung_pct", np.nan):.1f}/{o.get("EI_oow_lung_pred_pct", np.nan):.1f} '
                  f'in {o.get("EI_oow_in_pct", np.nan):.1f}/{o.get("EI_oow_in_pred_pct", np.nan):.1f} '
                  f'out {o.get("EI_oow_out_pct", np.nan):.1f}/{o.get("EI_oow_out_pred_pct", np.nan):.1f} '
                  f'snr_med {o.get("EI_snr_med", np.nan):.1f} resolvable {100*o.get("EI_resolvable_frac", np.nan):.0f}% '
                  f'L2 {o["apod_L2_fit"]:.0f}/{o["apod_L2_theory"]:.0f}', flush=True)
        except Exception as e:
            print(f'[trust] {k}: FAILED {e!r}', flush=True)
    t = pd.DataFrame(rows)
    gr = t.apply(grade, axis=1, result_type='expand')
    t['grade'], t['grade_parts'] = gr[0], gr[1]
    t.to_csv(OUT / ('split_trust.csv' if '--only' not in argv else 'split_trust_sample.csv'), index=False)
    if '--only' not in argv and len(t) > 5:
        cohort_fig(t)
    print(t[['key', 'phi_k0_deg', 'gain', 'EI_snr_med', 'EI_resolvable_frac', 'EI_oow_lung_pct', 'EI_oow_lung_pred_pct',
             'EI_oow_in_excess_pct', 'EI_oow_out_excess_pct', 'EI_oow_rim_excess_pct', 'EI_oow_excess_resolvable_pct',
             'EI_corr_rbcfrac_gmean', 'shift_dz_vox', 'rbcfrac_corr_fixed', 'rbcfrac_corr_shifted',
             'dismag_corr_fixed', 'dismag_corr_shifted', 'grade_parts', 'grade']].round(2).to_string(index=False))


if __name__ == '__main__':
    sys.exit(main())

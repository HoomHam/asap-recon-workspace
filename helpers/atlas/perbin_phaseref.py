#!/usr/bin/env python
"""Per-bin phase reference for the RBC/TM split — does it remove the static RBC imprint (F74)?

Idea (consortium-style, Kaushik 2016): the dissolved image of bin b sees the same per-bin phase
deviation (B0 drift, motion through the field, reference mismatch) as the gas image of bin b,
because both are acquired at the same TE on the same trajectory. calcb's b removes only the
cycle-average phase. The residual phase of the bin's complex gas image G_b = F_gas,b * b is
therefore a per-bin, per-voxel correction for the dissolved image of that bin.

Inputs (fork f5ac7c5 output.mrd): gas_phase_image real(F*b), gas_phase_complex F*b (16,Z,Y,X),
dissolved item = aRBC + i*aTP with meta rbc_tp_dphi_deg / rbc_tp_ph_rad / rbc_tp_target.
Steps per bin:
  z_rot = aRBC e^{i phi} + aTP                (stored basis, after the container's global rotation ph_b)
  z     = z_rot e^{-i ph_b}                   (the complex dissolved image as gridded)
  psi_b = angle( smooth(G_b, sigma) )         (per-bin residual gas phase, lung-masked, sigma vox)
  z'    = z e^{-i psi_b}                      (per-bin referenced)
  re-split: global sweep to the spectroscopic target (same rule as the container / resplit_rbctp)
Then the F74 static test on aRBC (fixed vs shifted EI/EE correlation) before vs after, plus the
out-of-wedge fractions and the imprint index. Also reports the "self" variant (psi from the
dissolved image itself, sigma large) as a bound: it removes ALL low-frequency phase, true or not.

Usage: perbin_phaseref.py <run_root> <key> [--sigma 3] [--fig]
  run_root: folder holding <key>/d/output.mrd (Ext tyger_gascplx_2026-09-25)
Outputs: outputs/perbin_phaseref_2026-09-25/<key>.json + fig/<key>_*.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mrd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import snr_calc as S  # noqa: E402
sys.path.insert(0, str(HERE.parent))
from resplit_rbctp import sweep_phase, split_maps  # noqa: E402

WS = HERE.parents[1]
OUT = WS / 'outputs' / 'perbin_phaseref_2026-09-25'
FIG = OUT / 'fig'


def read_mrd(path):
    gas = gcx = dis = None
    meta = {}
    with mrd.BinaryMrdReader(str(path)) as r:
        r.read_header()
        for it in r.read_data():
            if isinstance(it, mrd.StreamItem.NdArrayFloat) and it.value.meta.get('gas_phase_image'):
                gas = np.asarray(it.value.data, dtype=np.float64)
            elif isinstance(it, mrd.StreamItem.NdArrayComplexFloat):
                if it.value.meta.get('gas_phase_complex'):
                    gcx = np.asarray(it.value.data, dtype=np.complex128)
                elif it.value.meta.get('dissolved_phase_image'):
                    dis = np.asarray(it.value.data, dtype=np.complex128)
                    meta = {k: str(v[0].value) for k, v in it.value.meta.items()}
    return gas, gcx, dis, meta


def smooth_phase(c, mask, sigma):
    """Phase of the Gaussian-smoothed complex image, mask-weighted (no leakage from outside)."""
    w = mask.astype(float)
    num = ndimage.gaussian_filter(c.real * w, sigma) + 1j * ndimage.gaussian_filter(c.imag * w, sigma)
    den = ndimage.gaussian_filter(w, sigma) + 1e-9
    return np.angle(num / den)


def zshift(a, b, mask, rng=12):
    best, dz = -np.inf, 0
    for d in range(-rng, rng + 1):
        c = np.corrcoef(a[mask], np.roll(b, d, axis=0)[mask])[0, 1]
        if c > best:
            best, dz = c, d
    return dz, best


def static_test(A_ei, A_ee, m_ei, m_ee, dz):
    f = m_ei & m_ee
    s = m_ei & np.roll(m_ee, dz, axis=0)
    return (float(np.corrcoef(A_ei[f], A_ee[f])[0, 1]),
            float(np.corrcoef(A_ei[s], np.roll(A_ee, dz, axis=0)[s])[0, 1]))


def run(root, key, sigma, do_fig):
    gas, gcx, dis, meta = read_mrd(Path(root) / key / 'd' / 'output.mrd')
    assert gcx is not None, 'no gas_phase_complex item — wrong image?'
    phi = np.radians(float(meta['rbc_tp_dphi_deg']))
    ph = np.array([float(v) for v in meta['rbc_tp_ph_rad'].split()])
    target = float(meta['rbc_tp_target'])
    nb = gas.shape[0]
    corners = S.corner_box(gas.shape[1:])
    gm = gas.mean(0)
    bg = corners & ~ndimage.binary_dilation(S.lung_mask(gm, gm[corners].std()), iterations=S.EXCL_DILATE)
    masks = [S.lung_mask(gas[b], gas[b][bg].std()) for b in range(nb)]
    b_ei = int(np.argmax([m.sum() for m in masks])); b_ee = 0
    aRBC0, aTP0 = dis.real, dis.imag
    z_rot = aRBC0 * np.exp(1j * phi) + aTP0
    z = z_rot * np.exp(-1j * ph[:, None, None, None])           # un-rotate: gridded dissolved image
    variants = {'orig': (aRBC0, aTP0)}
    psi_store = {}
    for name, src, sig in (('gas', gcx, sigma), ('gas_wide', gcx, 2 * sigma), ('self', z, 3 * sigma)):
        R = np.zeros_like(aRBC0); T = np.zeros_like(aTP0); phs = []
        for b in range(nb):
            mk = masks[b]
            psi = smooth_phase(src[b], ndimage.binary_dilation(mk, iterations=2), sig)
            if name == 'gas' and b in (b_ei, b_ee):
                psi_store[b] = psi
            zc = z[b] * np.exp(-1j * psi)
            if not np.isfinite(ph[b]) or not mk.any():
                phs.append(np.nan); continue
            p, Rr, sR, sT, nv = sweep_phase(zc[mk].sum(), phi, target)
            phs.append(p)
            if np.isnan(p):
                continue
            r_, t_ = split_maps(zc * np.exp(1j * p), phi)
            R[b], T[b] = r_, t_
        variants[name] = (R, T)
    dz, gcorr = zshift(gas[b_ei], gas[b_ee], masks[b_ei])
    Zm = np.abs(z)
    res = {'key': key, 'sigma_vox': sigma, 'phi_k0_deg': float(np.degrees(phi)), 'target': target,
           'b_ei': b_ei, 'b_ee': b_ee, 'dz': dz, 'gas_corr_shifted': gcorr,
           'dismag_fixed_shifted': static_test(Zm[b_ei], Zm[b_ee], masks[b_ei], masks[b_ee], dz)}
    # how big is the per-bin residual gas phase inside the lung?
    for b in (b_ei, b_ee):
        p = psi_store[b][masks[b]]
        res[f'psi_gas_lung_deg_bin{b}'] = {'median_abs': float(np.degrees(np.median(np.abs(p)))),
                                           'p90_abs': float(np.degrees(np.percentile(np.abs(p), 90))),
                                           'mean': float(np.degrees(np.angle(np.exp(1j * p).mean())))}
    for name, (R, T) in variants.items():
        mk = masks[b_ei]
        oow = ((R[b_ei] < 0) | (T[b_ei] < 0))[mk].mean()
        fr = R / np.maximum(Zm, 1e-9)
        res[name] = {'rbc_fixed_shifted': static_test(R[b_ei], R[b_ee], masks[b_ei], masks[b_ee], dz),
                     'tm_fixed_shifted': static_test(T[b_ei], T[b_ee], masks[b_ei], masks[b_ee], dz),
                     'frac_fixed_shifted': static_test(fr[b_ei], fr[b_ee], masks[b_ei], masks[b_ee], dz),
                     'oow_lung_pct_ei': 100 * float(oow),
                     'negRBC_pct_ei': 100 * float((R[b_ei] < 0)[mk].mean()),
                     'negTM_pct_ei': 100 * float((T[b_ei] < 0)[mk].mean()),
                     'corr_rbc_vs_orig_ei': float(np.corrcoef(R[b_ei][mk], aRBC0[b_ei][mk])[0, 1])}
        rf, rs = res[name]['rbc_fixed_shifted']; df_, ds_ = res['dismag_fixed_shifted']
        res[name]['imprint_index'] = (rf - rs) - (df_ - ds_)
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
    json.dump(res, open(OUT / f'{key}_s{sigma}.json', 'w'), indent=1)
    if do_fig:
        make_fig(key, sigma, gas, gcx, masks, b_ei, b_ee, variants, psi_store, res)
    return res


def make_fig(key, sigma, gas, gcx, masks, b_ei, b_ee, variants, psi, res):
    ys = np.where(masks[b_ei].any(axis=(0, 2)))[0]; y = int(ys[len(ys) // 2])
    sl = lambda v: v[:, y, :][::-1]
    names = ['orig', 'gas', 'gas_wide', 'self']
    fig, axs = plt.subplots(2, 2 + len(names), figsize=(3.0 * (2 + len(names)), 6.6), facecolor='black')
    for r, b in enumerate((b_ei, b_ee)):
        mk = sl(masks[b])
        g = sl(gas[b]); ax = axs[r, 0]; ax.imshow(g, cmap='gray', vmin=0, vmax=np.percentile(g[mk], 99.5))
        ax.set_title(f'gas bin {b} ({"EI" if r == 0 else "EE"})', fontsize=7, color='white')
        ax = axs[r, 1]; ax.imshow(np.where(mk, np.degrees(sl(psi[b])), np.nan), cmap='twilight', vmin=-60, vmax=60)
        ax.set_title(f'residual gas phase ψ_b [deg] (σ={sigma})', fontsize=7, color='white')
        for j, n in enumerate(names):
            R = variants[n][0][b]; ax = axs[r, 2 + j]
            vmax = np.percentile(np.maximum(R, 0)[masks[b]], 99.5) + 1e-9
            ax.imshow(np.where(mk, sl(np.maximum(R, 0)), np.nan), cmap='gray', vmin=0, vmax=vmax)
            f, s = res[n]['rbc_fixed_shifted']
            ax.set_title(f'aRBC {n}\nfixed/shifted {f:.2f}/{s:.2f}  imprint {res[n]["imprint_index"]:+.2f}',
                         fontsize=7, color='white')
    for ax in axs.ravel():
        ax.set_facecolor('black'); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f'{key}  per-bin gas phase reference test · coronal y={y} · |dis| fixed/shifted '
                 f'{res["dismag_fixed_shifted"][0]:.2f}/{res["dismag_fixed_shifted"][1]:.2f}', fontsize=8, color='white')
    fig.tight_layout(); fig.savefig(FIG / f'{key}_s{sigma}.png', dpi=130, facecolor='black'); plt.close(fig)


if __name__ == '__main__':
    a = sys.argv[1:]
    sigma = float(a[a.index('--sigma') + 1]) if '--sigma' in a else 3.0
    r = run(a[0], a[1], sigma, '--fig' in a)
    for n in ('orig', 'gas', 'gas_wide', 'self'):
        v = r[n]
        print(f'{n:9s} aRBC fixed/shifted {v["rbc_fixed_shifted"][0]:.3f}/{v["rbc_fixed_shifted"][1]:.3f} '
              f'imprint {v["imprint_index"]:+.3f}  oow {v["oow_lung_pct_ei"]:.1f}%  corr_vs_orig {v["corr_rbc_vs_orig_ei"]:.3f}')
    print('|dis| fixed/shifted', np.round(r['dismag_fixed_shifted'], 3), 'dz', r['dz'])
    print('psi_gas lung:', {k: v for k, v in r.items() if k.startswith('psi_gas')})

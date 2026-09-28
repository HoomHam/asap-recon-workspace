"""In-lung error of P1 vs PF variants against a fully-sampled truth (follow-up to pf_test.py).

Hooman's eye (2026-09-26): the non-KLUGE PFs lung looks best, Tyger's worst, "structure and smoothness". pf_test
measured noise in the FAR background only; inside the lung the error is thermal noise + per-bin aliasing of the
lung's own signal (54 % k-space holes). This isolates the two:

  truth T   = real(F_full . b), all interleaves, no binning (fully sampled k-space on the grid)
  per bin   the real sampled-cell pattern S_bin and weights knorm_bin from the bin's gridding
  sim       K = K_full * S_bin  (+ optional complex noise on the sampled cells, level matched to the real far-bg sigma)
            -> P1 (oldway 1 + KLUGE) / PFs (oldway 0, Hermitian, knorm-weighted, no KLUGE) / PFs+k / PFc (POCS)
  error     inside the lung vs T: rms error / lung mean, correlation, and high-pass ("texture") error separately

usage: python pf_inlung.py [--key 2024-10-22_021CH] [--bins 0,4,8] [--noise-draws 3]
out:   outputs/pf_test/<key>/inlung_{report.md, metrics.csv, sim_bin*.png}
"""
import argparse
import csv
import os
import sys

import numpy as np
from scipy import ndimage as ndi
from scipy.io import loadmat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pf_test import (ASAP, RAW_TREE, PROD_TREE, OUT, grid_acc, renorm_oldway1, renorm_oldway0, to_image,  # noqa: E402
                     crop_slices, pocs_pf, lung_masks, log, read_input_mrd, exclude_mask, bin_sample_weight,
                     gvar, imgtype, graddir, traj, steve_raw_cls, KNORM_FLOOR)
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402


def recon_variants(k_acc, knorm, b, MS, IS, iters=20):
    cs = crop_slices(MS, IS)
    k1 = renorm_oldway1(k_acc, knorm)
    k0, _ = renorm_oldway0(k_acc, knorm, MS)
    k0k, _ = renorm_oldway0(k_acc, knorm, MS, kluge=True)
    out = {'P1': np.real(to_image(k1, MS)[cs, cs, cs] * b),
           'PFs': np.real(to_image(k0, MS)[cs, cs, cs] * b),
           'PFs+k': np.real(to_image(k0k, MS)[cs, cs, cs] * b)}
    xc, _ = pocs_pf(k1, knorm, b, MS, IS, iters=iters)
    out['PFc'] = np.real(xc * b)
    return out


def errors(img, T, lung, hp_sigma=2.0):
    d = img - T
    e = d[lung]
    hp = lambda v: v - ndi.gaussian_filter(v, hp_sigma)   # texture = what a 2-vox smooth removes
    dh = hp(img) - hp(T)
    return dict(rms_err=np.sqrt(np.mean(e ** 2)) / T[lung].mean(),
                bias=e.mean() / T[lung].mean(),
                corr=np.corrcoef(img[lung], T[lung])[0, 1],
                texture_err=np.sqrt(np.mean(dh[lung] ** 2)) / T[lung].mean(),
                lung_mean_ratio=img[lung].mean() / T[lung].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--key', default='2024-10-22_021CH')
    ap.add_argument('--bins', default='0,4,8')
    ap.add_argument('--noise-draws', type=int, default=3)
    ap.add_argument('--iters', type=int, default=20)
    ap.add_argument('--b-raw', action='store_true',
                    help='phase reference = our calcb b BEFORE the 4.4-vox low-pass (full-res measured phase in the mask)')
    ap.add_argument('--b-self', action='store_true',
                    help='phase reference = conj phase of the truth object itself (voxel-exact; F.b real by construction)')
    ap.add_argument('--synth-real', action='store_true',
                    help='replace the truth k-space by the exact k-space of the real, crop-supported image T (ideal PF object)')
    ap.add_argument('--lung-only', action='store_true',
                    help='with --synth-real: zero the object outside a 12-vox dilation of the lung (lung spectrum, no bg noise)')
    ap.add_argument('--truth-bin', type=int, default=None,
                    help='use this bin\'s own P1 k-space (single breathing phase, sharp) as the truth instead of all data')
    a = ap.parse_args()
    tag = ('' if a.truth_bin is None else f'_truthbin{a.truth_bin}') + ('_braw' if a.b_raw else '') + ('_bself' if a.b_self else '') + ('_synth' if a.synth_real else '') + ('_lungonly' if a.lung_only else '')
    out = os.path.join(OUT, a.key)
    bins_sel = [int(x) for x in a.bins.split(',')]
    rng = np.random.default_rng(0)

    header, arrs = read_input_mrd(os.path.join(RAW_TREE, a.key, 'd', 'input.mrd'))
    prod = loadmat(os.path.join(PROD_TREE, a.key, 'd', 'recon.mat'))
    ul = {p.name: p.value for p in header.user_parameters.user_parameter_long}
    ud = {p.name: p.value for p in header.user_parameters.user_parameter_double}
    g = gvar()
    g.MS, g.IS, g.nbins, g.gplb = int(ul.get('MS', 240)), int(ul.get('IS', 100)), int(ul.get('nbins', 16)), int(ul.get('gplb', 300))
    killpts = int(ul.get('killpts', 2))
    meta = {k: ud.get(k, 0.0) for k in ('TR', 'TE', 'DPoff', 'dtdyn', 'dtspec')}
    meta['numspec'] = int(ul.get('numspec', 0))
    g_traj = traj()
    g_traj.killpts = killpts
    g_traj.load_traj_from_array(arrs['gp'], arrs['dp'], int(ul.get('nusimg', 32)))
    dyn = arrs['dyn'][:, killpts:, :]
    ref = arrs['ref'][:, killpts:, :] if arrs['ref'] is not None else None
    g_raw = steve_raw_cls()
    g_raw.load_from_arr(g_traj, ref, dyn, arrs['pneumo'], 'mrd_siemens', meta)
    g_traj.rescale_to_MS(g.MS, g.IS)
    MS, IS, npts = g.MS, g.IS, g_raw.npts
    tg = np.stack([np.asarray(g_traj.gettraj(imgtype.GPDYN, d)) for d in (graddir.X, graddir.Y, graddir.Z)], 1)
    data = np.reshape(np.asarray(g_raw.getimg(imgtype.GPDYN))[:, 0, :], npts * g_raw.ntotalilvs, order='F') * exclude_mask(g_raw)
    b = np.asarray(prod['calcb_b'])[0].astype(np.complex128)
    b /= np.maximum(np.abs(b), 1e-12)
    if a.b_raw:
        from pf_test import calcb_numpy
        log('b: full-resolution calcb (no phase low-pass) ...')
        b_raw = calcb_numpy(tg, data, npts, MS, IS, g.gplb)
        dphi = np.angle(b_raw * np.conj(b))
        log(f'b: raw vs production phase diff rms {np.degrees(np.sqrt(np.mean(dphi**2))):.2f} deg over the crop')
        b = b_raw
    ilvvol = np.asarray(prod['nav_volume'], dtype=float).ravel()
    binarr = np.ascontiguousarray(np.asarray(g_raw.bin(ilvvol)) * g.nbins)
    cs = crop_slices(MS, IS)
    prod_gas = np.asarray(prod['gas_phase'])

    # ---- truth: all data, no binning
    if a.truth_bin is None:
        log('truth: gridding all interleaves ...')
        ka_full, kn_full = grid_acc(tg, data, npts, MS, g.gplb, None)
    else:
        log(f'truth: gridding bin {a.truth_bin} alone (single phase, its own holes are part of the object) ...')
        sw_t = bin_sample_weight(binarr, g.nbins, a.truth_bin, npts, len(data))
        ka_full, kn_full = grid_acc(tg, data, npts, MS, g.gplb, sw_t)
    K_full = renorm_oldway1(ka_full, kn_full)
    F_t = to_image(K_full, MS)[cs, cs, cs]
    if a.b_self:
        b = np.conj(F_t) / np.maximum(np.abs(F_t), 1e-12)
    T = np.real(F_t * b)
    m = lung_masks(np.abs(to_image(K_full, MS)[cs, cs, cs]))
    lung = m['lung']
    log(f'truth: lung {lung.sum()} vox, lung mean {T[lung].mean():.4g}, far-bg sd {T[m["corner"]].std():.4g} '
        f'(holes inside sphere {np.mean(kn_full[kn_full > 0] < KNORM_FLOOR):.3f} of touched cells)')
    resid = np.angle(F_t * b)[lung]
    log(f'truth: residual phase of F.b inside the lung: {np.degrees(np.sqrt(np.mean(resid**2))):.2f} deg rms, '
        f'|Im|/|F| {np.mean(np.abs(np.sin(resid))):.3f}')
    x_full = to_image(K_full, MS)
    e_in = np.sum(np.abs(x_full[cs, cs, cs]) ** 2); e_all = np.sum(np.abs(x_full) ** 2)
    log(f'truth: image energy outside the {IS}^3 crop = {1 - e_in / e_all:.4f} of total (240^3 grid)')
    if a.synth_real:
        if a.lung_only:
            T = T * ndi.binary_dilation(lung, iterations=12)
        xpad = np.zeros((MS, MS, MS), dtype=np.complex128)
        xpad[cs, cs, cs] = T
        K_full = np.reshape(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(xpad))), -1, order='F')
        b = np.ones_like(b)
        chk = np.real(to_image(K_full, MS)[cs, cs, cs])
        log(f'truth: synthetic real crop-supported object; round-trip error {np.linalg.norm(chk - T) / np.linalg.norm(T):.1e}')
    np.save(os.path.join(out, f'truth_T{tag}.npy'), T)

    rows, panels = [], {}
    for ibin in bins_sel:
        log(f'bin {ibin}: gridding for the sampling pattern ...')
        sw = bin_sample_weight(binarr, g.nbins, ibin, npts, len(data))
        ka_bin, kn_bin = grid_acc(tg, data, npts, MS, g.gplb, sw)
        real_far_sd = prod_gas[ibin][m['corner']].std()
        # --- noise-free: the bin's hole pattern applied to the truth k-space
        ka_sim = K_full * kn_bin                       # so cudarenorm-style weighting is the bin's own
        nf = recon_variants(ka_sim, kn_bin, b, MS, IS, a.iters)
        if ibin == bins_sel[0]:
            log('  pocs convergence (noise-free):')
            pocs_pf(renorm_oldway1(ka_sim, kn_bin), kn_bin, b, MS, IS, iters=a.iters, verbose=True)
        for name, img in nf.items():
            rows.append(dict(bin=ibin, noise='none', variant=name, **errors(img, T, lung)))
        log('  noise-free rms err / lung mean: ' + ', '.join(f'{k} {errors(v, T, lung)["rms_err"]:.4f}' for k, v in nf.items()))
        # --- calibrate complex noise on sampled cells so P1's far-bg sd matches the real bin
        S = kn_bin >= KNORM_FLOOR
        n0 = (rng.standard_normal(kn_bin.shape) + 1j * rng.standard_normal(kn_bin.shape)) * S
        far_sd_unit = np.real(to_image(n0, MS)[cs, cs, cs] * b)[m['corner']].std()
        sig = real_far_sd / far_sd_unit
        log(f'  noise: real far-bg sd {real_far_sd:.4g} -> per-cell sigma {sig:.4g}')
        acc = {k: [] for k in nf}
        for d in range(a.noise_draws):
            n = (rng.standard_normal(kn_bin.shape) + 1j * rng.standard_normal(kn_bin.shape)) * S * sig
            ka_n = (K_full + n) * kn_bin
            nv = recon_variants(ka_n, kn_bin, b, MS, IS, a.iters)
            for name, img in nv.items():
                acc[name].append(errors(img, T, lung))
                if d == 0:
                    panels[(ibin, name)] = img
        for name in nf:
            mean = {k: float(np.mean([e[k] for e in acc[name]])) for k in acc[name][0]}
            rows.append(dict(bin=ibin, noise='matched', variant=name, **mean))
        log('  with noise rms err / lung mean: ' + ', '.join(f'{k} {np.mean([e["rms_err"] for e in acc[k]]):.4f}' for k in nf))
        # also the REAL bin images (from pf_test arrays) vs T, for reference: motion + noise + aliasing
        try:
            arr = np.load(os.path.join(out, 'arrays.npz'))
            for name, key in (('P1', 'P1'), ('PFs', 'PFs'), ('PFs+k', 'PFsk'), ('PFc', 'PFc')):
                rows.append(dict(bin=ibin, noise='real data', variant=name, **errors(arr[f'gas_{ibin}_{key}'], T, lung)))
        except (FileNotFoundError, KeyError):
            pass

    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(out, f'inlung{tag}_metrics.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.5g}' if isinstance(v, float) else v) for k, v in r.items()})
    with open(os.path.join(out, f'inlung{tag}_report.md'), 'w') as f:
        f.write(f'# In-lung error vs fully-sampled truth — {a.key}\n\n'
                'rms err and texture err are inside the lung, relative to the truth lung mean. "none" = bin hole pattern only '
                '(pure aliasing); "matched" = plus complex noise on the sampled cells scaled so P1 far-bg σ equals the real bin; '
                '"real data" = the actual bin recon (motion + noise + aliasing) vs the all-data truth.\n\n')
        f.write('| bin | noise | variant | rms err | texture err | corr | bias | lung mean ratio |\n|---|---|---|---|---|---|---|---|\n')
        for r in rows:
            f.write(f"| {r['bin']} | {r['noise']} | {r['variant']} | {r['rms_err']:.4f} | {r['texture_err']:.4f} | {r['corr']:.4f} | "
                    f"{r['bias']:+.4f} | {r['lung_mean_ratio']:.4f} |\n")

    z = int(np.round(ndi.center_of_mass(lung)[0]))
    for ibin in bins_sel:
        vmax = np.percentile(T[lung], 99)
        fig, ax = plt.subplots(2, 5, figsize=(20, 8))
        ax[0, 0].imshow(T[z], cmap='gray', vmin=-0.15 * vmax, vmax=vmax)
        ax[0, 0].set_title('truth (all data)' if a.truth_bin is None else f'truth (real bin {a.truth_bin} P1)')
        ax[1, 0].axis('off')
        for j, name in enumerate(('P1', 'PFs', 'PFs+k', 'PFc'), 1):
            img = panels[(ibin, name)]
            ax[0, j].imshow(img[z], cmap='gray', vmin=-0.15 * vmax, vmax=vmax)
            ax[0, j].set_title(f'{name} (sim: bin {ibin} holes + matched noise)')
            ax[1, j].imshow((img - T)[z] * lung[z], cmap='RdBu_r', vmin=-0.3 * vmax, vmax=0.3 * vmax)
            ax[1, j].set_title('error inside lung (±30 %)')
        for a_ in ax.ravel():
            a_.axis('off')
        fig.suptitle(f'{a.key} bin {ibin}: truth resampled through the real bin hole pattern — what each recon does to the lung')
        fig.tight_layout()
        fig.savefig(os.path.join(out, f'inlung{tag}_bin{ibin}.png'), dpi=100)
        plt.close(fig)
    log('done ->', out)


if __name__ == '__main__':
    main()

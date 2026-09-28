"""Does Steve's partial-Fourier branch (recon.py::cudarenorm oldway = 0) have merit?

Offline, Tyger-exact rebuild of one session from input.mrd (XeCS numpy port of cudarecon), DIAPHRAGM bins
taken from the production recon.mat (nav_volume -> raw.bin), calcb b with the 04f445b phase low-pass.
From ONE gridded k-space per bin (k_acc, knorm before normalisation) four reconstructions are derived:

  P1      oldway = 1 : k/knorm, KLUGE                    -> real(F.b)         == Tyger production (checked)
  PFs     oldway = 0 as written (Aug-2024, no rephase) : Hermitian average of k and conj(k at -k), knorm-
          weighted, DC zeroed, no KLUGE                  -> Fh real            -> real(Fh.b)
  PFs+k   PFs with the KLUGE applied to knorm(k)+knorm(-k) (isolates the missing-KLUGE noise)
  PFc     the correct phase-constrained partial Fourier : POCS with the SAME phase reference b ---
          image realness after b inside the crop, zero outside, measured cells kept, holes filled

and the same for the dissolved k-space (dplb = 40, complex image). Reported per bin: lung mean, far-bg
sigma, near-lung shell sigma (aliasing), edge sharpness, hole fill fraction, identity test Fh == Re[F], and
the RBC/TP split (LS-mapped to the production aRBC/aTP) for P1 vs PF.

usage: python pf_test.py [--key 2024-10-22_021CH] [--bins 0,4,8] [--iters 20]
out:   workspace/outputs/pf_test/<key>/{report.md, metrics.csv, *.png, arrays.npz}
"""
import argparse
import csv
import os
import sys
import time

import numpy as np
from numpy.fft import fftn, fftshift, ifftn, ifftshift
from scipy import ndimage as ndi
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter

ASAP = os.path.expanduser('~/Hooman/Work/Codes/2026_ASAP_Recon')
XECS = os.path.expanduser('~/Hooman/Work/Codes/2026_XeCS_Recon')
sys.path.insert(0, ASAP)
sys.path.insert(0, os.path.join(ASAP, 'workspace/helpers/recon'))
sys.path.insert(0, os.path.join(XECS, 'recon'))
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from gtypes import gvar, imgtype, graddir  # noqa: E402
from raw import traj, raw as steve_raw_cls  # noqa: E402
from steve_kernel_numpy import EPS, KDIST0SQ, BXSZ, KNORM_FLOOR  # noqa: E402
from selftest_steve_tyger import read_input_mrd, exclude_mask, calcb_numpy, bin_sample_weight  # noqa: E402

RAW_TREE = '/Volumes/HoomHamExt/AIkill_Dynamic/'       # has input.mrd
PROD_TREE = '/Volumes/HoomHamExt/AIkill_Dynamic_b44/'  # production recon.mat (fork 04f445b)
OUT = os.path.join(ASAP, 'workspace/outputs/pf_test')
SIGMA_B = 4.4                                          # calcb_phase_sigma at 04f445b
C_SERIES = ['#2563EB', '#EA580C', '#6B7280']           # blue, orange, gray (fixed order)


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


# ---------------------------------------------------------------- gridding: k_acc / knorm (cudarecon)
def grid_acc(traj_grid, data, npts, MS, smoothing, sample_weight=None):
    """cudarecon accumulation only (copy of steve_kernel_numpy.steve_recon up to the renorm)."""
    data = np.asarray(data, dtype=np.complex128)
    idx = np.arange(len(data))
    rawval = data.copy()
    if smoothing > EPS:
        rawval = rawval * np.exp(-(((idx % npts) / smoothing) ** 2))
    keep = np.abs(rawval) >= EPS
    if sample_weight is not None:
        keep &= np.asarray(sample_weight) > 0
        sw = np.asarray(sample_weight, dtype=np.float64)[keep]
    rawval = rawval[keep]
    kidx = idx[keep] % traj_grid.shape[0]
    kx, ky, kz = (traj_grid[kidx, i] for i in range(3))
    cx, cy, cz = ((k + 0.5).astype(np.int64) for k in (kx, ky, kz))
    MS2 = MS * MS
    k_acc = np.zeros(MS ** 3, dtype=np.complex128)
    knorm = np.zeros(MS ** 3, dtype=np.float64)
    offs = range(-BXSZ, BXSZ)
    for dx in offs:
        ix = cx + dx
        vx = (ix >= np.maximum(0, cx - BXSZ)) & (ix < np.minimum(MS, cx + BXSZ)) & (ix >= 0) & (ix < MS)
        for dy in offs:
            iy = cy + dy
            vy = vx & (iy >= 0) & (iy < MS)
            for dz in offs:
                iz = cz + dz
                v = vy & (iz >= 0) & (iz < MS)
                if not v.any():
                    continue
                dsq = (kx[v] - ix[v]) ** 2 + (ky[v] - iy[v]) ** 2 + (kz[v] - iz[v]) ** 2
                wt = np.exp(-dsq / KDIST0SQ)
                if sample_weight is not None:
                    wt = wt * sw[v]
                flat = ix[v] * MS2 + iy[v] * MS + iz[v]
                k_acc += np.bincount(flat, weights=wt * rawval[v].real, minlength=MS ** 3) \
                    + 1j * np.bincount(flat, weights=wt * rawval[v].imag, minlength=MS ** 3)
                knorm += np.bincount(flat, weights=wt, minlength=MS ** 3)
    return k_acc, knorm


def renorm_oldway1(k_acc, knorm):
    k = np.zeros_like(k_acc)
    nz = knorm > 0
    k[nz] = k_acc[nz] / knorm[nz]
    k[knorm < KNORM_FLOOR] = 0
    return k


def renorm_oldway0(k_acc, knorm, MS, kluge=False):
    """cudarenorm oldway = 0, Aug-2024 form (no rephase), flat-index mirror exactly as his."""
    N = MS ** 3
    c = int(MS / 2 + 0.1) * (MS ** 2 + MS + 1)
    k = np.zeros_like(k_acc)
    idx = np.arange(0, c + 1)
    sym = 2 * c - idx
    ok = (sym < N) & (idx != c)                    # his: symmidx > len(k) or idx == cidx -> zeroed
    i, s = idx[ok], sym[ok]
    kn = knorm + 1e-16                             # cudarezero sets knorm = eps
    den = kn[i] + kn[s]
    th = (k_acc[i] + np.conj(k_acc[s])) / den
    if kluge:
        th[den < KNORM_FLOOR] = 0
    k[i] = th
    k[s] = np.conj(th)
    den_full = np.zeros(N)                         # summed weight per cell, both halves (zeroed cells stay 0)
    den_full[i] = den
    den_full[s] = den
    return k, den_full


def to_image(kflat, MS):
    k3 = np.reshape(kflat, (MS, MS, MS), order='F')
    return fftshift(fftn(ifftshift(k3)))


def crop_slices(MS, IS):
    MSc = int(MS / 2 + 0.1)
    ll = int(MSc - IS / 2 + 0.1)
    return slice(ll, ll + IS)


def pocs_pf(kflat, knorm, b, MS, IS, iters=20, verbose=False):
    """Phase-constrained partial Fourier with Steve's own b as the phase map: image = real after b inside the
    IS crop, zero outside (the FOV support Steve discards anyway), measured cells re-imposed each iteration."""
    K = np.reshape(kflat, (MS, MS, MS), order='F')
    S = np.reshape(knorm >= KNORM_FLOOR, (MS, MS, MS), order='F')
    Kmeas = K[S]
    cs = crop_slices(MS, IS)
    x = fftshift(fftn(ifftshift(K)))
    for it in range(iters):
        xc = x[cs, cs, cs] * b
        xc = np.real(xc) * np.conj(b)
        xn = np.zeros_like(x)
        xn[cs, cs, cs] = xc
        Kn = fftshift(ifftn(ifftshift(xn)))
        Kn[S] = Kmeas
        x_new = fftshift(fftn(ifftshift(Kn)))
        d = np.linalg.norm(x_new[cs, cs, cs] - x[cs, cs, cs]) / np.linalg.norm(x[cs, cs, cs])
        x = x_new
        if verbose:
            log(f'    pocs it {it}: rel change {d:.2e}')
    return x[cs, cs, cs], Kn


# ---------------------------------------------------------------- masks + metrics
def lung_masks(mean_img):
    lung = mean_img > 0.2 * np.percentile(mean_img, 99.9)
    lung = ndi.binary_opening(lung, iterations=1)
    lab, nl = ndi.label(lung)
    lung = lab == (1 + np.argmax(ndi.sum(lung, lab, range(1, nl + 1))))
    d8, d16 = ndi.binary_dilation(lung, iterations=8), ndi.binary_dilation(lung, iterations=16)
    shell = d16 & ~d8
    far = ~d16
    n = mean_img.shape[0]
    edge = np.zeros_like(far)
    s = max(4, n // 8)
    edge[:s] = edge[-s:] = True
    edge[:, :s] = edge[:, -s:] = True
    edge[:, :, :s] = edge[:, :, -s:] = True
    ring = ndi.binary_dilation(lung, iterations=1) & ~ndi.binary_erosion(lung, iterations=1)
    return dict(lung=lung, shell=shell, corner=far & edge, ring=ring)


def metrics(img, m, ref=None):
    g = np.sqrt(sum(np.gradient(img, axis=a) ** 2 for a in range(3)))
    out = dict(lung_mean=img[m['lung']].mean(), bg_sd=img[m['corner']].std(), bg_mean=img[m['corner']].mean(),
               shell_sd=img[m['shell']].std(), edge_grad=g[m['ring']].mean() / max(abs(img[m['lung']].mean()), 1e-30))
    out['snr'] = out['lung_mean'] / out['bg_sd']
    out['shell_over_bg'] = out['shell_sd'] / out['bg_sd']
    if ref is not None:
        out['corr_lung_vs_P1'] = np.corrcoef(img[m['lung']], ref[m['lung']])[0, 1]
        out['bg_noise_corr_vs_P1'] = np.corrcoef(img[m['corner']], ref[m['corner']])[0, 1]
    return out


def coverage(knorm, den, traj_grid, MS):
    """Fraction of grid cells inside the sampled k-sphere that have no data (knorm < floor), before and after
    the Hermitian fill, as a function of k radius."""
    ctr = MS / 2
    rmax = np.max(np.linalg.norm(traj_grid - ctr, axis=1))
    ax = np.arange(MS) - ctr
    r = np.sqrt(ax[:, None, None] ** 2 + ax[None, :, None] ** 2 + ax[None, None, :] ** 2)
    r = np.reshape(r, -1, order='F')                    # symmetric in x,y,z: layout irrelevant
    inside = r <= rmax
    hole1 = inside & (knorm < KNORM_FLOOR)
    hole0 = inside & (den < KNORM_FLOOR)
    edges = np.linspace(0, rmax, 13)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = inside & (r >= lo) & (r < hi)
        rows.append((0.5 * (lo + hi), hole1[sel].mean(), hole0[sel].mean()))
    return dict(frac_hole_P1=hole1.sum() / inside.sum(), frac_hole_PF=hole0.sum() / inside.sum(),
                filled=(hole1 & ~hole0).sum() / max(hole1.sum(), 1), radial=np.array(rows), rmax=rmax)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--key', default='2024-10-22_021CH')
    ap.add_argument('--bins', default='0,4,8')
    ap.add_argument('--iters', type=int, default=20)
    a = ap.parse_args()
    out = os.path.join(OUT, a.key)
    os.makedirs(out, exist_ok=True)
    bins_sel = [int(x) for x in a.bins.split(',')]

    # ---- load raw exactly like tyger_recon / selftest
    header, arrs = read_input_mrd(os.path.join(RAW_TREE, a.key, 'd', 'input.mrd'))
    prod = loadmat(os.path.join(PROD_TREE, a.key, 'd', 'recon.mat'))
    ul = {p.name: p.value for p in header.user_parameters.user_parameter_long}
    ud = {p.name: p.value for p in header.user_parameters.user_parameter_double}
    g = gvar()
    g.MS, g.IS, g.nbins, g.gplb = int(ul.get('MS', 240)), int(ul.get('IS', 100)), int(ul.get('nbins', 16)), int(ul.get('gplb', 300))
    g.dplb = int(ul.get('dplb', g.dplb))
    killpts = int(ul.get('killpts', 2))
    meta = {k: ud.get(k, 0.0) for k in ('TR', 'TE', 'DPoff', 'dtdyn', 'dtspec')}
    meta['numspec'] = int(ul.get('numspec', 0))
    g_traj = traj()
    g_traj.killpts = killpts
    g_traj.load_traj_from_array(arrs['gp'], arrs['dp'], int(ul.get('nusimg', 32)))
    dyn = arrs['dyn'][:, killpts:, :]
    ref = arrs['ref'][:, killpts:, :] if arrs['ref'] is not None else None
    assert dyn.shape[0] == 1, 'single channel only'
    g_raw = steve_raw_cls()
    g_raw.load_from_arr(g_traj, ref, dyn, arrs['pneumo'], 'mrd_siemens', meta)
    g_traj.rescale_to_MS(g.MS, g.IS)
    MS, IS, npts = g.MS, g.IS, g_raw.npts
    mask = exclude_mask(g_raw)
    tg = {it: np.stack([np.asarray(g_traj.gettraj(it, d)) for d in (graddir.X, graddir.Y, graddir.Z)], 1)
          for it in (imgtype.GPDYN, imgtype.DPDYN)}
    data = {it: np.reshape(np.asarray(g_raw.getimg(it))[:, 0, :], npts * g_raw.ntotalilvs, order='F') * mask
            for it in (imgtype.GPDYN, imgtype.DPDYN)}
    log(f'{a.key}: MS={MS} IS={IS} nbins={g.nbins} gplb={g.gplb} dplb={g.dplb} npts={npts} nilv={g_raw.ntotalilvs} '
        f'gas samples={len(data[imgtype.GPDYN])} dp samples={len(data[imgtype.DPDYN])}')

    # ---- DIAPHRAGM bins from the production navigator (tyger_recon: ilvbin = raw.bin(ilvvol))
    ilvvol = np.asarray(prod['nav_volume'], dtype=float).ravel()
    assert len(ilvvol) == g_raw.ntotalilvs, (len(ilvvol), g_raw.ntotalilvs)
    binarr = np.ascontiguousarray(np.asarray(g_raw.bin(ilvvol)) * g.nbins)
    occ = [int(np.sum(np.floor(binarr[binarr >= 0]).astype(int) == i)) for i in range(g.nbins)]
    log('bin occupancy', occ)

    # ---- b: calcb + 04f445b phase low-pass; check against production calcb_b
    log('calcb ...')
    b = calcb_numpy(tg[imgtype.GPDYN], data[imgtype.GPDYN], npts, MS, IS, g.gplb)
    bmask = (np.abs(b) > 0).astype(float)  # placeholder, recomputed below the way calcb does
    # calcb_numpy returns unit |b| already; rebuild bmask from the un-normalised image is not possible here,
    # so take bmask from the production b: inside the mask the production phase differs from the polynomial,
    # which is exactly where calcb kept the measured phase. Use the low-passed production b directly as the
    # phase reference instead (it IS what Tyger multiplied) and keep our calcb b for the identity check.
    b_prod = np.asarray(prod['calcb_b'])[0].astype(np.complex128)   # (IS,IS,IS) as multiplied by Tyger
    b_prod = b_prod / np.maximum(np.abs(b_prod), 1e-12)
    # layout: recon.mat arrays are Tyger's (z,y,x)=steve layout; our images are steve layout too
    dphi = np.angle(b_prod * np.conj(b))
    log(f'our calcb (pre low-pass) vs production b: phase diff rms {np.degrees(np.sqrt(np.mean(dphi**2))):.2f} deg '
        f'(expect ~ the 5.9 deg the low-pass removed inside the mask, ~0 outside)')
    b = b_prod
    np.save(os.path.join(out, 'b_prod.npy'), b)

    prod_gas = np.asarray(prod['gas_phase'])
    prod_re, prod_im = np.asarray(prod['dissolved_phase_real']), np.asarray(prod['dissolved_phase_imag'])
    cs = crop_slices(MS, IS)

    rows, figs = [], {}
    saved = {}
    for ibin in bins_sel:
        sw = bin_sample_weight(binarr, g.nbins, ibin, npts, len(data[imgtype.GPDYN]))
        # ================= GAS
        log(f'bin {ibin}: gridding gas ...')
        ka, kn = grid_acc(tg[imgtype.GPDYN], data[imgtype.GPDYN], npts, MS, g.gplb, sw)
        k1 = renorm_oldway1(ka, kn)
        k0, den = renorm_oldway0(ka, kn, MS)
        k0k, _ = renorm_oldway0(ka, kn, MS, kluge=True)
        F1 = to_image(k1, MS)[cs, cs, cs]
        Fh = to_image(k0, MS)[cs, cs, cs]
        Fhk = to_image(k0k, MS)[cs, cs, cs]
        P1 = np.real(F1 * b)
        PFs = np.real(Fh * b)
        PFsk = np.real(Fhk * b)
        relerr = np.linalg.norm(P1 - prod_gas[ibin]) / np.linalg.norm(prod_gas[ibin])
        log(f'  P1 vs Tyger production gas: relerr {relerr:.2e}  (< 1e-3 = exact rebuild)')
        log(f'  Hermitian output imag/real: {np.abs(Fh.imag).max() / np.abs(Fh.real).max():.1e} (should be ~0)')
        log(f'  pocs ({a.iters} it) ...')
        xc, _ = pocs_pf(k1, kn, b, MS, IS, iters=a.iters, verbose=False)
        PFc = np.real(xc * b)
        cov = coverage(kn, den, tg[imgtype.GPDYN], MS)
        log(f'  gas coverage: holes inside k-sphere P1 {cov["frac_hole_P1"]:.3f} -> after Hermitian {cov["frac_hole_PF"]:.3f} '
            f'(filled {cov["filled"]:.2f} of holes)')
        if ibin == bins_sel[0]:
            m = lung_masks(np.abs(F1))
            log(f'  masks: lung {m["lung"].sum()} shell {m["shell"].sum()} corner {m["corner"].sum()}')
        # identity test: Hermitian == Re[F] on doubly-sampled cells?
        both = np.reshape((kn >= KNORM_FLOOR), -1)
        ReF1 = np.real(F1)
        ident = np.corrcoef(Fh.real[m['lung']], ReF1[m['lung']])[0, 1]
        # predicted PFs = Re[F].Re[b] + Im... : Fh real -> PFs = Fh * Re(b)
        pred = Fh.real * np.real(b)
        log(f'  identity: corr(Fh, Re[F]) in lung {ident:.4f}; PFs == Fh*Re(b): {np.allclose(PFs, pred)}')
        cosb = np.real(b)[m['lung']]
        log(f'  cos(angle b) in lung: mean {cosb.mean():.3f}, 5th pct {np.percentile(cosb, 5):.3f}')
        for name, img in (('P1', P1), ('PFs', PFs), ('PFs+k', PFsk), ('PFc', PFc), ('|F| (mag)', np.abs(F1))):
            r = dict(bin=ibin, itype='gas', variant=name, **metrics(img, m, ref=P1))
            r.update(relerr_P1_vs_tyger=relerr, frac_hole_P1=cov['frac_hole_P1'], frac_hole_PF=cov['frac_hole_PF'])
            rows.append(r)
        saved[f'gas_{ibin}'] = dict(P1=P1, PFs=PFs, PFsk=PFsk, PFc=PFc, Fh=Fh.real, ReF=ReF1, cov=cov['radial'])

        # ================= DISSOLVED
        log(f'bin {ibin}: gridding dissolved (dplb {g.dplb}) ...')
        ka, kn = grid_acc(tg[imgtype.DPDYN], data[imgtype.DPDYN], npts, MS, g.dplb, sw)
        k1 = renorm_oldway1(ka, kn)
        k0, den = renorm_oldway0(ka, kn, MS)
        D1 = to_image(k1, MS)[cs, cs, cs] * b           # complex dissolved image, Tyger's before the split
        Dh = to_image(k0, MS)[cs, cs, cs] * b           # PF version
        # effective k support of the dissolved data: radius where fwt = e^-4 (sample index 2*dplb)
        r_eff_idx = min(2 * g.dplb, npts - 1)
        ctr = MS / 2
        r_eff = np.linalg.norm(tg[imgtype.DPDYN][r_eff_idx] - ctr)
        cov_d = coverage(kn, den, tg[imgtype.DPDYN][:r_eff_idx + 1], MS)
        log(f'  dissolved k support radius ~{r_eff:.1f} cells (sample {r_eff_idx}); holes inside it P1 {cov_d["frac_hole_P1"]:.4f}')
        # LS map [Re D1, Im D1] -> production (aRBC, aTP) in the lung (convention-free)
        A = np.stack([D1.real[m['lung']], D1.imag[m['lung']]], 1)
        T = np.stack([prod_re[ibin][m['lung']], prod_im[ibin][m['lung']]], 1)
        Mx, *_ = np.linalg.lstsq(A, T, rcond=None)     # (2,2): T ~ A @ Mx
        def split(D):
            v = np.stack([D.real.ravel(), D.imag.ravel()], 1) @ Mx
            return v[:, 0].reshape(D.shape), v[:, 1].reshape(D.shape)
        R1, T1 = split(D1)
        Rh, Th = split(Dh)
        c_r = np.corrcoef(R1[m['lung']], prod_re[ibin][m['lung']])[0, 1]
        c_t = np.corrcoef(T1[m['lung']], prod_im[ibin][m['lung']])[0, 1]
        log(f'  our split vs production: corr aRBC {c_r:.4f}, aTP {c_t:.4f}; map angle {np.degrees(np.arctan2(Mx[1,0], Mx[0,0])):.1f} deg')
        # how much of the dissolved lung signal lives in the imaginary channel (what realness would destroy)
        im_frac = np.abs(D1.imag[m['lung']]).mean() / np.abs(D1[m['lung']]).mean()
        log(f'  dissolved: mean |Im D|/|D| in lung {im_frac:.3f} (0 = image already real)')
        for name, (R, Tt, D) in (('P1', (R1, T1, D1)), ('PFs', (Rh, Th, Dh))):
            for ch, img in (('aRBC', R), ('aTP', Tt), ('|D|', np.abs(D))):
                r = dict(bin=ibin, itype='dis_' + ch, variant=name, **metrics(img, m, ref={'aRBC': R1, 'aTP': T1, '|D|': np.abs(D1)}[ch]))
                r.update(frac_hole_P1=cov_d['frac_hole_P1'], frac_hole_PF=cov_d['frac_hole_PF'], im_frac_lung=im_frac,
                         corr_split_vs_prod_rbc=c_r, corr_split_vs_prod_tp=c_t)
                rows.append(r)
        saved[f'dis_{ibin}'] = dict(R1=R1, T1=T1, Rh=Rh, Th=Th, absD1=np.abs(D1), absDh=np.abs(Dh),
                                    prodR=prod_re[ibin], prodT=prod_im[ibin])

    # ---- write metrics
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(out, 'metrics.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.5g}' if isinstance(v, float) else v) for k, v in r.items()})
    np.savez_compressed(os.path.join(out, 'arrays.npz'), lung=m['lung'], shell=m['shell'], corner=m['corner'],
                        **{f'{k}_{n}': v for k, d in saved.items() for n, v in d.items()})

    # ---- figures
    z = int(np.round(ndi.center_of_mass(m['lung'])[0]))
    y = int(np.round(ndi.center_of_mass(m['lung'])[1]))
    for ibin in bins_sel:
        s = saved[f'gas_{ibin}']
        vmax = np.percentile(s['P1'][m['lung']], 99)
        fig, ax = plt.subplots(2, 4, figsize=(16, 8))
        for j, (nm, key) in enumerate((('P1 (Tyger)', 'P1'), ('PFs (Steve oldway 0)', 'PFs'), ('PFs + KLUGE', 'PFsk'), ('PFc (POCS, same b)', 'PFc'))):
            ax[0, j].imshow(s[key][z], cmap='gray', vmin=-0.15 * vmax, vmax=vmax)
            ax[0, j].set_title(nm)
            d = s[key][z] - s['P1'][z]
            ax[1, j].imshow(d, cmap='RdBu_r', vmin=-0.3 * vmax, vmax=0.3 * vmax)
            ax[1, j].set_title('minus P1' if j else 'coronal slice z=%d' % z)
        for a_ in ax.ravel():
            a_.axis('off')
        fig.suptitle(f'{a.key} gas bin {ibin}: same gridded k-space, four renorm/recon variants (top), difference to P1 (bottom, ±30 % of lung P99)')
        fig.tight_layout()
        fig.savefig(os.path.join(out, f'gas_bin{ibin}.png'), dpi=110)
        plt.close(fig)

        d = saved[f'dis_{ibin}']
        vr = np.percentile(np.abs(d['prodR'][m['lung']]), 99)
        vt = np.percentile(np.abs(d['prodT'][m['lung']]), 99)
        fig, ax = plt.subplots(2, 3, figsize=(12, 8))
        for j, (nm, key, v) in enumerate((('production aRBC', 'prodR', vr), ('P1 aRBC (ours)', 'R1', vr), ('PFs aRBC', 'Rh', vr))):
            ax[0, j].imshow(d[key][z], cmap='gray', vmin=-0.2 * v, vmax=v)
            ax[0, j].set_title(nm)
        for j, (nm, key, v) in enumerate((('production aTP', 'prodT', vt), ('P1 aTP (ours)', 'T1', vt), ('PFs aTP', 'Th', vt))):
            ax[1, j].imshow(d[key][z], cmap='gray', vmin=-0.2 * v, vmax=v)
            ax[1, j].set_title(nm)
        for a_ in ax.ravel():
            a_.axis('off')
        fig.suptitle(f'{a.key} dissolved bin {ibin}: RBC / TP split of the complex image — P1 vs Hermitian (PF) k-space')
        fig.tight_layout()
        fig.savefig(os.path.join(out, f'dis_bin{ibin}.png'), dpi=110)
        plt.close(fig)

    # coverage vs k radius (gas, first bin) — one axis, two series, direct labels
    cov = saved[f'gas_{bins_sel[0]}']['cov']
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(cov[:, 0], cov[:, 1], color=C_SERIES[0], lw=2, marker='o', ms=5, label='P1 (measured cells only)')
    ax.plot(cov[:, 0], cov[:, 2], color=C_SERIES[1], lw=2, marker='o', ms=5, label='after Hermitian fill (oldway 0)')
    ax.set_xlabel('|k| (grid cells, MS = 240)')
    ax.set_ylabel('fraction of cells with no data')
    ax.set_title(f'{a.key} gas bin {bins_sel[0]}: per-bin k-space holes vs radius')
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'gas_coverage.png'), dpi=120)
    plt.close(fig)

    # ---- report
    with open(os.path.join(out, 'report.md'), 'w') as f:
        f.write(f'# PF test — {a.key} (bins {bins_sel}, POCS {a.iters} it)\n\n')
        f.write('| bin | type | variant | lung mean | bg σ | SNR | shell σ / bg σ | edge grad | corr lung vs P1 | bg noise corr vs P1 |\n|---|---|---|---|---|---|---|---|---|---|\n')
        for r in rows:
            f.write(f"| {r['bin']} | {r['itype']} | {r['variant']} | {r['lung_mean']:.4g} | {r['bg_sd']:.4g} | {r['snr']:.1f} | "
                    f"{r['shell_over_bg']:.2f} | {r['edge_grad']:.4f} | {r.get('corr_lung_vs_P1', float('nan')):.4f} | "
                    f"{r.get('bg_noise_corr_vs_P1', float('nan')):.3f} |\n")
        f.write('\nsee metrics.csv for hole fractions, split correlations, im_frac.\n')
    log('done ->', out)


if __name__ == '__main__':
    main()

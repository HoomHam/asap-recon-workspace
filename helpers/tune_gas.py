"""Tyger gas-image knob tuning against a known object (step 2 of notes/Two_Targets_Registration_vs_Truth.md).

Why: Steve's live knobs (MS, gplb) and hardcoded ones (kdist0sq 0.2, box -2..+1, KLUGE floor 1e-5) were set by eye.
pf_inlung.py resamples the truth ON the grid, so it cannot judge grid/kernel choices. Here the object is sampled at
the real trajectory points (finufft type 2), so every knob is judged by the same physics:

  object   T = real(F_all . b_prod), all data, production kernel, NO readout apodisation, zeroed outside a 12-vox
             lung dilation (a realistic lung spectrum; the same object for every setting)
  sampling the real unique trajectory, the real per-sample soft-bin weights of each bin (production DIAPHRAGM nav,
             bug-faithful ilvnum rounding) and the real exclude mask
  noise    complex white noise per acquired sample; sigma calibrated once so the production setting reproduces the
             real bin's far-background sigma (the object's own aliasing removed in quadrature)
  recon    Steve's normalized-convolution gridding collapsed exactly onto the unique samples (the trajectory repeats
             every nuniq samples, gridding is linear: sum over repeats of binwt*data and of binwt) -> one sparse
             kernel matrix per setting, many recons per matrix

Per setting and bin, inside the lung, relative to the object's lung mean:
  blur   |R_all - T|            all data, noise-free: resolution + rolloff + KLUGE bias of the setting itself
  alias  |R_bin - R_all|        the bin's own k-space holes (noise-free)
  noise  sd of the noise-only recon
  total  |R_bin + noise - T|    one draw, common random numbers across settings -- the T-image criterion
  PSF    FWHM (vox) and gain of delta probes at 4 lung voxels (all data and bin 0)
Real data: the production setting must rebuild Tyger's gas_phase to relerr < 1e-3 (asserted), and every setting's
real bin-0 image goes into a montage for the eye.

usage: python tune_gas.py --key 2024-11-13_025JC [--bins 0,4,8,12] [--only prod,gplb450]
out:   workspace/outputs/tune_gas/<key>/{metrics.csv, report.md, real_montage.png, tradeoff.png, run.log}
       image stacks montage_bin0.npz -> Ext /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tune_gas/<key>/
"""
import argparse
import csv
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy import fft as sfft
from scipy import ndimage as ndi
from scipy.io import loadmat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pf_test import (ASAP, RAW_TREE, PROD_TREE, read_input_mrd, exclude_mask, bin_sample_weight,  # noqa: E402
                     crop_slices, lung_masks, gvar, imgtype, graddir, traj, steve_raw_cls, KNORM_FLOOR)
from steve_kernel_numpy import EPS, steve_rolloff  # noqa: E402
import finufft  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

OUT = os.path.join(ASAP, 'workspace/outputs/tune_gas')             # small files (csv, md, png) stay on the laptop
EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tune_gas'    # image stacks (npz) -> Ext mirror (big-output rule)
OUTBOX = os.path.expanduser('~/Hooman/Outbox/Work/2026_ASAP_Recon/tune_gas')   # Ext unplugged -> porter relocates
STEVE_BOX = (-2, -1, 0, 1)                   # recon.py: range(cx-2, cx+2), asymmetric
SYM2, SYM3 = tuple(range(-2, 3)), tuple(range(-3, 4))

# name: (MS, kdist0sq [grid cells^2], box offsets, gplb [0 = off])
SETTINGS = {
    'prod':      (240, 0.2, STEVE_BOX, 300),
    'gplb0':     (240, 0.2, STEVE_BOX, 0),
    'gplb150':   (240, 0.2, STEVE_BOX, 150),
    'gplb200':   (240, 0.2, STEVE_BOX, 200),
    'gplb450':   (240, 0.2, STEVE_BOX, 450),
    'gplb700':   (240, 0.2, STEVE_BOX, 700),
    'k0.5':      (240, 0.5, SYM2, 300),
    'k1.0':      (240, 1.0, SYM3, 300),
    'k2.0':      (240, 2.0, SYM3, 300),
    'MS160':     (160, 0.2, STEVE_BOX, 300),
    'MS200':     (200, 0.2, STEVE_BOX, 300),
    'MS160k0.5': (160, 0.5, SYM2, 300),
    'MS200k0.5': (200, 0.5, SYM2, 300),
    'gplb100':   (240, 0.2, STEVE_BOX, 100),
    'MS160g450': (160, 0.2, STEVE_BOX, 450),
    'MS200g450': (200, 0.2, STEVE_BOX, 450),
    'MS160g700': (160, 0.2, STEVE_BOX, 700),
}
GROUPS = {   # --set: named subsets (default 'main' = the 13 settings of the first sweep)
    'main': ['prod', 'gplb0', 'gplb150', 'gplb200', 'gplb450', 'gplb700', 'k0.5', 'k1.0', 'k2.0', 'MS160', 'MS200',
             'MS160k0.5', 'MS200k0.5'],
    'stab': ['prod', 'gplb200', 'gplb450', 'MS160', 'MS200'],
    'snr':  ['gplb100', 'gplb150', 'gplb200', 'prod', 'gplb450', 'gplb700', 'MS160', 'MS200'],
    'reg':  ['prod', 'gplb450', 'gplb700', 'gplb0', 'MS160', 'MS200', 'MS160g450', 'MS200g450', 'MS160g700'],
}
# --dp: the complex dissolved image D = F.b (before the RBC/TP rotation). Only grid knobs: dplb is NOT swept because
# the static object has no RBC-TP drift along the readout (F92) -- a dplb sweep here would ignore the Dixon mixing.
DP_SETTINGS = {
    'prod':  (240, 0.2, STEVE_BOX, 40),
    'MS200': (200, 0.2, STEVE_BOX, 40),
    'MS160': (160, 0.2, STEVE_BOX, 40),
    'k0.5':  (240, 0.5, SYM2, 40),
}


def log(*a):
    msg = time.strftime('%H:%M:%S') + ' ' + ' '.join(str(x) for x in a)
    print(msg, flush=True)
    if LOGF:
        LOGF.write(msg + '\n')
        LOGF.flush()


LOGF = None


# ---------------------------------------------------------------- gridding, collapsed onto the unique samples
def bin_weight(binarr, nbins, ibin, npts, nsamples, bd0sq):
    """selftest bin_sample_weight with bindist0sq as a parameter (2.0 = recon.py); bug-faithful ilvnum rounding."""
    idx = np.arange(nsamples)
    ilvnum = np.minimum((idx + npts // 2) // npts, len(binarr) - 1)
    binctr = binarr[ilvnum]
    bindist = np.minimum.reduce([np.abs(ibin - binctr), np.abs(ibin + nbins - binctr), np.abs(ibin - nbins - binctr)])
    w = np.exp(-bindist * bindist / bd0sq)
    w[binctr < 0] = 0.0
    return w


def hp(v, s=2.0):
    return v - ndi.gaussian_filter(v, s)


def build_A(kdk, MS, IS, kd2, offs):
    """Sparse kernel matrix (MS^3 x nuniq), CSC, entries exp(-d^2/kd2) for Steve's cell neighbourhood.
    kdk: (nuniq, 3) trajectory in delta-k units, columns x, y, z. Flat cell index ix*MS^2 + iy*MS + iz (Steve)."""
    g = kdk * (MS / IS) + MS / 2                                  # raw.rescale_to_MS
    c = (g + 0.5).astype(np.int64)                                # cudarecon: int(k + 0.5)
    n, no = len(g), len(offs) ** 3
    rows = np.zeros((n, no), dtype=np.int32)
    vals = np.zeros((n, no), dtype=np.float64)
    j = 0
    for dx in offs:
        ix = c[:, 0] + dx
        for dy in offs:
            iy = c[:, 1] + dy
            for dz in offs:
                iz = c[:, 2] + dz
                v = (ix >= 0) & (ix < MS) & (iy >= 0) & (iy < MS) & (iz >= 0) & (iz < MS)
                dsq = (g[:, 0] - ix) ** 2 + (g[:, 1] - iy) ** 2 + (g[:, 2] - iz) ** 2
                rows[v, j] = (ix[v] * MS * MS + iy[v] * MS + iz[v]).astype(np.int32)
                vals[v, j] = np.exp(-dsq[v] / kd2)
                j += 1
    return sp.csc_matrix((vals.ravel(), rows.ravel(), np.arange(0, n * no + 1, no)), shape=(MS ** 3, n))


def grid_image(A, nums, den, MS, IS):
    """nums: list of complex (nuniq,) numerators (binwt*apod*data summed over repeats), den: (nuniq,) sum of binwt.
    Returns list of cropped complex images (Steve (z,y,x) layout) in object units (Steve's fftn image / MS^3, so
    images at different MS are directly comparable and a unit delta has gain ~1) and knorm."""
    kn = A @ den
    hole = kn < KNORM_FLOOR
    inv = np.zeros_like(kn)
    inv[kn > 0] = 1.0 / kn[kn > 0]
    inv[hole] = 0.0
    X = np.column_stack([v for z in nums for v in (z.real, z.imag)])
    acc = A @ X
    cs = crop_slices(MS, IS)
    out = []
    for i in range(len(nums)):
        k = (acc[:, 2 * i] + 1j * acc[:, 2 * i + 1]) * inv
        k3 = np.reshape(k, (MS, MS, MS), order='F')
        out.append(sfft.fftshift(sfft.fftn(sfft.ifftshift(k3), workers=-1))[cs, cs, cs] / MS ** 3)
    return out, kn


def hole_fraction(kn, kdk, MS, IS):
    rmax = np.max(np.linalg.norm(kdk, axis=1)) * MS / IS
    ax = np.arange(MS) - MS / 2
    r2 = ax[:, None, None] ** 2 + ax[None, :, None] ** 2 + ax[None, None, :] ** 2
    inside = np.reshape(r2, -1, order='F') <= rmax ** 2
    return float(np.mean(kn[inside] < KNORM_FLOOR))


# ---------------------------------------------------------------- PSF
def fwhm_1d(p, i0):
    half = p[i0] / 2
    if half <= 0:
        return np.nan
    j = i0
    while j + 1 < len(p) and p[j + 1] > half:
        j += 1
    if j + 1 >= len(p):
        return np.nan
    r = j + (p[j] - half) / (p[j] - p[j + 1])
    j = i0
    while j - 1 >= 0 and p[j - 1] > half:
        j -= 1
    if j - 1 < 0:
        return np.nan
    l_ = j - (p[j] - half) / (p[j] - p[j - 1])
    return r - l_


def psf_stats(img, pt):
    a, b, c = pt
    f = [fwhm_1d(img[:, b, c], a), fwhm_1d(img[a, :, c], b), fwhm_1d(img[a, b, :], c)]
    return float(np.nanmean(f)), float(img[a, b, c]), f


def probe_points(lung):
    idx = np.argwhere(lung)
    mz, mx = np.median(idx[:, 0]), np.median(idx[:, 2])
    pts = []
    for sz in (idx[:, 0] < mz, idx[:, 0] >= mz):
        for sx in (idx[:, 2] < mx, idx[:, 2] >= mx):
            q = idx[sz & sx]
            ctr = q.mean(0)
            pts.append(tuple(int(v) for v in q[np.argmin(np.sum((q - ctr) ** 2, 1))]))
    return pts


# ---------------------------------------------------------------- main
def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--key', default='2024-11-13_025JC')
    ap.add_argument('--bins', default='0,4,8,12')
    ap.add_argument('--only', default='', help='comma list of setting names (default: all)')
    ap.add_argument('--cg-control', action='store_true',
                    help='instead of the sweep: weighted least-squares CG (finufft, real-constrained) on the SAME simulated '
                         'samples -> is the blur floor the data or Steve\'s gridding? writes cg_control.csv')
    ap.add_argument('--set', default='main', help='setting group: ' + ', '.join(GROUPS))
    ap.add_argument('--bindist0sq', type=float, default=2.0, help='soft-bin width (recon.py: 2)')
    ap.add_argument('--nbins', type=int, default=0, help='override nbins (0 = header value, 16)')
    ap.add_argument('--sigma', type=float, default=0.0, help='fixed per-sample noise sigma (0 = calibrate on real bins)')
    ap.add_argument('--noise-mult', type=float, default=1.0, help='scale the noise (lower-SNR test)')
    ap.add_argument('--dp', action='store_true', help='dissolved image (complex D = F.b), grid knobs only (DP_SETTINGS)')
    ap.add_argument('--obj-gplb', type=int, default=0,
                    help='readout apodisation used to BUILD the object T (0 = none: T keeps ~1/2.5 of a bin\'s noise as fake '
                         'detail -> favours sharp settings; 300 = production: T loses real detail -> favours smooth ones)')
    a = ap.parse_args()
    tag = ('_dp' if a.dp else '') + (f'_obj{a.obj_gplb}' if a.obj_gplb else '') + (f'_{a.set}' if a.set != 'main' else '') \
        + (f'_bd{a.bindist0sq:g}' if a.bindist0sq != 2.0 else '') + (f'_nb{a.nbins}' if a.nbins else '') \
        + (f'_nm{a.noise_mult:g}' if a.noise_mult != 1.0 else '')
    out = os.path.join(OUT, a.key + tag)
    if a.only or len(a.bins.split(',')) < 4:                        # partial runs never overwrite a full sweep
        out = os.path.join(out, 'partial_' + (a.only.replace(',', '+') or 'all') + '_bins' + a.bins.replace(',', '-'))
    os.makedirs(out, exist_ok=True)
    LOGF = open(os.path.join(out, 'run.log'), 'a')
    bins = [int(x) for x in a.bins.split(',')]
    S = DP_SETTINGS if a.dp else SETTINGS
    names = [s for s in (a.only.split(',') if a.only else (S if a.dp else GROUPS[a.set]))]
    itype = imgtype.DPDYN if a.dp else imgtype.GPDYN
    gas = not a.dp

    # ---- load exactly like pf_test / tyger_recon
    header, arrs = read_input_mrd(os.path.join(RAW_TREE, a.key, 'd', 'input.mrd'))
    prod = loadmat(os.path.join(PROD_TREE, a.key, 'd', 'recon.mat'))
    ul = {p.name: p.value for p in header.user_parameters.user_parameter_long}
    ud = {p.name: p.value for p in header.user_parameters.user_parameter_double}
    g = gvar()
    g.MS, g.IS, g.nbins, g.gplb = int(ul.get('MS', 240)), int(ul.get('IS', 100)), int(ul.get('nbins', 16)), int(ul.get('gplb', 300))
    assert (g.MS, g.IS, g.gplb) == (SETTINGS['prod'][0], 100, SETTINGS['prod'][3]), (g.MS, g.IS, g.gplb)
    g.dplb = int(ul.get('dplb', 40))
    if a.nbins:
        g.nbins = a.nbins
    std_binning = (g.nbins == 16 and a.bindist0sq == 2.0)
    lb0 = g.dplb if a.dp else g.gplb
    assert lb0 == S['prod'][3], (lb0, S['prod'][3])
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
    MS0, IS, npts = g.MS, g.IS, g_raw.npts
    tg = np.stack([np.asarray(g_traj.gettraj(itype, d)) for d in (graddir.X, graddir.Y, graddir.Z)], 1)
    kdk = (tg - MS0 / 2) * IS / MS0                                  # delta-k units, columns x y z
    nuniq = len(kdk)
    data = np.reshape(np.asarray(g_raw.getimg(itype))[:, 0, :], npts * g_raw.ntotalilvs, order='F') \
        * exclude_mask(g_raw)
    nsamp = len(data)
    assert nuniq % npts == 0
    u = np.arange(nsamp) % nuniq
    nidx = np.arange(nsamp) % npts
    # cudarecon skips |raw * fwt| < EPS in BOTH k and knorm: with dplb 40 that drops every sample past n ~ 220
    keep = (np.abs(data) >= EPS) & (np.abs(data) * np.exp(-(nidx / lb0) ** 2 if lb0 > EPS else 0) >= EPS)
    log(f'{"DISSOLVED" if a.dp else "GAS"} {a.key}: MS={MS0} IS={IS} nbins={g.nbins} gplb={g.gplb} npts={npts} nilv={g_raw.ntotalilvs} '
        f'nuniq={nuniq} ({nuniq // npts} ilv, {nsamp / nuniq:.2f} repeats) dtdyn={meta["dtdyn"] * 1e6:.2f} us '
        f'TE={meta["TE"] * 1e3:.3f} ms kmax={np.max(np.linalg.norm(kdk, axis=1)):.2f} dk; kept {keep.mean():.4f}')
    amin = np.abs(data[keep]).min()
    lb_min = npts / np.sqrt(np.log(amin / EPS))
    log(f'min |data| kept {amin:.3g}: LB >= {lb_min:.0f} keeps every sample (EPS rule); last kept sample index {nidx[keep].max()}')
    for n_ in names:                                   # every swept LB must keep the same sample set as production
        assert S[n_][3] == lb0 or S[n_][3] == 0 or (S[n_][3] >= lb_min and lb0 >= lb_min), (n_, S[n_][3])

    # ---- bins (production DIAPHRAGM nav) and the collapsed sums per bin
    ilvvol = np.asarray(prod['nav_volume'], dtype=float).ravel()
    assert len(ilvvol) == g_raw.ntotalilvs, (len(ilvvol), g_raw.ntotalilvs)
    binarr = np.ascontiguousarray(np.asarray(g_raw.bin(ilvvol)) * g.nbins)
    col = {}
    for key in ['all'] + bins:
        sw = keep.astype(float) if key == 'all' else bin_weight(binarr, g.nbins, key, npts, nsamp, a.bindist0sq) * keep
        W = np.bincount(u, weights=sw, minlength=nuniq)
        W2 = np.bincount(u, weights=sw * sw, minlength=nuniq)
        D = np.bincount(u, weights=sw * data.real, minlength=nuniq) + 1j * np.bincount(u, weights=sw * data.imag, minlength=nuniq)
        col[key] = dict(W=W, W2=W2, D=D)
    log('collapsed sums done; bin weight totals', {k: round(float(v['W'].sum() / npts), 1) for k, v in col.items()})
    sidx = np.arange(nuniq) % npts                                  # sample index along the readout (after killpts)

    def apod(gplb):
        return np.exp(-(sidx / gplb) ** 2) if gplb > EPS else np.ones(nuniq)

    b = np.asarray(prod['calcb_b'])[0].astype(np.complex128)
    b /= np.maximum(np.abs(b), 1e-12)
    prod_gas = np.asarray(prod['gas_phase'])
    proj = np.real if gas else (lambda x: x)                   # gas T image = real(F.b); dissolved = complex F.b

    # ---- production kernel: exact-rebuild check, masks, object T
    ms, kd2, offs, gplb0 = S['prod']
    log('building production kernel ...')
    A = build_A(kdk, ms, IS, kd2, offs)
    ap0 = apod(gplb0)
    imgs, _ = grid_image(A, [ap0 * col[bins[0]]['D']], col[bins[0]]['W'], ms, IS)
    obj_lb = a.obj_gplb if (a.obj_gplb or gas) else lb0          # dissolved: an un-apodised object is pure noise
    (F_all_ap, F_all_noap), _ = grid_image(A, [ap0 * col['all']['D'], apod(obj_lb) * col['all']['D']], col['all']['W'], ms, IS)
    m = lung_masks(np.abs(F_all_ap))
    lung, corner = m['lung'], m['corner']
    if not std_binning:
        relerr = float('nan')
        log(f'non-production binning (nbins {g.nbins}, bindist0sq {a.bindist0sq:g}): no Tyger rebuild check')
    elif gas:
        P1 = np.real(imgs[0] * b) * ms ** 3                      # back to Tyger's units for the check
        relerr = np.linalg.norm(P1 - prod_gas[bins[0]]) / np.linalg.norm(prod_gas[bins[0]])
        log(f'CHECK production rebuild bin {bins[0]} vs Tyger gas_phase: relerr {relerr:.2e}')
        assert relerr < 1e-3, 'collapsed gridding does not reproduce Tyger'
    else:                                                        # production stores the RBC/TP rotation of D: LS-map
        D1 = imgs[0] * b
        Am = np.stack([D1.real[lung], D1.imag[lung]], 1)
        Tm = np.stack([np.asarray(prod['dissolved_phase_real'])[bins[0]][lung], np.asarray(prod['dissolved_phase_imag'])[bins[0]][lung]], 1)
        Mx, *_ = np.linalg.lstsq(Am, Tm, rcond=None)
        fit = Am @ Mx
        cr = [np.corrcoef(fit[:, i], Tm[:, i])[0, 1] for i in (0, 1)]
        relerr = 1 - min(cr)
        log(f'CHECK production rebuild bin {bins[0]} vs Tyger aRBC/aTP (LS 2x2 map): corr {cr[0]:.6f} / {cr[1]:.6f}')
        assert min(cr) > 0.9999, 'collapsed dissolved gridding does not reproduce Tyger'
    T = proj(F_all_noap * b) * ndi.binary_dilation(lung, iterations=12)
    muT = np.real(T[lung]).mean() if gas else np.abs(T[lung]).mean()
    log(f'object: lung {lung.sum()} vox, T lung mean {muT:.4g}; corner {corner.sum()} vox')

    # ---- forward model: s_u = sum_x T(x) exp(+2 pi i k_u . x / IS); image axes (z,y,x) <-> k columns (2,1,0)
    xyz = [np.ascontiguousarray(2 * np.pi * kdk[:, c] / IS) for c in (2, 1, 0)]
    s_T = finufft.nufft3d2(*xyz, T.astype(np.complex128), isign=1, eps=1e-10, modeord=0)
    pt = (IS // 2 + 7, IS // 2 - 5, IS // 2 + 3)
    dlt = np.zeros((IS,) * 3, complex)
    dlt[pt] = 1
    s_d = finufft.nufft3d2(*xyz, dlt, isign=1, eps=1e-10, modeord=0)
    x_pt = np.array([pt[2], pt[1], pt[0]]) - IS / 2                 # x, y, z pixel offsets
    s_exact = np.exp(2j * np.pi * (kdk @ x_pt) / IS)
    log(f'CHECK finufft axis/sign mapping vs analytic delta: max err {np.max(np.abs(s_d - s_exact)):.1e}')
    assert np.max(np.abs(s_d - s_exact)) < 1e-6
    pts = probe_points(lung)
    s_probe = [np.exp(2j * np.pi * (kdk @ (np.array([p[2], p[1], p[0]]) - IS / 2)) / IS) for p in pts]
    log(f'probe points (z,y,x): {pts}')

    # ---- noise calibration at the production setting (one sigma per acquired sample, all bins)
    rng = np.random.default_rng(12345)
    unit = {k: (rng.standard_normal(nuniq) + 1j * rng.standard_normal(nuniq)) for k in bins}   # common random numbers
    sig = []
    for ib in (bins if not a.sigma else []):
        c_ = col[ib]
        (Rr, Rt, Rn), _ = grid_image(A, [ap0 * c_['D'], ap0 * s_T * c_['W'], ap0 * unit[ib] * np.sqrt(c_['W2'])], c_['W'], ms, IS)
        sd_real = np.real(Rr * b)[corner].std()
        sd_al = np.real(Rt)[corner].std()
        sd_unit = np.real(Rn)[corner].std()
        sig.append(np.sqrt(max(sd_real ** 2 - sd_al ** 2, 0)) / sd_unit)
        log(f'  noise cal bin {ib}: real corner sd {sd_real:.4g}, object-alias sd {sd_al:.4g}, unit {sd_unit:.4g} -> sigma {sig[-1]:.4g}')
    sigma = float(np.median(sig)) if not a.sigma else a.sigma
    log(f'per-sample noise sigma (each of re/im) = {sigma:.4g}' + (f'  (spread {np.std(sig) / sigma:.1%})' if sig else ' (given)'))
    (Tn,), _ = grid_image(A, [apod(obj_lb) * unit[bins[0]] * sigma * np.sqrt(col['all']['W2'])], col['all']['W'], ms, IS)
    log(f'object T (obj LB {obj_lb}) carries noise sd {np.real(Tn)[lung].std() / muT:.3f} of its lung mean as fake detail')
    del A
    sig_run = sigma * a.noise_mult
    if a.noise_mult != 1.0:
        log(f'noise x{a.noise_mult:g}: sweep sigma {sig_run:.4g}')
    reg = ndi.binary_dilation(lung, iterations=3)                  # registration region: lung + boundary
    if a.cg_control:
        assert gas, 'cg control is gas-only'
        cg_control(a, out, kdk, xyz, col, bins, s_T, T, lung, muT, unit, sigma)
        return

    # ---- sweep
    rows, montage = [], {}
    for name in names:
        ms, kd2, offs, gp = S[name]
        t0 = time.time()
        A = build_A(kdk, ms, IS, kd2, offs)
        apv = apod(gp)
        roll = steve_rolloff(ms, IS, kd2)
        sig_dk = np.sqrt(kd2 / 2) * IS / ms
        ca = col['all']
        nums = [apv * s_T * ca['W']] + [apv * sp_ * ca['W'] for sp_ in s_probe]
        imgs, kn_all = grid_image(A, nums, ca['W'], ms, IS)
        R_all = proj(imgs[0])
        rms = lambda d: np.sqrt(np.mean(np.abs(d[lung]) ** 2)) / muT
        lmean = (lambda x: x[lung].mean()) if gas else (lambda x: np.abs(x[lung]).mean())
        for dea in (False, True):
            div = roll if dea else 1.0
            Ra = R_all / div
            psf_all = [psf_stats(np.real(im) / div, p) for im, p in zip(imgs[1:], pts)]
            blur = rms(Ra - T)
            bias_all = lmean(Ra) / muT
            for ib in bins:
                c_ = col[ib]
                if not dea:
                    nb = [apv * s_T * c_['W'], apv * unit[ib] * sig_run * np.sqrt(c_['W2']), apv * c_['D']]
                    if ib == bins[0]:
                        nb += [apv * sp_ * c_['W'] for sp_ in s_probe]
                    ims, kn = grid_image(A, nb, c_['W'], ms, IS)
                    cache = dict(R=proj(ims[0]), N=proj(ims[1]), real=proj(ims[2] * b), holes=hole_fraction(kn, kdk, ms, IS))
                    if ib == bins[0]:
                        cache['psf_imgs'] = [np.real(im) for im in ims[3:]]
                    col[ib]['cache'] = cache
                cache = col[ib]['cache']
                R, N, Rreal = cache['R'] / div, cache['N'] / div, cache['real'] / div
                row = dict(setting=name, deapod=int(dea), MS=ms, kdist0sq=kd2, sigma_kernel_dk=round(sig_dk, 4),
                           box=len(offs), gplb=gp, bin=ib, holes=cache['holes'],
                           blur=blur, bias_all=bias_all,
                           alias=rms(R - Ra), nf_err=rms(R - T), noise=N[lung].std() / muT, total=rms(R + N - T),
                           bias_bin=lmean(R) / muT,
                           real_snr=lmean(Rreal) / np.real(Rreal)[corner].std(),
                           fwhm_all=np.nanmean([p[0] for p in psf_all]),
                           gain_all=np.mean([p[1] for p in psf_all]))
                if 'psf_imgs' in cache:
                    psf0 = [psf_stats(im / div, p) for im, p in zip(cache['psf_imgs'], pts)]
                    row['fwhm_bin0'] = np.nanmean([p[0] for p in psf0])
                    row['gain_bin0'] = np.mean([p[1] for p in psf0])
                rows.append(row)
                if ib == bins[0] and not dea:
                    montage[name] = Rreal if gas else np.abs(Rreal)
            # registration proxies: anatomy is identical in every bin (static object), so fine texture that differs
            # between bins is what a registration would chase; texture shared with the all-data image is an anchor
            hp_all = hp(np.real(Ra))[reg]
            hpb = [hp(np.real(col[ib]['cache']['R'] + col[ib]['cache']['N']) / (roll if dea else 1.0))[reg] for ib in bins]
            pair = np.mean([np.corrcoef(hpb[i], hpb[j])[0, 1] for i in range(len(bins)) for j in range(i + 1, len(bins))])
            anat = np.mean([np.corrcoef(h, hp_all)[0, 1] for h in hpb])
            incons = np.mean([np.sqrt(np.mean((h - hp_all) ** 2)) for h in hpb]) / muT
            gmag = np.sqrt(sum(np.gradient(np.real(Ra), axis=ax_) ** 2 for ax_ in range(3)))
            edge = gmag[m['ring']].mean() / muT
            for r_ in rows:
                if r_['setting'] == name and r_['deapod'] == int(dea):
                    r_.update(hp_pair_corr=pair, hp_anat_corr=anat, hp_incons=incons, edge=edge, anchor_snr=edge / incons)
        r0 = [r for r in rows if r['setting'] == name and r['deapod'] == 0]
        log(f'{name:10s} MS={ms} kd2={kd2} (sigma {sig_dk:.3f} dk) gplb={gp}: FWHM all {r0[0]["fwhm_all"]:.2f} '
            f'bin0 {r0[0].get("fwhm_bin0", np.nan):.2f} vox | blur {r0[0]["blur"]:.3f} alias {np.mean([r["alias"] for r in r0]):.3f} '
            f'noise {np.mean([r["noise"] for r in r0]):.3f} total {np.mean([r["total"] for r in r0]):.3f} | holes '
            f'{np.mean([r["holes"] for r in r0]):.3f} | real SNR {np.mean([r["real_snr"] for r in r0]):.1f} ({time.time() - t0:.0f} s)')
        del A
        for ib in bins:
            col[ib].pop('cache', None)

    # ---- write
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(out, 'metrics.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.5g}' if isinstance(v, float) else v) for k, v in r.items()})
    big = os.path.join(EXT if os.path.isdir(os.path.dirname(EXT)) else OUTBOX, os.path.relpath(out, OUT))
    os.makedirs(big, exist_ok=True)
    np.savez_compressed(os.path.join(big, 'montage_bin%d.npz' % bins[0]), lung=lung, T=T.astype(np.float32),
                        **{k: v.astype(np.float32) for k, v in montage.items()})
    log('image stacks ->', big)
    summarize(rows, names, out, a.key + (' DISSOLVED (complex F.b, rms of |error|)' if a.dp else ' gas'), bins, montage, lung, sigma, relerr)
    log('done ->', out)


def cg_control(a, out, kdk, xyz, col, bins, s_T, T, lung, muT, unit, sigma, iters=(2, 5, 10, 20, 40, 80)):
    """Weighted LS: min sum_u W_u |A x - y|^2 over REAL x (IS^3), normal equations Re(A^H W A) x = Re(A^H D), plain CG
    from 0 (early stopping = the only regulariser). Noise-free all data / bin 0, and bin 0 with the same noise draw."""
    IS = T.shape[0]
    p2 = finufft.Plan(2, (IS, IS, IS), isign=1, eps=1e-7)
    p1 = finufft.Plan(1, (IS, IS, IS), isign=-1, eps=1e-7)
    p2.setpts(*xyz)
    p1.setpts(*xyz)

    def N(x, W):
        return np.real(p1.execute(W * p2.execute(x.astype(np.complex128))))

    b0 = bins[0]
    cases = [('all, noise-free', col['all']['W'], s_T * col['all']['W']),
             (f'bin {b0}, noise-free', col[b0]['W'], s_T * col[b0]['W']),
             (f'bin {b0}, matched noise', col[b0]['W'], s_T * col[b0]['W'] + unit[b0] * sigma * np.sqrt(col[b0]['W2']))]
    rows = []
    for name, W, D in cases:
        rhs = np.real(p1.execute(D))
        x = np.zeros_like(rhs)
        r = rhs.copy()
        pdir = r.copy()
        rs = np.vdot(r, r).real
        t0 = time.time()
        for it in range(1, max(iters) + 1):
            Ap = N(pdir, W)
            alpha = rs / np.vdot(pdir, Ap).real
            x += alpha * pdir
            r -= alpha * Ap
            rs_new = np.vdot(r, r).real
            pdir = r + (rs_new / rs) * pdir
            rs = rs_new
            if it in iters:
                e = np.sqrt(np.mean((x - T)[lung] ** 2)) / muT
                rows.append(dict(case=name, iters=it, rms_err=e, bias=x[lung].mean() / muT))
                log(f'  CG {name:24s} it {it:3d}: rms err {e:.4f}  lung mean ratio {x[lung].mean() / muT:.4f}  ({time.time() - t0:.0f} s)')
    with open(os.path.join(out, 'cg_control.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r_ in rows:
            w.writerow({k: (f'{v:.5g}' if isinstance(v, float) else v) for k, v in r_.items()})


def summarize(rows, names, out, key, bins, montage, lung, sigma, relerr):
    def agg(name, dea, k):
        v = [r[k] for r in rows if r['setting'] == name and r['deapod'] == dea and k in r]
        return float(np.mean(v)) if v else np.nan

    with open(os.path.join(out, 'report.md'), 'w') as f:
        f.write(f'# Tuning — {key} (bins {bins})\n\nProduction rebuild relerr {relerr:.1e}; per-sample noise sigma '
                f'{sigma:.4g}. Errors are rms inside the lung / object lung mean; mean over bins. deapod = image divided by '
                'the kernel rolloff (2steve/03).\n\n')
        f.write('| setting | deapod | MS | kd2 | σ_k (Δk) | gplb | holes | FWHM all | FWHM bin0 | gain bin0 | blur | alias | '
                'noise | **total** | bias bin | real SNR |\n|' + '---|' * 16 + '\n')
        for n in names:
            for dea in (0, 1):
                r = next(r for r in rows if r['setting'] == n and r['deapod'] == dea)
                f.write(f"| {n} | {dea} | {r['MS']} | {r['kdist0sq']} | {r['sigma_kernel_dk']} | {r['gplb']} | {agg(n, dea, 'holes'):.3f} | "
                        f"{agg(n, dea, 'fwhm_all'):.2f} | {agg(n, dea, 'fwhm_bin0'):.2f} | {agg(n, dea, 'gain_bin0'):.3f} | "
                        f"{agg(n, dea, 'blur'):.3f} | {agg(n, dea, 'alias'):.3f} | {agg(n, dea, 'noise'):.3f} | "
                        f"**{agg(n, dea, 'total'):.3f}** | {agg(n, dea, 'bias_bin'):.3f} | {agg(n, dea, 'real_snr'):.1f} |\n")

    # trade-off: noise vs FWHM (bin 0), point = setting; total error as label
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for n in names:
        x, y, t = agg(n, 0, 'fwhm_bin0'), agg(n, 0, 'noise'), agg(n, 0, 'total')
        c = '#EA580C' if n == 'prod' else ('#2563EB' if n.startswith('gplb') else ('#16A34A' if n.startswith('k') else '#6B7280'))
        ax[0].plot(x, y, 'o', color=c, ms=7)
        ax[0].annotate(n, (x, y), textcoords='offset points', xytext=(5, 4), fontsize=8)
    ax[0].set_xlabel('PSF FWHM, bin 0 (voxels, mean of 3 axes, 4 lung probes)')
    ax[0].set_ylabel('noise sd in lung / lung mean')
    ax[0].set_title('resolution vs noise (blue = gplb, green = kernel, gray = MS)')
    ax[0].grid(alpha=0.25)
    comps = ['blur', 'alias', 'noise']
    x = np.arange(len(names))
    bottom = np.zeros(len(names))
    for comp, c in zip(comps, ('#6B7280', '#2563EB', '#EA580C')):
        v = np.array([agg(n, 0, comp) ** 2 for n in names])
        ax[1].bar(x, v, bottom=bottom, color=c, label=comp + '²')
        bottom += v
    ax[1].plot(x, [agg(n, 0, 'total') ** 2 for n in names], 'k_', ms=16, label='total² (one draw)')
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax[1].set_ylabel('mean-square error / lung mean²')
    ax[1].set_title('error budget (not additive exactly: blur and alias correlate)')
    ax[1].legend(frameon=False, fontsize=8)
    for a_ in ax:
        for s in ('top', 'right'):
            a_.spines[s].set_visible(False)
    fig.suptitle(f'{key}: gas knobs against a known lung object sampled on the real trajectory')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'tradeoff.png'), dpi=120)
    plt.close(fig)

    # montage of the REAL bin-0 image per setting (coronal slice through the lung centroid)
    y = int(np.round(ndi.center_of_mass(lung)[1]))
    nc = 5
    nr = int(np.ceil(len(names) / nc))
    fig, ax = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.4 * nr))
    vmax = np.percentile(montage['prod'][lung], 99)
    for i, a_ in enumerate(ax.ravel()):
        a_.axis('off')
        if i < len(names):
            img = montage[names[i]]
            a_.imshow(img[:, y, :], cmap='gray', vmin=-0.15 * vmax, vmax=vmax * img[lung].mean() / montage['prod'][lung].mean())
            a_.set_title(names[i], fontsize=9)
    fig.suptitle(f'{key}: real data, bin {bins[0]}, coronal y={y} (window scaled to each lung mean)')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'real_montage.png'), dpi=110)
    plt.close(fig)


if __name__ == '__main__':
    main()

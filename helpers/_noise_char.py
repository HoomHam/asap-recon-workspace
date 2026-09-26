"""Background-noise character of Tyger 16-bin outputs (b44 production tree).
Answers: Gaussian? stationary? spatially correlated? per-bin sigma? inter-bin corr?"""
import sys
import numpy as np
import scipy.io as sio
from scipy import ndimage as ndi
from scipy.stats import skew, kurtosis

ROOT = '/Volumes/HoomHamExt/AIkill_Dynamic_b44'


def corners(n=100, c=12):
    m = np.zeros((n, n, n), bool)
    for x in (slice(0, c), slice(n - c, n)):
        for y in (slice(0, c), slice(n - c, n)):
            for z in (slice(0, c), slice(n - c, n)):
                m[x, y, z] = True
    return m


def acf(vol4, mask, lags=(1, 2, 3, 4, 6)):
    """normalized autocorrelation of noise inside mask, per axis, pooled over bins"""
    out = {}
    for ax in range(3):
        r = []
        for L in lags:
            a = np.roll(vol4, -L, axis=ax + 1)
            mk = mask & np.roll(mask, -L, axis=ax)
            x = vol4[:, mk]; y = a[:, mk]
            x = x - x.mean(1, keepdims=True); y = y - y.mean(1, keepdims=True)
            r.append(float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum())))
        out['xyz'[ax]] = r
    return out


def run(key):
    m = sio.loadmat(f'{ROOT}/{key}/d/recon.mat')
    g = m['gas_phase'].astype(np.float64)
    nb = g.shape[0]
    # lung/body support from the bin-mean magnitude, generously dilated
    gm = np.abs(g).mean(0)
    cm = corners()
    s0 = g[:, cm].std()
    sup = ndi.binary_dilation(gm > 5 * s0 / np.sqrt(nb) * 3, iterations=8)
    bg = ~sup
    # far background = corners minus support; "near" = shell 8-16 vox outside support
    far = cm & bg
    near = ndi.binary_dilation(sup, iterations=8) & ~sup
    print(f'\n=== {key}  (nch from b: {m["calcb_b"].shape[0]})  bg voxels far={far.sum()} near={near.sum()}')
    v = g[:, far].ravel()
    print(f'gas_phase (signed) far-bg: mean/sd={v.mean()/v.std():+.3f} skew={skew(v):+.3f} exkurt={kurtosis(v):+.3f}')
    if 'gas_phase_magnitude' in m:
        mg = m['gas_phase_magnitude'].astype(np.float64)[:, far].ravel()
        print(f'gas_phase_magnitude far-bg: mean/sd={mg.mean()/mg.std():.3f}  (Rayleigh 1.913, Gaussian-0 ~0)')
    sig = g[:, far].std(1)
    print('sigma per bin (far bg) rel. to median:', np.round(sig / np.median(sig), 2).tolist())
    # stationarity: sigma per corner, and near-lung shell vs corners (pooled over bins)
    cs = []
    n, c = 100, 12
    for x in (slice(0, c), slice(n - c, n)):
        for y in (slice(0, c), slice(n - c, n)):
            for z in (slice(0, c), slice(n - c, n)):
                mk = np.zeros_like(far); mk[x, y, z] = True; mk &= bg
                if mk.sum() > 200:
                    cs.append(g[:, mk].std())
    cs = np.array(cs)
    print(f'sigma across 8 corners: min/max rel = {cs.min()/np.median(cs):.2f}/{cs.max()/np.median(cs):.2f};'
          f' near-support shell / corners = {g[:, near].std()/np.median(cs):.2f}')
    # spatial autocorrelation
    print('gas spatial ACF (lags 1,2,3,4,6 vox):', {k: np.round(r, 2).tolist() for k, r in acf(g, far).items()})
    # inter-bin correlation of noise (circular)
    x = g[:, far]; x = x - x.mean(1, keepdims=True)
    ib = []
    for d in range(1, 6):
        y = np.roll(x, -d, 0)
        ib.append(float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum())))
    print('gas inter-bin noise corr Δ1..5:', np.round(ib, 2).tolist())
    # dissolved real/imag
    dr = m['dissolved_phase_real'].astype(np.float64); di = m['dissolved_phase_imag'].astype(np.float64)
    for nm, d in (('dis_real', dr), ('dis_imag', di)):
        vv = d[:, far].ravel()
        print(f'{nm} far-bg: mean/sd={vv.mean()/vv.std():+.3f} skew={skew(vv):+.3f} exkurt={kurtosis(vv):+.3f};'
              f' ACF x,y,z lag1/2/4 =', {k: np.round(np.array(r)[[0, 1, 3]], 2).tolist() for k, r in acf(d, far).items()})
    rr = dr[:, far].ravel(); ii = di[:, far].ravel()
    print(f'corr(dis_real, dis_imag) bg = {np.corrcoef(rr, ii)[0,1]:+.2f}; split={str(m["rbc_tp_separated"][0])}')


for k in sys.argv[1:]:
    run(k)

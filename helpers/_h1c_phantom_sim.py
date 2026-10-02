"""Scratch (h1 challenge): ground-truth phantom, step 2 = k-space data and the truth per phase.

Object and motion field from _h1c_phantom_build.py (phantom/object_v1.npz), parameters from _h1c_phantom_params_v1.json.
Per line: true amplitude a = normalised shared k0 surrogate + smooth random mismatch, quantised into amplitude states;
object warped along SI (lung diluted by the stretch); two virtual channels (1 and a perturbed smooth ratio);
NUFFT from the 1.75 mm grid at the real trajectory on a 4x finer time grid; receiver filter; decimation; unit thermal
noise; back to raw channels with the real whitening. Output has the layout of the XeCS cache.
Truth per phase b: object at the median true amplitude of the lines within 1/64 of cyclic coordinate b/16, band-limited
to the 3.5 mm grid and to the acquired k-sphere, cropped to the 100^3 Tyger frame.
"""
import os, sys, json, time
import numpy as np, finufft
from scipy import ndimage as ndi
from scipy.interpolate import CubicSpline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY

G = DY.G; G2 = tuple(2 * g for g in G); DXF = H.DX / 2
PD = os.path.join(H.EXT, 'phantom')
PAR = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_h1c_phantom_params_v1.json')))
VAR = sys.argv[1] if len(sys.argv) > 1 else 'v1'                        # v1b: lung parenchyma decays along the readout (T2*)
T2S = 1.5e-3 if VAR == 'v1b' else None
rng = np.random.default_rng(PAR['seed'] + 1)
O = np.load(os.path.join(PD, 'object_v1.npz'))
obj, phase, gf, dgdz, lungf = O['obj'], O['phase'], O['g'], O['dgdz'], O['lung']
parf = lungf & ~O['vessels']                                            # parenchyma voxels (v1b: these decay)
caud = int(O['caud']); A = PAR['displacement_mm']
D = H.load(); k, tr, vol = D['k'], D['tr'], D['vol']; L = vol.size
wh = np.load(os.path.join(H.EXT, 'whiten.npz')); W, U = wh['W'], wh['U']

# ---- true amplitude per line ----
a_sh = (vol - vol.min()) / (vol.max() - vol.min())
mm = ndi.gaussian_filter1d(rng.normal(size=L + 2000), PAR['surrogate_mismatch_corr_s'] / tr)[1000:-1000]
mm *= PAR['surrogate_mismatch_sd'] / mm.std()
a_true = np.clip(a_sh + mm, 0, None)
NS = PAR['amplitude_states']
lev = np.linspace(0, a_true.max(), NS)
st = np.rint(a_true / lev[1]).astype(int).clip(0, NS - 1)
print('true amplitude: max %.3f, corr with the shared surrogate %.4f, state step %.2f mm' % (a_true.max(), np.corrcoef(a_true, a_sh)[0, 1], lev[1] * A))


def warp(a, part=None):
    """part: None = whole object, 'par' = parenchyma only, 'rest' = everything else."""
    u = a * A * gf / DXF * caud
    zi = np.arange(G2[2], dtype=np.float32)[None, None, :] - u
    z0 = np.clip(np.floor(zi).astype(np.int32), 0, G2[2] - 2); f = (zi - z0).astype(np.float32)
    src = obj * np.where(lungf, 1 / (1 + a * A * dgdz), 1).astype(np.float32)
    if part == 'par':
        src = src * parf
    elif part == 'rest':
        src = src * ~parf
    return np.take_along_axis(src, z0, 2) * (1 - f) + np.take_along_axis(src, z0 + 1, 2) * f


# ---- second channel: perturbed smooth ratio on the fine grid ----
S = np.load(os.path.join(H.EXT, 'sens_g.npy')); R = S[1] / S[0]
pert = ndi.gaussian_filter(rng.normal(size=G), 8) + 1j * ndi.gaussian_filter(rng.normal(size=G), 8)
pert *= PAR['sensitivity_perturbation'] / np.sqrt((np.abs(pert) ** 2).mean())
Rt = R * (1 + pert)
coords = np.meshgrid(*[np.arange(n) / 2.0 for n in G2], indexing='ij')
Rf = (ndi.map_coordinates(Rt.real, coords, order=1, mode='nearest') + 1j * ndi.map_coordinates(Rt.imag, coords, order=1, mode='nearest')).astype(np.complex64)
del coords

# ---- sampling ----
OS = 4; sf = np.arange(0, H.NPTS - 1 + 1e-9, 1 / OS); nf = sf.size
kf = CubicSpline(np.arange(H.NPTS), k, axis=1)(sf)                      # (832, nf, 3) cycles/mm
z = np.load(os.path.join(H.EXT, 'shared', 'receiver_filter_v1.npz')); PADF = 128
ff = np.fft.fftfreq(nf + 2 * PADF, 1 / OS)
Hf = np.interp(ff, z['f'], z['H'].real) + 1j * np.interp(ff, z['f'], z['H'].imag)
keep = np.arange(H.NPTS) * OS
v = np.zeros((2, L, H.NPTS), np.complex64)
t0 = time.time()
for j in range(NS):
    lines = np.where(st == j)[0]
    if lines.size == 0:
        continue
    kap = (2 * np.pi * DXF * kf[lines % H.NILV].reshape(-1, 3)).astype(np.float32)
    xs = [np.ascontiguousarray(kap[:, i]) for i in range(3)]
    nu = lambda o_: finufft.nufft3d2(xs[0], xs[1], xs[2], o_, isign=-1, eps=1e-5, upsampfac=1.25).reshape(lines.size, nf)
    if T2S is None:
        comps = [((warp(lev[j]) * phase / 8).astype(np.complex64), None)]   # /8: fine voxels per 3.5 mm voxel
    else:
        dec = np.exp(-sf * H.DWELL / T2S).astype(np.float32)[None, :]       # from the first readout sample
        comps = [((warp(lev[j], 'rest') * phase / 8).astype(np.complex64), None), ((warp(lev[j], 'par') * phase / 8).astype(np.complex64), dec)]
    for c in range(2):
        y = 0
        for o, dk in comps:
            yc = nu(o if c == 0 else o * Rf)
            y = y + (yc if dk is None else yc * dk)
        yp = np.pad(y, ((0, 0), (PADF, PADF)), mode='edge')
        v[c, lines] = np.fft.ifft(np.fft.fft(yp, axis=1) * Hf[None, :], axis=1)[:, PADF + keep]
    if j % 8 == 7:
        print(f'  state {j+1}/{NS} ({time.time()-t0:.0f} s)', flush=True)
clean = v.copy()
v += ((rng.normal(size=v.shape) + 1j * rng.normal(size=v.shape)) / np.sqrt(2)).astype(np.complex64)
raw = np.einsum('ab,bls->als', np.linalg.inv(W) @ U, v).astype(np.complex64)
np.savez(os.path.join(PD, f'cache_phantom_{VAR}.npz'), raw=raw, k=k, tr=tr, vol=vol.astype(np.float32))
np.savez(os.path.join(PD, f'phantom_{VAR}_internal.npz'), a_true=a_true, a_shared=a_sh, state=st, levels=lev, v_clean_ch0=clean[0])

# ---- truth per phase ----
c = H.cyc_coord(vol)
kx = [np.fft.fftshift(np.fft.fftfreq(n, H.DX)) for n in G]
sph = (kx[0][:, None, None] ** 2 + kx[1][None, :, None] ** 2 + kx[2][None, None, :] ** 2) <= 0.110 ** 2
cc = [g // 2 for g in G]


def to_grid(of, sphere):
    X = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(of)))[cc[0]:cc[0] + G[0], cc[1]:cc[1] + G[1], cc[2]:cc[2] + G[2]] / 8
    if sphere:
        X = X * sph
    return np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X)))


tg, ts, ab = [], [], []
for b in range(16):
    d = np.abs(c - b / 16); d = np.minimum(d, 1 - d)
    a_b = float(np.median(a_true[(d < 1 / 64) & (np.arange(L) >= H.NILV)]))
    of = (warp(a_b) * phase).astype(np.complex64)
    tg.append(DY.crop(to_grid(of, False)).astype(np.complex64)); ts.append(DY.crop(to_grid(of, True)).astype(np.complex64)); ab.append(a_b)
    if b == 0:
        np.save(os.path.join(PD, f'truth_{VAR}_fine_phase0.npy'), np.abs(of).astype(np.float32))
np.savez(os.path.join(PD, f'truth_{VAR}.npz'), truth_grid=np.stack(tg), truth_ksphere=np.stack(ts), amplitude=np.array(ab), displacement_mm=np.array(ab) * A)
print('phase amplitudes (x %.0f mm = SI shift at the lung base):' % A, np.round(ab, 3))
print('done %.0f s' % (time.time() - t0))

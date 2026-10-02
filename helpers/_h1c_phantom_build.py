"""Scratch (h1 challenge): ground-truth phantom, step 1 = the object and its motion field (no k-space yet).

Reads the difficulty parameters from phantom_params_v1.json (agreed with XeCS before data are generated).
Fine grid: LR 352 x AP 240 x SI 320 at 1.75 mm (616 x 420 x 560 mm). Anatomy = magnitude of the real full-field
least-squares image (texture kept), sharp body edge on the fine grid; lungs replaced by a smooth parenchyma level
with an AP gravity gradient plus tapered tubes with graded intensity; smooth object phase from the real image.
Motion: displacement along SI only, u = a(line) * A * g(r); g ramps from 0 at the lung apex to 1 at the lung base,
stays 1 below it and tapers to 0 further down; it is 0 in the body wall and the arms (sliding at the pleura).
Lung parenchyma and vessels are diluted by the local stretch (1 / (1 + du/dz)).
Writes phantom/object_v1.npz on Ext and a review figure.
"""
import os, sys, json
import numpy as np
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY

G = DY.G; G2 = tuple(2 * g for g in G); DXF = H.DX / 2
PD = os.path.join(H.EXT, 'phantom'); os.makedirs(PD, exist_ok=True)
PAR = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '_h1c_phantom_params_v1.json')))
rng = np.random.default_rng(PAR['seed'])


def upsample(x):
    X = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x)))
    Xp = np.zeros(G2, X.dtype); c = [g // 2 for g in G]
    Xp[c[0]:c[0] + G[0], c[1]:c[1] + G[1], c[2]:c[2] + G[2]] = X
    return np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Xp))) * 8


sup = np.load(os.path.join(H.EXT, 'support_g.npy'))
fx = os.path.join(H.EXT, 'static_g_ls_ch0_v2.npy')
if os.path.exists(fx):
    x = np.load(fx).astype(np.complex128)
else:                                                                   # principal channel, receiver filter, thermal units
    ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)[0]; cnt = np.load(os.path.join(H.EXT, 'cnt_train.npy'))
    op = H.RxNufft(H.load()['k'], n=G, nt=1, single=False, eps=1e-5)
    dcf = np.load(os.path.join(H.EXT, 'dcf_g.npy')) * (cnt > 0)[:, None]; Wt = dcf / dcf.sum()
    rhs = sup * op.adj(Wt * ybar); n0 = lambda z: sup * op.adj(Wt * op.fwd(z)); sc = np.vdot(rhs, n0(rhs)).real / np.vdot(rhs, rhs).real
    x = H.cg(lambda z: n0(z) + 1e-3 * sc * z, rhs, it=25)
    np.save(fx, x.astype(np.complex64))
# ---- coarse-grid masks (3.5 mm) ----
mag = np.abs(x); ms = ndi.gaussian_filter(mag, 1.0)
ref = np.median(ms[sup & (ms > np.percentile(ms[sup], 60))])            # bright-tissue reference level
body = ndi.gaussian_filter(mag, 1.5) > 0.30 * ref
lab, n = ndi.label(body); sz = ndi.sum(body, lab, range(1, n + 1)); body = lab == (np.argmax(sz) + 1)
hull = body.copy()
hull = ndi.binary_closing(np.pad(hull, 14), iterations=12)[14:-14, 14:-14, 14:-14]
for i in range(G[2]): hull[:, :, i] = ndi.binary_fill_holes(hull[:, :, i])
dark = hull & (ndi.gaussian_filter(mag, 1.5) < 0.42 * ref) & ndi.binary_erosion(hull, iterations=4)
lab, n = ndi.label(dark); sz = ndi.sum(dark, lab, range(1, n + 1)); order = np.argsort(sz)[::-1]
lung = np.isin(lab, order[:2] + 1) if sz[order[1]] > 0.25 * sz[order[0]] else (lab == order[0] + 1)
lung = ndi.binary_opening(ndi.binary_closing(lung, iterations=2), iterations=1)
print('reference tissue level', ref, '| hull voxels', hull.sum(), '| lung voxels', lung.sum(), f'= {lung.sum()*H.DX**3/1e6:.2f} L')
# head / feet direction along SI (axis 2): the neck end has the smaller body cross-section
zs = np.where(lung.any((0, 1)))[0]; area = hull.sum((0, 1))
head_low = area[max(zs.min() - 25, 0):zs.min()].mean() < area[zs.max():zs.max() + 25].mean()
caud = +1 if head_low else -1                                           # +1: caudal = increasing SI index
z_apex, z_base = (zs.min(), zs.max()) if caud > 0 else (zs.max(), zs.min())
# posterior direction along AP (axis 1): the spine / back is the side where the lungs reach closest to the body surface in bulk
ys = np.arange(G[1]); lung_y = (lung.sum((0, 2)) * ys).sum() / lung.sum(); body_y = (hull.sum((0, 2)) * ys).sum() / hull.sum()
post = +1 if lung_y > body_y else -1                                    # lungs sit posterior of the body centroid
print(f'caudal = {"+" if caud > 0 else "-"}SI index, apex z {z_apex}, base z {z_base} | posterior = {"+" if post > 0 else "-"}AP index (lung centroid {lung_y:.1f}, body {body_y:.1f})')
# ---- motion field on the coarse grid: g(r) in [0, 1] ----
zc = np.arange(G[2], dtype=float)
t = (zc - z_apex) * caud / abs(z_base - z_apex)                          # 0 at apex, 1 at base, > 1 below
below = (t - 1) * abs(z_base - z_apex) * H.DX                            # mm below the lung base
gz = np.where(t < 0, 0, np.where(t <= 1, t, np.where(below < PAR['abdomen_plateau_mm'], 1.0,
     0.5 * (1 + np.cos(np.pi * np.clip((below - PAR['abdomen_plateau_mm']) / PAR['abdomen_taper_mm'], 0, 1))))))
edt = ndi.distance_transform_edt(hull) * H.DX
torso = edt > PAR['arm_radius_mm']                                       # opening: removes arms / neck
torso = ndi.distance_transform_edt(~torso) * H.DX < PAR['arm_radius_mm']
wall = np.clip((edt - PAR['wall_mm']) / PAR['wall_transition_mm'], 0, 1); wall = wall * wall * (3 - 2 * wall)
g = ndi.gaussian_filter(gz[None, None, :] * wall * (torso & hull), 1.5)
# ---- fine grid ----
xu = upsample(x)
magf = np.abs(xu).astype(np.float32)
ph = ndi.gaussian_filter((x / (np.abs(x) + 1e-6 * ref)).real * sup, 2.5) + 1j * ndi.gaussian_filter((x / (np.abs(x) + 1e-6 * ref)).imag * sup, 2.5)
phase = np.angle(ph).astype(np.float32)
z2 = lambda a, order=1: ndi.zoom(a.astype(np.float32), 2, order=order, mode='nearest', grid_mode=False, prefilter=order > 1)[:G2[0], :G2[1], :G2[2]]
up = lambda a, order=1: ndi.map_coordinates(a.astype(np.float32), np.meshgrid(*[np.arange(n) / 2.0 for n in G2], indexing='ij'), order=order, mode='nearest')
bodyf = up(ndi.gaussian_filter(mag, 1.0)) > 0.30 * ref
hullf = up(hull.astype(np.float32)) > 0.5
lungf = up(lung.astype(np.float32)) > 0.5
bodyf = ndi.binary_closing(bodyf, iterations=2) & hullf                 # sharp skin edge on the fine grid; pinholes closed, air gaps kept
bodyf |= lungf
gf = up(g); phasef_re = up(np.cos(phase)); phasef_im = up(np.sin(phase))
phf = np.exp(1j * np.arctan2(phasef_im, phasef_re)).astype(np.complex64)
obj = np.where(bodyf, magf, 0).astype(np.float32)
obj[bodyf] = np.maximum(obj[bodyf], 0.25 * ref)                          # no holes at the fold-in speckle minima
# lungs: parenchyma with gravity gradient (denser posterior)
yy = (np.arange(G2[1]) - G2[1] / 2) * DXF
ly = np.where(lungf.any((0, 2)))[0]; ymid = 0.5 * (yy[ly.min()] + yy[ly.max()]); yhalf = 0.5 * (yy[ly.max()] - yy[ly.min()])
grav = 1 + PAR['gravity_gradient'] * post * (yy - ymid) / yhalf
par = (PAR['parenchyma_level'] * ref * grav)[None, :, None] * np.ones(G2, np.float32)
# vessels: tapered tubes from a hilum point of each lung
ves = np.zeros(G2, np.float32)
labf, nl = ndi.label(lungf)
xmid = G2[0] / 2
for li in range(1, nl + 1):
    pts = np.argwhere(labf == li)
    if len(pts) < 5000: continue
    cen = pts.mean(0)
    med = pts[np.argsort(np.abs(pts[:, 0] - xmid))[:max(len(pts) // 50, 10)]]          # most medial 2 %
    hil = med[np.argsort(np.abs(med[:, 2] - cen[2]) + np.abs(med[:, 1] - cen[1]))[:20]].mean(0)
    for _ in range(PAR['vessels_per_lung']):
        end = pts[rng.integers(len(pts))].astype(float)
        ctl = 0.5 * (hil + end) + rng.normal(0, 0.15 * np.linalg.norm(end - hil), 3)
        r0 = rng.uniform(PAR['vessel_radius_mm'][0], PAR['vessel_radius_mm'][1])
        npt = int(np.linalg.norm(end - hil) * 1.5) + 2
        for s in np.linspace(0, 1, npt):
            p = (1 - s) ** 2 * hil + 2 * s * (1 - s) * ctl + s ** 2 * end
            r = (r0 * (1 - s) + PAR['vessel_end_radius_mm'] * s) / DXF                       # fine voxels
            amp = ref * (PAR['vessel_contrast'][0] + (PAR['vessel_contrast'][1] - PAR['vessel_contrast'][0]) * min(r * DXF / PAR['vessel_radius_mm'][1], 1))
            R = int(np.ceil(r + 1)); c = np.round(p).astype(int)
            sl = tuple(slice(max(c[a] - R, 0), min(c[a] + R + 1, G2[a])) for a in range(3))
            gr = np.meshgrid(*[np.arange(s_.start, s_.stop) - p[a] for a, s_ in enumerate(sl)], indexing='ij')
            dist = np.sqrt(sum(q ** 2 for q in gr))
            ves[sl] = np.maximum(ves[sl], amp * np.clip(r - dist + 0.5, 0, 1))
ves *= lungf
lung_img = np.maximum(par, ves)
obj = np.where(lungf, lung_img, obj).astype(np.float32)
dgdz = np.gradient(gf, DXF, axis=2).astype(np.float32) * caud          # d g / d(caudal mm)
np.savez(os.path.join(PD, 'object_v1.npz'), obj=obj, phase=phf, g=gf.astype(np.float32), dgdz=dgdz, lung=lungf, body=bodyf, vessels=(ves > 0),
         caud=caud, post=post, ref=ref, z_apex=z_apex, z_base=z_base)
json.dump(dict(ref_tissue_level=float(ref), caudal_sign=int(caud), posterior_sign=int(post), lung_litres=float(lungf.sum() * DXF ** 3 / 1e6),
               vessel_fraction_of_lung=float((ves > 0).sum() / lungf.sum()), fine_grid=G2, dx_mm=DXF), open(os.path.join(H.OUT, 'phantom_object_v1.json'), 'w'), indent=1)
print('lung %.2f L, vessels fill %.1f %% of the lung' % (lungf.sum() * DXF ** 3 / 1e6, 100 * (ves > 0).sum() / lungf.sum()))
# ---- review figure: end expiration (a = 0) and full inspiration (a = 1) ----
A = PAR['displacement_mm']
def state(a):
    u = a * A * gf / DXF * caud                                          # fine voxels along SI
    zi = np.arange(G2[2], dtype=np.float32)[None, None, :] - u
    z0 = np.clip(np.floor(zi).astype(np.int32), 0, G2[2] - 2); f = (zi - z0).astype(np.float32)
    dil = np.where(lungf, 1 / (1 + a * A * dgdz), 1).astype(np.float32)
    src = obj * dil
    return np.take_along_axis(src, z0, 2) * (1 - f) + np.take_along_axis(src, z0 + 1, 2) * f
s0, s1 = state(0.0), state(1.0)
cx, cy, cz = G2[0] // 2, int(np.where(lungf.any((0, 2)))[0].mean()), int(np.where(lungf.any((0, 1)))[0].mean())
xl = int(np.argwhere(lungf)[:, 0].min() + 0.25 * (np.argwhere(lungf)[:, 0].max() - np.argwhere(lungf)[:, 0].min()))
fig, ax = plt.subplots(2, 4, figsize=(22, 11)); vm = 1.3 * ref
for r, (s, nm) in enumerate(((s0, 'a = 0 (end expiration)'), (s1, f'a = 1 (deepest inspiration, {A} mm)'))):
    ax[r, 0].imshow(s[:, cy, :].T, cmap='gray', vmin=0, vmax=vm, origin='lower'); ax[r, 0].set_title(f'coronal, {nm}')
    ax[r, 1].imshow(s[xl, :, :].T, cmap='gray', vmin=0, vmax=vm, origin='lower'); ax[r, 1].set_title('sagittal through one lung')
    ax[r, 2].imshow(s[:, :, cz].T, cmap='gray', vmin=0, vmax=vm, origin='lower'); ax[r, 2].set_title('axial, mid lung')
ax[0, 3].imshow(gf[:, cy, :].T, cmap='viridis', vmin=0, vmax=1, origin='lower'); ax[0, 3].contour(lungf[:, cy, :].T, [0.5], colors='w', linewidths=.5); ax[0, 3].set_title('motion weight g (coronal)')
ax[1, 3].imshow((s1 - s0)[:, cy, :].T, cmap='RdBu', vmin=-ref, vmax=ref, origin='lower'); ax[1, 3].set_title('inspiration - expiration (coronal)')
for a_ in ax.ravel(): a_.axis('off')
fig.tight_layout(); fig.savefig(os.path.join(H.OUT, 'phantom_object_v1.png'), dpi=50)
